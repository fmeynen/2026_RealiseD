# ---- Contrasts and Wald-type tests ------------------------------------------------


# Helper functions

## Aanmaken van de contrast matrices (visit vs. mean; all pairwise or combined)

revision_contrasts <- function(K) {
  meanL<-cbind(matrix(0,K,2),diag(K)-matrix(1/K,K,K))
  vp<-combn(seq_len(K),2);pairL<-matrix(0,ncol(vp),K+2L)
  for(j in seq_len(ncol(vp))) {pairL[j,2L+vp[1,j]]<-1;pairL[j,2L+vp[2,j]]<- -1}
  rownames(meanL)<-paste0("Visit ",seq_len(K)," - coefficient mean")
  rownames(pairL)<-paste0("Visit ",vp[1,]," - Visit ",vp[2,])
  # Full-rank global basis. The maximum families retain ALL redundant rows.
  globalL<-cbind(matrix(0,K-1L,2),cbind(diag(K-1L),-1))
  list(treatment=list(Global=globalL,Holm=meanL,Single_step=meanL,
                      Tukey=pairL,Combined=rbind(meanL,pairL)),
       time=list(Global=matrix(c(1,-1,rep(0,K)),1),
                 Holm=matrix(c(1,-1,rep(0,K)),1),
                 Single_step=matrix(c(1,-1,rep(0,K)),1)))
}

## Berekenen van 1 p-waarde (maximum van de individuele test statistieken)
## die gebruik maakt van multivariate t of N verdeling 

revision_max_p <- function(z,R,df,cfg,seed) {
  m<-length(z);mx<-max(abs(z))
  if(m==1L)return(list(p=if(df==0)2*pnorm(-mx) else 2*pt(-mx,df),
                       integration_error=0,integration_message="",warning=""))
  if(min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)< -1e-8)
    stop("Indefinite contrast correlation; no covariance repair applied")
  old<-get(".Random.seed",envir=.GlobalEnv)
  on.exit(assign(".Random.seed",old,envir=.GlobalEnv),add=TRUE)
  set.seed(seed);warnings<-character()
  val<-withCallingHandlers({
    alg<-mvtnorm::GenzBretz(maxpts=cfg$maxpts,abseps=cfg$abseps,releps=0)
    if(df==0)mvtnorm::pmvnorm(lower=rep(-mx,m),upper=rep(mx,m),corr=R,algorithm=alg) else
      mvtnorm::pmvt(lower=rep(-mx,m),upper=rep(mx,m),df=df,corr=R,algorithm=alg)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  err<-attr(val,"error");if(is.null(err))err<-NA_real_
  msg<-paste(attr(val,"msg"),collapse="; ")
  if(is.finite(err)&&err>cfg$abseps)warnings<-c(warnings,"Integration error exceeds abseps")
  if(nzchar(msg)&&msg!="Normal Completion")warnings<-c(warnings,msg)
  p<-1-as.numeric(val)
  if(length(p)!=1||!is.finite(p)||p< -1e-10||p>1+1e-10)stop("Invalid integrated p-value")
  list(p=min(1,max(0,p)),integration_error=err,integration_message=msg,
       warning=paste(unique(warnings),collapse="; "))
}

## Berekenen van de p-waarden per model en methode
revision_test <- function(model,L,procedure,reference,N,cfg,seed) {
  tryCatch({
    if(!is.null(model$error))stop(model$error)
    b<-as.numeric(model$b);V<-(model$V+t(model$V))/2
    stopifnot(length(b)==ncol(L),all(is.finite(b)),all(is.finite(V)))
    delta<-drop(L%*%b);S<-L%*%V%*%t(L);S<-(S+t(S))/2
    if(any(diag(S)<=0)||any(!is.finite(S)))stop("Invalid contrast covariance")
    df<-if(reference=="Asymptotic")0 else N-1L
    z<-delta/sqrt(diag(S));ierr<-NA_real_;imsg<-warn<-""
    if(procedure=="Global") {
      e<-eigen(S,symmetric=TRUE,only.values=TRUE)$values
      if(min(e)<=1e-10*max(abs(e)))stop("Global contrast covariance not positive definite")
      Q<-sum(forwardsolve(t(chol(S)),delta)^2);r<-nrow(L)
      p<-if(df==0)pchisq(Q,r,lower.tail=FALSE) else pf(Q/r,r,df,lower.tail=FALSE)
      statistic<-if(df==0)Q else Q/r
    } else if(procedure=="Holm") {
      raw<-if(df==0)2*pnorm(-abs(z)) else 2*pt(-abs(z),df)
      p<-min(p.adjust(raw,"holm"));statistic<-max(abs(z))
    } else {
      ans<-revision_max_p(z,cov2cor(S),df,cfg,seed)
      p<-ans$p;ierr<-ans$integration_error;imsg<-ans$integration_message;warn<-ans$warning
      statistic<-max(abs(z))
    }
    data.frame(statistic=statistic,p_value=p,reject=p<cfg$alpha,
               df=if(df==0)Inf else df,status=if(nzchar(warn))"warning" else "ok",
               error=NA_character_,warning=warn,integration_error=ierr,integration_message=imsg)
  },error=function(e)data.frame(statistic=NA_real_,p_value=NA_real_,reject=NA,
                                df=if(reference=="Asymptotic")Inf else N-1L,status="failed",error=conditionMessage(e),
                                warning="",integration_error=NA_real_,integration_message=""))
}


# Toegepast op model, met data in dat. 

model <- list(
  b = as.numeric(coef(mod)),
  V = as.matrix(vcov(mod))
)

N <- length(unique(dat$id))    # Number of subjects
K <- length(unique(dat$time))  # Number of visits

cfg <- list(
  alpha = 0.05,
  maxpts = 100000L,
  abseps = 1e-5
)

# Also initializes .Random.seed, required by your pasted revision_max_p()
set.seed(123)

contrasts <- revision_contrasts(K)

# Asymptotic reference for the original/adjusted PIM comparison.
# Choose "Finite" for the finite-sample t/F analysis.
reference <- "Finite"

results <- do.call(rbind, lapply(names(contrasts), function(hypothesis) {
  
  do.call(rbind, lapply(
    names(contrasts[[hypothesis]]),
    function(procedure) {
      
      ans <- revision_test(
        model = model,
        L = contrasts[[hypothesis]][[procedure]],
        procedure = procedure,
        reference = reference,
        N = N,
        cfg = cfg,
        seed = 123L
      )
      
      cbind(
        data.frame(
          hypothesis = hypothesis,
          procedure = procedure,
          reference = reference
        ),
        ans
      )
    }
  ))
}))

print(
  results[, c(
    "hypothesis", "procedure", "reference",
    "statistic", "p_value", "status", "error"
  )],
  row.names = FALSE
)