# Plan: report comparing the two LSPIM engines (3GEE vs Logistic)

Agreed 2026-10-16 (grilling session). Branch `docs/lspim-engine-comparison-report` off `main`. The
user starts the implementation and opens the PR. Commit the `.qmd` and the rendered PDF.

Goal: a short PDF report showing that the original LSPIM engine (geessbin, three GEE fits) and the
improved engine (one penalized logistic fit plus three FW sandwiches) give the same results, and
explaining where and why they differ. Source run: `scripts/run_lspim_engine_comparison.R`
(12 scenarios, 5000 replicates, alpha = 0.05).

Ground rules:

- Multi-line `Rscript -e` segfaults on this machine: write R code to a file and run
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" file.R`.
- Quarto may not be on PATH. It is also at
  `C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe`.
- View PDF pages with `pdftools::pdf_convert()` (PNG) and `pdftools::pdf_text()`. Do not try `pdftoppm`.
- Do not add a `_quarto.yml` anywhere under `reports/`: `here::here()` would take it as the project
  root and then not find `results/`.

## Decisions

- **File.** `reports/drafts/2026-10-16-lspim-engine-comparison.qmd`, rendered to PDF next to it.
- **Header.** Copy the YAML header of `research_question/meeting_notes/update_261019.qmd`
  (A4, `geometry: margin=2cm`, `fig-width: 7`, `fig-height: 4`, `editor: visual`, `echo`,
  `warning` and `message` FALSE). Filter by relative path:
  `filters: - ../../research_question/meeting_notes/table-widths.lua` (no copy of the filter).
  Setup chunk: `knitr::opts_chunk$set(dev = "png", dpi = 300)`, data paths via `here::here()`.
- **Naming.** Long names "3GEE + FW" (three geessbin PGEE fits with Ford-Westgate SEs) and
  "Logistic + 3FW" (one penalized logistic fit with three Ford-Westgate sandwiches). Short labels
  **3GEE** and **Logistic** in tables and figures. Never "in-house"; do not use the engine codes
  `geessbin` / `pgee_fw` as labels (they may appear once in the setup, to link to the code).
- **Data.**
  - Comparison object:
    `results/data/lspim_engine_comparison/lspim_engine_comparison_c1feac084af6cb59_22612da3a9139167.rds`
    (same file also at `results/lspim_engine_comparison/`; read the `results/data/` copy).
    Elements: `side_by_side`, `decision_agreement`, `failure_summary`, `one_sided_failures`,
    `aggregation` (has `time_median_seconds`), `scenarios`, `analysis_run_hashes`
    (`pgee_fw = c1feac084af6cb59`, `geessbin = 22612da3a9139167`), `n_simulations`, `alpha`.
  - Per-replicate rows: `results/data/<hash>/analysis_combined_convenience.rds`, element `$results`
    (columns used: `scenario_id`, `sim_id`, `engine`, `status`, `convergence_status`,
    `interaction_rejected`, `warning_message`, `error_message`). Needed for the common-replicate
    rates and the warning table. Take the hashes from `analysis_run_hashes`, not hard-coded.
- **No Monte Carlo SEs** (user decision).

## Sections

1. **Setup.** Both methods in a few sentences:
   - Under the independence working correlation the GEE estimating equations are the logistic
     score equations, so geessbin's PGEE point estimate is a Firth-type penalized logistic
     regression on the 0 / 0.5 / 1 pseudo-scores, with the penalty weighted by the Pearson scale
     phi (re-estimated each iteration). The clustering does not enter the beta equation, so the
     three geessbin fits give the same beta; only the sandwiches (clustered by C1, C2 and the pair
     C3) differ. Logistic + 3FW fits beta once and computes the three FW sandwiches. Not a plain
     `glm()`: the penalty term is what makes it PGEE.
   - The test (Holm on the treatment-by-visit contrasts) is the same for both.
   - Scenarios: n = 10, 20, 50; three-observations-minimum dropout; linear beta3 = 0.22 and 0;
     log large crossing (beta1 = -2.42, beta3 = 1.95) and log null (identical arms); 12 scenarios,
     5000 replicates each, alpha = 0.05. The same generated datasets for both methods
     (generation run `e64245fe0a9d0284`).
2. **Decision agreement.** Table from `decision_agreement`: per scenario n_pairs, both reject,
   only 3GEE rejects, only Logistic rejects, neither, agreement. Expected: every off-diagonal 0,
   agreement 1 in all 12 scenarios.
3. **Rejection rates** (type I error for the nulls, power for the effects), table + figure. Per
   scenario and method, three columns:
   - **Own replicates:** each method's rate over its own replicates with a decision (as in the
     main simulation; from `side_by_side`).
   - **Common replicates:** each method's rate over the replicates where both gave a decision,
     **computed separately per method from the per-replicate rows** (user: show that they are
     actually identical, do not rely on construction). Both columns shown side by side.
   - **3GEE where only Logistic failed:** 3GEE's rate on those replicates (n in brackets).
   Figure: rate vs n, colour by method, facets by effect, own-replicate rates (dotted line at 0.05
   for the null facets). Expected pattern: identical at n = 20/50; at n = 10, 3GEE higher
   (0.153 vs 0.119 linear null, 0.146 vs 0.116 log null) because of the replicates Logistic
   refuses.
4. **Failures and warnings.**
   - Table per scenario and method: errors, non-convergence (from `failure_summary`), replicates
     where only one method failed.
   - Failure types (from `one_sided_failures`): Logistic leverage-1 FW error, clustering C1
     (1954) and C2 (1995); Logistic all-zero design columns (168, late `trt_visit` columns with no
     comparable pairs); 3GEE "infinite or missing values in 'x'" (62, all at n = 10), where
     Logistic returned a decision flagged not converged ("fitted probabilities numerically 0 or 1").
     Count them from the data, do not hard-code.
   - 3GEE's rejection rate on the replicates where Logistic failed, split by failure type
     (null scenarios: ~0.26-0.36 for leverage-1, ~0.41 for all-zero columns).
   - Brief warning table: share of replicates per method and scenario with the nearest-PSD warning
     ("Combined V has negative eigenvalues...") and with "fitted probabilities numerically 0 or 1".
     Match on substrings of `warning_message` (Logistic messages can combine both with " | ").
5. **Computation time.** Median seconds per replicate per method and n (from `aggregation`,
   `time_median_seconds`), table + log-scale figure. Expected: 3GEE ~0.8 / 2.7 / 41 s,
   Logistic ~0.009 / 0.014 / 0.06 s.
6. **Conclusion.** Short bullet points only:
   - Identical decisions on every replicate where both methods give one; identical rates at
     n = 20 and 50.
   - The n = 10 difference comes entirely from replicates where Logistic + 3FW stops (leverage-1
     cluster or all-zero design column) and 3GEE still returns a decision; on those, 3GEE rejects
     far more often than alpha under the null.
   - Earlier single check (2026-10-08 scratch verification, not this run): in such a case geessbin
     reported "converged" with a C2 variance of 3.4e-17 for the affected column. Flag it as one
     example.
   - Logistic + 3FW is several hundred times faster (~650x at n = 50).

## Verification

- Render with quarto, convert every page with `pdftools::pdf_convert()` and check layout (tables fit
  the page width, figures readable, no raw NA / code output).
- Cross-check two numbers in the PDF against the RDS (one agreement count, one median time).
- Commit `.qmd` + PDF with a conventional message (`docs: ...`), attribution trailer as usual.
