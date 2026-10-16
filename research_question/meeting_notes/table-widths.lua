-- Size table columns to their content instead of Pandoc's equal or dash-based widths.
--
-- Short columns get exactly the width of their longest cell (so they never wrap);
-- the remaining width is shared among the long columns in proportion to their length.
-- Tables that fit on one line are left at natural width. A table with an explicit
-- tbl-colwidths attribute is left alone.
--
-- In PDF output every table is also set in \small with tighter rows, so that a
-- section's table, figure and text fit on one page.
--
-- Metadata: table-width-chars (default 98, for \small text) is the number of characters
-- that fit in the text width; raise it if short columns get too much room, lower it if they wrap.

local line_chars = 98
local padding = 3

local function text_length(blocks)
  local n = 0
  pandoc.walk_block(pandoc.Div(blocks), {
    Str = function(s) n = n + utf8.len(s.text) end,
    Space = function() n = n + 1 end,
    Code = function(c) n = n + utf8.len(c.text) end,
    Math = function(m)
      -- Count a TeX command such as \beta as one glyph and drop markup characters.
      local t = m.text:gsub("\\[%a]+", "x"):gsub("[{}_^%s]", "")
      n = n + utf8.len(t)
    end,
  })
  return n
end

local function column_lengths(tbl)
  local lengths = {}
  for i = 1, #tbl.colspecs do lengths[i] = 0 end
  local function scan(rows)
    for _, row in ipairs(rows) do
      for j, cell in ipairs(row.cells) do
        if j <= #lengths then
          lengths[j] = math.max(lengths[j], text_length(cell.contents))
        end
      end
    end
  end
  scan(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do
    scan(body.head)
    scan(body.body)
  end
  scan(tbl.foot.rows)
  return lengths
end

-- Water-filling: fix every column that fits in its fair share of what is left.
local function widths_for(lengths)
  local needed, open = {}, {}
  for j, len in ipairs(lengths) do
    needed[j] = len + padding
    open[j] = true
  end
  local widths = {}
  local remaining = line_chars
  local changed = true
  while changed do
    changed = false
    local n_open = 0
    for _ in pairs(open) do n_open = n_open + 1 end
    if n_open == 0 then break end
    local share = remaining / n_open
    for j in pairs(open) do
      if needed[j] <= share then
        widths[j] = needed[j] / line_chars
        remaining = remaining - needed[j]
        open[j] = nil
        changed = true
      end
    end
  end
  local open_total = 0
  for j in pairs(open) do open_total = open_total + needed[j] end
  for j in pairs(open) do
    widths[j] = (remaining / line_chars) * needed[j] / open_total
  end
  return widths
end

local function Table(tbl)
  if tbl.attr.attributes["tbl-colwidths"] then return nil end
  local lengths = column_lengths(tbl)
  local total = 0
  for _, len in ipairs(lengths) do total = total + len + padding end
  local widths = nil
  if total > line_chars then widths = widths_for(lengths) end
  for j, spec in ipairs(tbl.colspecs) do
    tbl.colspecs[j] = { spec[1], widths and widths[j] or nil }
  end
  return tbl
end

local small_tables = [[
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\small\renewcommand{\arraystretch}{0.9}}
]]

local function Meta(meta)
  if meta["table-width-chars"] then
    line_chars = tonumber(pandoc.utils.stringify(meta["table-width-chars"])) or line_chars
  end
  if FORMAT:match("latex") then
    local includes = meta["header-includes"]
    if includes == nil then
      includes = pandoc.MetaList({})
    elseif includes.t ~= "MetaList" then
      includes = pandoc.MetaList({ includes })
    end
    includes:insert(pandoc.MetaBlocks({ pandoc.RawBlock("latex", small_tables) }))
    meta["header-includes"] = includes
    return meta
  end
end

return {
  { Meta = Meta },
  { Table = Table },
}
