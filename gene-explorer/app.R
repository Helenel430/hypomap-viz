# ============================================================
# Romanov_10x — Gene Expression Lookup
# A simple explorer (no UMAP/tSNE): type a gene, get its detection %,
# mean expression (among expressing cells), and expression percentile,
# globally and broken down by cell type.
#
# Loads only the precomputed gene-explorer/gene_stats.rds (built by
# ../compute_gene_stats.R). Run with:  shiny::runApp("gene-explorer")
# ============================================================
library(shiny)
library(DT)

stats <- readRDS("gene_stats.rds")
genes <- stats$genes
total <- stats$total_cells
cts   <- stats$celltypes
# genes detected in >=1 cell — the reference set for the global percentile
n_expressed <- sum(stats$global$n_expressing > 0)
# scale reference for the mean-expression bars: the 99th percentile of all
# gene x cell-type mean-when-expressed values (robust to outliers like Malat1).
# Values above it cap the bar and are flagged as over-scale.
mean_scale <- as.numeric(quantile(
  stats$byct$mean_expressing[stats$byct$n_expressing > 0], 0.99, na.rm = TRUE))

# Horizontal bar chart drawn as plain HTML/CSS — no R graphics device needed,
# so it renders reliably everywhere. Bars use an ABSOLUTE scale (0..maxval),
# with a header row of column titles and the raw expressing-cell count
# (#Expressing, in grey). Calculations are explained in each tab's caption.
htmlBars <- function(labels, values, counts, maxval, fmt = "%.2f", accent = "#3a6ea3",
                     valueTitle = "Value", countTitle = "#Expressing") {
  keep <- !is.na(values)
  labels <- labels[keep]; values <- values[keep]; counts <- counts[keep]
  if (length(values) == 0 || max(values) <= 0)
    return(tags$div(style = "color:#999;font-size:13px",
                    "Not expressed in any cell type."))
  ord <- order(values, decreasing = TRUE)
  labels <- labels[ord]; values <- values[ord]; counts <- counts[ord]
  valW <- "104px"; cntW <- "96px"
  hcell <- "white-space:nowrap"
  header <- tags$div(
    style = "display:flex;align-items:flex-end;font-size:11.5px;font-weight:700;color:#555;border-bottom:1px solid #e3e3e3;padding-bottom:5px;margin-bottom:4px",
    tags$div(style = "flex:1", "Cell type"),
    tags$div(style = sprintf("width:%s;text-align:right;padding-left:10px;%s", valW, hcell),
             valueTitle),
    tags$div(style = sprintf("width:%s;text-align:right;%s", cntW, hcell),
             countTitle))
  rows <- lapply(seq_along(values), function(i) {
    over <- isTRUE(values[i] > maxval)                 # value exceeds the scale?
    w    <- min(100, 100 * values[i] / maxval)
    fill <- tags$div(style = sprintf(
      "width:%.1f%%;height:100%%;border-radius:3px;background:%s", w, accent))
    # striped cap on the right end signals the bar runs off the scale
    cap  <- if (over) tags$div(style = sprintf(paste0(
      "position:absolute;top:0;right:0;height:100%%;width:16px;border-radius:0 3px 3px 0;",
      "background:repeating-linear-gradient(45deg,%s,%s 3px,#fff 3px,#fff 6px)"),
      accent, accent)) else NULL
    track <- tags$div(
      style = "flex:1;background:#eef0f2;border-radius:3px;height:16px;position:relative;overflow:hidden",
      fill, cap)
    valcell <- tags$div(
      style = sprintf("width:%s;text-align:right;padding-left:10px;font-size:13px;color:%s;%s",
                      valW, if (over) accent else "#333", if (over) "font-weight:600" else ""),
      if (over) HTML(paste0("&#9656; ", sprintf(fmt, values[i]))) else sprintf(fmt, values[i]))
    cntcell <- tags$div(
      style = sprintf("width:%s;text-align:right;font-size:12px;color:#999", cntW),
      format(counts[i], big.mark = ","))
    tags$div(style = "margin:9px 0",
      tags$div(style = "font-size:13px;color:#333;margin-bottom:2px", labels[i]),
      tags$div(style = "display:flex;align-items:center", track, valcell, cntcell))
  })
  tags$div(style = "margin-top:6px", header, rows)
}

# render "num / den" as a real stacked fraction (numerator over a rule over
# denominator) so caption formulas read like maths rather than inline text.
frac <- function(num, den) sprintf(paste0(
  "<span style='display:inline-block;vertical-align:middle;text-align:center;margin:0 4px'>",
  "<span style='display:block;padding:0 8px 1px;border-bottom:1.5px solid #888'>%s</span>",
  "<span style='display:block;padding:1px 8px 0'>%s</span></span>"), num, den)

ui <- fluidPage(
  titlePanel("Romanov_10x — Gene Expression Lookup"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      textInput("gene_in", "Gene symbol:", value = "Pomc",
                placeholder = "start typing, e.g. igf"),
      uiOutput("suggestions"),
      actionButton("go", "Look up", class = "btn-primary"),
      # let Enter in the text box trigger the Look up button
      tags$script(HTML(
        "document.addEventListener('keydown',function(e){if(e.key==='Enter'&&document.activeElement&&document.activeElement.id==='gene_in'){var b=document.getElementById('go');if(b)b.click();}});")),
      tags$br(), tags$br(),
      helpText("Type a gene symbol (case-insensitive) and click Look up or press Enter.",
               sprintf(paste("RNA assay, log-normalized; %s P23 cells across %d cell",
                             "types — the HypoMap RomanovDev10x subset (P23 only)."),
                       format(total, big.mark = ","), length(cts)),
               "Percentile ranks genes by mean expression among expressing cells.")
    ),
    mainPanel(
      width = 9,
      h2(textOutput("gtitle")),
      uiOutput("globalbox"),
      tags$hr(),
      h4("By cell type"),
      tabsetPanel(
        tabPanel("Table", br(), DTOutput("cttable"),
                 tags$p(style = "color:#777;font-size:12.5px;margin-top:10px;line-height:1.6",
                   HTML("How to read this table: each row is one cell type.
                         <b>#Cells</b> is the number of P23 cells belonging to that type, and
                         <b>#Expressing</b> (with <b>%Expressing</b>) is how many of those cells — and
                         what fraction of them — detectably express this gene.
                         <b>Mean expr. (expressing cells)</b> is the gene's average log-normalized
                         expression, calculated over only the cells that express it, so it reflects how
                         strongly the gene is expressed when it is on.
                         <b>Percentile (within cell type)</b> ranks the gene against every other gene
                         detected in that same cell type by that mean-when-expressed value, so a value of
                         95 means the gene is expressed more strongly than 95% of the genes active in that
                         cell type; <b>Top 5% (in cell type)</b> simply flags the genes that reach the
                         95th percentile or above.
                         <b>Percentile (global)</b> is the same kind of rank taken across the whole cohort
                         rather than within a single cell type, which is why it is identical in every row
                         and matches the summary above."))),
        tabPanel("% expressing", br(),
                 tags$p(style = "color:#555;font-size:12.5px;line-height:2;margin-bottom:12px",
                        HTML(paste0(
                          "<b>%Expressing</b> = ", frac("<b>#Expressing</b>", "<b>#Cells</b>"), " &times; 100",
                          "<br><b>#Expressing</b> = cells within the cell type that carry &ge;1 detected count of the gene.",
                          "<br><b>#Cells</b> = all cells of the type.",
                          "<br>Each bar is on an absolute <b>0&ndash;100%</b> scale."))),
                 uiOutput("pctbars")),
        tabPanel("Mean (expressing cells)", br(),
                 tags$p(style = "color:#555;font-size:12.5px;line-height:2;margin-bottom:12px",
                        HTML(sprintf(paste0(
                          "<b>Mean (expressing cells)</b> = ",
                          frac("<b>&Sigma; expression</b>", "<b>#Expressing</b>"),
                          "<br><b>&Sigma; expression</b> = the gene's log-normalized expression summed across those expressing cells.",
                          "<br><b>#Expressing</b> = cells within the cell type that carry &ge;1 detected count of the gene.",
                          "<br>The bar shows the gene's <b>absolute expression strength</b> when it's on, on one",
                          " fixed scale shared by all genes (so bars are comparable across genes): a full bar",
                          " marks the top 1%% of expression levels across every gene and cell type (<b>%.2f</b>).",
                          "<br>A near-full bar therefore means strong <i>absolute</i> expression. That is not the",
                          " same as ranking highly against other genes in this cell type; for that ranking, see",
                          " <b>Percentile (within cell type)</b> in the Table tab.",
                          "<br>A striped bar with a <b>&#9656;</b> exceeds the scale (refer to the Mean expr. value)."),
                          mean_scale))),
                 uiOutput("meanbars"))
      )
    )
  )
)

server <- function(input, output, session) {
  gene_lc  <- setNames(genes, tolower(genes))   # lowercased -> canonical
  genes_lc <- tolower(genes)

  # flexible matcher: exact > prefix > substring (all case-insensitive)
  match_genes <- function(q, n = 12) {
    q <- tolower(trimws(q)); if (!nzchar(q)) return(character(0))
    exact <- genes[genes_lc == q]
    pre   <- startsWith(genes_lc, q) & genes_lc != q
    sub   <- grepl(q, genes_lc, fixed = TRUE) & !startsWith(genes_lc, q)
    head(unique(c(exact, sort(genes[pre]), sort(genes[sub]))), n)
  }
  suggest <- function(q) {
    m <- match_genes(q, 6)
    if (length(m) == 0) return("No similar symbols found.")
    paste0("Did you mean: ", paste(m, collapse = ", "), "?")
  }

  current    <- reactiveVal("Pomc")   # canonical gene currently shown
  last_query <- reactiveVal("Pomc")

  # live suggestions as you type (clickable)
  output$suggestions <- renderUI({
    q <- trimws(input$gene_in)
    if (nchar(q) < 2) return(NULL)
    # already showing an exact valid gene? hide the list
    if (tolower(q) %in% genes_lc && !is.na(current()) &&
        tolower(current()) == tolower(q)) return(NULL)
    m <- match_genes(q, 12)
    if (length(m) == 0)
      return(tags$div(style = "color:#999;font-size:13px;margin:4px 0", "no matches"))
    tags$div(
      style = paste("margin:4px 0 8px; max-height:190px; overflow-y:auto;",
                    "border:1px solid #e3e3e3; border-radius:6px"),
      lapply(m, function(gn) tags$a(
        href = "#",
        onclick = sprintf(
          "Shiny.setInputValue('picked','%s',{priority:'event'});return false;", gn),
        style = paste("display:block; padding:5px 10px; text-decoration:none;",
                      "color:#1a4f8a; border-bottom:1px solid #f4f4f4"),
        gn)))
  })

  # resolve typed text on Look up / Enter
  observeEvent(input$go, {
    q <- trimws(input$gene_in); last_query(q)
    current(unname(gene_lc[tolower(q)]))   # canonical symbol, or NA if not found
  })
  # a suggestion was clicked (value is already canonical)
  observeEvent(input$picked, {
    current(input$picked); last_query(input$picked)
    updateTextInput(session, "gene_in", value = input$picked)
  })

  g <- reactive({
    cur <- current()
    validate(need(!is.null(cur) && !is.na(cur) && nzchar(cur), sprintf(
      "Gene '%s' not found — try one of the suggestions. %s",
      last_query(), suggest(last_query()))))
    cur
  })

  output$gtitle <- renderText(g())

  output$globalbox <- renderUI({
    gname <- g()
    row <- stats$global[stats$global$gene == gname, ]
    ncells <- format(total, big.mark = ",")
    if (row$n_expressing == 0)
      return(HTML(sprintf(
        "<ul style='font-size:15px;line-height:1.7'><li><i>%s</i> is not detected
         in any of the %s P23 cells in this dataset.</li></ul>", gname, ncells)))
    HTML(sprintf(
      "<ul style='font-size:15px;line-height:1.7'>
        <li><i>%s</i> is expressed in <b>%.2f%%</b> of cells — it is detected in
            %s of the %s P23 cells in this dataset.</li>
        <li>Averaged over only the cells that express it, its mean (log-normalized)
            expression is <b>%.3f</b>.</li>
        <li>Compared with the %s genes that are detectably expressed in this cohort,
            and ranked by that same mean-when-expressed level, <i>%s</i> sits in the
            <b>%.1fth percentile</b> — its expression, when on, is higher than about
            %.0f%% of those genes.</li>
      </ul>",
      gname, row$pct_expressing, format(row$n_expressing, big.mark = ","),
      ncells, row$mean_expressing,
      format(n_expressed, big.mark = ","), gname, row$percentile, row$percentile))
  })

  ctdf <- reactive({
    gg <- g()
    pct_within <- stats$byct$percentile[gg, ]
    pct_global <- stats$global$percentile[stats$global$gene == gg]
    data.frame(
      CellType       = cts,
      nCells         = as.integer(stats$cell_counts),
      nExpressing    = as.integer(stats$byct$n_expressing[gg, ]),
      pctExpressing  = stats$byct$pct_expressing[gg, ],
      meanExpressing = stats$byct$mean_expressing[gg, ],
      pctWithin      = pct_within,
      pctGlobal      = pct_global,
      top5           = ifelse(!is.na(pct_within) & pct_within >= 95, "Yes", ""),
      row.names = NULL, check.names = FALSE)
  })

  output$cttable <- renderDT({
    datatable(ctdf(), rownames = FALSE,
              colnames = c("Cell type", "#Cells", "#Expressing", "%Expressing",
                           "Mean expr. (expressing cells)",
                           "Percentile (within cell type)", "Percentile (global)",
                           "Top 5% (in cell type)"),
              options = list(pageLength = 18, dom = "t")) |>
      formatRound(c("pctExpressing", "meanExpressing", "pctWithin", "pctGlobal"), 2)
  })

  output$pctbars <- renderUI(
    htmlBars(ctdf()$CellType, ctdf()$pctExpressing, ctdf()$nExpressing,
             maxval = 100, fmt = "%.1f%%", accent = "#3aa37a",
             valueTitle = "%Expressing", countTitle = "#Expressing"))

  output$meanbars <- renderUI(
    htmlBars(ctdf()$CellType, ctdf()$meanExpressing, ctdf()$nExpressing,
             maxval = mean_scale, fmt = "%.2f", accent = "#3a6ea3",
             valueTitle = "Mean expr.", countTitle = "#Expressing"))
}

shinyApp(ui, server)
