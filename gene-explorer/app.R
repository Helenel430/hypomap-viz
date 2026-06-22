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
library(ggplot2)

stats <- readRDS("gene_stats.rds")
genes <- stats$genes
total <- stats$total_cells
cts   <- stats$celltypes

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
               "RNA assay, log-normalized; 51,245 cells across 18 cell types.",
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
                 tags$p(style = "color:#777;font-size:12.5px;margin-top:10px;line-height:1.55",
                   HTML("<b>%Expressing</b> / <b>#Expressing</b> are out of this cell type's <b>#Cells</b>.
                         <b>Mean expr.</b> averages log-normalized expression over only the expressing cells.
                         <b>Percentile (within cell type)</b> ranks this gene against all genes detected in
                         that same cell type, by mean-when-expressed; <b>Top 5% (in cell type)</b> = ≥95th
                         percentile within the cell type. <b>Percentile (global)</b> is the gene's
                         dataset-wide rank (identical in every row; matches the summary above)."))),
        tabPanel("% expressing", br(), plotOutput("pctplot", height = 430)),
        tabPanel("Mean (expressing cells)", br(), plotOutput("meanplot", height = 430))
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
    row <- stats$global[stats$global$gene == g(), ]
    HTML(sprintf(
      "<ul style='font-size:15px;line-height:1.7'>
        <li><b>%% of cells expressing (global):</b> %.2f%%
            &nbsp;<span style='color:#888'>(%s of %s cells)</span></li>
        <li><b>Mean expression (among expressing cells):</b> %.3f</li>
        <li><b>Expression percentile:</b> %.1f
            <span style='color:#888'>(ranked by mean-when-expressed)</span></li>
      </ul>",
      row$pct_expressing, format(row$n_expressing, big.mark = ","),
      format(total, big.mark = ","), row$mean_expressing, row$percentile))
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

  output$pctplot <- renderPlot({
    d <- ctdf(); d <- d[order(d$pctExpressing), ]
    d$CellType <- factor(d$CellType, levels = d$CellType)
    ggplot(d, aes(pctExpressing, CellType)) +
      geom_col(fill = "#3aa37a") +
      labs(x = "% of cells expressing", y = NULL,
           title = paste0(g(), " — detection by cell type")) +
      theme_minimal(base_size = 14)
  })

  output$meanplot <- renderPlot({
    d <- ctdf(); d <- d[order(d$meanExpressing), ]
    d$CellType <- factor(d$CellType, levels = d$CellType)
    ggplot(d, aes(meanExpressing, CellType)) +
      geom_col(fill = "#3a6ea3") +
      labs(x = "Mean expression (expressing cells)", y = NULL,
           title = paste0(g(), " — intensity by cell type")) +
      theme_minimal(base_size = 14)
  })
}

shinyApp(ui, server)
