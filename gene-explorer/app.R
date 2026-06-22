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
                placeholder = "e.g. Pomc"),
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
        tabPanel("Table", br(), DTOutput("cttable")),
        tabPanel("% expressing", br(), plotOutput("pctplot", height = 430)),
        tabPanel("Mean (expressing cells)", br(), plotOutput("meanplot", height = 430))
      )
    )
  )
)

server <- function(input, output, session) {
  # case-insensitive lookup map: lowercased symbol -> canonical symbol
  gene_lc <- setNames(genes, tolower(genes))

  # suggest canonical symbols that start with what was typed
  suggest <- function(q) {
    if (!nzchar(q)) return("")
    hits <- genes[startsWith(tolower(genes), tolower(q))]
    if (length(hits) == 0) return("No similar symbols found.")
    paste0("Did you mean: ", paste(head(hits, 6), collapse = ", "), "?")
  }

  # resolve only when "Look up" is clicked (or Enter); fires once on load too
  looked <- eventReactive(input$go, ignoreNULL = FALSE, {
    q <- trimws(input$gene_in)
    list(query = q, gene = unname(gene_lc[tolower(q)]))
  })

  g <- reactive({
    L <- looked()
    validate(need(nzchar(L$query), "Type a gene symbol and click Look up."))
    validate(need(!is.na(L$gene), sprintf(
      "Gene '%s' not found — symbols are case-sensitive (e.g. Pomc, Agrp). %s",
      L$query, suggest(L$query))))
    L$gene
  })

  output$gtitle <- renderText(g())

  output$globalbox <- renderUI({
    row <- stats$global[stats$global$gene == g(), ]
    flag <- if (isTRUE(row$top5pct))
      "<b style='color:#0a7a4f'>Yes — top 5%</b>" else "No"
    HTML(sprintf(
      "<ul style='font-size:15px;line-height:1.7'>
        <li><b>%% of cells expressing (global):</b> %.2f%%
            &nbsp;<span style='color:#888'>(%s of %s cells)</span></li>
        <li><b>Mean expression (among expressing cells):</b> %.3f</li>
        <li><b>Expression percentile:</b> %.1f
            <span style='color:#888'>(ranked by mean-when-expressed)</span></li>
        <li><b>Top 5%% of expressed genes?</b> %s</li>
      </ul>",
      row$pct_expressing, format(row$n_expressing, big.mark = ","),
      format(total, big.mark = ","), row$mean_expressing, row$percentile, flag))
  })

  ctdf <- reactive({
    gg <- g()
    data.frame(
      CellType       = cts,
      nCells         = as.integer(stats$cell_counts),
      nExpressing    = as.integer(stats$byct$n_expressing[gg, ]),
      pctExpressing  = stats$byct$pct_expressing[gg, ],
      meanExpressing = stats$byct$mean_expressing[gg, ],
      percentile     = stats$byct$percentile[gg, ],
      top5           = ifelse(stats$byct$percentile[gg, ] >= 95, "Yes", ""),
      row.names = NULL, check.names = FALSE)
  })

  output$cttable <- renderDT({
    datatable(ctdf(), rownames = FALSE,
              colnames = c("Cell type", "#Cells", "#Expressing", "%Expressing",
                           "Mean (expr. cells)", "Percentile", "Top 5%"),
              options = list(pageLength = 18, dom = "t")) |>
      formatRound(c("pctExpressing", "meanExpressing", "percentile"), 2)
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
