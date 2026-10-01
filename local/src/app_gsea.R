#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(ggplot2)
  library(fgsea)
})

parse_args <- function(args) {
  values <- list(
    host = "127.0.0.1",
    port = 3838L,
    rnk_dir = ".",
    gmt_dir = ".",
    results = "multi_GSEA.gz",
    gseaParam = 0
  )
  names_map <- c(
    "--host" = "host", "--port" = "port", "--rnk-dir" = "rnk_dir",
    "--gmt-dir" = "gmt_dir", "--results" = "results",
    "--gseaParam" = "gseaParam"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% names(names_map) || i == length(args)) {
      stop("Usage: app_gsea.R [--host HOST] [--port PORT] [--rnk-dir DIR] ",
           "[--gmt-dir DIR] [--results FILE] [--gseaParam WEIGHT]")
    }
    values[[names_map[[key]]]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values$port <- suppressWarnings(as.integer(values$port))
  if (is.na(values$port) || values$port < 1L || values$port > 65535L) {
    stop("--port must be an integer between 1 and 65535")
  }
  values$gseaParam <- suppressWarnings(as.numeric(values$gseaParam))
  if (!is.finite(values$gseaParam) || values$gseaParam < 0) {
    stop("--gseaParam must be a finite non-negative number")
  }
  values
}

read_ranking <- function(path) {
  # Same ranking preparation as shinySea/server.R; do not filter or collapse genes.
  ranking <- read.delim(path, header = FALSE, col.names = c("GeneID", "value"))
  ranking <- ranking[order(ranking$value, decreasing = TRUE), ]
  ranking$ranking_position <- seq_len(nrow(ranking))
  ranking
}

read_gsea_results <- function(path) {
  if (!file.exists(path)) stop("GSEA results not found: ", path)
  has_header <- grepl("header_added", basename(path), fixed = TRUE)
  results <- read.delim(path, header = has_header, stringsAsFactors = FALSE,
                        check.names = FALSE)
  if (!has_header) {
    expected <- c("contrast", "msigdb_type", "pathway", "pval", "padj",
                  "ES", "NES", "nMoreExtreme", "size", "leadingEdge")
    if (ncol(results) != length(expected)) {
      stop("Expected 10 columns in headerless GSEA results, found ", ncol(results))
    }
    names(results) <- expected
  }
  required <- c("contrast", "msigdb_type", "pathway", "pval", "padj",
                "NES", "size")
  missing <- setdiff(required, names(results))
  if (length(missing)) stop("Missing GSEA columns: ", paste(missing, collapse = ", "))
  for (column in c("pval", "padj", "ES", "NES", "size")) {
    if (column %in% names(results)) {
      results[[column]] <- suppressWarnings(as.numeric(results[[column]]))
    }
  }
  results
}

normalise_id <- function(x) sub("\\.rnk$", "", basename(as.character(x)))

options <- parse_args(commandArgs(trailingOnly = TRUE))
rnk_files <- list.files(options$rnk_dir, pattern = "\\.rnk$", full.names = TRUE)
gmt_files <- list.files(options$gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
if (!length(rnk_files)) stop("No .rnk files found in: ", options$rnk_dir)
if (!length(gmt_files)) stop("No .gmt files found in: ", options$gmt_dir)

rankings <- setNames(lapply(rnk_files, read_ranking), normalise_id(rnk_files))
pathways <- list()
for (gmt_file in gmt_files) {
  collection <- sub("\\.gmt$", "", basename(gmt_file))
  pathways[[collection]] <- fgsea::gmtPathways(gmt_file)
}
gsea_table <- read_gsea_results(options$results)
gsea_table$contrast_id <- normalise_id(gsea_table$contrast)
available_rankings <- intersect(names(rankings), unique(gsea_table$contrast_id))
if (!length(available_rankings)) {
  stop("No ranking name matches the contrasts in ", options$results,
       ". Rankings: ", paste(names(rankings), collapse = ", "),
       "; contrasts: ", paste(unique(gsea_table$contrast_id), collapse = ", "))
}
rankings <- rankings[available_rankings]

ui <- fluidPage(
  titlePanel("SnakeGSEA local explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("ranking", "Ranking", choices = names(rankings)),
      selectInput("collection", "MSigDB collection",
                  choices = c("ALL", sort(unique(gsea_table$msigdb_type)))),
      sliderInput("padj", "Maximum FDR (padj)", min = 0.001, max = 1,
                  value = 0.05, step = 0.01),
      radioButtons("direction", "Enrichment direction",
                   choices = c("All" = "all", "NES > 0" = "positive",
                               "NES < 0" = "negative")),
      helpText("Select a table row to display its enrichment curve and leading edge.")
    ),
    mainPanel(
      DTOutput("results_table"),
      helpText("The table above contains pipeline results. The selected pathway below is recalculated as in shinySea."),
      h4(textOutput("selection_title")),
      tableOutput("fgsea_results"),
      plotOutput("enrichment_plot", height = "380px"),
      h4("Gene-set genes and leading-edge membership"),
      downloadButton("downloadData", "Download leading_edge.tsv"),
      DTOutput("leading_edge_table")
    )
  )
)

server <- function(input, output, session) {
  filtered_results <- reactive({
    data <- gsea_table[gsea_table$contrast_id == input$ranking, , drop = FALSE]
    if (input$collection != "ALL") {
      data <- data[data$msigdb_type == input$collection, , drop = FALSE]
    }
    data <- data[is.finite(data$padj) & data$padj <= input$padj, , drop = FALSE]
    if (input$direction == "positive") data <- data[data$NES > 0, , drop = FALSE]
    if (input$direction == "negative") data <- data[data$NES < 0, , drop = FALSE]
    data
  })

  output$results_table <- renderDT({
    data <- filtered_results()
    shown <- intersect(c("msigdb_type", "pathway", "NES", "pval", "padj", "size"),
                       names(data))
    datatable(data[, shown, drop = FALSE], selection = "single", rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE))
  })

  selected_row <- reactive({
    data <- filtered_results()
    selected <- input$results_table_rows_selected
    if (!length(selected) && nrow(data)) selected <- 1L
    req(length(selected), nrow(data) >= selected)
    data[selected, , drop = FALSE]
  })

  custom_gsea <- reactive({
    row <- selected_row()
    pathway <- row$pathway[[1]]
    collection <- pathways[[row$msigdb_type[[1]]]]
    genes <- collection[[pathway]]
    validate(need(length(genes) > 0, "Selected pathway is missing from its local GMT collection."))
    rnk <- rankings[[input$ranking]]
    rnk_vector <- unlist(rnk$value)
    names(rnk_vector) <- rnk$GeneID
    gene_set <- as.data.frame(genes)
    colnames(gene_set) <- "GeneID"

    # Preserve the original shinySea single-pathway calculation and defaults.
    fgseaRes <- fgsea::fgseaSimple(
      pathways = gene_set, stats = rnk_vector, nperm = 1000,
      minSize = 5, maxSize = 5000, gseaParam = options$gseaParam
    )
    validate(need(nrow(fgseaRes) > 0,
                  "No result: the pathway may be outside shinySea's 5–5000 gene limits."))
    fgseaRes_tab <- subset(fgseaRes, select = -c(leadingEdge, pathway))
    fgseaPlot <- fgsea::plotEnrichment(
      pathway = gene_set[, 1], stats = rnk_vector, gseaParam = options$gseaParam
    )

    leadingEdge_tab <- as.data.frame(fgseaRes$leadingEdge[[1]])
    colnames(leadingEdge_tab) <- "Leading_Edge"
    leadingEdge_tab$hit <- rep("Yes", nrow(leadingEdge_tab))
    leadingEdge_tab <- merge(gene_set, leadingEdge_tab, all.x = TRUE,
                            by.y = "Leading_Edge", by.x = "GeneID")
    leadingEdge_tab[is.na(leadingEdge_tab)] <- "No"
    leadingEdge_tab <- merge(leadingEdge_tab, rnk)
    colnames(leadingEdge_tab)[3] <- "ranking_score"
    leadingEdge_tab <- leadingEdge_tab[order(leadingEdge_tab$ranking_position), ]
    list(title = gsub("_", " ", paste(row$msigdb_type[[1]], pathway, sep = "_")),
         plot = fgseaPlot, results = fgseaRes_tab, leading_edge = leadingEdge_tab)
  })

  output$selection_title <- renderText(custom_gsea()$title)
  output$fgsea_results <- renderTable(custom_gsea()$results, digits = 4)
  output$enrichment_plot <- renderPlot(plot(custom_gsea()$plot))
  output$leading_edge_table <- renderDT({
    datatable(custom_gsea()$leading_edge, rownames = FALSE,
              options = list(pageLength = 10, order = list()))
  })
  output$downloadData <- downloadHandler(
    filename = function() "leading_edge.tsv",
    content = function(file) {
      write.table(custom_gsea()$leading_edge, file, row.names = FALSE,
                  quote = FALSE, sep = "\t")
    }
  )
}

message("SnakeGSEA local explorer listening on http://", options$host, ":", options$port)
runApp(shinyApp(ui, server), host = options$host, port = options$port,
       launch.browser = FALSE)
