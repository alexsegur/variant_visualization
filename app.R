library(shiny)
library(shinydashboard)
library(DT)
library(VariantAnnotation)
library(GenomicRanges)
library(JBrowseR)
library(httr)
library(jsonlite)
source("functions.R", local = TRUE)


# ==================== UI ====================
ui <- dashboardPage(
  dashboardHeader(title = "Variant Browser"),
  
  dashboardSidebar(
                   
    h4("Búsqueda de variante", style = "padding-left: 10px; font-weight: bold;"),
    
    tags$div(style = "padding: 0px 12px;",
      radioButtons("variant_input_type", "Formato de entrada:",
        choices = list(
          "HGVSc" = "hgvsc",
          "Formato Genómico" = "genomic"
        ),
        selected = "hgvsc"
      ),
      
      conditionalPanel(
        condition = "input.variant_input_type == 'hgvsc'",
        textInput("hgvsc_input", "HGVSc:", 
                  placeholder = "Ej: NM_000059.4:c.-152T>C")
      ),
      
      conditionalPanel(
        condition = "input.variant_input_type == 'genomic'", 
        textInput("genomic_input", "Formato genómico:", 
                  placeholder = "Ej: 13 32889692 T/C")
      ),
      
      numericInput("window_size", "Entorno (pb):", 
                   value = 20, min = 20, max = 1000),
      
      actionButton("search_variant", "Buscar variante",
                   class = "btn-primary")
    ),
    
    hr(),
    h4("Carga de archivos", style = "padding-left: 10px; font-weight: bold;"),
    
    radioButtons("type_file_mode", "Selecciona tipo de VCF:",
                 choices = c(
                   "Subir VCF plano (.vcf)" = "plain",
                   "Subir VCF comprimido + índice (.vcf.gz + .vcf.gz.tbi)" = "bgzip"
                 ),
                 selected = "plain"),
    
    # Panel para opción 1: Plain VCF
    conditionalPanel(
      condition = "input.type_file_mode == 'plain'",
      fileInput("vcf_plain", "Cargar .vcf",
                accept = c(".vcf"))
    ),
    
    # Panel para opción 2: VCF comprimido + indice
    conditionalPanel(
      condition = "input.type_file_mode == 'bgzip'",
      fileInput("vcf_file", "Cargar .vcf.gz", 
                accept = ".vcf.gz"),
      fileInput("tbi_file", "Cargar .vcf.gz.tbi", 
                accept = ".tbi")
    ),
    
    # Botón para cargar
    actionButton("load_vcf_btn", "Visualizar", 
                 class = "btn-primary", icon = icon("upload"))
    
  ),
  
  dashboardBody(
    fluidRow(
      box(width = 12,
          title = "Lista de variantes clasificadas",
          DTOutput("variants_table"))
    ),
    fluidRow(
      box(width = 12, height = "100%",
          title = "Visualización JBrowseR - Contexto genómico",
          JBrowseROutput("jbrowse", height = "100%"))
    )
  )
)

# ==================== SERVER ====================
server <- function(input, output, session) {
  
  vcf_df <- reactiveVal(NULL)   # dataframe variante/es para tabla
  clinvar_df <- reactiveVal(NULL)   # variantes ClinVar
  jbrowse_region <- reactiveVal(NULL)   # región a mostrar
  jbrowse_server <- reactiveVal(NULL)
  
  observeEvent(input$search_variant, {
    req(input$search_variant)
    if (input$window_size < 20 || input$window_size > 1000) {
      showNotification(
        "El entorno debe estar entre 20 y 1000 pb",
        type = "error"
      )
      return(NULL)
    }
    
    parsed <- NULL
    
    if (input$variant_input_type == "genomic") {
      req(input$genomic_input, nzchar(trimws(input$genomic_input)))
      parsed <- parse_genomic_input(input$genomic_input)
      
      if (is.null(parsed)) {
        showNotification("Formato genómico inválido. Use: 'chr pos ref/alt'", 
                         type = "error", duration = 5)
      } else {
        showNotification(paste0("Formato genómico introducido correctamente.",
                                parsed[1]," ",parsed[2]," ",parsed[3],"/",parsed[4]), 
                         type = "warning", duration = 5)
      }
      
    } else if (input$variant_input_type == "hgvsc") {
      req(input$hgvsc_input, nzchar(trimws(input$hgvsc_input)))
      parsed <-parse_hgvsc(input$hgvsc_input)
      
      if (is.null(parsed)) {
        showNotification("Formato HGVSc inválido o no encontrado.", 
                         type = "error", duration = 5)
      } else {
        showNotification("Formato HGVSc introducido correctamente.", 
                         type = "warning", duration = 5)
      }
    }
    vcf_df(data.frame(parsed))
  })
  
  observeEvent(input$load_vcf_btn, {
    req(input$load_vcf_btn)
    
    has_plain <- !is.null(input$vcf_plain) &&
      !is.null(input$vcf_plain$datapath) &&
      nzchar(input$vcf_plain$datapath)
    
    has_gz <- !is.null(input$vcf_file) &&
      !is.null(input$vcf_file$datapath) &&
      nzchar(input$vcf_file$datapath)
    
    has_tbi <- !is.null(input$tbi_file) &&
      !is.null(input$tbi_file$datapath) &&
      nzchar(input$tbi_file$datapath)
    
    if (!has_plain && !(has_gz && has_tbi)) {
      showNotification("No se ha cargado ningun archivo.", type = "error")
      return(NULL)
    }
    
    tryCatch({
      temp_dir <- tempfile("jbrowse_data_")
      dir.create(temp_dir)
      showNotification(paste("Directorio temporal creado:", basename(temp_dir)), 
                       type = "message", duration = 5)
      gz_path <- NULL
      
      if (has_plain && input$type_file_mode == 'plain') {
        showNotification("Procesando VCF plano (ordenando + bgzip + tabix)...", 
                         type = "message", duration = 6)
        
        gz_path <- process_vcf_plain(input$vcf_plain$datapath, temp_dir)
        
      } else if (has_gz && has_tbi && input$type_file_mode == 'bgzip') {
        
        showNotification("Cargando archivos .vcf.gz + .tbi...", type = "message", duration = 6)
        
        file.copy(input$vcf_file$datapath, file.path(temp_dir, "variants.vcf.gz"), overwrite = TRUE)
        file.copy(input$tbi_file$datapath, file.path(temp_dir, "variants.vcf.gz.tbi"), overwrite = TRUE)
        gz_path <- file.path(temp_dir, "variants.vcf.gz")
      }
      
      # tabla
      vcf_df(vcf_to_df(gz_path))
      
      # JBrowse server
      jbrowse_server(serve_data(gz_path))
    }, error = function(e) {
      showNotification(e$message, type = "error")
    })
  })
  
  variants_table_data <- reactive({
    if (!is.null(vcf_df())) return(vcf_df())
    NULL
  })
  
  proxy <- dataTableProxy("variants_table")
  
  observeEvent(variants_table_data(), {
    req(variants_table_data())
    selectRows(proxy, 1)
  })
  
  output$variants_table <- renderDT({
    req(variants_table_data())
    datatable(variants_table_data(), selection = "single")
  })
  
  variant_select <- reactive({
    req(input$variants_table_rows_selected)
    
    row <- variants_table_data()[input$variants_table_rows_selected, ]
    
    paste0(as.character(GenomicRanges::seqnames(row)), ":", row$start)
  })
  
  # ==================== JBrowseR ====================
  
  observeEvent(input$variants_table_rows_selected, {
    row <- variants_table_data()[input$variants_table_rows_selected, ]
    jbrowse_region(make_region(row$CHR, row$POS, input$window_size))
  })
  
  observeEvent(jbrowse_region(), {
    output$jbrowse <- renderJBrowseR({
      config <- json_config("./config.json")
    
      JBrowseR(
        "JsonView",
        config = config,
        location = jbrowse_region()
      )
    })
  })

}

# Lanzar la app
shinyApp(ui = ui, server = server)