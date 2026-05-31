source("functions.R", local = TRUE)
library(jsonlite)
library(JBrowseR)


server <- function(input, output, session) {
  
  vcf_df <- reactiveVal(NULL) # dataframe variante/es para tabla
  jbrowse_region <- reactiveVal(NULL) # región a mostrar (chr:start..end)
  
  # Apartado busqueda:
  
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
      parsed <- parse_genomic(input$genomic_input)
      
      if (is.null(parsed)) {
        showNotification("Error. Algo fue mal.", 
                         type = "error", duration = 5)
      } else {
        showNotification(paste0("Formato genómico introducido correctamente.",
                                parsed[1]," ",parsed[2]," ",parsed[3],"/",parsed[4]), 
                         type = "warning", duration = 5)
      }
      
    } else if (input$variant_input_type == "hgvsc") {
      req(input$hgvsc_input, nzchar(trimws(input$hgvsc_input)))
      parsed <-parse_hgvsc(input$hgvsc_input)

    }
    
    vcf_df(data.frame(parsed))
  })
  
  
  # Apartado carga:
  
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
      gz_path <- NULL
      
      if (has_plain && input$type_file_mode == 'plain') {
        showNotification("Procesando VCF plano ...", 
                         type = "message", duration = 6)
        
        gz_path <- process_vcf_plain(input$vcf_plain$datapath, temp_dir)
        
      } else if (has_gz && has_tbi && input$type_file_mode == 'bgzip') {
        
        showNotification("Procesando VCF comprimido ...", type = "message", duration = 6)
        
        file.copy(input$vcf_file$datapath, file.path(temp_dir, "variants.vcf.gz"), overwrite = TRUE)
        file.copy(input$tbi_file$datapath, file.path(temp_dir, "variants.vcf.gz.tbi"), overwrite = TRUE)
        gz_path <- file.path(temp_dir, "variants.vcf.gz")
      }
      
      vcf_df(vcf_to_df(gz_path))
      
    }, error = function(e) {
      showNotification(e$message, type = "error")
    })
  })
  
  # Gestión de tabla de variantes:
  
  variants_table_data <- reactive({
    if (!is.null(vcf_df())) return(vcf_df())
    NULL
  })
  
  output$variants_table <- renderDT({
    req(variants_table_data())
    datatable(variants_table_data(),
              selection = "single",
              options = list(
                scrollX = TRUE,           
                scrollY = FALSE))
  })
  
  proxy <- dataTableProxy("variants_table")
  
  observeEvent(variants_table_data(), {
    req(variants_table_data())
    selectRows(proxy, 1)
  })
  
  observeEvent(input$variants_table_rows_selected, {
    row <- variants_table_data()[input$variants_table_rows_selected, ]
    jbrowse_region(make_region(row$CHR, row$POS, input$window_size))
  })
  
  #Gestion de navegador genómico JBrowse:
  
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
  
  