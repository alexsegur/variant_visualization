

ui <- dashboardPage(
  dashboardHeader(title = "Variant Browser"),
  
  dashboardSidebar( # BARRA IZQUIERDA
    
    # Apartado Busqueda
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
             # No sirve para limitar, se comprueba en server.R
             
             actionButton("search_variant", "Buscar variante",
                          class = "btn-primary")
    ),
    
    hr(), # separador
    
    # Apartado carga
    h4("Carga de archivos", style = "padding-left: 10px; font-weight: bold;"),
    
    tags$div(style = "padding: 0px 12px;",
             radioButtons("type_file_mode", "Selecciona tipo de VCF:",
                     choices = c(
                       "Subir VCF plano (.vcf)" = "plain",
                       "Subir VCF comprimido + índice (.vcf.gz + .vcf.gz.tbi)" = "bgzip"
                     ),
                     selected = "plain"),
        
    
            conditionalPanel(
              condition = "input.type_file_mode == 'plain'",
              fileInput("vcf_plain", "Cargar .vcf",
                        accept = c(".vcf"))
            ),
            
    
            conditionalPanel(
              condition = "input.type_file_mode == 'bgzip'",
              fileInput("vcf_file", "Cargar .vcf.gz", 
                        accept = ".vcf.gz"),
              fileInput("tbi_file", "Cargar .vcf.gz.tbi", 
                        accept = ".tbi")
            ),
        
            actionButton("load_vcf_btn", "Visualizar", 
                     class = "btn-primary", icon = icon("upload"))
    )
    
  ),

  dashboardBody( # ESPACIO CENTRAL
    
    fluidRow( # Tabla arriba
      box(width = 12,
          title = "Lista de variantes clasificadas",
          DTOutput("variants_table"))
    ),
    fluidRow( # JBrowse debajo
      box(width = 12, height = "100%",
          title = "Visualización JBrowseR - Contexto genómico",
          JBrowseROutput("jbrowse", height = "100%"))
    )
  )
)