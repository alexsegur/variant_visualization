library(VariantAnnotation)
library(GenomicRanges)
library(httr)
library(jsonlite)
library(dplyr)

## TODO: Falta separar validaciones de functions.R

`%||%` <- function(x, y) if (!is.null(x)) x else y


# De vcf + tabix a dataframe (para opción de carga de archivos)
vcf_to_df <- function(path) {
  
  vcf <- readVcf(path, genome = "hg38")
  gr <- rowRanges(vcf)
  n <- length(gr)
  df <- data.frame(
    ID = names(gr) %||% ".",
    CHR = as.character(seqnames(gr)),
    POS = as.character(start(gr)),
    REF = as.character(ref(vcf)),
    ALT = sapply(alt(vcf), function(x) paste(as.character(x), collapse = "")),
    QUAL = qual(vcf) %||% ".",
    FILTER = fixed(vcf)$FILTER %||% ".",
    INFO = sub("^INFO=", "", info(vcf)$INFO) %||% ".",
    stringsAsFactors = FALSE
  )
  return(df)
}




# De archivo vcf a vcf comprimido + tabix
process_vcf_plain <- function(path_in, temp_dir) {
  
  vcf <- readVcf(path_in, genome = "hg38")
  gr <- rowRanges(vcf)
  gr_sorted <- sort(gr, ignore.strand = TRUE)
  vcf_sorted <- vcf
  rowRanges(vcf_sorted) <- gr_sorted
  
  plain_file <- file.path(temp_dir, "variants.vcf")
  writeVcf(vcf_sorted, plain_file)

  gz_file <- file.path(temp_dir, "variants.vcf.gz")
  bgzip(plain_file, dest = gz_file, overwrite = TRUE)
  indexTabix(gz_file, format = "vcf")
  
  return(gz_file)
}




# Parsear HGVSc (NM_000059.4:c.7A>G o NC_000007.14:g.55249071T>A) en GRCh38
# Con anotación de Ensembl VEP API 
parse_hgvsc <- function(hgvsc) {
  
    server <- "https://rest.ensembl.org"
    ext <- paste0("/vep/human/hgvs/", hgvsc, "?",
                  "canonical=1&", #Selecciona los transcriptos canónicos
                  "mane=1&", #Selecciona el transcripto canónico, mostrando el ID
                  "AlphaMissense=1&", #Valora peligrosidad missense
                  "ClinPred=1&", #similar a alphamissense
                  "Enformer=1&", # Peligrosidad Variantes reguladoras
                  "CADD=1&", #CADDphred peligrosidad del 0 al 99
                  "REVEL=1&", #Peligrosidad del 0 al 1
                  "SIFT=1&") #Peligrosidad para funcion de proteina
    
    r <- tryCatch({
      httr::GET(
        paste0(server, ext),
        httr::content_type("application/json"),
        httr::timeout(20)
      )
    }, error = function(e) {
      showNotification( paste("Servicio de busqueda HGVS no disponible:",e), 
                       type = "error", duration = 5)
      return(NULL)
    })
    if (is.null(r)) return(NULL)
    
    codigo <- codigo <- status_code(r)
    if (codigo == 400) {
      showNotification("HGVS o parámetros inválidos.", 
                       type = "error", duration = 5)
      return (NULL)
    } else if (codigo != 200){
      showNotification(paste("Error inesperado con código:", codigo),
                       type = "error", duration = 5)
      return(NULL)
    }
    
    res <- httr::content(r, as = "parsed", simplifyVector = TRUE)

    if (length(res) == 0) return(NULL)
    
    chr <- res$seq_region_name %||% "."
    pos <- res$start %||% "."
    alleles <- res$allele_string %||% "."
    items <- strsplit(alleles, "/")[[1]] %||% "."
    ref <- toupper(items[1]) %||% "."
    alt <- toupper(items[2]) %||% "."
    
    most_severe <- res$most_severe_consequence %||% "."
    input_transcript <- sub(":.*$", "", hgvsc)
    
    
    consequences <- res$transcript_consequences
    selected <- NULL
    
    if (!is.null(consequences) && length(consequences) > 0) {
      
      df <- bind_rows(consequences)
      
      canonical_df <- df[df$canonical == 1, ]
      
      if (nrow(canonical_df) == 0){
        selected <- df[1,]
      } else{
        mane_rows <- canonical_df[!is.na(canonical_df$mane_select) & canonical_df$mane_select != "", ]
        if (nrow(mane_rows) > 0) {
          mane_match <- mane_rows[mane_rows$mane_select == input_transcript, ]
          if (nrow(mane_match) == 0){
            selected <- canonical_df[1,]
          } else{
            selected <- mane_match[1,]
            browser()
          }
        } else {
          selected <- canonical_df[1,]
        }
      } 
    }
    
    # Función para conseguir extraer el campo deseado
    get_value <- function(field, nested_field = NULL) {
      
      if (is.null(selected)) return("-")
      
      # Valor directo
      if (field %in% names(selected) && is.null(nested_field) &&
          !is.null(selected[[field]]) && !is.na(selected[[field]])) {
        return(selected[[field]])
      }
      
      # Valor anidado
      else if (!is.null(nested_field) &&
               field %in% names(selected)) {
        obj <- selected[[field]]
        
        if (is.data.frame(obj) &&
            nested_field %in% names(obj)) {
          value <- obj[[nested_field]]
          if (!is.null(value) && !is.na(value)) {
            return(value)
          }
        }
      }
      return("-")
    }
    
    gen <- get_value("gene_symbol")
    alpha_missense <- get_value("alphamissense","am_pathogenicity")
    alpha_pathogenicity<- get_value("alphamissense","am_class") 
    cadd_phred <- get_value("cadd_phred")
    clinpred_score <- get_value("clinpred")
    revel_score <- get_value("revel")
    sift_score <- get_value("sift_score")
    polyphen_score <- get_value("polyphen_score")
    polyphen_predict <- get_value("polyphen_prediction")
    enformer_score <- get_value("enformer_sar")
    
    return(list(CHR = chr, POS = pos, REF = ref, ALT = alt, Gen = gen,
                "Peor Consecuencia" = most_severe,
                Alphamissense = alpha_missense,
                Alpha.predict = alpha_pathogenicity,
                CADD = cadd_phred,
                Clinpred = ifelse(clinpred_score == "-", "-", round(as.numeric(clinpred_score),3)),
                REVEL = revel_score,
                SIFT = sift_score,
                Polyph. = polyphen_score,
                "Polyph predict" = polyphen_predict,
                Enformer = enformer_score
                ))
    
}




# Función para parsear formato genómico "13 32889692 T/C" usando VEP
# Únicamente funciona con formato: "Chr pos ref/alt", vacío = "-"
parse_genomic <- function(genomic_str) {
  
  # 1. Validar y parsear el input
  parts <- strsplit(trimws(genomic_str), "\\s+")[[1]]
  
  if (length(parts) != 3) {
    showNotification("Formato inválido. Use: 'Cromosoma Posición REF/ALT'", 
                     type = "error", duration = 5)
    return(NULL)
  }
  
  chr <- toupper(parts[1])
  pos <- suppressWarnings(as.numeric(parts[2]))
  alleles <- parts[3]
  items <- strsplit(alleles, "/")[[1]]
  ref <- toupper(items[1])
  alt <- toupper(items[2])
  
  # Validaciones básicas
  if (is.na(pos)) {
    showNotification("Posición debe ser numérica", type = "error", duration = 5)
    return(NULL)
  }
  
  # Validar cromosoma
  v_chr <- grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)$", chr, ignore.case = TRUE)
  if (!v_chr) {
    showNotification("Cromosoma inválido. Use 1-22, X, Y o MT", 
                     type = "error", duration = 5)
    return(NULL)
  }
  
  # Validar formato de alelos
  if (!validate_allele(ref) || !validate_allele(alt)) {
    showNotification("REF o ALT contienen caracteres inválidos", 
                     type = "error", duration = 5)
    return(NULL)
  } else if(ref == alt) {
    showNotification("REF y ALT deben ser diferentes", 
                     type = "error", duration = 5)
    return(NULL)
  }
  
  # 2. Construir la URL para VEP usando region
  server <- "https://rest.ensembl.org"
  
  # Determinar el formato de coordenadas según el tipo de variante
  if (ref == "-") {
    # Inserción: formato START-END/ALT con END = START - 1
    start_pos <- pos
    end_pos <- pos - 1
    variant_allele <- alt
  } else if (alt == "-") {
    # Deleción: formato START-END/-
    start_pos <- pos
    end_pos <- pos + nchar(ref) - 1
    variant_allele <- "-"
  } else {
    # SNP o substitución: formato START-END/ALT
    start_pos <- pos
    end_pos <- pos + nchar(ref) - 1
    variant_allele <- alt
  }
  browser()
  # Construir el endpoint
  region <- paste0(chr, ":", start_pos, "-", end_pos)
  ext <- paste0("/vep/human/region/", region, "/", variant_allele, "?",
                "canonical=1&",      # Selecciona transcriptos canónicos
                "mane=1&",           # Selecciona transcripto MANE
                "AlphaMissense=1&",  # Valora peligrosidad missense
                "ClinPred=1&",       # Similar a AlphaMissense
                "Enformer=1&",       # Peligrosidad variantes reguladoras
                "CADD=1&",           # CADD phred (0-99)
                "REVEL=1&",          # Peligrosidad (0-1)
                "SIFT=1&")           # Peligrosidad función proteica
  
  # 3. Realizar la llamada a la API
  r <- tryCatch({
    httr::GET(
      paste0(server, ext),
      httr::content_type("application/json"),
      httr::timeout(20)
    )
  }, error = function(e) {
    showNotification(paste("Servicio de búsqueda VEP no disponible:", e), 
                     type = "error", duration = 5)
    return(NULL)
  })
  
  if (is.null(r)) return(NULL)
  
  # 4. Manejar códigos de estado HTTP
  codigo <- status_code(r)
  if (codigo == 400) {
    showNotification("Coordenadas o parámetros inválidos.", 
                     type = "error", duration = 5)
    return(NULL)
  } else if (codigo == 404) {
    showNotification("Variante no encontrada en la región especificada.", 
                     type = "error", duration = 5)
    return(NULL)
  } else if (codigo != 200) {
    showNotification(paste("Error inesperado con código:", codigo),
                     type = "error", duration = 5)
    return(NULL)
  }
  
  # 5. Procesar la respuesta
  res <- httr::content(r, as = "parsed", simplifyVector = TRUE)
  
  if (length(res) == 0) {
    showNotification("No se encontraron consecuencias para esta variante", 
                     type = "warning", duration = 5)
    return(NULL)
  }
  
  # Extraer información básica
  chr <- res$seq_region_name %||% "."
  pos <- res$start %||% "."
  alleles <- res$allele_string %||% "."
  items <- strsplit(alleles, "/")[[1]] %||% "."
  ref <- toupper(items[1]) %||% "."
  alt <- toupper(items[2]) %||% "."
  
  most_severe <- res$most_severe_consequence %||% "."
  
  # 6. Procesar transcript consequences (misma lógica que parse_hgvsc)
  consequences <- res$transcript_consequences
  selected <- NULL
  
  if (!is.null(consequences) && length(consequences) > 0) {
    
    df <- bind_rows(consequences)
    
    # Verificar que exista la columna 'canonical'
    if (!"canonical" %in% names(df)) {
      df$canonical <- 0
    }
    
    canonical_df <- df[df$canonical == 1, ]
    
    if (nrow(canonical_df) == 0) {
      selected <- df[1, ]
    } else {
      # Verificar que exista 'mane_select'
      if (!"mane_select" %in% names(canonical_df)) {
        selected <- canonical_df[1, ]
      } else {
        mane_rows <- canonical_df[!is.na(canonical_df$mane_select) & 
                                    canonical_df$mane_select != "", ]
        if (nrow(mane_rows) > 0) {
          # Para genomic input, no tenemos un input_transcript específico
          # Así que tomamos el primero de MANE
          selected <- mane_rows[1, ]
        } else {
          selected <- canonical_df[1, ]
        }
      }
    }
  }
  
  # 7. Función auxiliar para extraer campos (idéntica a parse_hgvsc)
  get_value <- function(field, nested_field = NULL) {
    
    if (is.null(selected)) return("-")
    
    # Valor directo
    if (field %in% names(selected) && is.null(nested_field) &&
        !is.null(selected[[field]]) && !is.na(selected[[field]])) {
      return(selected[[field]])
    }
    
    # Valor anidado
    else if (!is.null(nested_field) && field %in% names(selected)) {
      obj <- selected[[field]]
      
      if (is.data.frame(obj) && nested_field %in% names(obj)) {
        value <- obj[[nested_field]]
        if (!is.null(value) && !is.na(value)) {
          return(value)
        }
      }
    }
    return("-")
  }
  
  # 8. Extraer todos los campos predictivos
  gen <- get_value("gene_symbol")
  alpha_missense <- get_value("alphamissense", "am_pathogenicity")
  alpha_pathogenicity <- get_value("alphamissense", "am_class") 
  cadd_phred <- get_value("cadd_phred")
  clinpred_score <- get_value("clinpred")
  revel_score <- get_value("revel")
  sift_score <- get_value("sift_score")
  polyphen_score <- get_value("polyphen_score")
  polyphen_predict <- get_value("polyphen_prediction")
  enformer_score <- get_value("enformer_sar")
  
  # 9. Retornar el resultado en el mismo formato que parse_hgvsc
  return(list(CHR = chr, 
              POS = pos, 
              REF = ref, 
              ALT = alt, 
              Gen = gen,
              "Peor Consecuencia" = most_severe,
              Alphamissense = alpha_missense,
              Alpha.predict = alpha_pathogenicity,
              CADD = cadd_phred,
              Clinpred = ifelse(clinpred_score == "-", "-", round(as.numeric(clinpred_score), 3)),
              REVEL = revel_score,
              SIFT = sift_score,
              Polyph. = polyphen_score,
              "Polyph predict" = polyphen_predict,
              Enformer = enformer_score))
}


#Funcion validacion formato alelo
validate_allele <- function(allele) {
  grepl("^[ACGT-]+$", allele)
}


#Funcion validacion con GRCh38
validate_ref_allele <- function(chr, pos, ref_input) {
  
  if (ref_input == "-") return(TRUE)
  
  tryCatch({
    
    server <- "https://rest.ensembl.org"
    end_pos <- pos + nchar(ref_input) - 1
    region <- paste0(chr, ":", pos, "..", end_pos)
    ext <- paste0("/sequence/region/human/", region, "?")
    
    r <- httr::GET(
      paste0(server, ext),
      httr::content_type("application/json"),
      httr::timeout(10)
    )
    
    httr::stop_for_status(r)
    
    res <- httr::content(r, as = "parsed")
    
    ref_genome <- toupper(res$seq)
    
    return(ref_genome == toupper(ref_input))
    
  }, error = function(e) {
    message("Error validando REF: ", e$message)
    return(NA)
  })
}




#Transforma chr, pos, window a string "10:1,110,692..1,110,710"
make_region <- function(chr, pos, window) {
  pos <- as.numeric(pos)
  
  start <- pos - (window/2)
  end <- pos + (window/2)
  
  start_str <- format(start, big.mark = ",", scientific = FALSE)
  end_str   <- format(end,   big.mark = ",", scientific = FALSE)
  
  chr <- sub("^chr", "",chr, ignore.case = TRUE)
  
  str <-paste0(chr, ":", start_str, "..", end_str)
  return (str)
}


