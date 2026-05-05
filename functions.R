
## TODO: Falta separar validaciones de functions.R


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




# Parsear HGVSc (NM_000059.4:c.-152T>C) en GRCh38
# Con anotación de Ensembl API
parse_hgvsc <- function(hgvsc) {
  
  tryCatch({
    
    server <- "https://rest.ensembl.org"
    ext <- paste0("/vep/human/hgvs/", hgvsc, "?")
    
    r <- httr::GET(
      paste0(server, ext),
      httr::content_type("application/json"),
      httr::timeout(10)
    )
    
    httr::stop_for_status(r)
    
    res <- httr::content(r, as = "parsed", simplifyVector = TRUE)
    
    if (length(res) == 0) return(NULL)
    
    chr <- res$seq_region_name
    pos <- res$start
    alleles <- identify_rest(res$allele_string)
    
    return(list(CHR = chr, POS = pos, REF = alleles[1], ALT = alleles[2]))
    
  }, error = function(e) {
    message("Error en parse_hgvsc: ", e$message)
    return(NULL)
  })
}




# Función para parsear formato genómico "13 32889692 T/C"
# Únicamente funciona con formato: "Chr pos ref/alt", vacío = "-"
parse_genomic_input <- function(genomic_str) {

  parts <- strsplit(trimws(genomic_str), "\\s+")[[1]]
  
  if (length(parts) != 3) {
    return(NULL)
  }
  
  chr <- toupper(parts[1])
  pos <- suppressWarnings(as.numeric(parts[2]))
  alleles <- parts[3]
  items <- strsplit(alleles, "/")[[1]]
  ref <- toupper(items[1])
  alt <- toupper(items[2])
  
  #Validacion posicion numerica
  if (is.na(pos)) {
    showNotification("Posición debe ser numérico",
                     type= "error", duration = 5)
    return(NULL)
    
  }
  
  #Validacion formato chr
  v_chr <- grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)$", chr, ignore.case = TRUE)
  
  if (!v_chr) {
    showNotification("Cromosoma inválido. Use 1-22, X, Y o MT",
                     type= "error", duration = 5)
    return(NULL)
  }
  
  #Validacion formato alelo
  if (!validate_allele(ref) || !validate_allele(alt)) {
    showNotification("REF o ALT contienen caracteres inválidos",
                     type = "error", duration = 5)
    return(NULL)
  } else if(ref == alt){
    showNotification("REF y ALT deben ser diferentes",
                     type = "error", duration = 5)
    return(NULL)
  }
  
  v <- validate_ref_allele(chr, pos, ref)
  
  #validacion con GRCh38
  if (!v || is.na(v)) {
    showNotification("El nucleótido/s referencia no coincide/n con la posición en GRCh38", 
                     type = "error", duration = 5)
    return(NULL)
  }

  return(list(CHR = chr, POS = pos, REF = ref, ALT = alt))
}


validate_allele <- function(allele) {
  grepl("^[ACGT-]+$", allele)
}


validate_ref_allele <- function(chr, pos, ref_input) {
  
  # Si inserción, no validar
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


