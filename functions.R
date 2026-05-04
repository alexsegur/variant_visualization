

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
# Únicamente funciona con formato: "Chr pos ref/alt"
parse_genomic_input <- function(genomic_str) {

  parts <- strsplit(trimws(genomic_str), "\\s+")[[1]]
  
  if (length(parts) != 3) {
    return(NULL)
  }
  
  chr <- parts[1]
  pos <- as.numeric(parts[2])
  alleles <- identify_rest(parts[3])
  return(list(CHR = chr, POS = pos, REF = alleles[1], ALT = alleles[2]))
}


identify_rest <- function(rest){
  
  items <- strsplit(rest, "/")[[1]]
  ref <- items[1]
  alt <- items[2]
  
  return(c(ref,alt))
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


