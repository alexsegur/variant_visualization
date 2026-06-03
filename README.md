# Variant Browser
 
Aplicación web desarrollada en R/Shiny que permite buscar y visualizar variantes genómicas en contexto genómico mediante JBrowseR, con soporte para entrada en formato HGVSc y genómico, y carga de archivos VCF.
 
---
 
## Requisitos previos
 
- **R 4.5.0** → https://cran.r-project.org/
- **RStudio 2026.01** → https://posit.co/download/rstudio-desktop/
---
 
## Instalación y puesta en marcha
 
### 1. Descargar el proyecto
 
Clona el repositorio o descárgalo como ZIP desde GitHub:
 
```bash
git clone https://github.com/alexsegur/variant_visualization.git
```
 
O bien: botón verde **Code → Download ZIP**, luego descomprimir.
 
### 2. Abrir como proyecto en RStudio
 
- En RStudio: `File → Open Project...`
- Navegar hasta la carpeta del proyecto y seleccionar el archivo `variant_visualization.Rproj`
> Es importante abrirlo como proyecto (no solo abrir los archivos `.R` sueltos) para que `renv` funcione correctamente.
 
### 3. Restaurar el entorno con `renv`
 
Al abrir el proyecto por primera vez, `renv` debería activarse automáticamente. En la consola de R ejecuta:
 
```r
renv::restore()
```
 
Esto instalará todas las dependencias en las versiones exactas registradas en `renv.lock`.
 
### 4. Paquetes que pueden requerir instalación manual
 
Algunos paquetes de Bioconductor no siempre se resuelven correctamente a través de `renv`. Si tras el paso anterior aparece algún error al lanzar la app, instálalos manualmente:
 
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
BiocManager::install(c(
  "VariantAnnotation",
  "GenomicRanges",
  "Rsamtools"
))
```
 
El paquete `JBrowseR` también puede requerir instalación manual si no se resuelve con `renv`:
 
```r
install.packages("JBrowseR")
```
 
### 5. Lanzar la aplicación
 
Desde la consola de R, con el proyecto abierto:
 
```r
shiny::runApp()
```
 
O bien abre `app.R` en RStudio y pulsa el botón **Run App** en la esquina superior derecha del editor.
 
---
 
## Uso de la aplicación
 
La interfaz tiene dos funcionalidades principales accesibles desde la barra lateral izquierda:
 
### Búsqueda de variante
 
- **HGVSc**: introduce por ejemplo `NM_000059.4:c.-152T>C`.
- **Formato genómico**: introduce cromosoma, posición y alelos separados por espacios, por ejemplo `13 32889692 T/C`.
- **Entorno (pb)**: número de pares de bases a mostrar alrededor de la variante en el navegador (entre 20 y 1000).

### Carga de archivos VCF
 
Permite cargar un VCF propio para visualizar todas sus variantes en la tabla y navegar entre ellas.
 
- **VCF plano (`.vcf`)**
- **VCF comprimido (`.vcf.gz` + `.vcf.gz.tbi`)**
 
---
 
## Notas
 
- La aplicación requiere conexión a internet.
- El genoma de referencia utilizado es **GRCh38/hg38**.
