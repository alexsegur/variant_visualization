library(shiny)
library(shinydashboard)
library(DT)
library(VariantAnnotation)
library(GenomicRanges)
library(JBrowseR)
library(httr)
library(jsonlite)
source("ui.R",local = TRUE)
source("server.R", local =TRUE)


# Lanzar la app
shinyApp(ui = ui, server = server)
