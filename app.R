library(shiny)
source("ui.R",local = TRUE)
source("server.R", local =TRUE)


# Lanzar la app
shinyApp(ui = ui, server = server)
