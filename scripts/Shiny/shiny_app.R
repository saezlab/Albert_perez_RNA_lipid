library(shiny)
library(shinyjs)
library(shinythemes)
library(shinybusy)
options(repos = BiocManager::repositories())


ui = fluidPage(
  useShinyjs(),
  theme = shinytheme("cerulean"),
  navbarPage("Albert P data", tabPanel("Home",sidebarPanel(
    textInput(inputId = "TF", label = "Select TF", value = "All"),
    numericInput(inputId = "rank", label = "TF rank", value = 1, min = 1, max=270),
    numericInput(inputId = "likelihood", label = "Target-Lipid likelihood threshold", value = 0.9, min = 0, max =1, step = 0.01),
    selectInput(inputId = "sort", 
                label = "Sort by:", 
                choices = c("Likelihood","TF rank")),
    actionButton("submit", "Generate data", class="btn btn-primary"),
    actionButton("reset", "Reset")
    
  ), mainPanel(add_busy_spinner(spin = "double-bounce", position = 'bottom-right',color = "red"),tabsetPanel(type = "tabs",
                           tabPanel("Download tab",htmlOutput('main')
                                    
                                   
                                    ,
                                    downloadButton(outputId = "sub_interactomes", "Download data", class = "btn-block"))
                           
                           
  
  ))),
  tabPanel("About", 
           div(includeMarkdown("about.Rmd")),
           h4("E-Mail: anis.mansouri.dz@gmail.com.")))
  
)#,






server = function(input, output, session) {
  all_data=reactive({
    data_mat=readRDS(url("https://zenodo.org/record/4736866/files/final_matrix.rds?download=1"))
    #data_mat=readRDS(url("https://zenodo.org/record/4737255/files/tst.rds?download=1"))
    
    return(list("my_data"=data_mat))
    
    })
    
    
  output$main=renderPrint({
    if (input$submit>0) {isolate(all_data())} else {return(h4("Server is ready to generate the data"))}
    shinyjs::show("sub_interactomes")
    })
  
  
  
  
  
  proc_data=reactive({
    
   
    ctrl_df=all_data()$my_data[, c(1,5:10)]
    OA_df=all_data()$my_data[, c(2,5:10)]
    PA_df=all_data()$my_data[, c(3,5:10)]
    PO_df=all_data()$my_data[, c(4,5:10)]
    
    if (input$TF=="All") {
    
    if (input$rank==1) {
      ctrl_df=ctrl_df[ctrl_df$Ctrl == 1,]
      OA_df=OA_df[OA_df$OA == 1,]
      PA_df=PA_df[PA_df$PA == 1,]
      PO_df=PO_df[PO_df$PO == 1,]
    }
    else {
    
    ctrl_df=ctrl_df[ctrl_df$Ctrl %in% c(1:input$rank),]
    OA_df=OA_df[OA_df$OA %in% c(1:input$rank),]
    PA_df=PA_df[PA_df$PA %in% c(1:input$rank),]
    PO_df=PO_df[PO_df$PO %in% c(1:input$rank),]}
    
    ctrl_df=ctrl_df[ctrl_df$likelihood >= input$likelihood,]
    OA_df=OA_df[OA_df$likelihood >= input$likelihood,]
    PA_df=PA_df[PA_df$likelihood >= input$likelihood,]
    PO_df=PO_df[PO_df$likelihood >= input$likelihood,]
    
    if (input$sort == "Likelihood") {
      ctrl_df=ctrl_df[order(ctrl_df$likelihood, decreasing=TRUE),]
      OA_df=OA_df[order(OA_df$likelihood, decreasing=TRUE),]
      PA_df=PA_df[order(PA_df$likelihood, decreasing=TRUE),]
      PO_df=PO_df[order(PO_df$likelihood, decreasing=TRUE),]}
    else {
      ctrl_df=ctrl_df[order(ctrl_df$Ctrl, decreasing=FALSE),]
      OA_df=OA_df[order(OA_df$OA, decreasing=FALSE),]
      PA_df=PA_df[order(PA_df$PA, decreasing=FALSE),]
      PO_df=PO_df[order(PO_df$PO, decreasing=FALSE),]
      
    }
    
    
    }
    
    else {
      ctrl_df=ctrl_df[ctrl_df$tf %in% input$TF,]
      OA_df=OA_df[OA_df$tf %in% input$TF,]
      PA_df=PA_df[PA_df$tf %in% input$TF,]
      PO_df=PO_df[PO_df$tf %in% input$TF,]
      
      ctrl_df=ctrl_df[order(ctrl_df$likelihood, decreasing=TRUE),]
      OA_df=OA_df[order(OA_df$likelihood, decreasing=TRUE),]
      PA_df=PA_df[order(PA_df$likelihood, decreasing=TRUE),]
      PO_df=PO_df[order(PO_df$likelihood, decreasing=TRUE),]
      
      
    }
    
    
    
    
    
    
    return(list("ctrl"=ctrl_df, "OA"=OA_df, "PA"= PA_df, "PO"=PO_df))
    
    })
  
  observeEvent(input$reset, {
    shinyjs::reset("TF")
    shinyjs::reset("rank")
    shinyjs::reset("likelihood")
    shinyjs::reset("sort")
    shinyjs::reset("conditions")
  })
  
  

   
    observe({
      shinyjs::hide("sub_interactomes")})
  
    
   output$sub_interactomes <- downloadHandler(
      filename = function() {
     paste("interactome", ".zip", sep = "")
      },
      content = function(file) {

          write.table(proc_data()$ctrl, file="Ctrl.txt", quote = FALSE, row.names = FALSE)
          write.table(proc_data()$OA, file="OA.txt", quote = FALSE, row.names = FALSE)
          write.table(proc_data()$PA, file="PA.txt", quote = FALSE, row.names = FALSE)
          write.table(proc_data()$PO, file="PO.txt", quote = FALSE, row.names = FALSE)
          files = c("Ctrl.txt","OA.txt","PA.txt","PO.txt")
        
        
        zip(file,files)
      
      }
    )

  
}
shinyApp(ui, server)
