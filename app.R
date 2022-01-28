#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

ui <-  fluidPage(
  h1("EasyFlowJo"),
  fileInput("ColDatafile", label = h3("ColData File input")),
  fileInput("csvFiles", label = h3(".csv File inputs"),multiple=TRUE),
  uiOutput("download"),
  tableOutput("testTable")
)

# Define server logic required to draw a histogram
server <- function( input, output, session ) {

  library(magrittr)

  mergeDataTables <- function(colData,csvList){
    AddFileNameAndGate <- function(csvFilename){
      csv <- csvList[[csvFilename]]
      # add column with filename
      csv$csvFilename <- csvFilename
      # add column with original .fcs filename
      fcsFilename <- gsub("([^_]+$)",".fcs",ignore.case=T,perl=T,csvFilename)
      csv$fcsFilename <- gsub("_.fcs",".fcs",ignore.case=T,perl=T,fcsFilename)
      csv$fcsFilename <- gsub("export_","",ignore.case=T,perl=T,csv$fcsFilename)
      # get gate from text after last underscore
      gate <- gsub(".*(?=_)","",perl=T,csvFilename)
      gate <- gsub(".csv","",gate)
      gate <- gsub("_","",gate)
      csv$gate <- gate
      return(csv)
    }
    ParentGate <- unique(colData$ParentGate)
    G1gate <- unique(colData$G1gate)
    Sgate <- unique(colData$Sgate)

    makeSubsetFilename <- function(filename,Gate){
      basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename))
      CSVFilename <- paste0("export_",basename,"_",Gate,".csv")
      return(CSVFilename)
    }

    AllFCSFileNames <- unique(c(colData$Control,colData$Treated))
    AllParentCSVFilenames <- sapply(AllFCSFileNames,makeSubsetFilename,ParentGate)
    AllG1CSVFilenames <- sapply(AllFCSFileNames,makeSubsetFilename,G1gate)
    AllSCSVFilenames <- sapply(AllFCSFileNames,makeSubsetFilename,Sgate)

    ParentCSVs <- lapply(AllParentCSVFilenames,AddFileNameAndGate)
    G1_times <- lapply(AllG1CSVFilenames,function(filename){csvList[[filename]] %>% .$Time})
    S_times <- lapply(AllSCSVFilenames,function(filename){csvList[[filename]] %>% .$Time})

    ParentCSVsWithGates <- lapply(ParentCSVs,function(csv){
      csvFilename <- unique(csv$fcsFilename)
      All_times <- csv$Time
      G1times <- G1_times[[csvFilename]]
      S_times <- S_times[[csvFilename]]
      csv$gate[na.omit(match(S_times,All_times))] <- "S"
      csv$gate[na.omit(match(G1times,All_times))] <- "G1"
      return(csv)
      #(na.omit(match(All_times,G1times)))
    })

    return(ParentCSVsWithGates)
  }

  # Your application server logic
  # read colData file into reactive
  coldata <- reactive({read.csv(input$ColDatafile$datapath)})
  # make a list of all csv data
  csvList <- reactive({
    req(input$csvFiles)
    csvList <- lapply(input$csvFiles$datapath,read.csv)
    names(csvList) <- input$csvFiles$name
    return(csvList)
  })
  csvList2 <- reactive({
    req(input$csvFiles)
    csvList <- lapply(csvList(),function(df){
      cdata <- coldata()
      DNAChannel <- gsub("-",".",cdata$DNAchanel[1])
      YChannel <- gsub("-",".",cdata$yaxischannel[1])
      names(df) <- gsub(DNAChannel,"x",names(df))
      names(df) <- gsub(YChannel,"y",names(df))
      return(df)
    })
    return(csvList)
  })


  # merge the data tables. this reactive makes a list of dataframes with all channel values and that csvFilename, fcsFilename, and gate columns
  mergedDataTables <- reactive({
    req(csvList())
    req(input$ColDatafile)
    mergeDataTables(coldata(),
                                csvList2())
  })
  #generate linear models from the control data. This reactive makes a list of linear models.
  linearModelsForBackground <- reactive({
    req(coldata())
    req(mergedDataTables())
    lmList <- lapply(coldata()$Control,function(cntrlName){
      df <- mergedDataTables()[[cntrlName]]
      lm(y~x,data=df)
    })
    names(lmList) <- coldata()$Treated
    return(lmList)
  })

  backgroundSubtracted <- reactive({
    req(coldata())
    req(mergedDataTables())
    req(linearModelsForBackground())
    backSubList <- lapply(coldata()$Treated,function(treatName){
      df <- mergedDataTables()[[treatName]]
      linMod <- linearModelsForBackground()[[treatName]]
      df$y2 <- df$y - (predict(linMod,data.frame(x=df$x)))
      return(df)
    })
    names(backSubList) <- coldata()$Treated
    return(backSubList)
  })

  output$downloadH <- downloadHandler(
    filename = function() {
      paste0("dataList", ".rds")
    },
    content = function(file) {
      saveRDS(backgroundSubtracted(), file)
    }
  )

  output$download<-renderUI({
    req(backgroundSubtracted())
    downloadButton("downloadH",label="download .rds")
  })

  output$testTable <- renderTable({
    req(backgroundSubtracted())
    head(backgroundSubtracted()[[1]])})
}

# Run the application
shinyApp(ui = ui, server = server)
