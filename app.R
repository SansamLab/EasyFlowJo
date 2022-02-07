#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinymeta)

############
ui <-  fluidPage(
  tabsetPanel(
   tabPanel("Background Subtract",
     h1("EasyFlowJo"),
  fileInput("ColDatafile", label = h3("ColData File input")),
  fileInput("csvFiles", label = h3(".csv File inputs"),multiple=TRUE),
  verbatimTextOutput("parentGate"),
  verbatimTextOutput("gates"),
  uiOutput("download"),
  uiOutput("downloadCode"),
  tableOutput("testTable")),
  tabPanel("Instructions",
    h1("EasyFlowJo Usage Instructions"),
    h2("Motivation"),
    p("The purpose of this shiny app is to easily calculate background subtracted values from bivariate cell cycle flow cytometry data. Once the conversions are made, you can save the data in the form of a .rds file, which can be loaded into R for statistical analysis and plotting."),
    h2("Test Data"),
    p("You may access test data from github."),
    p("https://github.com/SansamLab/EasyFlowJo"),
    h2("How To Use"),
    h3("Step 1:  Gate your data in flowjo."),
    img(src = "image001.jpg"),
    h3("Step 2:  Export channel values from flowjo to .csv files"),
    h4("2a.  Select the all parent, G1, and  S populations"),
    p("Select the parent gate (Single Cells), G1, and S.  Click select equivalent nodes to select those gates for all samples."),
    img(src = "image003.png"),
    h4("2b.  Export the channel values"),
    p("After parent, G1, and S gates are selected for all samples, click Export/Concatenate Populations."),
    img(src = "image005.png"),
    img(src = "image007.png"),
    h3("Step 3:  Create a colData file describing your data"),
    img(src = "image009.png"),
    h3("Step 4:  Load colData into shiny app"),
    img(src = "image011.png"),
    h3("Step 5:  Load channel value .csv files into shiny app"),
    img(src = "image013.png"),
    h3("Step 6:  Download .rds file"),
    img(src = "image015.png")),
  tabPanel("Rendered Code",verbatimTextOutput("code"))
  )
)
###########

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  library(magrittr)

  # read colData file
  coldata <- metaReactive2({
    req(input$ColDatafile$datapath)
    metaExpr({
      read.csv(..(input$ColDatafile$datapath))
    })
  })

  # make a list of all csv data
  csvList <- metaReactive2({
    req(input$csvFiles)
    req(coldata())
    metaExpr({
      cdata <- ..(coldata())
      # read csv's
      csvList <- lapply(..(input$csvFiles$datapath), read.csv) %>%
      # add ID column
      lapply(.,function(df){df$ID <- do.call(paste, c(df, sep = "")); df})
      # name list items
      names(csvList) <- ..(input$csvFiles$name)
      # add filename and gate columns
      csvList <- lapply(names(csvList),function(nme){
        df <- csvList[[nme]]
        df$gate <- gsub(".csv","",nme,ignore.case=T)%>%gsub("^(.*[_])","",.)
        df$csvFilename <- nme
        # add column with original .fcs filename
        fcsFilename <- gsub("([^_]+$)", ".fcs", ignore.case = T, perl = T, nme)
        df$fcsFilename <- gsub("_.fcs", ".fcs", ignore.case = T, perl = T, fcsFilename)
        df$fcsFilename <- gsub("export_", "", ignore.case = T, perl = T, df$fcsFilename)
        df
      })
      # rename DNA and y axis columns to x and y
      csvList <- lapply(csvList,function(df){
        DNAChannel <- gsub("-", ".", cdata$DNAchanel[1])
        YChannel <- gsub("-", ".", cdata$yaxischannel[1])
        names(df) <- gsub(DNAChannel, "x", names(df))
        names(df) <- gsub(YChannel, "y", names(df))
        df
      })
      # name list items
      names(csvList) <- ..(input$csvFiles$name)
      csvList
    })
  })

  parentGate <- metaReactive2({
    req(coldata())
    metaExpr({unique(..(coldata())$ParentGate)})
  })

  gateList <- metaReactive2({
    req(input$csvFiles)
    metaExpr(({
      gsub(".csv","",input$csvFiles$name,ignore.case=T)%>%
        gsub("^(.*[_])","",.)%>%unique()%>%.[-which(.==..(parentGate()))]
    }))
  })

  # addIDColumn <- function(csvList){
  #   lapply(csvList,function(df){
  #     df$ID <- do.call(paste, c(df, sep = ""))
  #     df
  #   })
  # }
  #
  # addGateColumn <- function(csvList){
  #   lapply(names(csvList),function(nme){
  #     df <- csvList[[nme]]
  #     df$gate <- gsub(".csv","",nme,ignore.case=T)%>%gsub("^(.*[_])","",.)
  #     df
  #   })
  # }


  # csvList2 <- metaReactive2({
  #   req(input$csvFiles)
  #   metaExpr({
  #     csvList <- lapply(..(csvList()), function(df) {
  #       cdata <- ..(coldata())
  #       # add ID column
  #       df$ID <- do.call(paste, c(df, sep = ""))
  #       DNAChannel <- gsub("-", ".", cdata$DNAchanel[1])
  #       YChannel <- gsub("-", ".", cdata$yaxischannel[1])
  #       names(df) <- gsub(DNAChannel, "x", names(df))
  #       names(df) <- gsub(YChannel, "y", names(df))
  #       df
  #     })
  #     csvList
  #   })
  # })


  # merge the data tables. this metaReactive2 makes a list of dataframes with all channel values and that csvFilename, fcsFilename, and gate columns
  mergedDataTables <- metaReactive2({
    req(csvList())
    req(coldata())
    req(parentGate())
    req(gateList())
    metaExpr({
      # split list of csv's by fcs filename
      fcsFilenames <- sapply(..(csvList()),function(df){df$fcsFilename[1]})
      splitCsvList <- split(..(csvList()),f=fcsFilenames)
      # apply merge over each fcsFilename set
      mergedCsvList <- lapply(splitCsvList,function(lst){
        # make vector of gates
        gates <- sapply(lst,function(df){df$gate[1]})
        # get parent gate data
        pGateCsv <- lst[[which(gates==..(parentGate()))]]
        # substitute gates for each of the child gates
        for (gte in ..(gateList())) {
          pGateCsv$gate[match(lst[[which(gates==gte)]]$ID, pGateCsv$ID)] <- gte
        }
        pGateCsv
      })
      mergedCsvList
      })
  })

  # generate linear models from the control data. This metaReactive2 makes a list of linear models.
  linearModelsForBackground <- metaReactive2({
    req(coldata())
    req(mergedDataTables())
    metaExpr({
      lmList <- lapply(..(coldata())$Control, function(cntrlName) {
        df <- ..(mergedDataTables())[[cntrlName]]
        lm(y ~ x, data = df)
      })
      names(lmList) <- ..(coldata())$Treated
      lmList
    })
  })


  backgroundSubtracted <- metaReactive2({
    req(coldata())
    req(mergedDataTables())
    req(linearModelsForBackground())
    metaExpr({
      backSubList <- lapply(..(coldata())$Treated, function(treatName) {
        df <- ..(mergedDataTables())[[treatName]]
        linMod <- ..(linearModelsForBackground())[[treatName]]
        df$y2 <- df$y - (predict(linMod, data.frame(x = df$x)))
        df
      })
      names(backSubList) <- ..(coldata())$Treated
      backSubList
    })
  })

  output$downloadH <- downloadHandler(
    filename = function() {
      paste0("dataList", ".rds")
    },
    content = function(file) {
      saveRDS(backgroundSubtracted(), file)
    }
  )

  output$download <- renderUI({
    req(backgroundSubtracted())
    downloadButton("downloadH", label = "download .rds")
  })

  output$testTable <- renderTable({
    req(backgroundSubtracted())
    head(backgroundSubtracted()[[1]])
  })

  code <- reactive({
    # swap out the implementation of coldata for ec
    ec <- newExpansionContext()
    ec$substituteMetaReactive(coldata, function() {
      metaExpr(read.csv(..(input$ColDatafile$name)))
    })
    ec$substituteMetaReactive(csvList, function() {
      metaExpr({csvList <- lapply(..(input$csvFiles$name), read.csv)
      names(csvList) <- ..(input$csvFiles$name)
      csvList})
    })
    formatCode(expandChain(
      .expansionContext = ec,
"#################################################################",
"##                    load colData file                        ##",
"#################################################################",
      invisible(coldata()),
"#################################################################",
"##                      load csv files                         ##",
"#################################################################",
      invisible(csvList()),
"#################################################################",
"##     make single dataframes with G1 and S events marked      ##",
"#################################################################",
      invisible(mergedDataTables()),
"#################################################################",
"##        calculate linear models for background values        ##",
"#################################################################",
      invisible(linearModelsForBackground()),
"#################################################################",
"##            calculate background subtracted values           ##",
"#################################################################",
      invisible(backgroundSubtracted()),
"#################################################################",
"##                     write .rds file                         ##",
"#################################################################",
      quote(saveRDS(backgroundSubtracted,"backGroundSubtracted.rds"))
    ))
  })

  output$parentGate <- renderPrint({
    req(parentGate())
    paste("Parent gate:  ",parentGate(),sep="")})

  output$gates <- renderPrint({
    req(gateList())
    paste("Gate detected:  ",gateList(),sep="")})

  output$code <- renderPrint({code()})

  output$downloadCode <- renderUI({
    req(backgroundSubtracted())
    downloadButton("downloadCodeH", label = "download code")
  })


  output$downloadCodeH <- downloadHandler(
    filename = function() {
      paste("script-", Sys.Date(), ".zip", sep="")
    },
    content = function(file) {
      buildScriptBundle(
        code(), file, render=FALSE
      #,
      #  render_args = list(output_format = "pdf_document")
      )
    }
  )

}


# Run the application
shinyApp(ui = ui, server = server)
