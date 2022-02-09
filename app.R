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
ui <- fluidPage(
  tags$button(
    id = "close",
    type = "button",
    class = "btn action-button",
    onclick = "setTimeout(function(){window.close();},500);", # close browser
    "Close window"
  ),
  tabsetPanel(
    tabPanel(
      "Background Subtract",
      h1("EasyFlowJo"),
      sidebarLayout(
        sidebarPanel(fileInput("ColDatafile", label = h3("ColData File input"))),
        mainPanel(wellPanel(
          h4("the colData .csv should have the following columns:"),
          h5("Control,Treated,yaxischannel,DNAchannel,ParentGate"),
          img(src = "colDataExample.png")
        ))
      ),
      sidebarLayout(
        sidebarPanel(
          fileInput("csvFiles", label = h3(".csv File inputs"), multiple = TRUE)
        ),
        mainPanel(
          wellPanel(
            h4(".csv filenames must have the following format:"),
            uiOutput(outputId = "csvFileText"),
            h4("For example, the \"G1\" gate for \"myFile.fcs\" would be:"),
            uiOutput(outputId = "csvFileTextExample")
          )
        )
      ),
      verbatimTextOutput("parentGate"),
      verbatimTextOutput("gates"),
      uiOutput("download"),
      uiOutput("downloadCode"),
      uiOutput("columnSelection"),
      tableOutput("testTable")
    ),
    tabPanel(
      "Instructions",
      h1("EasyFlowJo Usage Instructions"),
      h2("Motivation"),
      p("The purpose of this shiny app is to easily calculate background subtracted values from bivariate cell cycle flow cytometry data. Once the conversions are made, you can save the data in the form of a .rds file, which can be loaded into R for statistical analysis and plotting."),
      h2("Test Data"),
      p("You may access test data from github."),
      p("https://github.com/SansamLab/EasyFlowJo"),
      h2("How To Use"),
      h3("Step 1:  Gate your data in flowjo."),
      img(src = "image001.jpg"),
      h3("Step 2:  Export scale values from flowjo to .csv files"),
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
      img(src = "image015.png")
    ),
    tabPanel("Rendered Code", verbatimTextOutput("code"))
  )
)
###########

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  library(magrittr)

  observe({
    if (input$close > 0) stopApp()                             # stop shiny
  })

  # read colData file
  coldata <- metaReactive2({
    req(input$ColDatafile$datapath)
    metaExpr({
      read.csv(..(input$ColDatafile$datapath))
    })
  })

  # make a list of all csv data
  csvListLoad <- metaReactive2({
    req(input$csvFiles)
    req(coldata())
    metaExpr({
      cdata <- ..(coldata())
      # read csv's
      csvList <- lapply(..(input$csvFiles$datapath), read.csv)})
    csvList
  })

  csvList <- metaReactive2({
    req(input$csvFiles)
    req(coldata())
    req(csvListLoad())
    metaExpr({
      cdata <- ..(coldata())
      # read csv's
      csvList <- ..(csvListLoad()) %>%
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
        DNAChannel2 <- gsub("-", ".", cdata$DNAchannel[1])
        DNAChannel <- paste0(DNAChannel,DNAChannel2)
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
      gsub(".csv","",..(input$csvFiles$name),ignore.case=T)%>%
        gsub("^(.*[_])","",.)%>%unique()%>%.[-which(.==..(parentGate()))]
    }))
  })

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
      saveRDS(backgroundSubtractedSubset(), file)
    }
  )

  output$download <- renderUI({
    req(backgroundSubtracted())
    downloadButton("downloadH", label = "download .rds")
  })

  output$testTable <- renderTable({
    req(backgroundSubtractedSubset())
    head(backgroundSubtractedSubset()[[1]],5)
  })

  code <- reactive({
    # swap out the implementation of coldata for ec
    ec <- newExpansionContext()
    ec$substituteMetaReactive(coldata, function() {
      metaExpr(read.csv(..(input$ColDatafile$name)))
    })
    ec$substituteMetaReactive(csvListLoad, function() {
      metaExpr({csvList <- lapply(..(input$csvFiles$name), read.csv)
      names(csvList) <- ..(input$csvFiles$name)
      csvList})
    })
    formatCode(expandChain(
      .expansionContext = ec,
      quote(library(magrittr)),
"#################################################################",
"##                    load colData file                        ##",
"#################################################################",
      invisible(coldata()),
"#################################################################",
"##                      load csv files                         ##",
"#################################################################",
      invisible(csvListLoad()),
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
      invisible(backgroundSubtractedSubset()),
      quote(saveRDS(backgroundSubtractedSubset,"backGroundSubtracted.rds"))
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
    req(backgroundSubtractedSubset())
    downloadButton("downloadCodeH", label = "download code")
  })


  output$downloadCodeH <- downloadHandler(
    filename = function() {
      paste("script-", Sys.Date(), ".zip", sep="")
    },
    content = function(file) {
      buildScriptBundle(
        code(), file, render=FALSE
      )
    }
  )

  output$csvFileText <- renderText({
    HTML(paste0("<b style='color:red'>","export","</b>", "_fcsFileBaseName_","<b style='color:green'>","gate","</b>",".csv"))
  })
  output$csvFileTextExample <- renderText({
    HTML(paste0("<b style='color:red'>","export","</b>", "_myFile_","<b style='color:green'>","G1","</b>",".csv"))
  })
  output$columnSelection <- renderUI({
    req(backgroundSubtracted())
    selectizeInput('columnSelect',
                   label="Select the columns to export",
                   choices=names(backgroundSubtracted()[[1]]),
                   selected = c("csvFilename","fcsFilename","gate","x","y","y2"),
                   multiple = TRUE,
                   options = NULL)
  })
  backgroundSubtractedSubset <- metaReactive2({
    req(output$columnSelection)
    metaExpr({
      lapply(..(backgroundSubtracted()),function(df){df[..(input$columnSelect)]})
    })
  })


}


# Run the application
shinyApp(ui = ui, server = server)
