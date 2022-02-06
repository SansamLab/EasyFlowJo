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
  uiOutput("download"),
  uiOutput("downloadCode"),
  tableOutput("testTable")),
  tabPanel("Instructions",
    h1("EasyFlowJo Usage Instructions"),
    h2("Motivation"),
    p("The purpose of this shiny app is to easily calculate background subtracted values from flow cytometry data. Once the conversions are made, you can save the data in the form of a .rds file, which can be loaded into R for statistical analysis and plotting."),
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
    metaExpr({
      csvList <- lapply(..(input$csvFiles$datapath), read.csv)
      names(csvList) <- ..(input$csvFiles$name)
      csvList
    })
  })

  # swap out the implementation of csvList for ec2
  # ec2 <- newExpansionContext()
  # ec2$substituteMetaReactive(csvList, function() {
  #   pname <- input$ColDatafile$datapath
  #   metaExpr(read.csv("coldata.csv"))
  # })

  csvList2 <- metaReactive2({
    req(input$csvFiles)
    metaExpr({
      csvList <- lapply(..(csvList()), function(df) {
        cdata <- ..(coldata())
        DNAChannel <- gsub("-", ".", cdata$DNAchanel[1])
        YChannel <- gsub("-", ".", cdata$yaxischannel[1])
        names(df) <- gsub(DNAChannel, "x", names(df))
        names(df) <- gsub(YChannel, "y", names(df))
        df
      })
      csvList
    })
  })


  # merge the data tables. this metaReactive2 makes a list of dataframes with all channel values and that csvFilename, fcsFilename, and gate columns
  mergedDataTables <- metaReactive2({
    req(csvList())
    req(input$ColDatafile)
    metaExpr({
      mergeDataTables <- function(colData, csvList) {
        AddFileNameAndGate <- function(csvFilename) {
          csv <- csvList[[csvFilename]]
          # add column with filename
          csv$csvFilename <- csvFilename
          # add column with original .fcs filename
          fcsFilename <- gsub("([^_]+$)", ".fcs", ignore.case = T, perl = T, csvFilename)
          csv$fcsFilename <- gsub("_.fcs", ".fcs", ignore.case = T, perl = T, fcsFilename)
          csv$fcsFilename <- gsub("export_", "", ignore.case = T, perl = T, csv$fcsFilename)
          # get gate from text after last underscore
          gate <- gsub(".*(?=_)", "", perl = T, csvFilename)
          gate <- gsub(".csv", "", gate)
          gate <- gsub("_", "", gate)
          csv$gate <- gate
          csv
        }
        ParentGate <- unique(colData$ParentGate)
        G1gate <- unique(colData$G1gate)
        Sgate <- unique(colData$Sgate)

        makeSubsetFilename <- function(filename, Gate) {
          basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename))
          CSVFilename <- paste0("export_", basename, "_", Gate, ".csv")
          CSVFilename
        }

        AllFCSFileNames <- unique(c(colData$Control, colData$Treated))
        AllParentCSVFilenames <- sapply(AllFCSFileNames, makeSubsetFilename, ParentGate)
        AllG1CSVFilenames <- sapply(AllFCSFileNames, makeSubsetFilename, G1gate)
        AllSCSVFilenames <- sapply(AllFCSFileNames, makeSubsetFilename, Sgate)

        ParentCSVs <- lapply(AllParentCSVFilenames, AddFileNameAndGate)
        G1_times <- lapply(AllG1CSVFilenames, function(filename) {
          csvList[[filename]] %>% .$Time
        })
        S_times <- lapply(AllSCSVFilenames, function(filename) {
          csvList[[filename]] %>% .$Time
        })

        ParentCSVsWithGates <- lapply(ParentCSVs, function(csv) {
          csvFilename <- unique(csv$fcsFilename)
          All_times <- csv$Time
          G1times <- G1_times[[csvFilename]]
          S_times <- S_times[[csvFilename]]
          csv$gate[na.omit(match(S_times, All_times))] <- "S"
          csv$gate[na.omit(match(G1times, All_times))] <- "G1"
          csv
          # (na.omit(match(All_times,G1times)))
        })

        return(ParentCSVsWithGates)
      }
      mergeDataTables(
        ..(coldata()),
        ..(csvList2())
      )
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
      invisible(csvList2()),
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

  output$code <- renderPrint({class(code())})

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
