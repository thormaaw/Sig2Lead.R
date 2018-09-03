#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(visNetwork)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
  
  # Application title
  titlePanel("Sig2Lead"),
      fluidPage(
        tabsetPanel(
          tabPanel("Search", value ="search", icon = NULL,
                   # fluidRow(
                     br(),
                            selectInput("Signature", NULL, list("Input a Gene", "Upload a Signature"), selected = "Input a Gene"),
                     br(),
                     
                     conditionalPanel(
                       condition = "input.Signature == 'Input a Gene'",
                            textInput("gene_knockdown", "Input a Gene", value="bcl2a1")
                     ),
                     conditionalPanel(
                        condition = "input.Signature == 'Upload a Signature'",
                            fileInput("UploadSignature", "Upload Signature", multiple = FALSE, accept = c(".csv"), width = NULL, 
                                         buttonLabel = "Browse...", placeholder = "No file selected")
                     ),
                            br(),
                     fileInput("AddedCompounds", "Add compounds in SMILES (Optional)", multiple = TRUE, accept = c(".smi", ".csv"), width = NULL,
                               buttonLabel = "Browse...", placeholder = "No file selected"),
                               textInput("AddedLabel", "Label for added compounds", value="Added"),
                            br(),
                   selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
                            actionButton(label="Go!", "Go")
           # )
        ),
        
#        mainPanel(
#        textOutput("ItWorked"),
          tabPanel("LINCS Compounds", value = "LINCS Compounds", icon =NULL,
          DT::dataTableOutput("SMILES"),
          downloadButton("SMILESDownload", label = "Download SMILES"),
          downloadButton("SDFDownload", label = "Download SDF")
        ),
      tabPanel("Heatmap", value="heatmap", icon = NULL, 
#               actionButton(label="GenerateHeatmap", "GenerateHeatmap"),
               textInput("CutHeight", "Tanimoto Similarity", value="0.75"),
               textInput("ClusterSize", "Minimum Cluster Size", value="3"),
               actionButton(label="Get MDS Plot", "GetRepresentatives"),
               plotOutput("distPlot", width = "850px", height = "800px")
               
              
        ),
      tabPanel("MDS Plot", value="MDS", icon=NULL,
               plotOutput("MDSPlot", width = "900px", height = "800px"),
               DT::dataTableOutput("Representatives"),
               downloadButton("RepDownload", label="Download Representatives"),
               downloadButton("ClusterDownload", label = "Download Clusters"),
               actionButton(label="Get Related NCI Compounds", "GetNCI")
               ),

      tabPanel("Similar Compounds", value="NCI", icon=NULL,
              DT::dataTableOutput("Similar_NCI"),
              downloadButton("NCIDownload", label = "Download Similar NCI")),

      tabPanel("STITCH Network", value = "STITCH", icon = NULL,
               visNetworkOutput("STITCHPlot", width = "900px", height = "800px"),
               textInput("ClusterNumber", "Cluster Number"),
               textInput("Connections", "Number of Connections"),
               textInput("Gene", "Gene of Interest", value = "bcl2a1"),
               actionButton(label="View with STITCH", "GetSTITCH")
               )
        )
      )
  )
)


