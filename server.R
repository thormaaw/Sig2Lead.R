#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

#########Load all required libraries here, they cannot be sourced
library(shiny)
library(httr)
library(jsonlite)
library(DT)
#library(heatmaply)
#library(shinyHeatmaply)
library(visNetwork)
library(ChemmineOB)
library(ChemmineR)
library(plyr)
library("dendextend")
library("colorspace")
library(ggforce)
library(rlist)
library(scatterpie)
library(ggrepel)
library(bazar)
library(XML)
library(RCurl)
library(bitops)
library(scrapeR)
library(igraph)


######Discordant or Concordant
######Benchmarking -- Behrouz


shinyServer(function(input, output) {
 observeEvent(input$Go,
#  output$Go <- eventReactive(input$Go,  
    {
      withProgress(message = "Running Sig2Lead", value = 0, {
        source("lib/GeneKD.R", local = TRUE)
        source("lib/ConDisconSignatures.R", local = TRUE)
      #1. Query iLINCS for user input gene
        if (input$Signature == "Input a Gene"){
          define_knockdown(input$gene_knockdown)
        print("Found your knockdowns") 
        incProgress(1/7, message = "Finding Concordant Signatures")
        ConDisconSignatures(sigid)
        }
        else if (input$Signature == "Upload a Signature"){
          source("lib/UploadSignature.R")
          uploadsig <- input$UploadSignature
          sigUpload(paste(getwd(), uploadsig$name, sep="/"))
          print("Found your signature")
          incProgress(1/7, message = "Finding Concordant Signatures")
          UploadConDiscon(rjson)
          
        }
        
        #incProgress(1/7, detail = "Finding Concordant Signatures")
        
      #2. Get a list of all concordant compounds
        
        
        print("Found concordant signatures")
        incProgress(1/7, message = "Obtaining SMILES")
        source("lib/GetSMILES.R", local = TRUE)
      #3. Generate SMILES of compounds
        GetSMILES()
        print("SMILES acquired")
        incProgress(1/7, message = "Converting SMILES to SDF...Slow")
      #4. Generate Data Table
        display <- final_compounds
        colnames(display) <- c("LSM_ID", "SMILES")
        #output$SMILES <<- DT::renderDataTable(display)
        #output$SMILES <<- DT::renderDataTable({DT::datatable(display, selection = list(selection = "single", target = "cell"))}, escape = FALSE)
        
        makeLink <- function(val) {
          #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
          paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
        }
        output$SMILES <- DT::renderDataTable({
          
          #display$link <- makeLink(display$LSM_ID)
          display$LSM_ID <- makeLink(display$LSM_ID)
          return(display)
        }, 
        escape = FALSE, rownames = FALSE)
        
        #http://lincsportal.ccs.miami.edu/SmallMolecules/view/LSM-5467
        display <<- display
        
        source("lib/ChemmineOB.R", local = TRUE)
        
        ###This is the slow step, see if we can go straight from smiles to fingerprints and seperately convert to SDF
        print("Converting SMILES to SDF")
        Get_SDF()
        
        ###
        
        print("SDF conversion complete")
        incProgress(3/7, message = "Clustering Related Compounds")
      #5. Cluster compounds identified through LINCS, generate heatmaps
        source("lib/chemmineR.R", local = TRUE)
        
      #5.1. Include Added Compounds if applicable
         if (!is.null(input$AddedCompounds))      {
                      adds<-input$AddedCompounds
                      adds_SMI<-read.SMIset(paste(getwd(), adds$name, sep="/"))
                      adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
                      adds_csv <- adds_csv[c(2,1)]
                      colnames(adds_csv) <- c("LSM_ID", "SMILES")
                      adds_csv <<- adds_csv
                      test2 <<- rbind(adds_csv, display)
                      colnames(test2) <<- c("Compound_ID", "SMILES")
                      adds_SDF<-smiles2sdf(adds_SMI)
                      sdf_smiles <<- c(sdf_smiles, adds_SDF)
                      cluster_compounds()
                      print("Your compounds have been clustered with LINCS compounds")
                      source("lib/ColorMap.R", local = TRUE)
                      color_dend()
                      #color_dend() currently not rendering
                      #output$distPlot<<- renderPlot(heatmap.2(dist_mat, Rowv=dend, Colv=dend, colRow = heatmap_colors, colCol = heatmap_colors, col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
                }#,
      #5.2. Don't add compounds
        else
                {
                  cluster_compounds()
                  print("LINCS clustering complete!")
                  output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
                }
        #)
        
        
        
        
#        cluster_compounds()
        
#        output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
        incProgress(7/7, message = "Complete!")
    })
    })
  output$SMILESDownload <- downloadHandler(
    filename = function() {
      paste(input$gene_knockdown, '_compounds.csv', sep='')
    },
    content = function(filename){
      write.csv(final_compounds, filename)
    }
    )
  
  output$SDFDownload <- downloadHandler(
    filename = function() {
      paste(input$gene_knockdown, '.sdf', sep='')
    },
    content = function(filename){
      write.SDF(sdf_smiles, filename)
    }
  )
  
#  output$Representatives <- eventReactive(input$GetRepresentatives,
  observeEvent(input$GetRepresentatives,
        {
          Added_Label <<- input$AddedLabel
          withProgress(message = "Plotting Representatives", value=0, {
            #source("lib/Centroid.R", local=TRUE)
            source("lib/Centroid2.R", local=TRUE)
            #Add a threshold option for cut_tree and Cluster size
          #6. Cut tree at user specified similarity score, identify representative based on minimum distance from all other members
            cut_tree((1-as.numeric(input$CutHeight)))
          #6.1. Include Added compounds
            if (!is.null(input$AddedCompounds)){
              find_centroid_adds(ClusterMembers, (as.numeric(input$ClusterSize)-1))
            }
          #6.2. No compounds were added
            else{
            find_centroid(ClusterMembers, (as.numeric(input$ClusterSize)-1))
            }
            df_centroid <- as.data.frame(unlist(centroid), stringsAsFactors = FALSE)
            colnames(df_centroid)="Representative"
            df_centroid<-df_centroid
            centroid_clusters(centroid, ClusterMembers)
            df_centroid <- cbind(centroid = df_centroid, Cluster = unlist(ClusterCentroid))
            
            
            makeLink <- function(val) {
              #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
              paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
            }
            
            makeLink2 <- function(val) {
              #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
              paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
            }
      

            if (!is.null(input$AddedCompounds)){
              test5 <- data.frame("Compound_ID"=0, "SMILES"=0)
            for (i in 1:length(df_centroid$Representative)){    
                          if (length(which(test2$Compound_ID %in% df_centroid$Representative[[i]]))>0){
                            test3 <- which(test2$Compound_ID %in% df_centroid$Representative[[i]])
                                 for (j in 1:length(test3)){
                                   test4 <- test2[test3[[j]], c("Compound_ID", "SMILES")]
                                   if (test5$Compound_ID[[1]]!=0){
                                   test5 <- rbind(test5, test4)
                                   }
                                   else{
                                     test5 <- test4
                                   }
                                   }
                          }
            }
            }
            else {
              test5 <- data.frame("LSM_ID"=0, "SMILES"=0)
              for (i in 1:length(df_centroid$Representative)){    
                if (length(which(display$LSM_ID %in% df_centroid$Representative[[i]]))>0){
                  test3 <- which(display$LSM_ID %in% df_centroid$Representative[[i]])
                  for (j in 1:length(test3)){
                    test4 <- display[test3[[j]], c("LSM_ID", "SMILES")]
                    if (test5$LSM_ID[[1]]!=0){
                    test5 <- rbind(test5, test4)
                    }
                    else{
                      test5 <- test4
                    }
                  }
                }
              }
            }
            colnames(test5) <- c("Representative", "SMILES")
            test5 <<- test5
              
            
            df_centroid <- merge(test5, df_centroid, stringsAsFactors = FALSE)
            df_centroid <- data.frame(lapply(df_centroid, as.character), stringsAsFactors = FALSE)
            df_centroid$Cluster <- as.integer(df_centroid$Cluster)
            df_centroid <- df_centroid[order(df_centroid$Cluster),]
            df_centroid<<-df_centroid
        
            output$Representatives <<- DT::renderDataTable({
            for (i in 1:length(df_centroid$Representative)){
              if (!is.na((grep("LSM", df_centroid$Representative[[i]])) && (grep("LSM", df_centroid$Representative[[i]])==1))){
                df_centroid$Representative[[i]] <- makeLink(df_centroid$Representative[[i]])
              }
              else if(!is.null(grep("LSM", df_centroid$Representative[[i]]))){
                df_centroid$Representative[[i]] <- gsub("ZINC0", "ZINC", df_centroid$Representative[[i]])              
                df_centroid$Representative[[i]] <- makeLink2(df_centroid$Representative[[i]])
              }
            }
              return(df_centroid)
              df_centroid <<- df_centroid
            }, escape = FALSE, rownames = FALSE)
            #print(df_centroid$Representative)
            
            
            print(paste("Got the centroids as ", centroid[[1]], sep=""))
            incProgress(1/2, detail = "Plotting MDS")
            source("lib/MDS.R")
          #7. Generate MDS plot from representatives distance to one another, pie chart radius corresponds to cluster size
            MDS_plot(centroid)
#            centroid_clusters(centroid, ClusterMembers)
            #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
            #                             + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
            output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
                                         + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", Added_Label)))
            #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid), hjust=1, vjust=2)) 
              #                           + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
            incProgress(2/2, detail = "Complete!")                                                         
            })
        })
            #p <- (ggplot() + geom_text_repel(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid))) 
            #    + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added")))
            
            #move_layers(p, geom_text_repel(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid)), position = "top"))
            
            #output$MDSPlot<<- renderPlot(p)
                                                                                                      
            #                            })
                                         
            
  #####Set it up so you can download Clusters the Representatives are part of
  
  output$RepDownload <- downloadHandler(
    filename = "representatives.csv",
    content = function(filename){
      write.csv(df_centroid, filename)
    }
  
)
  output$ClusterDownload <- downloadHandler(
    filename = "clusters.csv",
    content = function(filename){
      write.csv(ClusterMembers, filename)
    }
  )
#  output$Similar_NCI <- eventReactive(input$GetNCI, {
  observeEvent(input$GetNCI, {
    print("Finding similar NCI compounds")
    source("./lib/SimilaritySearch.R")
    nci_list <- list()
    for(i in 1:length(centroid)){
      similarity_search(centroid[i])
      nci_list[[i]] <- output_nci}
    #nci_df <- as.data.frame(cbind(unlist(nci_list, recursive = FALSE), unlist(ClusterCentroid)))
    nci_df <- mapply(c,ClusterCentroid, nci_list)
    nci_df2 <- list()
    for(i in 1:length(nci_df)){
      if(is.character(nci_df[[i]])){ 
        nci_df2 <- list.append(nci_df2, nci_df[[i]])}
    }
    
    nci_df3 <- data.frame()
    for (j in 1:length(nci_df2)){
      for (k in 2:length(nci_df2[[j]])){
        cluster <- nci_df2[[j]][1]
        nsc <- nci_df2[[j]][k]
        cmpd_row <- cbind(nsc, cluster)
        nci_df3 <- rbind(nci_df3, cmpd_row)
      }
    } 
    #cluster_number <- unlist(ClusterCentroid)
    #nci_unlist <- unlist(nci_list)
    
    makeLink2 <- function(val) {
      #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
      paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
    }
    output$Similar_NCI <- DT::renderDataTable({
      
      #nci_df3$link <- makeLink2(nci_df3$nsc)
      nci_df3$nsc <- makeLink2(nci_df3$nsc)
      return(nci_df3)
      
    }, escape = FALSE)
    
    
    
#    output$Similar_NCI <- renderDataTable(nci_df3)
    #print(nci_list)})
  })

######Not yet working, df2 is not a global variable yet.
    observeEvent(input$GetSTITCH, {
    print("Getting Script")
    #source("./lib/StitchNetFn3.R")
    #source("./lib/StitchNetFn7.R")
    source("./lib/StitchNetFn8.1.R")
    print("sourced")
    #stitchNet("sirolimus", "mtor")
    SelectedCluster <- input$ClusterNumber
    h <- 1
    CompoundLSM <- list()
    for (i in 1:length(ClusterMembers)){
#      if (ClusterMembers[[i]] %in% LINCSCompounds)
      if (ClusterMembers[[i]] == SelectedCluster) {
        CompoundLSM[[h]] <- labels(ClusterMembers)[[1]][[i]]
        h <- h+1
      }
    }
    CompoundLSM <<- CompoundLSM
    print("IDs obtained")
    h <- 1

    CompoundName <- list()
    for (j in 1:length(CompoundLSM)){
      if (length(grep("LSM-", CompoundLSM[[j]]))!=0 && grep("LSM-", CompoundLSM[[j]])==1){
        #for (f in 1:length(df2$lincspertid)){
          if (CompoundLSM[[j]] %in% df2$lincspertid){
            index <- which(df2$lincspertid==CompoundLSM[[j]])
            if (length(index)>1){
              index <- index[[1]]
            }
            #CompoundName[[h]] <- df2$stitchID[[index]]
            CompoundName[[h]] <- df2$compound[[index]]
            h <- h+1
        #  }
        }
      }
      else if (!("LSM" %in% CompoundLSM[[j]])){
        CompoundName[[h]] <- CompoundLSM[[j]]
        h <- h+1
      }
      
    }
    CompoundName <<- CompoundName
    
    chemInput <- unlist(CompoundName)
    chemInput <<- chemInput
    stitchNet(chem = chemInput, gene = input$Gene, limit = input$Connections)
    print("function ran")
    output$STITCHPlot <- renderVisNetwork(visNetwork(molnodes, interScore) %>% 
                                            visNodes(physics = FALSE) %>%
                                            visEdges(physics= FALSE, smooth = FALSE, color = list(color = "gray", highlight = "blue")) %>%
                                            visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
                                                      color = list(background = "ffa7a7", border = "black", 
                                                                   highlight = list(background = "ffa7a7", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Gene KD", shape = "circle",
                                                      color = list(background = "fdff00", border = "black",
                                                                   highlight = list(background = "feff8f", 
                                                                                    border = "red")) ) %>%
                                            visGroups(groupname = "STITCH connectors", shape = "dot", size = 10,
                                                      color = list(background = "ddd", border = "black", 
                                                                   highlight = list(background = "caff44", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                                                      color = list(background = "4e91fc", border = "253085", 
                                                                   highlight = list(background = "4e91fc", 
                                                                                    border = "red"))) %>%
                                            visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = FALSE) %>%
                                            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                                            visInteraction(navigationButtons = TRUE) %>%
                                            visLegend(width = 0.25) %>%
                                            visExport(type="png", name="export-network", float="left")
                                          )
    print("finished")

  })
########
    
  })
