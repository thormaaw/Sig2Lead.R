library(XML)
library(RCurl)
library(bitops)
library(scrapeR)
library(visNetwork)

# Creating the stitchNet function
stitchNet <- function(chem = 'NULL', 
                    gene = 'NULL', 
                    confidence = '400', 
                    species = '9606',
                    limit = '10'){

# Get API url to pull the chemical-target interaction data from STITCH (stitch.embl.de)
url_stitch <- 'http://stitch.embl.de/api/psi-mi-tab/interactionsList?'
#chemInput <<- chem
#if(grepl("\\(\\(", chemInput)){
#  print(paste("Filtered out:", " ", chemInput[grep("\\(\\(", chemInput)]))
#}
#chemInput <- chemInput[!grepl("\\(\\(", chemInput)] 
#chemInput <- chemInput[!grepl(",", chemInput)] 

# Clean chemInput to solve an api compatibility issue
chemInput_clean <- chemInput
for (i in 1:length(chemInput_clean)){
  if(grepl("\\(\\(", chemInput_clean[i])){ 
    print(paste("Filtered out to prevent errors:", " ",chemInput_clean[grep("\\(\\(", chemInput_clean[i])]))
    chemInput_clean <- chemInput_clean[!grepl("\\(\\(", chemInput_clean)] 
    } else if(grepl(",", chemInput_clean[i])){
      print(paste("Filtered out to prevent errors:", " ", chemInput_clean[i]))
      chemInput_clean <-  chemInput_clean[!grepl(",", chemInput_clean)]
      }
}

chems <- paste(chemInput_clean, collapse = '%0D')
chems <- gsub(" ", "%20", chems)
genes <<- gene 
score <- '400' #score 400 = moderate confident (set as a default value), score 700 = strong confident, score >900 = very strong confident, score range = 0-999
speciesID <- '9606' #By default, species ID = 9606 = H.sapians; Mus Musculus txid = 10090, Zebra fish (Danio rerio) txid= 7955
limitInt <- '10' #By default, maximum numbers of connected nodes= 10; a range of 0 - 25 is suggestive

api <- paste0(url_stitch, 'identifiers=', chems,'%0D', genes,
              '&required_score=', score, 
              '&species=', speciesID, 
              '&limit=', limitInt
              )
print(api)
# Get interScore (edge) for visNetwork
doc <<- readLines(api)
doc <- data.frame(id = doc)
newCol1 <- strsplit(as.character(doc$id), '\t',fixed=TRUE)
doc <- data.frame(doc, do.call(rbind, newCol1))
newCol2 <-strsplit(as.character(doc$X15),'|',fixed=TRUE)
doc <- data.frame(doc, do.call(rbind, newCol2))
newCol3 <- strsplit(as.character(doc$X1.1),':',fixed=TRUE)
doc <- data.frame(doc,do.call(rbind, newCol3))
#interScore <- data.frame(doc[, c(4, 5, 22)])
interScore <- data.frame(doc$X3, doc$X4, doc$X2.2)
names(interScore) <- c('from', 'to', 'score')
width <- (as.numeric(as.character(interScore$score))*2)^2
interScore <- data.frame(interScore, width)
names(interScore) <- c('from', 'to', 'score', 'width')
interScore <- data.frame(interScore, title = paste0("<p>", "Score:", "<br>", interScore$score, "</p>"))
interScore <<- interScore

# Get molnodes (node) for visNetwork
mol <- c(as.character(interScore$from), as.character(interScore$to))
mol <- unique(mol)
mol <<- mol
group <- replace(mol, toupper(mol) %in% toupper(chemInput), "Input chemicals present in STITCH")
group <- replace(group, toupper(group) %in% toupper(genes), "Gene KD")
group <- replace(group, !group %in% c("Input chemicals present in STITCH", "Gene KD"), "STITCH connectors")
molnodes <- data.frame(label = mol, id = mol, group = group, 
                       title = paste0('<a target= "_blank"', '', 'href=', '', 
                                      'https://en.wikipedia.org/wiki/', mol, '>', mol, '</a>')
                      )

# create HTML links for molnodes
link_modnodes <- character()
for (i in seq_along(mol)){
  if(grepl(toupper(genes), mol[i]) == 1){
    link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                               'https://www.genecards.org/cgi-bin/carddisp.pl?gene=', 
                               mol[i], '>', mol[i], '</a>')
  } else if(grepl("CHEMBL", mol[i]) == 1){
    link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                      'https://www.ebi.ac.uk/chembl/beta/compound_report_card/', 
                      mol[i], '>', mol[i], '</a>')
  } else if(grepl("ZINC", mol[i]) == 1){
    link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                      'http://zinc15.docking.org/substances/', 
                      mol[i], '>', mol[i], '</a>')
  } else if(grepl("LSM", mol[i]) == 1){
    link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                               'http://lincsportal.ccs.miami.edu/SmallMolecules/view/', 
                               mol[i], '>', mol[i], '</a>')
  } else if(grepl("-", mol[i]) == 1){
    link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                               ' https://pubchem.ncbi.nlm.nih.gov/search/#query=', 
                               mol[i], '>', mol[i], '</a>')
  } else {
    link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                      #'https://pubchem.ncbi.nlm.nih.gov/compound/', 
                      'https://en.wikipedia.org/wiki/', 
                      mol[i], '>', mol[i], '</a>')
  }
}

molnodes$title <- link_modnodes


## get information of LSM compounds (undetected by stitch) by hyperlink to iLINCS
inputChem <- toupper(chemInput) %in% toupper(mol)
inputChem <<- inputChem
stitch_undetect <- chemInput[!inputChem]

## When stitch_undetect produces any missing value, replacing missing value with NULL
stitch_undetect <- if(as.numeric(length(stitch_undetect))==0){"NULL"} else{stitch_undetect} 
undetect_Node <- data.frame(label = stitch_undetect, id = stitch_undetect, group = "Undetected by STITCH", 
                           # title = paste0('<a target= "_blank"', '', 'href=', '', 
                             #      'http://lincsportal.ccs.miami.edu/SmallMolecules/view/', 
                              #      stitch_undetect, '>', stitch_undetect, '</a>')
                            title = paste0('<a target= "_blank"', '', 'href=', '', 
                                           'https://en.wikipedia.org/wiki/', 
                                           stitch_undetect, '>', stitch_undetect, '</a>')
                             )

# create HTML links for undetect_Node
link_undetect_Node <- character()
for (i in seq_along(stitch_undetect)){
  if(grepl("CHEMBL", stitch_undetect[i]) == 1){
    link_undetect_Node[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                   'https://www.ebi.ac.uk/chembl/beta/compound_report_card/', 
                   stitch_undetect[i], '>', stitch_undetect[i], '</a>')
  } else if(grepl("ZINC", stitch_undetect[i]) == 1){
    link_undetect_Node[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                   'http://zinc15.docking.org/substances/', 
                   stitch_undetect[i], '>', stitch_undetect[i], '</a>')
  } else if(grepl("LSM", stitch_undetect[i]) == 1){
    link_undetect_Node[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                    'http://lincsportal.ccs.miami.edu/SmallMolecules/view/', 
                    stitch_undetect[i], '>', stitch_undetect[i], '</a>')
  } else {
    link_undetect_Node[i] <- paste0('<a target= "_blank"', '', 'href=', '', 
                   'https://pubchem.ncbi.nlm.nih.gov/search/#query=', 
                   stitch_undetect[i], '>', stitch_undetect[i], '</a>')
  }
}

undetect_Node$title <- link_undetect_Node

## Adding undetected compounds into molnodes
molnodes <- rbind(molnodes, undetect_Node)
molnodes <<- molnodes

# Generate a network!
visNetwork(molnodes, interScore) %>%
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
  #visExport only works under Shiny environment, so enable a code below when running in Shiny environment
  visExport(type="png", name="export-network", float="left") 

}

