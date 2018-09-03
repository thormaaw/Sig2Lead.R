########################################
# Find Con/Dis/cordant Signtures
########################################

ConDisconSignatures <- function(signatureID)
  {
  con_comma <<- list("character")
  print(input$ConOrDiscon)
  for (i in 1:length(signatureID))
  {
#url2 <- paste("http://www.ilincs.org/api/SignatureMeta/findConcordantSignatures?sigID=%22",signatureID,"%22&lib=%22LIB_5%22", sep="")
    url2 <- paste("http://www.ilincs.org/api/SignatureMeta/findConcordantSignatures?sigID=%22",signatureID[[i]],"%22&lib=%22LIB_5%22", sep="")
df <- fromJSON(url2)
if (length(df)>0)
df2<<- df[order(df$similarity), ]


con <- df2[df2$similarity > 0.321,]
discon <- df2[df2$similarity < -0.235, ]

if (input$ConOrDiscon=="Concordant"){
  concordant_list <- list("character")
      for (j in 1:length(con))
        {
       concordant_list[j] <- con[j]
      }

  con_comma[[i]] <- concordant_list
}

else if (input$ConOrDiscon=="Discordant"){
  discordant_list <- list("character")
  for (j in 1:length(discon))
  {
    discordant_list[j] <- discon[j]
  }
  
  con_comma[[i]] <- discordant_list
}

#con_comma <<- paste(con$perturbagenID, collapse="%22,%22")
#print(con_comma)
  }
  con_comma <<- con_comma
}

UploadConDiscon <- function(uploadresult){
  con_comma <- list("character")
  
  df2 <<- uploadresult
  con <<- df2$lincsPertID[df2$similarity > 0.321]
  discon <<- df2$lincsPertID[df2$similarity < -0.235]

  TrueCon <- list()
  if (input$ConOrDiscon=="Concordant"){    
  NoNullIndex <- list()
  v=1
  for (o in 1:length(con)){
    if (!is.null(con[[o]])){
      NoNullIndex[[v]] <- o
      v <- v+1
      }
  }
  
  for (p in 1:length(NoNullIndex)){
  TrueCon[[p]] <- con[[NoNullIndex[[p]]]]
  }
  }
  else if (input$ConOrDiscon=="Discordant"){
    NoNullIndex <- list()
    v=1
    for (o in 1:length(discon)){
      if (!is.null(discon[[o]])){
        NoNullIndex[[v]] <- o
        v <- v+1
      }
    }
    
    for (p in 1:length(NoNullIndex)){
      TrueCon[[p]] <- discon[[NoNullIndex[[p]]]]
    }
  }
  
  con_comma <- TrueCon
  con_comma <- as.list(con_comma)
  con_comma <<- as.list(con_comma)
  }
