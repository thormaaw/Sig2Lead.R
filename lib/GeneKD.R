####################################
# User Input Variable (Gene KD)
####################################


define_knockdown <- function(gene)
{
  usr <- gene
  
  url1 <- paste("http://www.ilincs.org/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%2C%20%22libraryid%22%3A%22LIB_6%22%7D%7D", sep="")
  
  raw.result <- GET(url = url1)
  raw.result$status_code
  sigid <<- list("character")
  for (i in 1:length(content(raw.result)))
  {sigid[i] <- content(raw.result)[[i]]$signatureid
  }
  sigid <<- sigid
  #sigid <<- content(raw.result)[[1]]$signatureid
#  print(sigid)
}
