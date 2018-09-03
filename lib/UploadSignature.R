sigUpload <- function(filename){
  #Do stuff here
  #signature <- read.table(file = paste("%22", filename, "%22", sep = ""))
  #ftemp = tempfile(pattern="filename", fileext=".xls", tmpdir = tempdir())
  #write.csv(filename, ftemp)
  r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(filename)))
  #r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file("TestSig1.txt")))
  
  rjson <- toJSON(content(r)$status$concordanceTable)
  rjson <<- fromJSON(rjson)
  
  #####Down the road status$data will likely become something else (status$concordanceTable)
  #l <- lapply(content(r)$status$data, function(x)unlist(x))
  #l <- lapply(content(r)$status$concordanceTable, function(x)unlist(x))
  #####
  
  
  #ilincsresult <<- data.frame(t(sapply(l, c)), stringsAsFactors = FALSE)
  #ilincsresult <- data.frame(fromJSON(l))
  #colnames(ilincsresult) <- c("similarity", "pValue", "nGenes", "compound", "lincsPertID", "concentration", "time", "_row", "signatureid", "cellline")
  #ilincsresult <<- ilincsresult
  #print(content(r))
  #r$status_code
}

