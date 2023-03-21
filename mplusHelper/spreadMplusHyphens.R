### Define a function function for string splitting mplus variables into column names.
#Note that this does not reinsert inherently in order....

spreadMplusHyphens <- function(namesare){
  
  #Split the full string by whitespace
  spnm <- strsplit(namesare,split = " ")[[1]]
  
  #If any elements are completely blank (created by spaces in the document), drop them
  if(any(spnm=="")) spnm<- spnm[-which(spnm=="")]
  
  #Identify whether any of the elements are collapsed using a hyphen
  hashyphen<- grep("-",spnm)
  if(length(hashyphen)>0){
    #Identify and split patterns with hyphenation
    m<- gregexpr("[A-Z]+[0-9]+(-)[A-Z]+[0-9]+",namesare)
    hyph<- regmatches(namesare,m)
    spl_hyph<- strsplit(hyph[[1]],split = "-")
    
    #Extract the range of values and convert to a sequence
    ranges <- lapply(spl_hyph,gsub,pattern="[A-Z]+", replacement="")
    
    ranges <- lapply(ranges, dorange <- function(x){
      num<- as.numeric(x)
      num[1]:num[2]
    })
    
  #Extract the strings
  prefix <- lapply(spl_hyph, gsub,pattern="[0-9]+", replacement="")
  
  #Recombine the sequence of numbers with the strings
  rejoin<- mapply(paste0, prefix, ranges, SIMPLIFY = FALSE)
  
  #Make a list to store all 'spread' column names. Start by adding the first element of rejoin
  all <- list(c(spnm[1:(hashyphen[1]-1)], rejoin[[1]])) #Add in the first row of rejoined elements
  
  #If there were more elements to the rejoin list. then loop across each remaining element
  if(length(rejoin)>1){
    for(i in 2:length(hashyphen)){all[[i]] <- c(spnm[(hashyphen[i-1]+1):(hashyphen[i]-1)],rejoin[[i]])}
    
    #If values continue past the current spread hyphen string, append more
    if(!is.na(spnm[(hashyphen[i]+1)])){all[[i+1]] <- spnm[(hashyphen[i]+1):length(spnm)]}
  }
  
  #Then, if any values continue past the final spread hyphen string, append the remaining names
  if(!is.na(spnm[(hashyphen[length(hashyphen)]+1)])){all[[length(all)+1]] <- spnm[(hashyphen[length(hashyphen)]+1):length(spnm)]}

  #Unlist the 'all' list
  fullrejoin <- unlist(all)

  #Check to see if any hyphens have slipped into full-rejoin (occurs when elements of hashyphen are adjacent)
  #If this happens, remove the indexes indicated by slip
  slip <- grep("-",fullrejoin)
  if(length(slip)>0){fullrejoin <- fullrejoin[-slip]}
    
  } else {
    fullrejoin <- spnm
  }
  return(fullrejoin)
}
