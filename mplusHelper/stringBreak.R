#---
#Date: 10/04/2022
#Author: Thomas P Spargo <thomas.spargo@kcl.ac.uk>
#
#---

#Simple function for adding line breaks to character strings longer than the value defined in the `buffer` argument
#This was originally written to streamline data prep when using the MplusAutomation package,
  #and ensure character strings to be parsed by Mplus as input files do not exceed the 90 character limit
  #semicolons can optionally be appended to character strings using the `add` argument, as Mplus expects these at the end of each command.

#Arguments:
#string - the character string to adjust

#sep - the character to supply to regex for breaking strings
  #e.g. "/" for a long directory filepath; " " to break across empty space; "\t" for a tab delimiter
  #Note that as regex is used, some metacharacters may need to be escaped.
#buffer - the max number of characters to retain before adding a linebreak. Defaults to 80, but should never exceed 90 if string is to be passed to mplus.
#add - Logical, defaults to `FALSE`. if `TRUE`, append semicolon at the end of character string. A semicolon must feature at the end of every command.
#accomHash - Logical, defaults to `FALSE`. Specify whether or not to leave space in the final line of string to add a 32-character md5 hash (as produced by MplusAutomation::mplusModeler if `hashfilename=TRUE`)
#breakchr - define the character used in breaking strings, defaults to "\n", which will be a line-break

stringBreak <- function(string,sep, buffer=80, add=FALSE, accomHash = FALSE, breakchr="\n"){
  
  if(buffer>90){ #Warn when output strings will certainly be too long for Mplus
    warning(paste0("All lines in the Mplus .inp file must contain less than 90 characters. ",
    "Decrease the value of the buffer argument to avoid issues with Mplus parsing strings."))
  } else if (buffer>89 && add==TRUE){ #Warn when output strings will certainly be too long for Mplus
    warning(paste0("All lines in the Mplus .inp file must contain less than 90 characters.",
                   "Appending the semicolon at the end of this strings may exceed this limit when buffer=90",
                   "Decrease the value of the buffer argument to avoid issues with Mplus parsing strings."))
  }
  
  if(!is(string,"character")){ #Return error if supplied an object not of character class
    stop("The object supplied to the string argument does not have character class, please convert to character and try again.")
  }
  
  #Create pattern for string subsetting, broadly: Match to last instance of `sep` within string
  pattern <- paste0(sep,"[^",sep,"]+$")
  
  
  #Check to see if the string exceeds the allowed buffer, and adjust if so
  if(nchar(string)>=buffer){

    #if the string is shorter than 2* buffer, subset only once
    if(nchar(string)<buffer*2){
      
      stringstart<- sub(pattern,"", substr(string,1,buffer))  #Subset string to the first 1:buffer characters and break this substring as specified in pattern
      string <- sub(stringstart,"", string)                   #Then, drop the stringstart portion from string, which may be the end of the string
      
    } else if(nchar(string)>=buffer*2){     #Subset using a while loop if the string >= buffer*2

      stringlist<- vector(mode="list")      #create empty list to store string subsets
      while(nchar(string)>=buffer){         #loop until stringend no longer exceeds buffer length
        subs<- sub(pattern,"", substr(string,1,buffer))              #Subset string to the first 1:buffer characters and break as specified in pattern
        stringlist[[length(stringlist)+1]] <- subs                   #Append subs object as next element of stringlist
        string <- sub(subs,"", string)                               #Drop the current subset portion from string
      } # end loop
      
      stringstart <- paste0(unlist(stringlist),collapse=breakchr) #paste unlisted stringlist and final string portion, collapsing with line breaks
      
    }
    
    #If there is not enough room for addition of a 32-character hash to the string, subset string once more
    if(accomHash==TRUE && buffer-nchar(string)<=33){
      pattern <- paste0(sep,"[^",sep,"]+$")
      
      fitHash<- sub(pattern,"", substr(string,1,buffer))  #Subset string to the first 1:buffer characters and break this substring as specified in pattern
      string <- sub(fitHash,"", string)                   #Then, drop the fitHash portion from string
      
      string <- paste0(fitHash,breakchr,string)
    }
    
    stringout<- paste0(stringstart,breakchr,string) #Paste string with linebreak between start and end
    
  } else {
    # #If there is not enough room for addition of a 32-character hash to the string, subset string
    if(accomHash==TRUE && buffer-nchar(string)<=33){
      fitHash<- sub(pattern,"", substr(string,1,buffer))  #Subset string to the first 1:buffer characters and break this substring as specified in pattern
      string <- sub(fitHash,"", string)                   #Then, drop the fitHash portion from string

      stringout <- paste0(fitHash,breakchr,string)
      message("The string was adjusted to accommodate the 32-character MD5 hash")
    } else {
      
    #If the string characters are already < buffer (taking into account accomHash), do nothing and return message
    message("The string already contains fewer characters than the buffer value defined, no changes have been made")
    stringout <- string
    }
  }
  if(add==TRUE) stringout <- paste0(stringout,";") #Optionally append semi-colon indicating end of mplus line
  return(stringout)
}
