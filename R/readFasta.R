#' readFasta
#' 
#' @export
#' @docType methods
#' @rdname readFasta-methods
if (!isGeneric("readFasta"))
    setGeneric("readFasta",
               function(object,file,...)
                   standardGeneric("readFasta"))

infogrep <- function(x)
{
    return(sub("^>([a-zA-Z0-9]+) .+","\\1",x,perl=TRUE))
}

seqgrep <- function(x)
{
    return(gsub("\\*","",x))
}

#' @rdname readFasta-methods
#' @aliases readFasta,AASequenceList-method
setMethod("readFasta"
          ,signature(object="AASequenceList")
          ,function(object
                    ,file
                    ,grepinfo=infogrep
                    ,grepseq=seqgrep)
          {
              con <-file(file,"r" )
              all <- readLines(con,n=-1)
              pos <-  grep(">",all)
              dat <- vector("list",(length(pos)-1))
              nam <- character(length(pos)-1)
              
              if(length(pos)>1)
              {
                  for(x in 1:(length(pos)-1))
                  {
                      info <- grepinfo(all[pos[x]]) # get the info
                      #cat("x = ",x, "  | pos[x] = ",pos[x], " | ", info , "\n" )
                      seq <- paste(all[(pos[x]+1):(pos[x+1]-1)],collapse="")
                      seq <- grepseq(seq)
                      tmp <- AASequence(info,seq)
                      nam[x] <- info
                      dat[[x]] <- tmp
                  }
                  names(dat)<-nam
              }
              else
              {	
                  
              }
              as(object,"list") <- dat
              return(object)
          }
)

#promptMethods("readFasta")
