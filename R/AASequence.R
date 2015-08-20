##
# AASequence
##

setClass("AASequence"
         ,representation(
             info="character"
         )
         ,contains="character"
)


#promptClass("AASequence")
AASequence <- function(info,sequence)
{
    if(!missing(info) && !missing(sequence))
        new("AASequence",sequence,info)
    else if(!missing(sequence))
        new("AASequence",sequence)
    else
        new("AASequnece")
}

setMethod("initialize"
          ,signature(.Object="AASequence")
          ,function(.Object,sequence,info,alphabet=new("AAAlphabet"))
          {
              if(!missing(sequence))
              {
                  sequence<-toupper(sequence)
                  #first check if the string is constructed out of alphabet letters.
                  seqlevel <- as.character(levels(as.factor(unlist(strsplit(sequence , "" )))))
                  alphlevel <- levels(alphabet)
                  if(sum(is.element(seqlevel,alphlevel)) != length(seqlevel))
                  {
                      print(sequence)
                      stop("This chars are not in the alphabet: ",paste(seqlevel[is.element(seqlevel,alphlevel)==FALSE]
                                                                        ,collapse=" "
                                                                        ,"!!!\n"))
                  }
                  .Object@.Data <- sequence
              }
              if(!missing(info))
                  .Object@info<-info
              .Object
          }
)


setMethod("show",signature(object="AASequence")
          ,function(object)
          {
              cat("info : ", object@info,"\n")
              if(nchar(object)>23)
                  cat("sequence:\n",substr(object,1,10),"...",substr(object,nchar(object)-10,nchar(object)),"\n",sep="")
              else
                  cat("sequence:\n",as(object,"character"),"\n")
          }
)



setMethod("frequency",signature(x="AASequence"),
          function(x,alphabet=new("AAAlphabet"),...)
          {
              res<-table(strsplit(x,""))
              res <- res[alphabet]
              res[is.na(res)]<-0
              res<-as.numeric(res)
              names(res)<-alphabet
              return(res)
          }
)


