##Copyright 2001, W. Wolski, all rights reserved.
##
# Class Submatrix
##
library("methods") 

#' Class "Submatrix" Substitution Matrix
#' Amino Acid Substitution Matrices
#' @rdname Submatrix
setClass("Submatrix",
         representation(copyright="character",info="character",head="character",alphabet="character")
         ,contains="matrix",prototype(copyright="GNU GENERAL PUBLIC LICENSE Version 2, June 1991")
)

#'
#'@export
#'@docType methods
#'@rdname subFromEmboss-methods
if (!isGeneric("subFromEmboss"))
    setGeneric("subFromEmboss",
               function(object,path,...)
                   standardGeneric("subFromEmboss"))

#' create substitution matrix from emboss file
#' @rdname subFromEmboss-methods
#' @aliases subFromEmboss,Submatrix,character,ANY-method
setMethod("subFromEmboss",signature(object="Submatrix",path="character"),
          function(object,path,...)
          {
              con<-file(path,"r")
              res <- readLines(con=con,n=-1)
              close(con)
              path<-unlist(strsplit(path,"/"))
              object@info<-path[length(path)]
              head <- res[grep("#",res)]
              head<-paste(head,collapse="\n")
              head<-paste(head,"\n",sep="")
              object@head <- head
              res <- res[-grep("#",res)]
              alphabet <- unlist(strsplit(res[1]," +"))
              object@alphabet <- alphabet[2:length(alphabet)]
              res2<-NULL
              for(x in 2:length(res))
              {
                  if(nchar(res[x])!=0)
                  {
                      tmp<-unlist(strsplit(res[x]," +"))
                      res2<-rbind(res2,as.numeric(tmp[2:length(tmp)]))
                  }
              }
              object@.Data <- res2
              colnames(object)<-object@alphabet
              rownames(object)<-object@alphabet
              object
          })


#' show submatix
#' @export
setMethod("show","Submatrix",function(object)
{
    cat("info : ",object@info,"\n")
    cat("copyright : ",object@copyright,"\n")
    cat("head :\n" ,object@head)
    cat("alphabet : ",paste(object@alphabet,collapse=" "),"\n")
    print(as(object,"matrix"))
})

#' construct Submatrix
#' @export
Submatrix<-function()
{
    new("Submatrix")
}
