setClass("AAAlignment"
         ,representation(
             info1="character"
             ,info2="character"
             ,selfs1="numeric"
             ,selfs2="numeric"
             ,score="numeric"
             ,identity="numeric"
             ,alignsimilarity="numeric"
             ,lch1="numeric"
             ,lch2="numeric"
             ,alig1="character"
             ,alig2="character"
             ,beautify="character"
         )
)

#promptClass("AAAlignment")
##
##v \item{similarity}{
##v This is a count of the number of positions over the length of the alignment where >= 51% of the residues or bases at that position are similar. 
##v Any two residues or bases are defined as similar when they have positive comparisons (as defined by the comparison matrix being used in the alignment algorithm). 
##v  }
##v \item{identity}{
##v This is a count of the number of positions over the length of the alignment where all of the residues or bases at that position are identical. 
##v }

setMethod("show"
          ,signature(object="AAAlignment")
          ,function(object)
          {
              lalign<-nchar(object@alig1)
              cat("selfscore 1: ",object@selfs1,"; seq length 1 :",object@lch1,"\n")
              cat("selfscore 2: ",object@selfs2,"; seq length 2 :",object@lch2,"\n")
              cat("alig lenght: ",lalign,"\n")
              cat("score      : ",object@score,"\n")
              cat("FM         : ",object@score/(sqrt(object@selfs1*object@selfs2)),"\n")
              cat("identity   : ",object@identity,"/",min(object@lch1,object@lch2),"\n")
              cat("similarity : ",object@alignsimilarity,"/",min(object@lch1,object@lch2),"\n")
          }
)

setMethod("summary"
          ,signature(object="AAAlignment")
          ,function(object)
          {
              identity <- sum(unlist(strsplit(object@alig1,""))==unlist(strsplit(object@alig2,"")))
              nammax <- max(nchar(object@info1),nchar(object@info2))
              lalign<-nchar(object@alig1)
              cat("selfscore 1: ",object@selfs1,"\n")
              cat("selfscore 2: ",object@selfs2,"\n")
              cat("alig lenght: ",lalign,"\n")
              cat("score      : ",object@score,"\n")
              cat("FM(score)  : ",object@score/(sqrt(object@selfs1*object@selfs2)),"\n")
              cat("identity   : ",object@identity,"/",min(object@lch1,object@lch2),"\n")
              cat("similarity : ",object@alignsimilarity,"/",min(object@lch1,object@lch2),"\n")
              tmp<-c(seq(1,lalign,60),lalign)
              for(x in 1:(length(tmp)-1))
              {
                  s1 <- substr(object@alig1 , tmp[x] , tmp[x+1] )
                  s2 <- substr(object@alig2 , tmp[x] , tmp[x+1] )
                  beauti <- substr(object@beautify ,tmp[x],tmp[x+1])
                  cat( format( object@info1 , width=nammax ) , s1 , "\n" )
                  cat( format( " ", width =nammax), beauti, "\n")
                  cat( format( object@info2 , width=nammax ) , s2 , "\n" )
                  cat("\n")
              }
          }
)
