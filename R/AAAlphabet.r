#
#class AAAlphabet
#
setClass("AAAlphabet",representation(info="character")
         ,contains="character",
         prototype(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","Z","X","*")
                   ,info="AminoAcid"
         )
)

#promptClass("AAAlphabet")

setMethod("levels",signature(x="AAAlphabet"),
          function(x)
          {
              return(levels(as.factor(unlist(strsplit(x,"")))))
          }
)

setMethod("show",signature(object="AAAlphabet"),
          function(object)
          {
              cat("info : ",object@info,"\n")
              cat("Alphabet :\n",paste(as(object,"character"),collapse=" "),"\n")
          }
)
