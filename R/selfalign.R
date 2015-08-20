##
## computes the score of the sequence with itself.
##

if (!isGeneric("selfalign"))
    setGeneric("selfalign",
               function(object,sub,...)
                   standardGeneric("selfalign"))

setMethod("selfalign"
          ,signature(object="AASequence",sub="Submatrix")
          ,function(object,sub)
          {
              sdiag<-sub[diag(dim(sub)[1])==1]
              names(sdiag)<-colnames(sub)
              fre<-frequency(object)
              res<-NULL
              res<-0
              for(x in names(fre))
              {
                  res<-res+ sdiag[x]*fre[x]
              }
              return(res)
          }
)

#promptMethods("selfalign")
