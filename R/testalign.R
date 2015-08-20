#promptMethods("align")

##
##Read write functions for Fasta format
##
if (!isGeneric("testalign"))
    setGeneric("testalign",
               function(obj1,obj2,...)
                   standardGeneric("testalign")
    )

setMethod("testalign"
          ,signature(obj1="AASequence",obj2="AASequence")
          ,function(
              obj1
              ,obj2
              ,sub
              ,delta=-4
              ,gapext = delta
              ,alignment="global"
              ,scoring="score"
          )
          {
              res<-.Call("alignScoreSEXP"
                         ,obj1
                         ,obj2
                         ,sub
                         ,delta
                         ,gapext
                         ,alignment
                         ,scoring
              )
              return(res)
          }
)

