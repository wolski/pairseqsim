##
##Class : AASequenceList
##

setClass("AASequenceList"
         ,representation(info="character",names="character"
                         #    ,list="list" #change
         )
         ,contains="list"
)


setMethod("show","AASequenceList"
          ,function(object)
          {
              cat("info:    ", object@info,"\n")
              cat("length : ", length(as(object,"list")),"\n")
          }
)	


setReplaceMethod("[[", "AASequenceList"
                 , function(x, i, j,..., value)
                 {
                     if( !extends(class(value),"AASequence") )
                     {
                         stop(paste("This is an AASequenceList!"
                                    ,"so dont try to assing a object of class:\",the object is class"
                                    ,class(value)
                                    ,"\n"
                                    ,sep=" ")
                         )
                     }
                     tt<-as(x,"list")
                     tt[[i]]<-value
                     names(tt)[i]<-value@info
                     as(x,"list")<-tt
                     x
                 })



setMethod("[",
          "AASequenceList",
          def = function(x, i, j, ..., drop = F)
          {
              y <- as(x,"list")
              names(y)<-names(x)
              as(x,"list") <- y[i]
              return(x)
          }
)

setReplaceMethod("[","AASequenceList"
                 ,function(x,i,j,...,value)
                 {
                     if( !extends(class(value),"AASequenceList") )
                     {
                         stop(paste("This is an AASequenceList!"
                                    ,"so dont try to assing a object of class:\n"
                                    ,class(value)
                                    ,"\n"
                                    ,"Only Objects of class AASequenceList can be assigned.\n"
                                    ,sep=" ")
                         )
                     }
                     y<-as(x,"list")
                     y[i] <- value
                     as(x,"list") <- y
                     x
                 }
)

##
##Alignment Methods
##
