if (!isGeneric("listdist"))
    setGeneric("listdist",
               function(object,...)
                   standardGeneric("listdist"))

setMethod("listdist"
          ,signature(object="list")
          ,function(object,FUN,diag=FALSE,...)
          {
              lo<-length(object)
              if(length(object)==1)
                  return(dist(1))
              res<-numeric(lo*(lo-1)/2)
              aa <- 1
              for(rr in 1:lo)
              {
                  tt<-(rr+1)
                  if(tt <= lo)
                  {
                      SL <- object[tt:lo] #sublist
                      ee <- (aa-1) + length(SL)
                      tmp <- object[[rr]]
                      res[aa:ee] <- unlist(lapply(SL,FUN,tmp,...))
                  }
                  aa <- ee + 1
              }
              ans <- res
              attributes(ans) <- NULL
              attr(ans,"Labels") <- names(object)
              attr(ans,"Size") <- length(object)
              attr(ans, "call") <- match.call()
              class(ans) <- "dist"
              attr(ans,"Diag") <- diag
              attr(ans,"Upper") <- TRUE
              return(ans)
          }
)

