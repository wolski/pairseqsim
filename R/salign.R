if (!isGeneric("salign"))
    setGeneric("salign",
               function(obj1,obj2,...)
                   standardGeneric("salign")
    )

setMethod("salign"
          ,signature(obj1="AASequenceList",obj2="NULL")
          ,function(obj1
                    ,obj2
                    ,sub
                    ,delta= -4
                    ,gapext = delta
                    ,alignment = "global"
                    ,scoring = "identity"
                    ,diag = FALSE
          )
          {
              res<-listdist(obj1
                            ,testalign
                            ,diag=diag
                            ,sub
                            ,delta=delta
                            ,gapext=gapext
                            ,alignment=alignment
                            ,scoring=scoring
              )
              if(scoring=="identity"||scoring=="similarity"||scoring=="scoreN")
              {
                  res[1:length(res)] <- (1-as.numeric(res))
              }
              else if(scoring=="score")
              {
                  tmp <- as.numeric(res)
                  tmp <- (tmp - mean(tmp))/sqrt(var(tmp)) 
                  res[1:length(res)] <- pnorm(tmp,mean=0,sd=1,lower.tail=FALSE)
              }
              else if(scoring=="pozitive")
              {
                  res[1:length(res)] <- pnorm(res,mean=0,sd=1,lower.tail=FALSE)
              }
              return(res)
          }
)



setMethod("salign"
          ,signature(obj1="AASequence",obj2="AASequenceList")
          ,function(obj1,obj2,sub,delta=-4,gapext = delta, alignment="global",scoring="pozitive")
          {
              res<-salign(obj2,obj1,sub, delta=delta , gapext=gapext , alignment=alignment,scoring = scoring)
          }
)

setMethod("salign"
          ,signature(obj1="AASequenceList",obj2="AASequence")
          ,function(obj1,obj2,sub,delta=-4,gapext = delta, alignment="global",scoring="pozitive")
          {
              res<-lapply(obj1,testalign,obj2,sub,delta=delta,gapext=gapext,alignment=alignment,scoring=scoring)
              res <- unlist(res)
          }
)


##t Pairwise sequence Aligment of Amino Acid Sequence
##- Pairwise sequence Aligment of Amino Acid Sequence. Function can compute
##- global, local or overlap alignment of two amino acid sequences.
##+ ret : what to return? AAAlignment = object of class AAAlignement, identity/alignment length, similarity, score.
##+ delta : gap opening penalty
##+ gepext : gap extension penalty
##+ type :  type of alignemnt. (e.g. global,local,overlap)
##+ sub : similarity matrix BLOSUM; PAM similarity matrix.
##+ obj1 : object of class AASequence
##+ obj2 : object of class AASequence
##e 
##e 

setMethod("salign",signature(
    obj1="AASequence"
    ,obj2="AASequence"
)
,function(
    obj1
    ,obj2
    ,sub
    ,delta=-4
    ,gapext = delta
    ,alignment="global"
    ,scoring="AAAlignment"
)
{
    mret<-c("AAAlignment","identity","similarity","score","scoreN")
    if(!(scoring %in% mret))
    {
        stop("scoring argument can be either ", paste(mret,collapse=" "),"\n")
    }
    
    ####real computing.
    res<-.Call("alignSEXP"
               ,obj1
               ,obj2
               ,sub
               ,delta
               ,gapext
               ,alignment
    )
    
    if(scoring=="AAAlignment")
    {
        if(length(grep("   ",res[["errmsg"]])) == 0)
        {
            stop("ERROR:", errmsg, "\n")
        }
        
        ss1<-paste(res[["alig1"]],"*",sep="")
        ss2<-paste(res[["alig2"]],"*",sep="")
        lalign<-nchar(ss1)
        tmp<-c(seq(1,lalign,40),lalign)
        wmatch<-rep(" ",lalign)
        beauti<-""
        for(x in 1:(length(tmp)-1))
        {
            vs1 <- substr(ss1 , tmp[x] , tmp[x+1]-1 )
            vs2 <- substr(ss2 , tmp[x] , tmp[x+1]-1 )
            vs1<-unlist(strsplit(vs1,""))
            vs2<-unlist(strsplit(vs2,""))
            match <- wmatch[tmp[x]:(tmp[x+1]-1)]
            ##find similarities.
            ## what exactly means the star in the blosum matrix.
            ## you anyway will count only positive values
            ## * AA are always smaller than 0 so it does what we are looking for.
            vs1t<-vs1
            vs2t<-vs2
            vs1t[vs1=="-"] <- "*"
            vs2t[vs2=="-"] <- "*"
            tmpsim <- sub[vs1t,vs2t] # get the values of the diagonal
            
            tmpsim <- tmpsim[diag(dim(tmpsim)[1])==1] # get the diagonal
            names(tmpsim)<-NULL
            match[tmpsim>0]<-":" #mark similarities
            match[vs1==vs2]<-"|" #mark identities
            beauti<-paste(beauti,paste(match,collapse=""),sep="")
        }
        res<- new("AAAlignment"
                  ,info1 = obj1@info
                  ,info2 = obj2@info
                  ,selfs1 = res[["selfscore1"]]
                  ,selfs2 = res[["selfscore2"]]
                  ,score=res[["score"]]
                  ,identity = res[["identity"]]
                  ,alignsimilarity = res[["alignsimilarity"]]
                  ,alig1 = res[["alig1"]]                      
                  ,alig2 = res[["alig2"]]
                  ,lch1= nchar(obj1)
                  ,lch2= nchar(obj2)
                  ,beautify = beauti
        )
        return(res)
    }
    else if(scoring=="similarity")
    {
        return(res[["alignsimilarity"]]/min(nchar(obj1),nchar(obj2)))
    }
    else if(scoring=="identity")
    {
        return(res[["identity"]]/min(nchar(obj1),nchar(obj2)))
    }
    else if(scoring=="scoreN")
    {
        sc <- ifelse(res[["score"]]>0,res[["score"]],0)
        return(sc/min(res[["selfscore1"]],res[["selfscore2"]]))
    }
    else if(scoring=="score")
    {
        return(res[["score"]])
    }
}
)

