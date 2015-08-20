#setwd("/home/witek/devel/pairseqsim2/inst/doc/")
#source("/home/witek/devel/pairseqsim2/R/seqsim.R")
#dyn.load("/home/witek/devel/pairseqsim2/src/pairseqsim2.so")
#load("/home/witek/devel/pairseqsim2/data/EPAM110.rda")
#load("/home/witek/devel/pairseqsim2/data/EBLOSUM62.rda")

library(pairseqsim)
data(EPAM110)
data(EBLOSUM62)
mySequlist <-new("AASequenceList",info="my sequence list")
mySequlist<-readFasta(mySequlist,"ex.fasta",grepinfo=infogrep,grepseq=seqgrep)

seq1<-new("AASequence","MEDQVGFGFRPNDEEL",info="seq1")
seq1<-new("AASequence","VAISEVNICSYDPWNL",info="seq1")
seq1<-new("AASequence","VAISEVNICSY",info="seq1")
seq1<-new("AASequence","MEDQVGFGFRPNDEELVGHYLRNKIEGNTSRDVEVAISEVNICS",info="seq1")
seq2<-new("AASequence","MAASEHRCVGCGFRVKSLFIQYSPGNIRLMKCGNCKEVADEYIECERMIIFIDLILHRPKVYRHVLYNAINPATVNIQHLLWKLVFAYLLLDCYRSLLLRKSDEESSFSDSPVLLSIKVRSFLFNGLN",info="seq2")
seq2<-new("AASequence","MAASEHRCVGCGFRV",info="seq2")
res<-testalign(seq1,seq2,EBLOSUM62,delta=-10,gapext=-1,alignment="global")
res<-salign(mySequlist[[1]],mySequlist[[2]],EPAM110,delta=-4,gapext=-1,alignment="global",scoring="score")
res<-testalign(mySequlist[[1]],mySequlist[[2]],EPAM110,delta=-4,gapext=-1,alignment="global",scoring="pozitive")

#files<-dir()
#files<-files[2:length(files)]
#for(x in files)
#  {
#    assign(x,subFromEmboss(new("Submatrix"),x))
#    save(list=x,file=paste(x,".rda",sep=""))
#  }

