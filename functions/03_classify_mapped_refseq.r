# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-M", "--mapped"), type="character", default="08_bwamem_onlyflavoreads/83527_ID2400_1-20220112-1-SANO_S115_highqual.txt",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-R", "--reference"), type="character", default="02_reference/Ref_plus_biselli.fasta", 
              help="Config file [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="08_bwamem_onlyflavoreads/83527_ID2400_1-20220112-1-SANO_S115_highqual_name.txt", 
              help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$mapped)) {
  stop("WARNING: No mapped specified with '-I' flag.")
} else {  cat ("mapped is ", opt$mapped, "\n")
  mapped <- opt$mapped  
  }

if (is.null(opt$reference)) {
  stop("WARNING: No reference specified with '-I' flag.")
} else {  cat ("reference is ", opt$reference, "\n")
  reference <- opt$reference  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-I' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

classify<-function() 
{
    library(data.table)
	ref<-fread(cmd=paste("grep '^>'",reference, "| sed -e 's/>//g'"),sep="\t",header=F,data.table=F)
	ref<-unlist(ref)
	refID<-unlist(lapply(strsplit(ref," "),"[",1))
	refgen<-unlist(lapply(strsplit(ref," "),"[",2))
	refspe<-unlist(lapply(strsplit(ref," "),"[",3))
	refbin<-paste(refgen,refspe)
	myref<-data.frame(ID=refID,species=refbin)
	mapres<-fread(mapped,data.table=F)
	names(mapres)<-c("Count","ID")
	finalres<-merge(myref,mapres,sort=F)
	finalres<-finalres[order(finalres$Count,decreasing=T),]
	write.table(finalres,outfile,sep="\t",quote=F,row.names=F)
}
classify()
