###############################
# PURPOSE: Check if a command line argument is provided
.my.checkArg <- function(idx, args, usage, required=TRUE){
###############################
  if (length(idx) == 0 && required){
    cat(file=stderr(), usage)
    stop(paste("Required arguments are missing.", sep=""))
  }
  else {
    if (length(idx) == 0){
      return(NULL)  # Argument is not required, so return NULL
    }
    if (length(idx) > 1){
      stop(paste("Argument \'", args[idx],"\' is defined more than once.", sep=""))
    }
    return(args[idx[1]])
  }
}

###############################
# PURPOSE: Read a GTF or GFF file and convert it to a transcriptDb object
.my.gtf2txdb <- function(x=NULL){
###############################

  suppressPackageStartupMessages(require("GenomicFeatures"))

  if(is.null(x)){ stop(paste("Please provide input file for argument \'x\'\n", sep="")) }
  if(!file.exists(x)){ stop(paste("File \'", x,"\' does not exist!\n", sep="")) }
  
  txdb <- NULL
  if(length(grep("\\.gtf$", x, perl=TRUE))==1){
    txdb <- makeTranscriptDbFromGFF(file=x, format="gtf", species="")
  }else if(length(grep("\\.gff$", x, perl=TRUE))==1){
    txdb <- makeTranscriptDbFromGFF(file=x, format="gff", species="")
  }
  return(txdb)
}

###############################
# PURPOSE: Check the metadata column of a GTF file (column 9) and add a unique prefix to transcript and gene identifiers of novel transcripts.
# RETURNS: The input data frame with replaced transcript and gene identifiers in column 9.
.my.checkAndReplaceMetadata <- function(x=NULL, prefix=NULL, list.pattern=NULL){
###############################

  # check arguments
  if(is.null(x)){ stop(paste("Please provide a value for argument \'x\'\n", sep="")) }
  if(is.null(list.pattern)){ stop(paste("Please provide a list of patterns for argument \'list.pattern\'\n", sep="")) }
  if(ncol(x) != 9){ stop(paste("Number of columns exceeds 9.\n", sep=""))  }

  # check format of column 9
  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Check format in column 9 -- metadata... ", "\n", sep=""))
  if(length(grep(list.pattern[["pattern.all"]], x$V9, perl=TRUE)) != nrow(x)){
    stop(paste("Not all lines match regular expression: \'", list.pattern[["pattern.all"]],"\'.\n", sep=""))
  }
  if(length(grep(paste(list.pattern[["pattern.txID.novel"]],"|",list.pattern[["pattern.txID.known"]], sep=""), x$V9, perl=TRUE)) != nrow(x)){
    stop(paste("Not all lines match either \'", list.pattern[["pattern.txID.known"]],"\' or \'", list.pattern[["pattern.txID.novel"]],"\'.\n", sep=""))
  }
  if(length(grep(paste(list.pattern[["pattern.geneID.novel"]],"|",list.pattern[["pattern.geneID.known"]], sep=""), x$V9, perl=TRUE)) != nrow(x)){
    stop(paste("Not all lines match either \'", list.pattern[["pattern.geneID.known"]],"\' or \'", list.pattern[["pattern.geneID.novel"]],"\'.\n", sep=""))
  }
  if(length(grep(prefix, x$V9)) > 0){
    stop(paste("Please provide a different prefix, hence the current one \'", prefix,"\' matches a substring.\n", sep=""))
  }
  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " OK... ", "\n", sep=""))

  # create additional column with original txIDs -- required to extract exons later on
  x$txID <- sub(list.pattern[["pattern.all"]],"\\1", x$V9, perl=TRUE)

  # replace transcript and gene identifiers
  if(!is.null(prefix)){
    cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Add prefix to novel transcript and gene IDs... ", "\n", sep=""))
    x$V9 <- sub(list.pattern[["pattern.txID.novel"]], paste("transcript_id ", prefix, "_\\1; ",sep=""), x$V9, perl=TRUE)
    x$V9 <- sub(list.pattern[["pattern.geneID.novel"]], paste("gene_id ", prefix, "_\\1; ",sep=""), x$V9, perl=TRUE)
  }

  # in case ";" are not followed by a space, add space again
  x$V9 <- gsub(";(\\S+)", "; \\1", x$V9, perl=TRUE)
  
  return(x)
}

# Global variables
MODE_NO_OV_TRANSCRIPT="no_transcript_overlap"

############################################
# Command line arguments processing
############################################
# load command line arguments, if any are provided
options(echo=FALSE) # set to TRUE if you want to see all commands in standard output
args <- commandArgs(trailingOnly = FALSE)   # returns also the script name as '--file=doCorrelation.R'
script.name <- args[0]
if (length(args)==0){ stop("No command line arguments? ") }

### Retrieve script name from command line arguments
idx <- grep("\\-\\-file=", args, perl=TRUE)
script.name <- gsub("\\-\\-file=","",args[idx])
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Retrieve command line arguments... ", "\n", sep=""))
  
usage <- paste("\n#############\nusage: Rscript ", script.name, " --inputA <file> --inputB <file> --outputGTF <file> --outputSQLite <file> --prefixA <string> --prefixB <string> [--mode <", MODE_NO_OV_TRANSCRIPT,">]\n\n",
               "Rscript merging two GTF files according to a selected mode. Currently only mode \'", MODE_NO_OV_TRANSCRIPT,"\' is implemented,\n",
               "such that assembly in file A is extended by all novel transcripts in file B which do NOT overlap any novel transcript in file A.\n",
               "The overlap of entries in B with entries in A is calculated based on transcript annotation only.\n",
               "Final assembly is written to a GTF file. A user defined prefix is added to novel gene and transcript identifiers. Currently only cufflinks\n",
	       "identifiers are supported to detect novel transcripts in A and B (\'XLOC\' and \'TCONS\').\n\n",
               "--mode <no_transcript_overlap>\tMode of invocation, currently all those novel transcripts in file B overlapping transcripts in file A are discarded.\n",
               "--inputA <file>\t\t\tGTF input file A\n",
               "--inputB <file>\t\t\tGTF input file B\n",
               "--output <file>\t\t\tGTF output file\n",
               "--prefixA <string>\t\tPrefix added to gene and transcript identifier for all novel transcripts in file A when written to output file.\n",
               "--prefixB <string>\t\tPrefix added to gene and transcript identifier for all novel transcripts in file B when written to output file.\n",
               "--help\t\t\t\tPrint this usage.\n\n",
               "Author: Kristin Reiche\n#############\n",
               sep="")

ARG.mode <- MODE_NO_OV_TRANSCRIPT
ARG.inputA <- NULL
ARG.inputB <- NULL
ARG.output <- NULL
ARG.prefixA <- NULL
ARG.prefixB <- NULL

### Check input parameters
idx <- grep("\\-h$|\\-\\-help$", args, perl=TRUE)
if (length(idx) != 0){
  cat(file=stderr(), usage)
  quit(save="no", status=1)
}

idx <- grep("\\-m$|\\-\\-mode$", args, perl=TRUE)
ARG.mode <- .my.checkArg(idx+1, args, usage, required=TRUE)
if(ARG.mode != MODE_NO_OV_TRANSCRIPT){ cat(file=stderr(), usage); stop(paste("Unknown mode \'", ARG.mode,"\' for merging.", sep=""))  }

idx <- grep("\\-a$|\\-\\-inputA$", args, perl=TRUE)   
ARG.inputA <- .my.checkArg(idx+1, args, usage, required=TRUE)
if(!file.exists(ARG.inputA)){ stop(paste("Input file \'", ARG.inputA,"\' does not exist!\n", sep="")) }
if(file.info(ARG.inputA)$isdir == 1){ stop(paste(ARG.inputA,"\' is a directory!\n", sep="")) }
if(length(grep("\\.gtf$|\\.gff$", ARG.inputA, perl=TRUE))==0){ cat(file=stderr(), usage); stop(paste("Unknown file type for input file: \'", ARG.inputA,"\'.", sep=""))  }

idx <- grep("\\-b$|\\-\\-inputB$", args, perl=TRUE)   
ARG.inputB <- .my.checkArg(idx+1, args, usage, required=TRUE)
if(!file.exists(ARG.inputB)){ stop(paste("Input file \'", ARG.inputB,"\' does not exist!\n", sep="")) }
if(file.info(ARG.inputB)$isdir == 1){ stop(paste(ARG.inputB,"\' is a directory!\n", sep="")) }
if(length(grep("\\.gtf$|\\.gff$", ARG.inputB, perl=TRUE))==0){ cat(file=stderr(), usage); stop(paste("Unknown file type for input file: \'", ARG.inputB,"\'.", sep=""))  }

idx <- grep("\\-g$|\\-\\-output$", args, perl=TRUE)      
ARG.output <- .my.checkArg(idx+1, args, usage, required=TRUE)
if(file.exists(ARG.output)){ stop(paste("Output file \'", ARG.output,"\' does already exist!\n", sep="")) }

idx <- grep("\\-p$|\\-\\-prefixA$", args, perl=TRUE)   
ARG.prefixA <- .my.checkArg(idx+1, args, usage, required=TRUE)

idx <- grep("\\-q$|\\-\\-prefixB$", args, perl=TRUE)   
ARG.prefixB <- .my.checkArg(idx+1, args, usage, required=TRUE)

if(ARG.prefixA == ARG.prefixB){ stop(paste("Prefix for file A is identical to prefix for file B, but they must be distinct!\n", sep=""))  }

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]",
      " Script arguments: ",
      "\n\tARG.mode=",ARG.mode,
      "\n\tARG.inputA=",ARG.inputA,
      "\n\tARG.inputB=",ARG.inputB,
      "\n\tARG.prefixA=",ARG.prefixA,
      "\n\tARG.prefixB=",ARG.prefixB,
      "\n\tARG.output=",ARG.output,
      "\n", sep=""))

############################################
# The script
############################################
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Load R libraries... ", "\n", sep=""))
suppressPackageStartupMessages(library("GenomicFeatures"))

#########################
# Create transcript database according to GTF file
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Call \'makeTranscriptDbFromGFF\' for input file A... ", "\n", sep=""))
txdb.A <- .my.gtf2txdb(x=ARG.inputA)
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Call \'makeTranscriptDbFromGFF\' for input file B... ", "\n", sep=""))
txdb.B <- .my.gtf2txdb(x=ARG.inputB)

#########################
# read input GTF files as data frame --> easier export to final file
list.pattern <- list(pattern.all = ".+transcript_id\\s+(\\S+)\\;.+",
                     pattern.txID.novel = "transcript_id\\s+(TCONS\\S+)\\;",
                     pattern.txID.known = "transcript_id\\s+(ENS\\S+)\\;",
                     pattern.geneID.novel = "gene_id\\s+(XLOC\\S+)\\;",
                     pattern.geneID.known = "gene_id\\s+(ENS\\S+)\\;"
                     )
list.pattern.novelIDs <- list(txID="TCONS\\S+")

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Read GTF file A as a data frame... ", "\n", sep=""))
gtf.A <- read.delim(ARG.inputA, header=FALSE, stringsAsFactors=FALSE, comment.char = "#")
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Read GTF file B as a data frame... ", "\n", sep=""))
gtf.B <- read.delim(ARG.inputB, header=FALSE, stringsAsFactors=FALSE, comment.char = "#")

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Check metadata column of file A... ", "\n", sep=""))
gtf.A.new <- .my.checkAndReplaceMetadata(x=gtf.A, prefix=ARG.prefixA, list.pattern=list.pattern)
if(nrow(gtf.A) != nrow(gtf.A.new)){ stop(paste("Something went wrong while checking metadata column for file A !\n", sep=""))   }
gtf.A <- gtf.A.new
rm(gtf.A.new)

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Check metadata column of file B... ", "\n", sep=""))
gtf.B.new <- .my.checkAndReplaceMetadata(x=gtf.B, prefix=ARG.prefixB, list.pattern=list.pattern)
if(nrow(gtf.B) != nrow(gtf.B.new)){ stop(paste("Something went wrong while checking metadata column for file B!\n", sep=""))   }
gtf.B <- gtf.B.new
rm(gtf.B.new)

#########################
# Calculate overlap of transcripts 
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Calculate overlap according to mode \'", ARG.mode,"\'... ", "\n", sep=""))
gtf.final <- NULL
if(ARG.mode == MODE_NO_OV_TRANSCRIPT){

  # extract novel transcripts in A and B
  transcripts.A <- transcripts(txdb.A, columns=c("tx_id", "tx_name"))
  txID.novel.A <- grep(list.pattern.novelIDs[["txID"]], unique(transcripts.A$tx_name), value=TRUE)
  transcripts.A.novel <- transcripts(txdb.A, columns=c("tx_id", "tx_name"), vals=list(tx_name=txID.novel.A))
  
  transcripts.B <- transcripts(txdb.B, columns=c("tx_id", "tx_name"))
  txID.novel.B <- grep(list.pattern.novelIDs[["txID"]], unique(transcripts.B$tx_name), value=TRUE)
  transcripts.B.novel <- transcripts(txdb.B, columns=c("tx_id", "tx_name"), vals=list(tx_name=txID.novel.B))
  
  # which novel transcripts in B overlap with novel transcripts in A by at least one nucleotide?
  ov <- countOverlaps(transcripts.B.novel, transcripts.A.novel)

  # test overlap of transcripts
  if(length(ov) != length(transcripts.B.novel)){ stop(paste("Length of overlap vector does not equal number of novel transcripts in B!\n", sep=""))  }
  transcripts.B.noOv <- transcripts.B.novel[which(ov==0),]
  ov.test <- countOverlaps(transcripts.B.noOv, transcripts.A.novel, type=c("any"))
  if(length(which(ov.test>0))>0){ stop(paste("Something went wrong while identifying non-overlapping novel transcripts in B (test of transcripts)!\n", sep=""))  }

  # create final set of transcripts
  txID.B.noOv <- unique(transcripts.B.noOv$tx_name)
  gtf.B.noOv <- gtf.B[gtf.B$txID %in% txID.B.noOv,]

  if(length(which(!txID.B.noOv %in% gtf.B.noOv$txID)) > 0){
    stop(paste("Some selected transcripts cannot be found in original data frame!\n", sep="")) 
  }
  if(!identical(colnames(gtf.A), colnames(gtf.B.noOv))){
    stop(paste("Column names differ in final data frames!\n", sep="")) 
  }

  # Create final data frame and delete helper column $txID == column 10
  gtf.final <- rbind(gtf.A, gtf.B.noOv)
  gtf.final <- gtf.final[,1:9]

}


#########################
# Save GTF file
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Save transcript database as GTF file... ", "\n", sep=""))
if(is.null(gtf.final) || length(gtf.final)==0){
  stop(paste("Final data frame is empty!\n", sep="")) 
}
write.table(gtf.final, file=ARG.output, append=FALSE, quote=FALSE, sep="\t", dec = ".", row.names=FALSE, col.names=FALSE)


#########################
# Clean up and quit
if(length(warnings())>0){ 
        cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Warnings:\n", sep=""))
        warnings() 
}
sessionInfo()
quit(save="no", status=0)
  

