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
  
usage <- paste("\n#############\nusage: Rscript ", script.name, " --input <file> --output <file>\n\n",
               "Rscript transforming GTF files to a TxDb object (in SQLite format).\n\n",
               "--input <file>\t\tGTF input file\n",
               "--output <file>\t\tSQLite output file\n",
               "--help\t\t\tPrint this usage.\n\n",
               "Author: Kristin Reiche\n#############\n",
               sep="")

ARG.input <- NULL
ARG.output <- NULL

### Check input parameters
idx <- grep("\\-h$|\\-\\-help$", args, perl=TRUE)
if (length(idx) != 0){
  cat(file=stderr(), usage)
  quit(save="no", status=1)
}

idx <- grep("\\-i$|\\-\\-input$", args, perl=TRUE)   
ARG.input <- .my.checkArg(idx+1, args, usage, required=TRUE)
if(!file.exists(ARG.input)){ stop(paste("Input file \'", ARG.input,"\' does not exist!\n", sep="")) }
if(file.info(ARG.input)$isdir == 1){ stop(paste(ARG.input,"\' is a directory!\n", sep="")) }
if(length(grep("\\.gtf$|\\.gff$", ARG.input, perl=TRUE))==0){ cat(file=stderr(), usage); stop(paste("Unknown file type for input file: \'", ARG.input,"\'.", sep=""))  }

idx <- grep("\\-o$|\\-\\-output$", args, perl=TRUE)      
ARG.output <- .my.checkArg(idx+1, args, usage, required=TRUE)
if(file.exists(ARG.output)){ stop(paste("Output file \'", ARG.output,"\' does already exist!\n", sep="")) }

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]",
      " Script arguments: ",
      "\n\tARG.input=",ARG.input,
      "\n\tARG.output=",ARG.output,
      "\n", sep=""))


############################################
# The script
############################################
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Load R libraries... ", "\n", sep=""))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("GenomicFeatures"))

#########################
# Create transcript database according to GTF file
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Call \'makeTranscriptDbFromGFF\'... ", "\n", sep=""))
txdb <- NULL
if(length(grep("\\.gtf$", ARG.input, perl=TRUE))==1){
  txdb <- makeTranscriptDbFromGFF(file=ARG.input, format="gtf", species="")
}else if(length(grep("\\.gff$", ARG.input, perl=TRUE))==1){
  txdb <- makeTranscriptDbFromGFF(file=ARG.input, format="gff", species="")
}
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Save TxDb object... ", "\n", sep=""))
saveDb(txdb, file=ARG.output)

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " ALL DONE\n", sep=""))
if(length(warnings())>0){ 
        cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Warnings:\n", sep=""))
        warnings() 
}
quit(save="no", status=0)
  

