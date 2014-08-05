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

.my.getNewTranscriptID <- function(txdb){

  if(!is(txdb,"TranscriptDb"))
    stop("'txdb' must be a TranscriptDb object")
  SQL <- paste("SELECT * ",
               "FROM (SELECT transcript._tx_id AS tx_id, (tx_name || \"::\" || gene.gene_id) AS tx_name, tx_chrom, tx_strand, tx_start, tx_end ",             
					  "FROM gene JOIN transcript ",
                				    "ON gene._tx_id=transcript._tx_id ",
					  ") AS A ",
                "JOIN (SELECT _tx_id AS tx_id, exon_rank, exon_start, exon_end ",
                	  "FROM splicing  JOIN exon ON splicing._exon_id=exon._exon_id) AS B ",
                "ON A.tx_id = B.tx_id"
            	)
  data <- GenomicFeatures:::queryAnnotationDb(txdb, SQL);
  trans <- unique(data.frame(data)[c(1,2,3,4,5,6)]);
  splice <- data.frame(data)[c(1,8,9,10)];
  newTxdb <- makeTranscriptDb(trans,splice)
  
 return(newTxdb)
}

LIST.SPECIES <- list(hg19 = list(name="Homo sapiens", BSgenome="BSgenome.Hsapiens.UCSC.hg19"),
                     mm10 = list(name="Mus musculus", BSgenome="BSgenome.Mmusculus.UCSC.mm10")
                     )

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

## Help section
if("--help" %in% args) {
  cat("
      gtf2fasta.R
 
      Arguments:
      --[gtf/sql] <FILE>                Specify file type and name 
				        * gtf : GTF format with UNIQUE transcripts/exons
				        * sql : SQLite TranscriptDB (Bioconductor GenomicFeatures) 

      --species <mm10|hg19>             An abbreviation for the species, this is required to assign the correct genome and metadata.

      [--five_prime_flank <integer>]    Optional parameter to add 5' flanking region to transcript

      [--three_prime_flank <integer>]   Optional parameter to add 3' flanking region to transcript

      [--bed <FILE>]                    Optional file with genomic region, which are 
                                        use to filer transcripts, e.g. intergenic regions

      [--fasta <FILENAME>]              Name of output file, which contains DNA sequences

      [--help]                          print this text
 
      Example:
      Rscript gtf2Fasta.R --gtf input.gtf

      Authors: Karolin Wiedemann, Kristin Reiche\n\n")
  
  q(save="no")
}

ARG.gtf=NULL
ARG.sql=NULL
ARG.species=NULL
ARG.fivePrimeFlank=NULL
ARG.threePrimeFlank=NULL
ARG.fasta=NULL
ARG.bed=NULL

mode<-0;

idx <- grep("\\-\\-gtf$", args, perl=TRUE)      
ARG.gtf <- .my.checkArg(idx+1, args, usage, required=FALSE)
if(!is.null(ARG.gtf)){ 
  message(paste("GTF:", ARG.gtf, sep="\t")); 
  if(file.exists(ARG.gtf)) {mode <- mode+1};
}

idx <- grep("\\-\\-sql$", args, perl=TRUE)      
ARG.sql <- .my.checkArg(idx+1, args, usage, required=FALSE)
if(!is.null(ARG.sql)){ 
  message(paste("SQL:", ARG.sql, sep="\t"));
  if(file.exists(ARG.sql)){mode <- mode+2}; 
}

print(mode);

if(mode==0){
  message("GTF/SQL file do not exist!");
  q(save="no");
}
if(mode==3){
  message("2 input files (GTF and SQL) given. Please choose only one!");
  q(save="no");
}

idx <- grep("\\-\\-fasta$", args, perl=TRUE)      
ARG.fasta <- .my.checkArg(idx+1, args, usage, required=FALSE)
if(!is.null(ARG.fasta)){message(paste("FASTA:", ARG.fasta, sep="\t"));}

idx <- grep("\\-\\-bed$", args, perl=TRUE)      
ARG.bed <- .my.checkArg(idx+1, args, usage, required=FALSE)
if(!is.null(ARG.bed)){message(paste("BED:", ARG.bed, sep="\t"));}

idx <- grep("\\-\\-species$", args, perl=TRUE)      
ARG.species <- .my.checkArg(idx+1, args, usage, required=TRUE)
if(!ARG.species %in% names(LIST.SPECIES)){ stop(paste("Species ", ARG.species," is currently not supported!\n", sep="")) }

idx <- grep("\\-\\-five_prime_flank$", args, perl=TRUE)      
ARG.fivePrimeFlank <- .my.checkArg(idx+1, args, usage, required=FALSE)
if(!is.null(ARG.fivePrimeFlank) && length(grep("^\\d+$", ARG.fivePrimeFlank, perl=TRUE)) == 0 ){ stop(paste("Length of 5' flanking region must be a number!\n", sep="")) }
if(!is.null(ARG.fivePrimeFlank)){ ARG.fivePrimeFlank <- as.integer(ARG.fivePrimeFlank) }

idx <- grep("\\-\\-three_prime_flank$", args, perl=TRUE)      
ARG.threePrimeFlank <- .my.checkArg(idx+1, args, usage, required=FALSE)
if(!is.null(ARG.threePrimeFlank) && length(grep("^\\d+$", ARG.fivePrimeFlank, perl=TRUE)) == 0 ){ stop(paste("Length of 3' flanking region must be a number!\n", sep="")) }
if(!is.null(ARG.threePrimeFlank)){ ARG.threePrimeFlank <- as.integer(ARG.threePrimeFlank) }

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]",
      " Script arguments: ",
      "\n\tARG.gtf=", ARG.gtf,
      "\n\tARG.sql=", ARG.sql,
      "\n\tARG.species=", ARG.species,
      "\n\tARG.fivePrimeFlank=", ARG.fivePrimeFlank,
      "\n\tARG.threePrimeFlank=", ARG.threePrimeFlank,
      "\n\tARG.fasta=", ARG.fasta,
      "\n\tARG.bed=", ARG.bed,
      "\n", sep=""))


############################################
# The script
############################################
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Load R libraries... ", "\n", sep=""))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("plyr"))

# Read GTF file if specified
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Read input file... ", "\n", sep=""))
if(mode==1){
 message("LOAD GTF FILE");
 myTab <- read.table(file=ARG.gtf, sep="\t");
 myTabLines <- apply(myTab, 1, paste, collapse="_");
 
 # test, if transcripts unique in gtf file
 if(length(unique(myTabLines))!=length(myTabLines)){
   message("Transcripts in GTF file have to be unique!");
   q(save="no");
 }
}

txdb <- NULL
if(mode==2){  # read SQL file if specified
  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Load sqlite file... ", "\n", sep=""))
  txdb <- loadDb(ARG.sql)
} else{  
  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Create and save sqlite transcriptDb... ", "\n", sep=""))
  txdb <- makeTranscriptDbFromGFF(file=ARG.gtf, format="gtf", species=LIST.SPECIES[[ARG.species]]$name)
  sqliteDBFileName <- paste(ARG.gtf,"sqlite",sep=".")
  saveDb(txdb, sqliteDBFileName)
  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Sslite transcript database saved in \'", sqliteDBFileName,"\'... ", "\n", sep=""))
}

cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Change TRANSCRIPT_NAME into ENSTxxxxx.x::ENSGxxxxx.x... ", "\n", sep=""))
txdb <- .my.getNewTranscriptID(txdb)

if(!is.null(ARG.bed)){
  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Extract transcripts from predefined regions... ", "\n", sep=""))
  if(!file.exists(ARG.bed)){"BED file does not exist!"; q(save="no");}

  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Read region file... ", "\n", sep=""))
  intergenic <- import(ARG.bed, format="BED") 

  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Select transcripts... ", "\n", sep=""))
  tx <- transcripts(txdb)
  itx <- subsetByOverlaps(tx,intergenic,maxgap=0L,minoverlap=1L,type="within", ignore.strand=TRUE)
  names <- itx$tx_name
  iex <- select(txdb,keys=names, columns=columns(txdb), keytype="TXNAME")

  trans <- unique(data.frame(iex)[c(15,1,17,18,19,20)])
  colnames(trans) <- c("tx_id", "tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end")
  splice <- data.frame(iex)[c(15,16,12,13)]
  colnames(splice) <- c("tx_id", "exon_rank", "exon_start", "exon_end")

  txdb <- makeTranscriptDb(trans,splice);
}

### Flanking regions upstream and/or downstream?
if(!is.null(ARG.fivePrimeFlank) || !is.null(ARG.threePrimeFlank)){
  cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Add flanking regions... ", "\n", sep=""))
  if(is.null(ARG.fivePrimeFlank)){ ARG.fivePrimeFlank=0 }
  if(is.null(ARG.threePrimeFlank)){ ARG.threePrimeFlank=0 }

  # get exons - ordered by the splicing.exon_rank field in internal database
  # use.names=FALSE ==> use internal names, i.e. tx_id, required below for makeTranscriptDb
  exons.grl <- exonsBy(txdb, by=c("tx"), use.names=FALSE)
  exons.grl <- lapply(exons.grl, function(x){
    # order exons according to start position
    x <- x[order(start(x))]
    startPos <- start(x)
    endPos <- end(x)
    n <- length(startPos)
    if(unique(strand(x))=="+"){
      startPos[1] <- startPos[1] - ARG.fivePrimeFlank
      endPos[n] <- endPos[n] + ARG.threePrimeFlank
    }
    else if(unique(strand(x))=="-"){
      startPos[1] <- startPos[1] - ARG.threePrimeFlank
      endPos[n] <- endPos[n] + ARG.fivePrimeFlank
    }
    else{
      stop("Unknown strand.")
    }
    start(x) <- startPos
    end(x) <- endPos
    x
  })

  exons.df <- ldply(names(exons.grl), function(aName){
    x <- exons.grl[[aName]]
    data.frame("tx_id" = as.integer(rep(aName, length(x))),
               "exon_rank" = elementMetadata(x)$exon_rank,
               "exon_start" = start(x),
               "exon_end" = end(x),
               stringsAsFactors=FALSE
               )
  })

  transcripts <- transcripts(txdb)
  idx.posStrand <- which(strand(transcripts)=="+")
  idx.negStrand <- which(strand(transcripts)=="-")
  start(transcripts[idx.posStrand]) <- start(transcripts[idx.posStrand]) - ARG.fivePrimeFlank
  end(transcripts[idx.posStrand]) <- end(transcripts[idx.posStrand]) + ARG.threePrimeFlank
  start(transcripts[idx.negStrand]) <- start(transcripts[idx.negStrand]) - ARG.threePrimeFlank
  end(transcripts[idx.negStrand]) <- end(transcripts[idx.negStrand]) + ARG.fivePrimeFlank

  transcripts.df <- data.frame("tx_id" = as.integer(elementMetadata(transcripts)$tx_id),
                               "tx_name" = elementMetadata(transcripts)$tx_name,
                               "tx_chrom" = as.vector(seqnames(transcripts)),
                               "tx_strand" = as.vector(strand(transcripts)),
                               "tx_start" = start(transcripts),
                               "tx_end" = end(transcripts),
                               stringsAsFactors=FALSE
                               )
  
  txdb <- makeTranscriptDb(transcripts=transcripts.df, splicings=exons.df)

}

### Get FASTA
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Extract sequences... ", "\n", sep=""))
require(LIST.SPECIES[[ARG.species]]$BSgenome, character.only=TRUE)
genome <- eval(parse(text = LIST.SPECIES[[ARG.species]]$BSgenome)); #string into variable name
tx_seqs <- extractTranscriptsFromGenome(genome, txdb, use.names=TRUE);

### Export
cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Export fasta... ", "\n", sep=""))
if(is.null(ARG.fasta)){
  myGtfName <- basename(ARG.gtf);
  myName <- unlist(strsplit(myGtfName,"\\.gtf"));
  ARG.fasta <- paste(paste(myName,"_transcripts",sep=""), "fa", sep=".");
}
export(tx_seqs, ARG.fasta, format="fasta");


cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " ALL DONE\n", sep=""))
if(length(warnings())>0){ 
        cat(file=stderr(), paste("[", script.name,"][",Sys.time(),"]", " Warnings:\n", sep=""))
        warnings() 
}
sessionInfo()
quit(save="no", status=0)
