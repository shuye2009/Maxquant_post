#library(proteus)
library(ggpubr)
library(magrittr)
library(dplyr)
library(data.table)
library(gplots)
library(pheatmap)
library(PerformanceAnalytics)
library(missForest)

## report spectral count from evidence table
EvidenceToContrast <- function(evidence) {
  #evidence <- evi_table
  #Remove potential contaminants
  columns <- colnames(evidence)
  Con.value <- columns %like% "contaminant"
  NonContaminants <- evidence[, Con.value] == ""
  evidence <- evidence[NonContaminants,]

  #Remove Reverse Peptides
  NotReverse <- evidence[, "Reverse"] == ""|is.na(evidence[, "Reverse"])
  evidence <- evidence[NotReverse,]

  #Simplify evidence file
  ColumnList <- c("Leading.razor.protein", "Experiment", "MS.MS.count")
  SimplifiedEvidence <- evidence[, ColumnList]
  colnames(SimplifiedEvidence)[1] <- "Protein"
  SimplifiedEvidence$Protein <- as.character(SimplifiedEvidence$Protein)
  SimplifiedEvidence$Experiment <- as.character(SimplifiedEvidence$Experiment)

  #Aggregate proteins in evidence file
  EvidenceList <- aggregate(. ~Protein + Experiment, SimplifiedEvidence, FUN=sum)

  #Get List of Experiments and Preys
  ProteinNames <- c(unique(EvidenceList[, "Protein"]))
  ExperimentNames <- c(unique(EvidenceList[,"Experiment"]))
  NamesList <- list(ProteinNames, ExperimentNames)

  #Make matrix with experiments and preys
  Contrast <- matrix(nrow = length(ProteinNames), ncol = length(ExperimentNames), dimnames = NamesList)

  for (i in (1:nrow(EvidenceList))) {
    Experiment <- EvidenceList[i,"Experiment"]
    Protein <- EvidenceList[i,"Protein"]
    Spectral <- EvidenceList[i,"MS.MS.count"]
    Contrast[Protein, Experiment] <- Spectral
  }

  Contrast[is.na(Contrast)] <- 0

  #return contrast of spectral counts
  return(Contrast)
}

## report peptide count from evidence table
EvidenceToPeptideCount <- function(evidence) {
  #evidence <- evi_table
  #Remove potential contaminants
  columns <- colnames(evidence)
  Con.value <- columns %like% "contaminant"
  NonContaminants <- evidence[, Con.value] == ""
  evidence <- evidence[NonContaminants,]

  #Remove Reverse Peptides
  NotReverse <- evidence[, "Reverse"] == ""|is.na(evidence[, "Reverse"])
  evidence <- evidence[NotReverse,]

  #Simplify evidence file
  ColumnList <- c("Leading.razor.protein", "Experiment", "Peptide.ID")
  SimplifiedEvidence <- evidence[, ColumnList]
  colnames(SimplifiedEvidence)[1] <- "Protein"
  SimplifiedEvidence$Protein <- as.character(SimplifiedEvidence$Protein)
  SimplifiedEvidence$Experiment <- as.character(SimplifiedEvidence$Experiment)

  #Aggregate proteins in evidence file
  EvidenceList <- aggregate(. ~Protein + Experiment, SimplifiedEvidence, FUN=length)

  #Get List of Experiments and Preys
  ProteinNames <- c(unique(EvidenceList[, "Protein"]))
  ExperimentNames <- c(unique(EvidenceList[,"Experiment"]))
  NamesList <- list(ProteinNames, ExperimentNames)

  #Make matrix with experiments and preys
  PeptideCount <- matrix(nrow = length(ProteinNames), ncol = length(ExperimentNames), dimnames = NamesList)

  for (i in (1:nrow(EvidenceList))) {
    Experiment <- EvidenceList[i,"Experiment"]
    Protein <- EvidenceList[i,"Protein"]
    TotalCount <- EvidenceList[i,"Peptide.ID"]
    PeptideCount[Protein, Experiment] <- TotalCount
  }

  PeptideCount[is.na(PeptideCount)] <- 0

  #return PeptideCount of spectral counts
  return(PeptideCount)
}

## report either total, Razor + unique or Unique peptides from proteinGroup, default Razor+unique peptide count
## use 'other' as type for single sample
getPeptideCount <- function(proteinGroup, type="Razor + unique") {
  #proteinGroup <- pg_table
  #type <- "Razor + unique"
  #Remove potential contaminants
  columns <- colnames(proteinGroup)
  
  Con.value <- columns[columns %like% "contaminant"]
  Rev.value <- columns[columns %like% "Reverse"]
  
  proteinGroup[is.na(proteinGroup[, Rev.value]), Rev.value] <- ""
  proteinGroup[is.na(proteinGroup[, Con.value]), Con.value] <- ""
  ## use .data[[x]] to convert variable x into a column name
  proteinGroup %>% 
     filter(.data[[Con.value]] != "+") %>%
     filter(.data[[Rev.value]] != "+") -> proteinGroup
  
  
  #Simplify proteinGroup file
  print(type)

  if(type == "Unique"){
    prefix <- "Unique.peptides."
    expColumn <- columns[columns %like% "Unique.peptides."]

  }else if(type == "Razor + unique"){
    prefix <- "Razor...unique.peptides."
    expColumn <- columns[columns %like% "Razor...unique.peptides."]
  }else if(type =="Peptide"){
    prefix <- "Peptides."
    expColumn <- columns[columns %like% "Peptides."]
  }else{
    prefix <- " " ## for single sample
    expColumn <- columns[columns %like% "Peptides"]
  }

  proteinColumn <- "Majority.protein.IDs" # Maxquant recommmends this over "Protein.IDs"
  proteinIDs <- as.character(proteinGroup[,proteinColumn])
  protein <- unlist(lapply(proteinIDs, function(x){unlist(strsplit(x, ";", fixed=T))[1]}))

  PeptideCount <- as.matrix(proteinGroup[,expColumn])
  cns <- sub(prefix, "", expColumn, fixed=T)
  row.names(PeptideCount) <- protein
  colnames(PeptideCount) <- cns


  PeptideCount <- PeptideCount[!grepl("REV__", row.names(PeptideCount)), ]
  PeptideCount[is.na(PeptideCount)] <- 0

  #return PeptideCount of spectral counts
  PeptideCount <- as.matrix(PeptideCount)
  #if(dim(PeptideCount)[2] == 1){colnames(PeptideCount) <- "PeptideCount"}
  return(PeptideCount)
}


## report spectral counts from proteinGroup
getSpectralCount <- function(proteinGroup) {
  #proteinGroup <- pg_table
  #Remove potential contaminants and reverse matches
  columns <- colnames(proteinGroup)
  Con.value <- columns[columns %like% "contaminant"]
  Rev.value <- columns[columns %like% "Reverse"]
  proteinGroup[is.na(proteinGroup[, Rev.value]), Rev.value] <- ""
  proteinGroup[is.na(proteinGroup[, Con.value]), Con.value] <- ""
  
  proteinGroup %>% 
     filter(.data[[Con.value]] != "+") %>%
     filter(.data[[Rev.value]] != "+") -> proteinGroup
  
  
  #Simplify proteinGroup file
  prefix <- "MS.MS.count."
  expColumn <- columns[columns %like% "MS.MS.count."]
  if(length(expColumn) == 0){
    prefix <- " "
    expColumn <- columns[columns %like% "MS.MS.count"] ## for single sample
  }
  
  proteinColumn <- "Majority.protein.IDs" # Maxquant recommmends this over "Protein.IDs"
  proteinIDs <- as.character(proteinGroup[,proteinColumn])
  protein <- unlist(lapply(proteinIDs, function(x){unlist(strsplit(x, ";", fixed=T))[1]}))

  SpectralCount <- as.matrix(proteinGroup[,expColumn])
  cns <- sub(prefix, "", expColumn, fixed=T)
  row.names(SpectralCount) <- protein
  colnames(SpectralCount) <- cns

  SpectralCount <- SpectralCount[!grepl("REV__", row.names(SpectralCount)), ]
  SpectralCount[is.na(SpectralCount)] <- 0

  #return PeptideCount of spectral counts
  SpectralCount <- as.matrix(SpectralCount)
  #if(dim(SpectralCount)[2] == 1){colnames(SpectralCount) <- "SpectralCount"}
  return(SpectralCount)
}

## report LFQ intensity from proteinGroup
getLFQintensity <- function(proteinGroup, LFQ=TRUE) {
  #proteinGroup <- pg_table
  #LFQ <- FALSE
  #Remove potential contaminants
  columns <- colnames(proteinGroup)
  Con.value <- columns[columns %like% "contaminant"]
  Rev.value <- columns[columns %like% "Reverse"]
  proteinGroup[is.na(proteinGroup[, Rev.value]), Rev.value] <- ""
  proteinGroup[is.na(proteinGroup[, Con.value]), Con.value] <- ""
  
  proteinGroup %>% 
     filter(.data[[Con.value]] != "+") %>%
     filter(.data[[Rev.value]] != "+") -> proteinGroup
  

  #print(dim(proteinGroup))
  #Simplify proteinGroup file
  prefix <- "LFQ.intensity."
  expColumn <- columns[columns %like% "LFQ.intensity."]
  if(!LFQ){
    prefix <- "Intensity."
    expColumn <- columns[columns %like% "Intensity."]
  }
  if(length(expColumn) == 0){
    prefix <- " "
    expColumn <- columns[columns %like% "Intensity"] ## for single sample
  }

  #print(expColumn)
  
  proteinColumn <- "Majority.protein.IDs" # Maxquant recommmends this over "Protein.IDs"
  proteinIDs <- as.character(proteinGroup[,proteinColumn])
 
  protein <- unlist(lapply(proteinIDs, function(x){unlist(strsplit(x, ";", fixed=T))[1]}))

  LFQintensity <- as.matrix(proteinGroup[,expColumn])
  cns <- sub(prefix, "", expColumn, fixed=T)
  row.names(LFQintensity) <- protein
  colnames(LFQintensity) <- cns
  
  LFQintensity <- LFQintensity[!grepl("REV__", row.names(LFQintensity)), ]
  LFQintensity[is.na(LFQintensity)] <- 0

  #return PeptideCount of spectral counts
  LFQintensity <- as.matrix(LFQintensity)

  return(LFQintensity)
}

getMolecularWeight <- function(weightList, seq){
  H <- 1.007276 # hydrogen mass
  mw <- 0
  aa <- unlist(strsplit(seq, split=""))
  for(a in aa){
    mw <- mw + weightList[a] - 2*H ## every aa lose 2 hydrogen
  }
  mw < mw + 2*H  ## fist and last aa only lose one hydrogen
  return(mw)
}
## extract information from uniprot fasta header, to get a dataframe of Accession, geneName, description
build_id_map <- function(wList, fasta, organism){
  
  uniprot_fasta <- readLines(fasta)

  ## get sequence length
  seq_len <- NULL
  seq_mw <- NULL
  headers <- NULL
  sequences <- NULL
  header <- NULL
  aseq <- NULL
  for(seq in uniprot_fasta){
    if(grepl(">", seq)){
      if(is.null(header)){
        header <- seq
      }else{
        headers <- c(headers, header)
        sequences <- c(sequences, aseq)
        seq_len <- c(seq_len, nchar(aseq))
        seq_mw <- c(seq_mw, getMolecularWeight(wList, aseq))
        header <- seq
        aseq <- NULL
      }

    }else{
      aseq <- paste(aseq, seq, sep="")
    }
  }
  headers <- c(headers, header)
  sequences <- c(sequences, aseq)
  seq_len <- c(seq_len, nchar(aseq))
  seq_mw <- c(seq_mw, getMolecularWeight(wList, aseq))

  seq_info <- data.frame(headers, sequences, seq_len, seq_mw)

  head(seq_info)


  uniprot_genes <- headers
  if(is.na(organism) || organism == "NA" || is.null(organism) || organism == ""){
    uniprot_acc <- uniprot_name <- gsub(">", "", uniprot_genes)
    uniprot_description <- rep(NA, length(uniprot_genes))
    uniprot_organism <- rep(NA, length(uniprot_genes))
  }else{
    uniprot_acc <- unlist(lapply(uniprot_genes, function(x){unlist(strsplit(x, split="|", fixed=T))[2]}))
    pattern <- "GN=.+PE="
    uniprot_name <- unlist(lapply(uniprot_genes, function(x){y <- regmatches(x, regexpr(pattern, x)); if(length(y) == 0) y <- "-"; y}))
    uniprot_name <- unlist(lapply(uniprot_name, function(x){gsub("GN=|PE=", "", x)}))
    uniprot_name <- unlist(lapply(uniprot_name, function(x){trimws(x)}))
    pattern <- paste("_", organism, ".+OS=", sep="")
    uniprot_description <- unlist(lapply(uniprot_genes, function(x){y <- regmatches(x, regexpr(pattern, x)); if(length(y) == 0) y <- "-"; y}))
    uniprot_description <- unlist(lapply(uniprot_description, function(x){gsub(paste("_", organism, "|OS=", sep=""), "", x)}))
    uniprot_description <- unlist(lapply(uniprot_description, function(x){trimws(x)}))
    uniprot_organism <- rep(organism, length(uniprot_acc))
  }

  geneNameDesc <- data.frame(uniprot_acc, uniprot_name, seq_len, uniprot_description, uniprot_organism, seq_mw)

  noname <- geneNameDesc$uniprot_name == "-"
  geneNameDesc[noname, "uniprot_name"] <- geneNameDesc[noname, "uniprot_acc"]

  geneNameDesc
}

## extract information from refseq fasta header, to get a dataframe of Accession, geneName, description
build_id_map_sars2 <- function(wList, fasta, organism){
  refseq_fasta <- readLines(fasta)

  ## get sequence length
  seq_mw <- NULL
  seq_len <- NULL
  headers <- NULL
  sequences <- NULL
  header <- NULL
  aseq <- NULL
  for(seq in refseq_fasta){
    if(grepl(">", seq)){
      if(is.null(header)){
        header <- seq
      }else{
        headers <- c(headers, header)
        sequences <- c(sequences, aseq)
        seq_len <- c(seq_len, nchar(aseq))
        seq_mw <- c(seq_mw, getMolecularWeight(wList, aseq))
        header <- seq
        aseq <- NULL
      }

    }else{
      aseq <- paste(aseq, seq, sep="")
    }
  }
  headers <- c(headers, header)
  sequences <- c(sequences, aseq)
  seq_len <- c(seq_len, nchar(aseq))
  seq_mw <- c(seq_mw, getMolecularWeight(wList, aseq))

  seq_info <- data.frame(headers, sequences, seq_len, seq_mw)

  head(seq_info)

  refseq_genes <- headers
  refseq_acc <- unlist(lapply(refseq_genes, function(x){unlist(strsplit(x, split="|", fixed=T))[1]}))
  refseq_acc <- unlist(lapply(refseq_acc, function(x){gsub(">", "", x)}))
  refseq_name <- unlist(lapply(refseq_genes, function(x){unlist(strsplit(x, split="|", fixed=T))[2]}))

  pattern <- paste("_", organism, ".+\\[organism=", sep="")

  refseq_description <- unlist(lapply(refseq_genes, function(x){y <- regmatches(x, regexpr(pattern, x)); if(length(y) == 0) y <- "-"; y}))
  refseq_description <- unlist(lapply(refseq_description, function(x){gsub(paste("_", organism, "|\\[organism=", sep=""), "", x)}))
  refseq_description <- unlist(lapply(refseq_description, function(x){trimws(x)}))
  refseq_organism <- rep("SARS2", length(refseq_acc))

  #geneNameDesc <- data.frame(refseq_acc, refseq_name, seq_len, refseq_description, refseq_organism)
  # proteinGroups.txt use protein name in the place of protein accession for SARS2
  geneNameDesc <- data.frame(refseq_name, refseq_name, seq_len, refseq_description, refseq_organism, seq_mw)
  geneNameDesc
}


#annotate, start by getting uniprot_acc ID
AnnotateWithSymbol <- function(Contrast, idMap){
  #Contrast <- Spectral_Count
  #idMap <- geneNameDesc

  idMap$uniprot_acc <- toupper(idMap$uniprot_acc)
  rowName <- as.character(row.names(Contrast))
  colName <- as.character(colnames(Contrast))

  ## extract uniprot acc if the row name is in the format of "tr|A0A024RCN7|A0A024RCN7_HUMAN"
  luni <- lapply(rowName, function(prot) {
    uniprot <- prot
    if(grepl("\\|", prot)) {
      uniprot <- unlist(strsplit(prot, "|", fixed=TRUE))[2]
    }
    uniprot <- toupper(uniprot)
    c(prot, uniprot)
  })

  ids <- as.data.frame(do.call(rbind, luni))
  colnames(ids) <- c("protein", "uniprot_acc")
  head(ids)
  tail(ids)

  ## create a uniprot_acc column to merge with idMap which also has a uniprot_acc column
  Contrastm <- as.data.frame(Contrast) %>%
                mutate(uniprot_acc=ids$uniprot_acc)

  dim(Contrastm)
  head(Contrastm)
  ContrastOut <- merge(idMap, Contrastm, by="uniprot_acc", all.y=T)
  noname <- is.na(ContrastOut$uniprot_name)
  ContrastOut[noname, "uniprot_name"] <- ContrastOut[noname, "uniprot_acc"]
  rownames(ContrastOut) <- ContrastOut$uniprot_acc
  ContrastOut <- ContrastOut[, c("uniprot_name", "seq_len", "uniprot_organism", "uniprot_description", colName)]

  colnames(ContrastOut) <- c("Symbol", "Length", "Organism", "Description", colName)
 
  return(ContrastOut)
}

formatForSAINT <- function(interMat, design, dataName, type="spc", useProhits=FALSE, gfp_from_prohits){
 #interMat <- data_mat
 
 #dataName <- "Control"

  design <- design %>% filter(SUBSET == dataName)
  ## create bait table from design
  IPs <- design$Experiment[design$BAIT != ""]
  baits <- design$BAIT[design$BAIT != ""]
  baits <- toupper(baits)
  groups <- design$GROUP[design$BAIT != ""]
  groups <- toupper(groups)

  baitDat <- data.frame(IPs, baits, groups)

  ## create prey table from interMat
  #selected_acc <- unique(row.names(interMat))
  #preyDat <- geneMap[geneMap$uniprot_acc %in% selected_acc, c("uniprot_acc",  "seq_len", "uniprot_name")]
  preyDat <- interMat[, c("X", "Length", "Symbol")]
  colnames(preyDat) <- c("protein_name", "length", "gene_name")
  preyDat <- preyDat[!duplicated(preyDat$protein_name),]

  preyDat$protein_name <- toupper(preyDat$protein_name)
  preyDat$gene_name <- toupper(preyDat$gene_name)

  # only consider selected IPs
  rownames(interMat) <- interMat[,"X"]
  interMat <- interMat[, IPs]
  ## create inter table from interMat
  interDat <- NULL
  for(i in 1:dim(interMat)[2]){
    #i <- 5
    IP <- colnames(interMat)[i]
    #bait <- unlist(strsplit(IP, "_"))[1]
    bait <- design[design$Experiment == IP, "BAIT"]
    bait <- toupper(bait)
    for(j in 1:dim(interMat)[1]){

      prey <- row.names(interMat)[j]
      prey <- toupper(prey)
      inter <- 0
      if(!is.na(interMat[j,i])){
        inter <- interMat[j,i]
      }

      if((prey %in% preyDat$protein_name) && (inter > 0)){
        interDat <- rbind(interDat, c(IP, bait, prey, inter))
      }
    }
  }

      ## incorporate GFP control from Prohits

  if(useProhits){
    gfp_selected <- c("GFP", "GFP1", "GFP2", "GFP3", "GFP4", "GFP_1", "GFP_2", "GFP_3", "GFP_4")
    gfp_from_prohits <- gfp_from_prohits[gfp_from_prohits$Bait.Gene.ID %in% gfp_selected, ]

    gfp_prey <- unique(gfp_from_prohits$Protein.ID)
    gfp_bait_id <- unique(gfp_from_prohits$Bait.ID)
    length(gfp_bait_id)
    IPs <- c(IPs, gfp_bait_id)
    baits <- c(baits, rep("GFP", length(gfp_bait_id)))
    controlList <- c(controlList, gfp_bait_id)

    groups <- rep("T", length(IPs))
    groups[IPs %in% controlList] <- "C"
    baitDat <- data.frame(IPs, baits, groups)

    ## create prey table
    selected_acc <- unique(c(selected_acc, gfp_prey))
    preyDat <- geneMap[geneMap$uniprot_acc %in% selected_acc, c("uniprot_acc",  "seq_len", "uniprot_name")]
    colnames(preyDat) <- c("protein_name", "length", "gene_name")
    preyDat <- preyDat[!duplicated(preyDat$protein_name),]

    preyDat$protein_name <- toupper(preyDat$protein_name)
    preyDat$gene_name <- toupper(preyDat$gene_name)

    ## create interaction table

    for(i in 1:dim(gfp_from_prohits)[1]){
      IP <- gfp_from_prohits$Bait.ID[i]
      bait <- "GFP"
      prey <- gfp_from_prohits$Protein.ID[i]
      inter <- gfp_from_prohits$Total.Number.Peptide[i]
      if(prey %in% preyDat$protein_name){
        interDat <- rbind(interDat, c(IP, bait, prey, inter))
      }
    }
  }

  colnames(interDat) <- c("IP", "bait", "prey", "interaction")
  unique(interDat[,2])

  write.table(baitDat, "bait.dat", sep="\t", quote=F, row.names=F, col.names=F)
  write.table(preyDat, "prey.dat", sep="\t", quote=F, row.names=F, col.names=F)
  write.table(interDat, paste(type, "interaction.dat", sep="_"), sep="\t", quote=F, row.names=F, col.names=F)

  return(list(baitDat, preyDat, interDat))
}

formatForMiSTR <- function(interMat, idMap, design, dataName, baitAsPrey=TRUE){
  #interMat <- Spectral_Count
  #idMap <- geneNameDesc
  #controlList <- controlIP
  #print(dim(interMat))
  #print(dim(idMap))
  #print(dataName)


  design <- design[design$SUBSET == dataName, ]
  #print(design)
  ## create bait table from design
  IPs <- design$Experiment[design$BAIT != ""]
  baits <- design$BAIT[design$BAIT != ""]
  baits <- toupper(baits)


  baitDat <- data.frame(IPs, baits)

  ## caution: only works for viral baits, and viral-viral interaction is not considered
  if(!baitAsPrey){
    interMat <- interMat[!interMat$Symbol %in% baits, ]
  }
  ## END caution

  interDat <- NULL
  for(IP in IPs){
    #IP <- IPs[3]
    bait <- baitDat[baitDat$IPs == IP, "baits"]

    for(preyAcc in row.names(interMat)){
      #preyAcc <- row.names(interMat)[1]
      preyAcc <- toupper(preyAcc)
      preyGene <- interMat[preyAcc, 1]
      preyMW <- idMap[idMap$uniprot_acc==preyAcc, "seq_mw"]
      SPC <- 0
      if(!is.na(interMat[preyAcc, IP])){
        SPC <- interMat[preyAcc, IP]
      }

      if(SPC > 0){
        interDat <- rbind(interDat, c(IP, preyAcc, preyGene, preyMW, SPC, bait))
      }
    }
  }
  interDat <- data.frame(interDat)
  colnames(interDat) <- c("IP", "preyAcc", "preyGene","preyMW", "SPC", "bait")
  interDat <- interDat[interDat$preyGene != interDat$bait, ]

  write.table(baitDat, paste(dataName,"mist_keys.txt",sep="_"), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(interDat, paste(dataName,"mist_data.txt",sep="_"), sep="\t", quote=F, row.names=F, col.names=T)

  return(list(baitDat, interDat))
}

## for online version of MiST at https://modbase.compbio.ucsf.edu/mist/
formatForMiSTweb <- function(interMat, design, dataName, baitAsPrey=TRUE){
  #interMat <- Spectral_Count

  design <- design[design$SUBSET == dataName, ]
  ## create bait table from design
  IPs <- design$Experiment[design$BAIT != ""]
  baits <- design$BAIT[design$BAIT != ""]
  baits <- toupper(baits)

  ## caution: only works for viral baits, and viral-viral interaction is not considered
  if(!baitAsPrey){
    interMat <- interMat[!interMat$Symbol %in% baits, ]
  }
  ## END caution

  firstRow <- c(rep("#", 3), "Exps", IPs)
  secondRow <- c(rep("#", 3), "Baits", baits)
  thirdRow <- c("Prey", "Description", "Length", "BaitSims", baits)
  dataF <- data.frame(interMat$Symbol, interMat$Description,  interMat$Length,  c("#"), interMat[, IPs])
  dataF <- rbind(firstRow, secondRow, thirdRow, dataF)
  dataF <- na.omit(dataF)

  return(dataF)
}


evaluateSAINTscore <- function(saintResults, threshold = 0.05){
  #saintResults <- "Clone12_peptideCount_saint_results.tab"
  #setwd("C:/data/raw/EDYTA/PROTEIN/210820_Zuyao/maxquant_processed")

  saint_table <- read.delim(saintResults)

  qp <- ecdf(saint_table$BFDR)

  ## redefine threshold in case of intensity based saint score
  ## use 5% quantile of BFDR as threshold instead of a hard one
  ## seems not better than hard threshold, not used

  #if(grepl("_intensity_", saintResults)){
    #q_df <- data.frame(bfdr=saint_table$BFDR, qt_bfdr=qp(saint_table$BFDR))
    #q_df <- q_df[order(q_df$bfdr),]
    #q_df_thr <- q_df[q_df$qt_bfdr < 0.05, ]
    #threshold <- max(q_df_thr$bfdr)
  #}

  outfile <- paste("threshold", threshold, saintResults, sep="_")
  re <- hist(type.convert(saint_table$BFDR), plot=F)
  outfig <- gsub("\\.tab", "\\.pdf", outfile)
  pdf(outfig)
  print(plot(re))
  print(plot(qp, xlab="BFDR", ylab="Quantile"))
  dev.off()

  saint_top <- saint_table[saint_table$BFDR <= threshold, ]
  write.table(saint_top, outfile, quote=F, row.name=F, sep="\t")

}
## translate uniprot sp|id|name to gene name
TranslateWithSymbol <- function(idList, idMap){
  #idList <- top
  #idMap <- geneName

  luni <- lapply(idList, function(prot) {
    uniprot <- prot
    if(grepl("sp\\|", prot)) {
      uniprot <- unlist(strsplit(prot, "|", fixed=TRUE))[2]
    }
    c(prot, uniprot)
  })

  ids <- as.data.frame(do.call(rbind, luni))
  names(ids) <- c("protein", "uniprot")
  head(ids)


  annotations.id <- merge(ids, idMap, by.x="uniprot", by.y="uniprot_acc", all.x=T)
  annotations.id <- unique(annotations.id)
  row.names(annotations.id) <- annotations.id$protein
  annotations.id <- annotations.id[idList,]
  annotations.id[is.na(annotations.id$uniprot_name), "uniprot_name"] <- annotations.id[is.na(annotations.id$uniprot_name), "uniprot"]

  return(annotations.id$uniprot_name)
}

pValueAgaintT0 <- function(df){
  #df <- compsMatt
  lv <- levels(as.factor(df[,1]))
  cn <- colnames(df)
  pMat <- matrix(1, nrow=length(lv), ncol=length(cn)-1, dimnames=list(lv, cn[2:length(cn)]))
  for(colm in 2:dim(df)[2]){
    #colm <- 2
    v0 <- as.numeric(df[df[,1] == lv[1], colm])
    print(v0)
    for(ro in 1:length(lv)){
      #ro <- 1
      v1 <- as.numeric(df[df[,1] == lv[ro], colm])

      if(sd(v1) == 0){v1[1] <- v1[1] + 0.0001}

      resT <- t.test(v0, v1)
      p <- round(resT$p.value, digits=4)
      p[is.na(p)] <- 1
      pMat[ro, colm-1] <- p

    }
  }
  pMat
}

plotMatrix <- function(data_name, my_data, timep, reps, plotLog="", legendPos="topleft", ylabel="", yrange=c(0, 50)){
  if(0){
    data_name <- "Viral proteins intensity"
    my_data <- subMat
    timep <- tp
    reps <- 3
    plotLog=""
    legendPos="topleft"
    ylabel=paste("Log2(", data_name, ")", sep="")
  }

  ntp <- length(timep)
  comps <- colnames(my_data)
  extreme <- length(comps)

  timePoints <- rep(timep, each=reps)
  geneCol <- rainbow(extreme)
  geneCol[1] <- "black"
  xx <- rep.int(seq(1:ntp), extreme)
  repcols <- rep(geneCol[1:extreme], each=ntp)

  compsMatt <- as.data.frame(my_data) %>%
    mutate(timePoints=timePoints)
  col_order <- c("timePoints", colnames(compsMatt)[!grepl("timePoints", colnames(compsMatt))])
  compsMatt <- compsMatt[, col_order]
  #compsMatt <- as.data.frame(cbind(timePoints, compsMat))

  pmat <- pValueAgaintT0(compsMatt)
  #write.table(pmat, paste(data_name, "pValue against time0.tab"), col.names=NA, row.names=TRUE, sep="\t", quote=F)

  compsMattA <- aggregate(. ~ timePoints, data=compsMatt, FUN=mean)
  compsMattSD <- aggregate(. ~ timePoints, data=compsMatt, FUN=sd)
  low <- as.matrix(compsMattA[,2:ncol(compsMattA)] - compsMattSD[,2:ncol(compsMattSD)])
  high <- as.matrix(compsMattA[,2:ncol(compsMattA)] + compsMattSD[,2:ncol(compsMattSD)])

  sigmat <- matrix("", nrow=dim(pmat)[1], ncol=dim(pmat)[2], dimnames=dimnames(pmat))
  sigmat[pmat < 0.05] <- "*"
  matplot(seq(1:ntp), compsMattA[,2:(extreme+1)], xaxt="n", log=plotLog, type="b", lwd=3, xlab="Time (hrs)", ylab=ylabel, main=data_name, ylim=yrange, col=geneCol, pch=20)
  axis(1, at=seq(1:ntp), labels=timep)
  arrows(xx, unlist(low[,1:extreme]), xx, unlist(high[,1:extreme]), col = repcols, angle = 90, length = 0.03, code = 3)
  text(xx, unlist(high[,1:extreme]), label=unlist(sigmat), col = repcols, cex=3, pos=3)
  legend(legendPos, legend=comps, col=geneCol, lty=1, pch=20, lwd=3)
}


## normalize data frame by log(sum) of columns
normalize_intensity <- function(df, start_col, end_col){
  #setwd("C:/data/raw/EDYTA/PROTEIN/210820_Zuyao/maxquant_processed")
  #LFQintensity <- read.delim("LFQintensity_from_210820_Zuyao.tab", header=T, sep="\t")
  #df <- LFQintensity
  #start_col <- 6
  #end_col <- ncol(LFQintensity)

  tobe <- df[,start_col:end_col]
  tobe[is.na(tobe)] <- 0

  c_sum <- log(apply(tobe, 2, sum))
  #constant_factor <- 1e+12 ## a rough average of sum of intensity in a sample
  scale_factor <- c_sum/median(c_sum) #constant_factor ## using total intensity
  scale_factor[is.na(scale_factor)] <- 0 # turn NA NaN into 0

  #print(scale_factor)
  ## a sample is considered as an outlier, if its median is 10 times greater than the median
  ## or less than 10% of the median, then it will not be scaled, original values are used.

  for(i in 1:length(scale_factor)){
    if(scale_factor[i] > 10 | scale_factor[i] < 0.1){
      scale_factor[i] <- 1
    }
  }

  scaled <- t(t(tobe)/scale_factor)

  scaled_sum <- log(apply(scaled, 2, sum))
  #print(scaled_sum)


  out <- cbind(df[,1:start_col-1], scaled)

  return(out)
}


plot_stats_treat_dose_group <- function(intensity, groupDesign, data_name, comparisons){
  #print(groupDesign)
  #intensity <- imputed_viral
  #data_name <- modification
  
  message("plotting STATS")
  #print(head(intensity))
  
  categs <- unlist(strsplit(comparisons, split=","))
  categ <- unlist(lapply(categs, function(x) unlist(strsplit(x, split="|", fixed=T))[1]))
  categRef <- unlist(lapply(categs, function(x) unlist(strsplit(x, split="|", fixed=T))[2]))
  names(categRef) <- categ

  colnames(intensity) <- gsub("Intensity.", "", colnames(intensity), fixed=T)
  #pheatmap(intensity)

  ## to avoid long processing time and large pdf files
  if(nrow(intensity) > 10000){
    sampling <- sample(rownames(intensity), 10000)
    intensity <- intensity[sampling, ]
    print("downsized to 10000 randomly selected genes")
  }

  ## turn headers into variables, so the header can be changed as needed
  NAME <- colnames(groupDesign)[1] 
  EXPERIMENT <- colnames(groupDesign)[2]
  TREATMENT <- colnames(groupDesign)[3] 
  DOSE <- colnames(groupDesign)[4] 
  GROUP <- colnames(groupDesign)[5]
  CELL <- colnames(groupDesign)[6]
  
  gene_list <- lapply(rownames(intensity), function(mod){
    exp_list <- lapply(colnames(intensity), function(exp){
      ## one exp might be assigned to several treatment,
      treat <- groupDesign[groupDesign[, EXPERIMENT] == exp, TREATMENT]
      Treatment <- treat
      Dose <- groupDesign[groupDesign[, EXPERIMENT] == exp, DOSE]
      Group <- groupDesign[groupDesign[, EXPERIMENT] == exp, GROUP]
      Cell <- groupDesign[groupDesign[, EXPERIMENT] == exp, CELL]

      Log2Intensity <- rep(intensity[mod, exp], length(treat))
      Gene <- rep(mod, length(treat))
      Experiment <- rep(exp, length(treat))
      data.frame(Gene, Experiment, Treatment, Dose=as.factor(Dose), 
            Group, Cell, Log2Intensity)
    })
    bind_rows(exp_list)
  })

  stats_df <- bind_rows(gene_list)
  colnames(stats_df) <- c("Gene", EXPERIMENT, TREATMENT, DOSE, 
                         GROUP, CELL, "Log2Intensity")
  #colnames(stats_df) <- c("Gene", "Experiment", "Treatment", "Dose", "Group", "Log2Intensity")
  stats_df <- stats_df[order(stats_df[, GROUP], decreasing=T), ]
  
  print(summary(stats_df))

  if(nrow(intensity) <= 0){
    message("The intensity table is empty, exiting")
    return()
  }
  
  for(cell in unique(groupDesign[, CELL])){
    cell_df <- stats_df %>%
      dplyr::filter(Cell == cell)
    if(nrow(cell_df) <= 0) next
    
    print(unique(cell_df$Treatment))
    print(cell)
    print(head(cell_df))
    
    pdf(paste(data_name, cell, "allGene_boxplot.pdf", sep="_"), height=10, width=8)
    #ggexport(p, filename=figname)

    p1a <- ggboxplot(cell_df, x=DOSE, y="Log2Intensity",
                   color = "black", palette = "jco", title=data_name, fill=DOSE,
                   xlab=FALSE, add="jitter", facet.by=c(TREATMENT, GROUP)) +
      geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, 
                 color="black") + # Add horizontal line at base mean
      stat_compare_means(label = "p.signif", ref.group=categRef[DOSE], 
                         method="t.test")   # Pairwise comparison against dose 0
    print(p1a)

 
    p1b <- ggboxplot(cell_df, x=TREATMENT, y="Log2Intensity", fill=TREATMENT,
                   color = TREATMENT, palette = "jco", title=data_name,
                   xlab=FALSE, add="jitter", facet.by=c(GROUP, DOSE)) +
      geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, color="black") + # Add horizontal line at base mean
      stat_compare_means(method = "anova", label.y = max(cell_df$Log2Intensity)+1)+        # Add global annova p-value
      stat_compare_means(label = "p.signif", ref.group=categRef[TREATMENT], method="t.test") # Pairwise comparison against NT
    print(p1b)
    
    p1c <- ggboxplot(cell_df, x=GROUP, y="Log2Intensity", fill=GROUP,
                     color = GROUP, palette = "jco", title=data_name,
                     xlab=FALSE, add="jitter", facet.by=c(TREATMENT, DOSE)) +
      geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, color="black") + # Add horizontal line at base mean
      stat_compare_means(method = "anova", label.y = max(cell_df$Log2Intensity)+1)+        # Add global annova p-value
      stat_compare_means(label = "p.signif", ref.group=categRef[GROUP], method="t.test") # Pairwise comparison against MOCK
    print(p1c)
    

    dev.off()
    pdf(paste(data_name, cell, "allGene_violinplot.pdf", sep="_"), height=10, width=8)
    #ggexport(p, filename=figname)

    p3a <- ggviolin(cell_df, x=DOSE, y="Log2Intensity", fill=DOSE,
                   color = "black", palette = "jco", title=data_name,
                   xlab=FALSE, add="boxplot", add.params = list(fill="white", color="black"), 
                   facet.by=c(TREATMENT, GROUP)) +
      geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, 
                 color="black") + # Add horizontal line at base mean
      stat_compare_means(label = "p.signif", ref.group=categRef[DOSE], 
                         method="t.test")   # Pairwise comparison against dose 0
    print(p3a)

  
    p3b <- ggviolin(cell_df, x=TREATMENT, y="Log2Intensity", fill=TREATMENT,
                   color = TREATMENT, palette = "jco", title = data_name,
                   xlab=FALSE, add="boxplot", add.params = list(fill="white", color="black"), 
                   facet.by=c(GROUP, DOSE)) +
      geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, color="black") + # Add horizontal line at base mean
      stat_compare_means(method = "anova", label.y = max(cell_df$Log2Intensity)+1) +      # Add global annova p-value
      stat_compare_means(label = "p.signif", ref.group=categRef[TREATMENT], method="t.test") # Pairwise comparison against NT
    print(p3b)
    
    p3c <- ggviolin(cell_df, x=GROUP, y="Log2Intensity", fill=GROUP,
                    color = GROUP, palette = "jco", title = data_name,
                    xlab=FALSE, add="boxplot", add.params = list(fill="white", color="black"), 
                    facet.by=c(TREATMENT, DOSE)) +
      geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, color="black") + # Add horizontal line at base mean
      stat_compare_means(method = "anova", label.y = max(cell_df$Log2Intensity)+1) +      # Add global annova p-value
      stat_compare_means(label = "p.signif", ref.group=categRef[GROUP], method="t.test") # Pairwise comparison against MOCK
    print(p3c)
    
    dev.off()


  if(length(unique(cell_df$Gene)) < 50){
    pdf(paste(data_name, cell, "perGene_boxplot.pdf", sep="_"), height=10, width=8)

    for(mod in unique(cell_df$Gene)){
      #mod <- "GPR89B|203"
      sub_df <- cell_df[cell_df$Gene == mod,]
      p5a <- ggboxplot(sub_df, x=DOSE, y="Log2Intensity",
                   color = DOSE, palette = "ucscgb", title=data_name, subtitle=mod,
                   xlab=FALSE, add="jitter", facet.by=c(TREATMENT, GROUP)) +
        stat_compare_means(label = "p.signif", ref.group=categRef[DOSE], method="t.test")
      print(p5a)

      p5b <- ggboxplot(sub_df, x=TREATMENT, y="Log2Intensity",
                     color = TREATMENT, palette = "ucscgb", title=data_name, subtitle=mod,
                     xlab=FALSE, add="jitter", facet.by=c(GROUP, DOSE)) + 
        stat_compare_means(label = "p.signif", ref.group=categRef[TREATMENT], method="t.test")
      print(p5b)
      p5c <- ggboxplot(sub_df, x=GROUP, y="Log2Intensity",
                       color = GROUP, palette = "ucscgb", title=data_name, subtitle=mod,
                       xlab=FALSE, add="jitter", facet.by=c(TREATMENT, DOSE)) + 
        stat_compare_means(label = "p.signif", ref.group=categRef[GROUP], method="t.test")
      print(p5c)
    }
    dev.off()
  }

  pdf(paste(data_name, cell, "perGene_violinplot.pdf", sep="_"), height=10, width=8)
  for(mod in unique(cell_df$Gene)){
    #mod <- "GPR89B|203"
    sub_df <- cell_df[cell_df$Gene == mod,]

    p7a <- ggviolin(sub_df, x=GROUP, y="Log2Intensity", fill=GROUP,
                   color = GROUP, palette = "jco", title=data_name, subtitle=mod,
                   xlab=FALSE, add="boxplot", add.params = list(fill="white", color="black"), 
                   facet.by=c(TREATMENT, DOSE)) +
      stat_compare_means(label = "p.signif", ref.group=categRef[GROUP], method="t.test")   # Pairwise comparison against MOCK
    print(p7a)

    p7b <- ggviolin(sub_df, x=TREATMENT, y="Log2Intensity", fill=TREATMENT,
                   color = TREATMENT, palette = "jco", title = data_name, subtitle=mod,
                   xlab=FALSE, add="boxplot", add.params = list(fill="white", color="black"),
                   facet.by=c(GROUP, DOSE)) +
      stat_compare_means(method = "anova", label.y = max(sub_df$Log2Intensity)+1) +      # Add global annova p-value
      stat_compare_means(label = "p.signif", ref.group=categRef[TREATMENT], method="t.test") # Pairwise comparison against NT
    print(p7b)
    
    p7c <- ggviolin(sub_df, x=DOSE, y="Log2Intensity", fill=DOSE,
                    color = DOSE, palette = "jco", title = data_name, subtitle=mod,
                    xlab=FALSE, add="boxplot", add.params = list(fill="white", color="black"), 
                    facet.by=c(TREATMENT, GROUP)) +
      stat_compare_means(method = "anova", label.y = max(sub_df$Log2Intensity)+1) +      # Add global annova p-value
      stat_compare_means(label = "p.signif", ref.group=categRef[DOSE], method="t.test") # Pairwise comparison against dose 0
    print(p7c)

  }
  dev.off()
  }
}

plot_stats_treat_dose_group_wrProteo <- function(path1, groupDesign, data_name, comparisons, mainSpec="HUMAN"){
  #data_name <- "SARS2"

  library("knitr")
  library("wrMisc")
  library("wrProteo")
  library("wrGraph")
  library("fdrtool")
  library("EnhancedVolcano")
   

  source("C:/GREENBLATT/Rscripts/RNAseq/Function_analysis_for_RNAseq_lib.R")

  cwd <- getwd()
  outd <- file.path(cwd, "Differential_analysis")
  if(!dir.exists(outd)){
    dir.create(outd)
  }
  setwd(outd)
  
  print(groupDesign)
  ## turn headers into variables, so the header can be changed as needed
  NAME <- colnames(groupDesign)[1] 
  EXPERIMENT <- colnames(groupDesign)[2]
  TREATMENT <- colnames(groupDesign)[3] 
  DOSE <- colnames(groupDesign)[4] 
  GROUP <- colnames(groupDesign)[5]
  CELL <- colnames(groupDesign)[6]
  
  categs <- unlist(strsplit(comparisons, split=","))
  
  categ <- unlist(lapply(categs, function(x) 
    unlist(strsplit(x, split="|", fixed=T))[1]))
  categRef <- unlist(lapply(categs, function(x) 
    unlist(strsplit(x, split="|", fixed=T))[2]))
  names(categRef) <- categ
  
  conditions <- apply(groupDesign[, categ], 1, function(x) paste(x, collapse="_"))
  conditions <- gsub(" ", "", conditions)
  names(conditions) <- groupDesign[, EXPERIMENT]
  

  pdf("median_normalization.pdf", width=16, height=8)
  specPr <- c(conta="CON__", mainSpecies=mainSpec, spike=data_name)
  dataMQ <- readMaxQuantFile(path1, specPref=specPr, normalizeMeth="median",
                             sampleNames = unique(groupDesign[, EXPERIMENT]))
  dev.off()

  dim(dataMQ$quant)
  summary(dataMQ$quant)
  head(dataMQ$quant)

  grp9 <- conditions[unique(names(conditions))] ## make sure sample names in the data matrix and those in design

  pdf("missing_value_inspection.pdf", width=8, height=8)
  matrixNAinspect(dataMQ$quant, gr=grp9, retnNA = TRUE, tit="Missing value inspection")
  dev.off()

  #dataMQimp <- matrixNAneighbourImpute(dataMQ$quant, gr=grp9, tit="Imputation results")
  
  pdf("missing_value_imputation.pdf", width=8, height=8)
  testMQall <- testRobustToNAimputation(dataMQ, gr=grp9, plotHist=T)
  dev.off()
  
  head(testMQall$datImp)
  print("contrasts:")
  print(testMQall$contrasts)

  pdf("PCA_plot.pdf", width=8, height=10)
  plotPCAw(testMQall$datImp, sampleGrp=grp9, tit="PCA on MaxQuant (NAs imputed)", rowTyName="proteins", useSymb2=0)
  dev.off()

  #MAplotW(testMQ)

  pdf("volcano_plots.pdf")
  for(i in 1:ncol(testMQall$contrasts)){
    # i <- 6
    contrast <- colnames(testMQall$contrasts)[i]
    factors <- unlist(strsplit(contrast, split="-",fixed=T))
    treatments <- unlist(lapply(factors, function(x)unlist(strsplit(x, split="_",fixed=T))[1]))
    doses <- unlist(lapply(factors, function(x)unlist(strsplit(x, split="_",fixed=T))[2]))
    groups <- unlist(lapply(factors, function(x)unlist(strsplit(x, split="_",fixed=T))[3]))
    
    grp_order <- 0
    contrastN <- NULL
    
    ## reorder group contrast such that background group is fixed as the second element of the pair
    if(treatments[1] == treatments[2] && doses[1] == doses[2] && categRef[GROUP] %in% groups){
       
       if(groups[1] == categRef[GROUP]){
          grp_order <- 1
          contrastN <- paste0(factors[2], "_vs_", factors[1])
       }else{
          grp_order <- -1
          contrastN <- paste0(factors[1], "_vs_", factors[2])
       }
       
    }
    
    if(groups[1] == groups[2] && doses[1] == doses[2] && categRef[TREATMENT] %in% treatments){
       
       if(treatments[1] == categRef[TREATMENT]){
          grp_order <- 1
          contrastN <- paste0(factors[2], "_vs_", factors[1])
       }else{
          grp_order <- -1
          contrastN <- paste0(factors[1], "_vs_", factors[2])
       }
    }
    
    if(groups[1] == groups[2] && treatments[1] == treatments[2] && categRef[DOSE] %in% doses){
      
      if(doses[1] == categRef[DOSE]){
        grp_order <- 1
        contrastN <- paste0(factors[2], "_vs_", factors[1])
      }else{
        grp_order <- -1
        contrastN <- paste0(factors[1], "_vs_", factors[2])
      }
    }
    
    if(grp_order != 0){
       
       res_all <- extractTestingResults(testMQall, compNo=i, statTy="BH", addTy=NULL, nSign=3, thrsh=1, FCthrs=0)
       res_all <- res_all[!grepl("CON__", res_all$Accession),]
       res_all[,5] <- grp_order*res_all[,5]
       selected <- which(colnames(res_all) %in% paste0("av.", factors))
       res_all <- res_all[, c(1:5, selected)]
       res_all <- res_all[order(res_all[, "FDR"]), ]
       write.table(res_all, paste(contrastN, "results_ALL.tab", sep="_"), sep="\t", row.names=F, quote=F)
       
       fc_list <- res_all[!is.na(res_all$GeneName),5]
       names(fc_list) <- res_all$GeneName[!is.na(res_all$GeneName)]
       
       try(run_gseGO_simpleList(fc_list, contrastN))
       
       lfc_cutoff <- 1
       p_cutoff <- 0.05
       p <- EnhancedVolcano(res_all,
                            lab = res_all$GeneName,
                            x = colnames(res_all)[5],
                            y = colnames(res_all)[4],
                            title = paste(contrastN),
                            subtitle = paste("N =", nrow(res_all)),
                            caption = paste("log2 fold change cutoff = +/-", lfc_cutoff,  " adjusted pvalue cutoff = ", p_cutoff, sep=""),
                            pCutoff = p_cutoff,
                            FCcutoff = lfc_cutoff,
                            labSize = 4,
                            pointSize = 3.0,
                            #drawConnectors = TRUE,
                            widthConnectors = 0.75)
       print(p)
    }
  }
  dev.off()

  return(testMQall$datImp)
}

extract_modifiction_intensity <- function(modification_table, design, modification){
  columns <- colnames(modification_table)
  Con.value <- columns %like% "contaminant"
  NonContaminants <- modification_table[, Con.value] == ""
  modification_table <- modification_table[NonContaminants,]

  #Remove Reverse Peptides
  NotReverse <- modification_table[, "Reverse"] == ""|is.na(modification_table[, "Reverse"])
  modification_table <- modification_table[NotReverse,]

  intensity_col <- paste("Intensity", design$Experiment, sep=".")
  intensity_table <- modification_table[, c(colnames(modification_table)[1:7], intensity_col)]
  dim(intensity_table)
  proteinColumn <- "Gene.names" # Maxquant recommmends this over "Protein.IDs"
  proteinIDs <- as.character(intensity_table[,proteinColumn])
  protein <- unlist(lapply(proteinIDs, function(x){unlist(strsplit(x, ";", fixed=T))[1]}))
  positionColumn <- "Positions.within.proteins" # Maxquant recommmends this over "Protein.IDs"
  positions <- as.character(intensity_table[,positionColumn])
  position <- unlist(lapply(positions, function(x){unlist(strsplit(x, ";", fixed=T))[1]}))

  rownames(intensity_table ) <- paste(protein, position, sep="_")
  write.table(intensity_table, paste(modification,"_intensity.tab",sep=""), sep="\t", col.names=NA, quote=F)
  return(as.matrix(intensity_table[,8:ncol(intensity_table)]))
}

## replace missing value with the minimum non-missing value of the matrix (LOD, limit of detection)
impute_intensity <- function(intensity_mat, rm_empty=FALSE, method="RF"){
  intensity_mat[intensity_mat==0] <- NaN
  na_per_row <- apply(intensity_mat, 1, function(x)sum(is.na(x)))
  na_per_col <- apply(intensity_mat, 2, function(x)sum(is.na(x)))
  # remove rows with NaN only
  intensity_mat <- intensity_mat[na_per_row < ncol(intensity_mat),]


  if(rm_empty){
    # remove columns with NaN only
    intensity_mat <- intensity_mat[, na_per_col < nrow(intensity_mat)]
  }
  NAs <- sum(na_per_row)/(nrow(intensity_mat)*ncol(intensity_mat))
  if(NAs > 0.5){
    # if fraction of missing value is greater than 50%, force use of LOD
    method <- "LOD"
  }

  imputed_mat <- intensity_mat
  if(method == "LOD"){
    min_intensity <- min(intensity_mat, na.rm = T)
    imputed_mat[is.na(imputed_mat)] <- min_intensity
  }

  if(method == "RF"){
    imputed_mat <- missForest(intensity_mat)$ximp
  }

  return(imputed_mat)
}

combine_inputs_for_SAINT <- function(saintDirs, outDir){

  if(!dir.exists(outDir)) dir.create(outDir, recursive=T)
  dataNames <- names(saintDirs)
  types <- c("spc", "peptideCount", "intensity")
  for(atype in types){
    #atype <- types[1]
    bait_out <- NULL
    prey_out <- NULL
    int_out <- NULL

    for(aName in dataNames){
      #aName <- dataNames[1]
      dir <- saintDirs[aName]
      bait_df <- read.delim(file.path(dir, atype, "bait.dat"), sep="\t", header=F)
      prey_df <- read.delim(file.path(dir, atype, "prey.dat"), sep="\t", header=F)
      int_df <- read.delim(file.path(dir, atype, paste(atype, "interaction.dat", sep="_")), sep="\t", header=F)

      bait_df[,1] <- paste(aName, bait_df[,1], sep="_")
      int_df[,1] <- paste(aName, int_df[,1], sep="_")

      bait_out <- rbind(bait_out, bait_df)
      prey_out <- rbind(prey_out, prey_df)
      int_out <- rbind(int_out, int_df)
    }

    prey_out <- unique(prey_out)

    if(!dir.exists(file.path(outDir, atype))) dir.create(file.path(outDir,atype), recursive=T)
    setwd(file.path(outDir, atype))

    write.table(bait_out, "bait.dat", col.names=F, row.name=F, sep="\t", quote=F)
    write.table(prey_out, "prey.dat", col.names=F, row.name=F, sep="\t", quote=F)
    write.table(int_out, "interaction.dat", col.names=F, row.name=F, sep="\t", quote=F)

    print("Running SAINT ...")
    SAINT_path <- "C:/GREENBLATT/SAINT/SAINTexpress_v3.6.3__2018-03-09/Precompiled_binaries/Windows64"
    nControl <- length(bait_out[bait_out[,3] == "C",3])
    output <- paste(atype, "saint_results.tab", sep="_")

    saintCommand <- "SAINTexpress-spc"
    if(atype == "intensity"){
      saintCommand <- "SAINTexpress-int"
    }

    ex <- system(paste(file.path(SAINT_path, saintCommand), paste("-L",nControl,sep=""), "interaction.dat", "prey.dat", "bait.dat"))

    if(ex == 0){
      print("SAINT is run successfully!")
      system(paste("mv", "list.txt", output))
      evaluateSAINTscore(output, 0.01)
    }else{
      print(paste("SAINT exited with error code", ex))
    }
  }
}


combine_processed_results <- function(inDirs, outDir){
  
  if(!dir.exists(outDir)) dir.create(outDir, recursive=T)
  dataNames <- names(inDirs)
    
  for(prefix in c("Total_spectral_count_from_", "Peptide_count_from_","LFQintensity_from_")){
    
    a_list <- list()
    for(aName in dataNames){
      #aName <- dataNames[1]
      dir <- inDirs[aName]
      a_list[[aName]] <- read.delim(file.path(dir, paste0(prefix, aName, ".tab")), sep="\t", header=T)
    }
  
    a_df <- Reduce(full_join, a_list)
    samps <- colnames(a_df)[5:ncol(a_df)]
    if(prefix == "Total_spectral_count_from_") {
      writeLines(samps, file.path(outDir,"experimentalDesignTemplate_S.txt"))
    }
    
    write.table(a_df, file.path(outDir, paste0(prefix, paste(dataNames,collapse="_"), "_combined.tab")), sep="\t", col.names=T, row.names=F, quote=F)
  }
  
}

