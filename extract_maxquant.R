#library(proteus)
#library(proteusLabelFree)


library(limma)
library(plotrix)


extract_spc <- function(pg_table, folder, geneNameDesc){

  Spectral_Count <- getSpectralCount(pg_table)

  if(dim(Spectral_Count)[1] > 5 & dim(Spectral_Count)[2] > 5){
    pdf(paste("Total_spectral_count_from_proteinGroup_heatmap.pdf", sep=""),  width=12, height=12)
    pheatmap(Spectral_Count, margins=c(10,10))
    dev.off()
  }

  Spectral_Count <- AnnotateWithSymbol(Spectral_Count, geneNameDesc)

  write.table(Spectral_Count, paste("Total_spectral_count_from_",folder,".tab", sep=""), col.names=NA, row.names=TRUE, sep="\t", quote=F)

  return(Spectral_Count)
}

extract_peptideCount <- function(pg_table, folder, geneNameDesc){

  Peptide_Count <- getPeptideCount(pg_table)
  
  longName <- colnames(Peptide_Count)
  shortName <- gsub("Peptides.", "", longName)
  colnames(Peptide_Count) <- shortName

  #Peptide_Count <- as.matrix(Peptide_Count[sort(rownames(Peptide_Count)),])

  if(dim(Peptide_Count)[1] > 5 & dim(Peptide_Count)[2] > 5){
    pdf(paste("Peptide_count_from_proteinGroup_heatmap.pdf", sep=""),  width=12, height=12)
    pheatmap(Peptide_Count, margins=c(10,10))
    dev.off()
  }

  Peptide_Count <- AnnotateWithSymbol(Peptide_Count, geneNameDesc)

  write.table(Peptide_Count, paste("Peptide_count_from_",folder,".tab", sep=""), col.names=NA, row.names=TRUE, sep="\t", quote=F)
  
  return(Peptide_Count)
}

extract_intensity <- function(pg_table, folder, geneNameDesc, design, isLFQ){

  LFQintensity <- NULL
  experiments <- design$Experiment
  
  experiments <- gsub("#", ".", experiments, fixed=T)
  
  LFQintensity <- getLFQintensity(pg_table, LFQ=isLFQ)
  
  
  if(!is.null(LFQintensity)){

    LFQintensity <- AnnotateWithSymbol(LFQintensity, geneNameDesc)
    
    write.table(LFQintensity, paste("LFQintensity_from_",folder,".tab", sep=""), col.names=NA, row.names=TRUE, sep="\t", quote=F)
    #half_minimum <- min(LFQintensity[LFQintensity > 0])/2

    logLFQintensity <- LFQintensity  ## log transformed for subsequent analysis
    logLFQintensity[, experiments] <- log2(LFQintensity[, experiments] + 1)
    if(dim(LFQintensity)[1] > 5 & dim(LFQintensity)[2] > 9){
       pdf(paste("LFQintensity_from_proteinGroup_heatmap.pdf", sep=""),  width=12, height=12)
       pheatmap(logLFQintensity[, experiments], margins=c(10,10))
       dev.off()
    }
    
  }

  return(LFQintensity)
}

## type=c("spc", "peptideCount", "intensity")
run_saint <- function(datafile, designfile, type="spc", SAINT_path, dataName, useProhits, gfp_from_prohits, outdir){
  #data_mat <- adata
  #type <- atype
  #dataName <- "Clone1"
  data_mat <- read.delim(datafile, sep="\t", header=TRUE)
  design <- read.delim(designfile, sep="\t", header=TRUE)

  print(paste("Running SAINT on ...", datafile))
  nControl <- nrow(filter(design, GROUP == "C" & SUBSET == dataName))
  
  setwd(outdir)
  output <- paste(dataName, type, "saint_results.tab", sep="_")

  saintCommand <- "SAINTexpress-spc"
  if(type == "intensity"){
    saintCommand <- "SAINTexpress-int"
  }

  if(type == "peptideCount"){
    formatForSAINT(interMat=data_mat, design=design, dataName, type=type, useProhits, gfp_from_prohits)
  }else{
    formatForSAINT(interMat=data_mat, design=design, dataName, type=type, FALSE, NULL)
  }

  ex <- system(paste(file.path(SAINT_path, saintCommand), paste("-L",nControl,sep=""), paste(type, "interaction.dat", sep="_"), "prey.dat", "bait.dat"))

  if(ex == 0){
    print("SAINT is run successfully!")
    system(paste("mv", "list.txt", output))
    evaluateSAINTscore(output, 0.01)
  }else{
    print(paste("SAINT exited with error code", ex))
  }
}

prepare_for_MiST <- function(data_mat, geneNameDesc, expOrder, design, dataName, baitAsPrey, type){
  reOrder <- c(colnames(data_mat)[1:4], expOrder)
  if(type == "intensity"){
    data_mat[, 5:ncol(data_mat)] <- log2( data_mat[, 5:ncol(data_mat)] + 1)
  }
  formatForMiSTR(data_mat[, reOrder], geneNameDesc, design, dataName, baitAsPrey)
  MiSTinput <- formatForMiSTweb(data_mat, design, dataName, baitAsPrey)
  write.table(MiSTinput, paste(dataName, type,"forMiST.tab", sep="_"), sep="\t", quote=F, row.names=F, col.names=F)
}

plot_viral_protein <- function(logNormLFQintensity, groupDesign, viral_name="SARS2", folder){
  #groupDesign <- read.delim(paste("C:/data/raw/EDYTA/PROTEIN", folder, "combined", "experimentalDesignTemplate_TS.txt", sep="/"))
  data_name <- "LFQintensity"

  experiments <- groupDesign$Experiment
  ngroups <- nrow(unique(groupDesign[, c("Time", "Group", "Cell")]))
  nreplicates <- length(experiments)/ngroups
  CellLines <- unique(groupDesign$Cell)

  intensity_data <- logNormLFQintensity

  vpName <- rownames(intensity_data[intensity_data$Organism==viral_name,])
  print(data_name)
  print(dim(intensity_data))

  data_name_copy <- data_name

  for(cellLine in CellLines){
    #cellLine <- CellLines[1]
    exps <- groupDesign$Group[groupDesign$Cell == cellLine]
    tps <- groupDesign$Time[groupDesign$Cell == cellLine]
    conditions <- paste(exps, tps, sep="_")
    exps <- unique(exps)
    exps <- exps[exps != ""]
    tp <- sort(unique(tps))


    #data_name <- paste(cellLine, data_name_copy, sep="_")
    dataMat <- (intensity_data[, grepl(cellLine, colnames(intensity_data))])
    #rownames(dataMat) <- as.character(intensity_data$Symbol)

    #rowsum <- apply(dataMat, 1, sum)
    #rowsum0 <- rowsum[which(rowsum == 0)]
    #dataMat <- dataMat[!rownames(dataMat) %in% names(rowsum0), ]
    #dataMat["N", "Huh7_0hrs_1"] <- 0 # outlier
    col_labels <- colnames(dataMat)

    head(dataMat)
    plotMDS(dataMat, top=10000, labels=col_labels)

    corMatrix <- cor(t(dataMat))

    #pheatmap(dataMat)
    DTC <- as.matrix(corMatrix[row.names(corMatrix) %in% vpName, !colnames(corMatrix) %in% vpName])
    DTC[is.na(DTC)] <- 0
    corTable <- as.data.frame(t(DTC))

    #corTable <- corTable[order(corTable$P0DTC9),]
    write.table(corTable, paste(cellLine, data_name, folder, "correlation_with_viral_proteins.tab", sep="_"), row.names = T, col.names = NA, sep="\t", quote=F)

    if(ncol(DTC) > 3 && nrow(DTC) > 3){
      pdf(paste(cellLine, data_name, folder, "viral_correlation_heatmap.pdf", sep="_"), width=8, height=8)
      pheatmap(t(DTC), main="correlation between viral and human proteins")
      dev.off()
    }

    pdf(paste(cellLine, data_name, folder, "viral_protein_expression.pdf", sep="_"), width=12, height=6)
    opar <- par(mfrow=c(1,length(exps)))

    vpsMat <- t(dataMat[vpName, col_labels])
    head(vpsMat)

    yrange <- c(min(vpsMat)*0.75, max(vpsMat)*1.25)
    yrange

    for(exp in exps){
      #exps <- "INF"
      subgroup <- groupDesign[groupDesign$Group == exp & groupDesign$Cell == cellLine,]
      subgroup <- subgroup[order(subgroup$Time),]
      subMat <- as.data.frame(vpsMat[subgroup$Experiment, ])
      colnames(subMat) <- colnames(vpsMat)

      plotMatrix(exp, subMat, timep=tp, reps=nreplicates, ylabel=paste("Log2(", data_name, ")", sep=""), yrange=yrange)
    }

    par(opar)
    dev.off()

    pdf(paste(cellLine, data_name, folder, "correlation.pdf", sep="_"), width=12, height=8)

    opar <- par(mfrow=c(2,length(exps)))
    nSelect  <- 5
    for(vp in vpName){
      #vp <- vpName[1]
      DTCv <- corMatrix[vp,]
      DTCv <- DTCv[order(DTCv, decreasing=T)]
      DTCv <- DTCv[!is.na(DTCv)]
      #DTCv_0_5 <- DTCv[abs(DTCv) > 0.5]
      #DTCv_0_5df <- as.data.frame(DTCv_0_5, row.names=names(DTCv_0_5))
      #colnames(DTCv_0_5df) <- "cor"
      #DTCv_0_5A <- AnnotateWithSymbol(DTCv_0_5df, geneName)
      #DTCv_0_5A <- DTCv_0_5A[order(as.numeric(DTCv_0_5A[, "cor"]), decreasing=T),]
      DTCv <- DTCv[!names(DTCv) %in% vpName]

      ## plot only when top and bottom correlated genes exist
      if(length(DTCv) > nSelect*2){
        top <- c(vp, names(head(DTCv, nSelect)))
        DTCv[top]
        bottom <- c(vp, names(tail(DTCv, nSelect)))
        DTCv[bottom]


        matList <- list()
        matList[["top"]] <- t(dataMat[top, col_labels])
        matList[["bottom"]] <- t(dataMat[bottom, col_labels])
        for(Mat in c("top", "bottom")){
          vpsMat <- matList[[Mat]]
          for(exp in exps){
            subgroup <- groupDesign[groupDesign$Group == exp & groupDesign$Cell == cellLine,]
            subgroup <- subgroup[order(subgroup$Time),]
            subMat <- as.data.frame(vpsMat[subgroup$Experiment, ])
            colnames(subMat) <- colnames(vpsMat)
            plotMatrix(paste(vp, exp, Mat, sep=":"), subMat, timep=tp, reps=nreplicates, ylabel=paste("Log2(", data_name, ")", sep=""), yrange=yrange)
          }
        }
      }
    }

    par(opar)
    dev.off()

  }
}

collect_protein_info <- function(aaInfo, fasta, organism, viral, viral_fasta, viral_name, id_map){
  if(0){
     aaInfo <- "resource/Amino_acids_info.tab"
     fasta <- "resource/mouse_uniprot_reference_protoeme_10090.fasta/UP000000589_10090.fasta"
     organism <- "MOUSE"
     id_map <- NULL
     viral <- F
     viral_fasta <- "resource/sars2_uniprot_reference_proteome_2697049.fasta/REFSEQ_Wuhan_Hu_1_nr.fasta"
     viral_name <- "SARS2"
  }
     
  aaTable <- read.delim(aaInfo, header=T)
  wList <- aaTable$Molecular_Weight
  names(wList) <- aaTable$Code

  if(!is.null(id_map)){
    geneNameDesc <- read.delim(id_map)
  }else{
    geneNameDesc <- build_id_map(wList, fasta, organism)
  }

  if(viral){

    #if(viral_name == "SARS2"){
    # viral_map <- build_id_map_sars2(wList, viral_fasta, "SARS2")
    #}else{
      viral_map <- build_id_map(wList, viral_fasta, viral_name)
    #}
    colnames(viral_map) <- colnames(geneNameDesc)
    geneNameDesc  <- rbind(geneNameDesc, viral_map)

  }
  return(geneNameDesc)
}
## this is the work-horse of the program
run_extract <- function(hwd, folder, geneNameDesc, isLFQ, outdir, TS, service){
   if(0){
     hwd <- "C:/data/raw/EDYTA/PROTEIN"
     folder <- "220316_Shamira_TH"
     outdir <- "maxquant_processed"
     isLFQ <- FALSE
     TS <- FALSE
   }

  wd <- file.path(hwd, folder, outdir, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd)
  }
  setwd(wd)

  design <- read.delim(file.path(hwd, folder, "combined", "experimentalDesignTemplate.txt"), colClasses = c("character", "numeric", "character", "character"))
  print(design)
  protein_info <- c("Symbol", "Length", "Organism", "Description")

  proteinGroup_file <- file.path(hwd, folder,"combined/txt/proteinGroups.txt")
  pg_table <- read.delim(proteinGroup_file, header=T, stringsAsFactors = T)
  print(hwd)
  print(folder)
  #print(pg_table)

  ## spectral count ####
  print("processing spectral count")
  try(Spectral_Count <- extract_spc(pg_table, folder, geneNameDesc))

## peptide count ####
  print("processing peptide count")
  try(Peptide_Count <- extract_peptideCount(pg_table, folder, geneNameDesc))

## LFQ intensity ####
  print("processing LFQ intensity")
  try(LFQintensity <- extract_intensity(pg_table, folder, geneNameDesc, design, isLFQ))

  logLFQintensity <- logNormLFQintensity <- LFQintensity
  norm_data <- normalize_intensity(LFQintensity, 5, ncol(LFQintensity))
  logLFQintensity[, 5:ncol(LFQintensity)] <- log2(LFQintensity[, 5:ncol(LFQintensity)] + 1)
  logNormLFQintensity[, 5:ncol(norm_data)] <- log2(norm_data[, 5:ncol(norm_data)] + 1)
  
  if(!service){
    write.table(norm_data, paste("TotalIntensity_normalized_LFQintensity_from_",folder,".tab", sep=""), col.names=NA, row.names=TRUE, sep="\t", quote=F)
    write.table(logNormLFQintensity, paste("TotalIntensity_normalized_log2LFQintensity_from_",folder,".tsv", sep=""), col.names=NA, row.names=TRUE, sep="\t", quote=F)
  }


  if(TS){
    groupDesign <- read.delim(file.path(hwd, folder, "combined", "experimentalDesignTemplate_TS.txt"))
    ngroups <- nrow(unique(groupDesign[, c("Time", "Group", "Cell")]))
    nreplicates <- length(experiments)/ngroups
  }else{
    groupDesign <- design
    ngroups <- nreplicates <- 1
  }

  experiments <- groupDesign$Experiment
  experiments <- gsub("#", ".", experiments, fixed=T)
  lfq <- colSums(LFQintensity[, experiments])
  loglfq <- colSums(logLFQintensity[, experiments])
  cs <- colSums(norm_data[, experiments])
  lcs <- colSums(logNormLFQintensity[, experiments])


  if(!service){
    pdf(paste("Total_intensity", folder, "per_sample.pdf", sep="_"))
    cl <- rep(c(1:ngroups),each=nreplicates)
    opar <- par(mar=c(2, 20, 3, 1)+0.1)
    barplot2(lfq, main="Total intensity", yaxt = "n", horiz=T, col=cl,  width=0.85, xpd=T)
    staxlab(2,1:length(cs), experiments, srt=0, line.spacing=2, cex=1, col=cl, top.line=0, adj=1)
    barplot2(cs, main="Total normalized intensity", yaxt = "n", horiz=T, col=cl,  width=0.85, xpd=T)
    staxlab(2,1:length(cs), experiments, srt=0, line.spacing=2, cex=1, col=cl, top.line=0, adj=1)
    barplot2(loglfq, main="Total log intensity", yaxt = "n", horiz=T, col=cl,  width=0.85, xpd=T)
    staxlab(2,1:length(lcs), experiments, srt=0, line.spacing=2, cex=1, col=cl, top.line=0, adj=1)
    barplot2(lcs, main="Total log normalized intensity", yaxt = "n", horiz=T, col=cl,  width=0.85, xpd=T)
    staxlab(2,1:length(lcs), experiments, srt=0, line.spacing=2, cex=1, col=cl, top.line=0, adj=1)
    par(opar)
    dev.off()
  }



  ## show relationship between quantities ####
  print("processing counts correlation")
  pdf("relationship_between_quantity_measures.pdf")
  opar <- par(mfrow=c(2,3))
  LFQintensity <- LFQintensity[rownames(Peptide_Count),]
  Spectral_Count <- Spectral_Count[rownames(Peptide_Count),]
  for(i in 5:dim(Peptide_Count)[2]){
    if(ncol(Spectral_Count) == ncol(Peptide_Count)){plot(Spectral_Count[, i]+1, Peptide_Count[, i]+1, log="xy", main=colnames(Peptide_Count)[i])}
    if(ncol(Spectral_Count) == ncol(LFQintensity)){plot(Spectral_Count[, i]+1, LFQintensity[, i]+1, log="xy", main=colnames(Peptide_Count)[i])}
    if(ncol(Peptide_Count) == ncol(LFQintensity)){plot(Peptide_Count[, i]+1, LFQintensity[, i]+1, log="xy", main=colnames(Peptide_Count)[i])}

  }
  par(opar)
  dev.off()

  return(list("spc"=Spectral_Count, "peptideCount"=Peptide_Count, "intensity"=LFQintensity, "normLFQintensity"=norm_data, "logNormLFQintensity"=logNormLFQintensity))
}


compute_SAINT <- function(hwd, folder, design, prohits, useProhits, baitAsPrey=FALSE, SAINT_path, outdir, dataName){
  
  if(0){
    hwd <- "C:/data/raw/EDYTA/PROTEIN"
    folder <- "220104_Nabeel_GZ"
    design <- file.path(hwd, folder, "combined", "experimentalDesignTemplate_S.txt")
    outdir <- "maxquant_processed"
    useProhits <- FALSE
    SAINT_path <- "C:/GREENBLATT/SAINT/SAINTexpress_v3.6.3__2018-03-09/Precompiled_binaries/Windows64"
    dataName <- "GROUP1"
  }

  gfp_from_prohits <- NULL
  if(useProhits){gfp_from_prohits <- read.delim(prohits, header=T)}

  setwd(file.path(hwd, folder, outdir))


  if(0){
    ## get MS run order for carry over inference, used for mist
    p <- file.path(hwd, folder)
    rawFiles <- list.files(path=p, pattern="raw")
    fileInfo <- file.info(paste(p, rawFiles, sep="/"))
    fileInfo <- fileInfo[order(fileInfo$mtime),]
    expOrder <- basename(rownames(fileInfo))
    expOrder <- gsub(".raw", "", expOrder)
    #expOrder <- setdiff(expOrder, "ORF10WU_GFP_1")
  }

  spc_file <- paste("Total_spectral_count_from_",folder,".tab",sep="")
  peptideCount_file <- paste("Peptide_count_from_",folder,".tab",sep="")
  intensity_file <- paste("TotalIntensity_normalized_LFQintensity_from_",folder,".tab",sep="")
  
  files <- c(spc_file, peptideCount_file, intensity_file)
  names(files) <- c("spc", "peptideCount", "intensity")

  #if(!dir.exists("SAINT")) dir.create("SAINT")
  #setwd("./SAINT")
  for(i in 1:length(files)){
    if(0) i <- 1
    atype <- names(files)[i]
    afile <- files[[i]]

    #prepare_for_MiST(adata, geneNameDesc, expOrder, design, dataName, baitAsPrey, type=atype)
    wd <- file.path(hwd, folder, outdir, paste("SAINT", dataName, sep="_"), atype)
    if(!dir.exists(wd)){
      dir.create(wd, recursive=TRUE)
    }
    
    run_saint(file.path(hwd, folder, outdir, afile), design, type=atype, SAINT_path, dataName, useProhits, gfp_from_prohits, outdir=wd)
  }

}

process_TimeSeries <- function(hwd, folder, logNormLFQintensity, viral_name, outdir){
  if(0){
    hwd <- "C:\\data\\raw\\EDYTA\\PROTEIN"
    folder <- "210617_SARS2_DRUGS_EXP1"
    viral_name <- "SARS2"
    outdir <- "maxquant_processed"
  }

  wd <- file.path(hwd, folder, outdir)
  if(!dir.exists(wd)){
    dir.create(wd)
  }
  setwd(wd)

  if(0){
    intensity_file <- file.path(wd, paste("TotalIntensity_normalized_log2LFQintensity_from_",folder,".tab",sep=""))
    logNormLFQintensity <- read.delim(intensity_file, sep="\t", header=TRUE)
    rownames(logNormLFQintensity) <- logNormLFQintensity[,1]
    logNormLFQintensity <- logNormLFQintensity[, -1]
  }


  print("procesing time series")
  #print(head(LFQintensity))
  groupDesign <- read.delim(file.path(hwd, folder, "combined", "experimentalDesignTemplate_TS.txt"))
  # print(groupDesign)
  # design format: Name	Fraction	Experiment	Time	Group	Cell
  plot_viral_protein(logNormLFQintensity, groupDesign, viral_name, folder)
}


do_stats <- function(hwd, folder, data_name, outdir, comparisons){
  if(0){
    hwd <- "C:/data/raw/EDYTA/PROTEIN"
    folder <- "220222_INF5_DIA"
    data_name <- "SARS2"
    outdir <- "maxquant_processed"
    comparisons <- "TREATMENT|ND,DOSE,GROUP|MOCK"
  }

  wd <- file.path(hwd, folder, outdir)
  if(!dir.exists(wd)){
    dir.create(wd)
  }
  setwd(wd)

  groupDesign <- read.delim(file.path(hwd, folder, "combined", "experimentalDesignTemplate_STAT.txt"))
  # print(groupDesign)
  # design format: Name	Fraction	Experiment	Time	Group	Cell

  ## plot group comparison
  path1 <- file.path(hwd, folder, "combined/txt")

  imputed_logLFQ <- plot_stats_treat_dose_group_wrProteo(path1, groupDesign, data_name, comparisons)

  ## process viral genes only
  if(!is.na(data_name)){
    intensity_file <- file.path(wd, paste("TotalIntensity_normalized_LFQintensity_from_",folder,".tab",sep=""))
    LFQintensity <- read.delim(intensity_file, sep="\t", header=TRUE)

    viral_prot <- LFQintensity[LFQintensity$Organism == data_name, "X"]

    imputed_viral <- imputed_logLFQ[viral_prot,]
    plot_stats_treat_dose_group(imputed_viral, groupDesign, data_name=paste(data_name, "_protein", sep=""), comparisons)

  }

}

process_modifications <- function(hwd, folder, outdir, modification){
  if(0){
    hwd <- "C:\\data\\raw\\EDYTA\\PROTEIN"
    folder <- "210617_SARS2_DRUGS_EXP1"
    outdir <- "maxquant_processed"
    modification <- "GlyGly (K)"
  }
  print(modification)
  wd <- file.path(hwd, folder, outdir, sep="/")
  if(!dir.exists(wd)){
    dir.create(wd)
  }
  setwd(wd)

  ## get signal intensity table: row=gene_position, col=experiment
  design <- read.delim(file.path(hwd, folder, "combined", "experimentalDesignTemplate.txt"))
  modification_file <- file.path(hwd, folder,"combined/txt", paste(modification, "Sites.txt", sep=""))
  modification_table <- read.delim(modification_file, header=T, stringsAsFactors = T)
  colnames(modification_table)
  modification_intensity <- extract_modifiction_intensity(modification_table, design, modification)

  dim(modification_intensity)
  summary(modification_intensity)
  #x <- na.omit(modification_intensity)
  #m.s = imputeLCMD::model.Selector(modification_intensity)
  imputed_intensity <- impute_intensity(modification_intensity, rm_empty=FALSE, method="LOD")
  while(!is.null(dev.list())){
    dev.off()
  }
  pdf(paste(modification, "heatmap.pdf", sep="_"), height=10, width=8)
  p <- pheatmap(imputed_intensity)
  print(p)
  dev.off()
  #pheatmap(modification_intensity)
  ## plot group comparison
  groupDesign <- read.delim(file.path(hwd, folder, "combined", "experimentalDesignTemplate_STAT.txt"))
  plot_stats_treat_dose_group(imputed_intensity, groupDesign, modification)
}


