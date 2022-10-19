#! /usr/bin/Rscript --vanilla

cat(">> LOADING DEPENDENCIES...\n")
suppressMessages(library(getopt))
suppressMessages(library(yaml))

#########################
## MAIN FILE #######

spec = matrix(c(
  'verbose', 'v', 2, "integer", "",
  'help'   , 'h', 0, "logical", "available arguments (this screen)",
  'config'  , 'c', 1, "character", "configuration file in YAML format"),
  byrow=TRUE, ncol=5)

opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

PIPELINE=T

# set source directory
args <- commandArgs(trailingOnly = F)
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))  # get the path to where main.R is located

## load all externeal files
source(paste(scriptPath,"/src/extract_maxquant_lib.R",sep=""))
source(paste(scriptPath,"/src/extract_maxquant.R",sep=""))


# check if libraries are installed
checkForLibraries <- function(required_libraries){
  pkgs = data.frame(installed.packages(), stringsAsFactors=F)
  #if the packages haven't been installed on the computer yet
  if(!all(required_libraries %in% pkgs$Package)){
    missing_packages = required_libraries[ !required_libraries %in% pkgs$Package ]
    stop( paste("The following packages need to be installed on your computer in order to run MIST:\n"), paste('\t',missing_packages,'\n', sep=""), call.=FALSE )

  }else{
    cat(">> All required R libraries accounted for.\n")
  }
}

# check for all packages reqiured
required_libraries = c('yaml','getopt','optparse','PerformanceAnalytics','plotrix','limma','data.table','stats','grDevices','graphics','pheatmap','RColorBrewer','ggplot2','gridExtra','gplots')
checkForLibraries(required_libraries)

# some qc to make sure config file exists and is well formatted
getConfig <- function(config_file){
  if( !file.exists( config_file )){
    stop( cat(paste("ERROR!!! The yml configuration file:\n\t",paste(getwd(),config_file,sep='/'), '\ndoes not exist. Please enter a correct path/filename.\n', sep='')), call.=F)
  }

  x = readLines(config_file)
  x = x[x!=""]  #remove \n\n cases (blank Lines)
  x = gsub(':  ',': ', gsub(":", ': ',x) )   # make sure there is a space between ':' and any character
  x = gsub('\t', '  ', x)
  config = paste(x, collapse='\n')


  config = tryCatch(yaml.load(config), error = function(e) { print("!!! Error loading the config file. Please make sure the file follows YAML format."); stop()} )

  # check formatting of config file
  config = formatConfig(config)
  return(config)
}

# chec out YAML config file and make sure it's setup properly
formatConfig = function(config){
  # add : to C/.
  config$SAINT$excutable_path <- gsub('C/','C:/',config$SAINT$excutable_path)
  config$Files$data_dir <- gsub('C/','C:/',config$Files$data_dir)

  return(config)
}


main <- function(opt){
  # read and qc config file
  config = getConfig(opt$config)
  print(config)
  ##  create an outputdir if it doesn't exist
  if(is.null(config$Extract$out_dir) || config$Extract$out_dir == '') config$Extract$out_dir <- "maxquant_processed"
  outpath <- file.path(config$Files$data_dir, config$Files$folder, config$Extract$out_dir)


  if(!dir.exists(outpath)){
    dir.create(outpath, showWarnings = T)
  }


  ## main switches between parts of the pipeline

  cat(">> COLLECTING PROTEIN INFORMATION\n")
  geneNameDesc = collect_protein_info(aaInfo=config$Files$aaInfo, fasta=config$Files$fasta, organism=config$Files$organism, viral=config$Virus$enabled, viral_fasta=config$Virus$viral_fasta, viral_name=config$Virus$viral_name, id_map=config$Files$id_map)
  print(head(geneNameDesc))
  if(config$Extract$enabled){

    data_extracted <- run_extract(hwd=config$Files$data_dir, folder=config$Files$folder, geneNameDesc, isLFQ=config$Extract$isLFQ, outdir=config$Extract$out_dir, TS=config$TimeSeries$enabled, service=config$Files$service)
    LFQintensity <- data_extracted[["logNormLFQintensity"]]
  }else{ ## use previous data matrix instead of the one from pre-processing call
    intensity_file <- file.path(outpath, paste("TotalIntensity_normalized_log2LFQintensity_from_",config$Files$folder,".tsv",sep=""))
    LFQintensity <- read.delim(intensity_file, sep="\t", header=TRUE)
    rownames(LFQintensity) <- LFQintensity[,1]
    LFQintensity <- LFQintensity[, -1]
  }

  if(config$SAINT$enabled){
    cat(">> SAINT SCORING\n")

    if(is.null(config$SAINT$data_name)){
      stop("data_name is not defined in the yml file!")
    }

    for(data_name in unlist(strsplit(config$SAINT$data_name, split=",", fixed=T))){
            compute_SAINT(hwd=config$Files$data_dir,
                    folder=config$Files$folder,
                    design=file.path(config$Files$data_dir, config$Files$folder, "combined", "experimentalDesignTemplate_S.txt"),
                    #geneNameDesc,
                    #data_list,
                    prohits=config$SAINT$prohits,
                    useProhits=config$SAINT$useProhits,
                    baitAsPrey=config$SAINT$baitAsPrey,
                    SAINT_path=config$SAINT$excutable_path,
                    outdir=config$Extract$out_dir,
                    dataName=data_name)
    }
  }


  if(config$TimeSeries$enabled){
    cat(">> PLOT TIME SERIES\n")

    process_TimeSeries(config$Files$data_dir, folder=config$Files$folder, LFQintensity, viral_name=config$Virus$viral_name, outdir=config$Extract$out_dir)
  }
  if(!is.null(config$STATS)){
    if(config$STATS$enabled){
      cat(">> PERFORMING STAT ANALYSIS\n")
      do_stats(config$Files$data_dir, folder=config$Files$folder, data_name=config$Virus$viral_name, outdir=config$Extract$out_dir, comparisons=config$STATS$comparison)
    }
  }
  if(!is.null(config$Modifications)){
    if(config$Modifications$enabled){
      cat(">> PROCESSING MODIFICATIONS\n")
      modifications <- unlist(strsplit(config$Modifications$types, split=","))
      for(modification in modifications){
        process_modifications(config$Files$data_dir, folder=config$Files$folder, outdir=config$Extract$out_dir, modification)
      }
    }
  }


  cat(">> MAXQUANT POSTPROCESSING FINISHED!\n")
}

if(!exists('DEBUG') || DEBUG==F) main(opt)




#setwd("C:\\GREENBLATT\\Rscripts\\Maxquant_postprocess")
#system("Rscript maxquant_analyzer.R -c yml/220720_MaryChang.yml")
#system("Rscript maxquant_analyzer.R -c yml/220520_Meena_SG.yml")
#system("Rscript maxquant_analyzer.R -c yml/220524_Shamira_TH.yml")
#system("Rscript maxquant_analyzer.R -c yml/220525_Sabrina_MM.yml")
#system("Rscript maxquant_analyzer.R -c yml/220527_MaryChang1.yml")
#system("Rscript maxquant_analyzer.R -c yml/220527_MaryChang2.yml")

#system("Rscript maxquant_analyzer.R -c yml/220328_HZ_Bonin.yml")

#system("Rscript maxquant_analyzer.R -c yml/220316_Shamira_TH.yml")
#system("Rscript maxquant_analyzer.R -c yml/220222_INF5_DIA.yml")

#system("Rscript maxquant_analyzer.R -c yml/220104_Nabeel_GZ.yml")
#system("Rscript maxquant_analyzer.R -c yml/220104_Umamam_LF.yml")
#system("Rscript maxquant_analyzer.R -c yml/211221_INF5_2.yml")
#system("Rscript maxquant_analyzer.R -c yml/211206_SHAMIRA_TH.yml")
#system("Rscript maxquant_analyzer.R -c yml/211206_inf5_SARS2.yml")
#system("Rscript maxquant_analyzer.R -c yml/211115_BABAK_JM.yml")
#system("Rscript maxquant_analyzer.R -c yml/210324_Nabeel.yml")
#system("Rscript maxquant_analyzer.R -c yml/210629_Nabeel_GZ.yml")
#system("Rscript maxquant_analyzer.R -c yml/211025_Babak_JM.yml")
#system("Rscript maxquant_analyzer.R -c yml/211025_umama_LF.yml")
#system("Rscript maxquant_analyzer.R -c yml/210929_SARS2_INF4.yml")
#system("Rscript maxquant_analyzer.R -c yml/210820_Zuyao.yml")
#system("Rscript maxquant_analyzer.R -c yml/210811_Nishanth_CM.yml")
#system("Rscript maxquant_analyzer.R -c yml/210803_Calu3_drugs_Inf3.yml")
#system("Rscript maxquant_analyzer.R -c yml/210727_Calu3_drugs_INF2_regular.yml")
#system("Rscript maxquant_analyzer.R -c yml/210727_Calu3_drugs_INF2_targetted.yml")


#system("Rscript maxquant_analyzer.R -c yml/210722_Babak_JM_BIoID.yml")

#system("Rscript maxquant_analyzer.R -c yml/210621_SERVICE_SG_MED.yml")
#system("Rscript maxquant_analyzer.R -c yml/210617_SARS2_DRUGS_EXP1.yml")
#system("Rscript maxquant_analyzer.R -c yml/201027_sars2_infected.yml")
#system("Rscript maxquant_analyzer.R -c yml/210212_SaRS2_Caco2_full.yml")
#system("Rscript maxquant_analyzer.R -c yml/201102_sars2_infected_calu3_caco2.yml")
#system("Rscript maxquant_analyzer.R -c yml/200825_Vero.yml")
#system("Rscript maxquant_analyzer.R -c yml/200827_caco2_huh7.yml")
#system("Rscript maxquant_analyzer.R -c yml/210510_AB_TH.yml")
