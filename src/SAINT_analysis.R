scriptPath <- "C:\\GREENBLATT\\Rscripts\\Maxquant_postprocess"

source(paste(scriptPath,"/src/extract_maxquant_lib.R",sep=""))
source(paste(scriptPath,"/src/extract_maxquant.R",sep=""))
if(0){
   saintDirs <- c("C:/data/raw/EDYTA/PROTEIN/210629_Nabeel_GZ/maxquant_processed/SAINT_210629_Nabeel_Dorothy",
                  "C:/data/raw/EDYTA/PROTEIN/210324_Nabeel/maxquant_processed/SAINT_210324_Nabeel_Dorothy")
   names(saintDirs) <- c("210629_Nabeel", "210324_Nabeel")
   outDir <- "C:/data/raw/EDYTA/PROTEIN/210629_Nabeel_GZ/maxquant_processed/SAINT_210324_210629_Nabeel_Dorothy"
   
   combine_inputs_for_SAINT(saintDirs, outDir)
}

if(0){
   inDirs <- c("C:/data/raw/EDYTA/PROTEIN/211115_BABAK_JM/maxquant_processed",
               "C:/data/raw/EDYTA/PROTEIN/211025_Babak_JM/maxquant_processed")
   names(inDirs) <- c("211115_BABAK_JM", "211025_Babak_JM")
   outDir <- "C:/data/raw/EDYTA/PROTEIN/211115_BABAK_JM/maxquant_processed/combined"
   
   combine_processed_results(inDirs, outDir)
   
   SAINT_path <- "C:/GREENBLATT/SAINT/SAINTexpress_v3.6.3__2018-03-09/Precompiled_binaries/Windows64"
   designfile <- "C:/data/raw/EDYTA/PROTEIN/211115_BABAK_JM/maxquant_processed/combined/experimentalDesignTemplate_S.txt"
   datafile <- "Total_spectral_count_from_211115_BABAK_JM_211025_Babak_JM_combined.tab"
   outdir <- "C:/data/raw/EDYTA/PROTEIN/211115_BABAK_JM/maxquant_processed/combined"
   for(dataName in c("Control", "Treatment")){
      run_saint(datafile, designfile, type="spc", SAINT_path, dataName, useProhits, gfp_from_prohits, outdir)
   }
   
   
}
