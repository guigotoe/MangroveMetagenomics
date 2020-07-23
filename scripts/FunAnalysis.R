argv <- commandArgs(trailingOnly = FALSE)
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
toolbox <- paste0(mypath,'/toolbox.R')
source(toolbox)

requirements <- c('tidyverse','tictoc','phyloseq','vegan','metagenomeSeq','Rlof',"RColorBrewer","RAM",'curatedMetagenomicData','ape','matrixStats','DESeq2','pheatmap','Maaslin2')
packages(requirements)
theme_set(theme_bw())
##


results <- 'results/Functional/'
if(!dir.exists(results)) dir.create(results,showWarnings=F)
if(!file.exists(paste0(results,'merged_abundance_table.txt'))){
  ##copy from server
  filepath <- '~/hamburg/mangroveCol/functprofile/'
  pat <- 'ko_relab.tsv'
  fun <- c('ko','ecs','eggnog','pfam')
  pats <- c(paste0(fun,'_relab.tsv'),paste0(fun,'_cpm.tsv'),'pathabundance.tsv')
  for pat in pats{
    humannfiles <- list.files(filepath,pattern=pat,recursive = T,full.names = T)
    ffiles <- list()
    for (i in 1:length(humannfiles)){
      ffiles[[i]] <- read_tsv(humannfiles[i])
    }
    fdf <- plyr::join_all(ffiles,by='# Gene Family',type='full')
    write_tsv(fdf,paste0(results,'merged_',pat))
  }
  
}
raw_metadata <- read_tsv(paste0('docs/metadata.txt'))
raw_metadata$ID
#raw_metadata <- raw_metadata%>%mutate(libname=map(raw_metadata$original.ID, function(x) x%>%str_split('\\.')%>%unlist%>%.[1])%>%unlist)%>%transmute(libname,!!!.)

mtax_raw <- read.table(paste0(results,'/merged_abundance_table.txt'),quote='',sep="\t",skip='#',header=T,stringsAsFactors=F)%>%as_tibble()
mtax <- mtax_raw
colnames(mtax) <- c('fullTax',colnames(mtax)[-1]%>%str_replace('_taxonomic_profile','')%>%str_replace('X',''))
taxlev <- map(mtax$fullTax,function(x) x%>%str_split('\\|')%>%unlist%>%length)%>%unlist #1=k,2=p,3=c,4=o,5=f,6=g,7=s,8=t
taxid <- map(mtax$fullTax,function(x) x%>%str_split('\\|')%>%unlist%>%rev(.)%>%.[1])%>%unlist
mtax <- mtax%>%mutate(taxlevel=taxlev,tax=taxid)%>%transmute(taxlevel,tax, !!!.)
libs <- colnames(mtax)[-c(1:3)]%>%map(.,function(x) x%>%str_split("_|\\.")%>%unlist%>%.[1])%>%unlist
samples <- tibble(libname=libs,sample=colnames(mtax)[-c(1:3)])
samples_info <- samples%>%left_join(raw_metadata,by=c('libname'='ID'))%>%
  mutate(libname=factor(libname,levels=c("MmtII1","MmtII2","MmtII3","3A","3B","3C","2A","2B","2C","4A","4B","4C",
                                         "1","4","8","10","13","16","17","18","19","20","21","22","23","24","25","26","27","28")))#%>%dplyr::select(-original.ID,-databank_used)
samples_info <- as.data.frame(lapply(samples_info,function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
