#!/usr/bin/env Rscript
################     *_*     #######################
# By Guillermo Torres PhD                          #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
#############  About this script  ##################
# Major update: August 2018
# Created: 12 August 2019
#
# This is written as part of microbiome analysis,
# but could be splitted to serve different purposes.
#
# toolbox.R contains a set of functions
# commonly used by the bees' microbiome pipeline.
#==================================================#
########---- load/install packages ----#############
#==================================================#
#cat test.fa | sed -e '1!{s/^>.*/NNNNNNNNNNNNNNNNNNNN/g;}' | sed ':a;N;$!ba;s/\n//2g' | sed '1!s/.\{80\}/&\n/g' > beeG.fa
packages <- function(requirements,quiet=FALSE){
  has   <- requirements %in% rownames(installed.packages())
  if(any(!has)){
    message("Installing packages...")
    setRepositories(ind=c(1:7))
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.uni-muenster.de/"
    #options(install.packages.check.source = "no")
    install.packages(requirements[!has],repos=r)
  }
  if(quiet){
    for(r in requirements){suppressMessages(require(r,character.only=TRUE))}
  }else for(r in requirements){message(paste(r,suppressMessages(require(r,character.only=TRUE)),sep=': '))}
}
#==========================================================================================#
########---- Executer- execute with programs with system and check output ----#############
#========================================================================================#
executer <- function(checkout,exec,succ,fail,logfile=NULL,checkin=NULL,pass=FALSE){
  if(!is.null(checkin)){
    if(!file.exists(checkin) | (file.exists(checkin) & file.size(checkin)==0)){
      if(!is.null(logfile)){cat(paste0(checkin,' does not exist or is empty!, check and try again :)'),file=logfile,append=TRUE)}
      if(pass){warning(fail);return()
      }else{stop(paste0(checkin,' does not exist or is empty!, check and try again :)'))}
    }
  }
  if(file.exists(checkout) & file.size(checkout)!=0){message(succ)
  }else{
    sys <- try(system(exec,intern=TRUE))
    if(is.null(attributes(sys))){message(succ);return(sys)
    }else{
      if(!is.null(logfile)){cat(paste0(sys),file=logfile,append=TRUE)}
      if(pass){warning(fail);return()
      }else{stop(fail)}
    }
    
  }
}
#=============================================================================#
########---- Phyloseq merge_sample function fixed ----#############
#==============================================================================#
mergesamplesps_mean <- function(ps,group){
  library(phyloseq);library(tidyverse)
  mps <- merge_samples(ps,group)
  replicates <- sample_data(ps)$libname%>%sort()%>%table%>%as.data.frame#as_tibble()
  otu_table(mps) <- otu_table(mps)/replicates$Freq
  newmetdat <- sample_data(ps)%>%as_tibble()%>%distinct(libname, .keep_all = TRUE)%>%arrange(libname)%>%as.data.frame()
  rownames(newmetdat) <- rownames(otu_table(mps))
  sample_data(mps) <- newmetdat
  return(mps)
}
#=============================================================================#
########---- TaxOverview function ----#############
#==============================================================================#
taxoverview <- function(ps,tl,results,top10=F){
  #change the type, color and shape for the plot according to the needs.
  #tl <- 'Family'
  #otu_table(ps_f)%>%colSums()
  #ps <- ps_f
  ps.ord <- ordinate(ps, "NMDS", "bray")
  p1 <- plot_ordination(ps, ps.ord, type="samples", color="location", shape="type")+ geom_point(size=5) + ggtitle("Samples")
  ggsave(paste0(results,tl,'_NMDS_all_1.pdf'),p1,width=9, height=7)
  pdf(paste0(results,tl,'_BarPlot.pdf'),width=7, height=7,useDingbats=F)
  print(plot_bar(ps, fill = tl))
  dev.off()
  
  mps <- mergesamplesps_mean(ps,"libname")
  pdf(paste0(results,tl,'_BarPlot_merged.pdf'),width=7, height=7,useDingbats=F)
  print(plot_bar(mps, fill = tl))
  dev.off()
  mps.ord <- ordinate(mps, "NMDS", "bray")
  p3 = plot_ordination(mps, mps.ord, type="samples", color="location", shape="salinity")+ geom_point(size=5)+ggtitle(paste0("Ordinal plot based on ",tl))
  ggsave(paste0(results,tl,'_NMDS_merged.pdf'),p3,width=9, height=7)
  p4 = plot_ordination(mps, mps.ord, type="split", color=tl, shape="location",)+
    ggtitle(paste0("Ordinal split plot based on top 10 ",tl))+geom_point(size = 4)+geom_text(aes(label=libname),hjust=-0.3,vjust=1,size=4)
  ggsave(paste0(results,tl,'_NMDS_merged_split.pdf'),p4,width=9, height=7)
  if(top10){
    tax.sum = tapply(taxa_sums(ps), tax_table(ps)[, tl], sum, na.rm=TRUE)
    top10 = names(sort(tax.sum, TRUE))[1:10]
    psx = prune_taxa((tax_table(ps)[, tl] %in% top10), ps)
    
    mpst10 <- mergesamplesps_mean(psx,"libname")
    pdf(paste0(results,tl,'_BarPlot_Top10merged.pdf'),width=7, height=7,useDingbats=F)
    print(plot_bar(mpst10, fill = tl))
    dev.off()
    mpst10.ord <- ordinate(mpst10, "NMDS", "bray")
    p3 = plot_ordination(mpst10, mpst10.ord, type="samples", color="location", shape="salinity")+ geom_point(size=5) +  ggtitle(paste0("Ordinal plot based on top 10 ",tl))
    ggsave(paste0(results,tl,'_NMDS_merged.pdf'),p3,width=9, height=7)
    p4 = plot_ordination(mpst10, mpst10.ord, type="split", color=tl, shape="location")+
      ggtitle(paste0("Ordinal split plot based on top 10 ",tl))+geom_point(size = 4)+geom_text(aes(label=libname),hjust=-0.3,vjust=1,size=4)
    p4
    ggsave(paste0(results,tl,'_NMDS_merged_2.pdf'),p4,width=9, height=7)
  }
  
}
#=============================================================================#
########---- from metaplhlan tsv to metagenomeSeq and Phyloseq ----#############
#==============================================================================#
MRdata <- function(mtaxlv,samples_info){
  counts <- mtaxlv[,-c(1:3)]%>%as.data.frame()
  taxonomy <- mtaxlv[,c(1:3)]
  xnames = taxonomy$fullTax
  shortnames = gsub(".+\\|", "", xnames)
  x2 = strsplit(xnames, split="|", fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) =  taxonomy$fullTax
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxtable_df <- left_join(taxonomy,taxmat%>%as_tibble()%>%mutate(fullTax=rownames(taxmat)),by='fullTax')%>%dplyr::select(-1,-3)
  taxtable <- taxtable_df%>%dplyr::select(-1)%>%as.data.frame()
  rownames(taxtable) <- taxtable_df$tax;rownames(counts) <- taxtable_df$tax
  metadata <- samples_info%>%as.data.frame()
  rownames(metadata) <- metadata[,2]
  ord = match(colnames(counts),rownames(metadata))
  metadata = metadata[ord,]
  length(colnames(counts));length(union(colnames(counts),rownames(metadata)))
  data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxtable))
  return(data)
}
pseqGen <- function(mrdata,tree=F){
  getMetaphlanTree <- function(removeGCF=TRUE, simplify=TRUE){
    if (!requireNamespace("ape")) {
      stop("Please install the ape package to read Newick trees")
    }
    nwkfile <- bzfile(system.file("extdata/metaphlan2_selected.tree.reroot.nwk.bz2",
                                  package="curatedMetagenomicData"))
    tree <- ape::read.tree(nwkfile)
    close(nwkfile)
    if(removeGCF)
      tree$tip.label <- sub("\\|GCF_[0-9]+$", "", tree$tip.label)
    if(simplify)
      tree$tip.label <- gsub(".+\\|", "", tree$tip.label)
    return(tree)
  }
  taxmat <- phyloseq::tax_table(fData(mrdata)%>%as.matrix())
  otutab <- phyloseq::otu_table(MRcounts(mrdata), taxa_are_rows=TRUE)
  if(tree){
    tree <- getMetaphlanTree()
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(pData(mrdata)),tree)
  }else res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(pData(mrdata)))
  return(res)
}
##
metaphlanToPhyloseq <- function(metaphlandir,pat='taxonomic_profile.tsv',metadat=NULL,simplify=TRUE){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplify=TRUE, use only the most detailed level of taxa names in the final object
  ## metaphlanToPhyloseq("~/Downloads/metaphlan_bugs_list")
  .getMetaphlanTree <- function(removeGCF=TRUE, simplify=TRUE){
    if (!requireNamespace("ape")) {
      stop("Please install the ape package to read Newick trees")
    }
    nwkfile <- bzfile(system.file("extdata/metaphlan2_selected.tree.reroot.nwk.bz2",
                                  package="curatedMetagenomicData"))
    tree <- ape::read.tree(nwkfile)
    close(nwkfile)
    if(removeGCF) tree$tip.label <- sub("\\|GCF_[0-9]+$", "", tree$tip.label)
    if(simplify) tree$tip.label <- gsub(".+\\|", "", tree$tip.label)
    return(tree)
  }
  .joinListOfMatrices <- function(obj) {
    rnames <- Reduce(union, lapply(obj, rownames))
    cnames <- names(obj)
    if (!all(isUnique(cnames))) {
      stop("Column names are not unique.")
    }
    output <- matrix(0,
                     nrow = length(rnames),
                     ncol = length(cnames),
                     dimnames = list(rnames, cnames)
    )
    for (i in seq_along(obj)) {
      output[match(rownames(obj[[i]]), rownames(output)), i] <- obj[[i]][, 1]
    }
    return(output)
  }
  fnames <- list.files(metaphlandir,pattern=pat,recursive=T,full.names=T)
  bug.l <- lapply(fnames, function(x){
    res <- read.delim(file.path(metaphlandir, x), stringsAsFactors = FALSE, row.names = 1)
    colnames(res) <- x
    return(res)
  })
  names(bug.l) <- fnames
  tax <- .joinListOfMatrices(bug.l)
  tax <- mtax_raw%>%as.data.frame()
  xnames = rownames(tax)
  shortnames = gsub(".+\\|", "", xnames)
  if(simplify){
    rownames(tax) = shortnames
  }
  x2 = strsplit(xnames, split="|", fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  tree <- .getMetaphlanTree(simplify=simplify)
  if(is.null(metadat)){
    metadat <- data.frame(file=fnames, row.names=fnames, stringsAsFactors = FALSE)
  }
  res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat), tree)
  return(res)
}
##

#========================================================#
########---- get symlinks from seq libraries ----#############
#======================================================#
local_rename <- function(name,namesdf,n=2){
  if(!name%in%namesdf) {return(name)
  }else{local_rename(paste0(name,'_',n),namesdf,n=n+1)}
}
genSymLink <- function(from,to,pattern){
  kbepattern <- "^Kbe[[:alnum:]]*_[[:alpha:]]{1,4}$"
  from <- from #''
  to <- to #~/hamburg/BeesMicrobiome/rawdata/Beeslibs/'
  kbe <- dir(from,pattern=pattern,full.names=T,ignore.case=T,recursive=T,include.dirs=T)
  names <- rep(0,length(kbe))
  for(i in 1:length(kbe)){
    name <- tail(unlist(strsplit(kbe[i],'/')),1)
    name <- rename(name,names)
    names[i] <- name
    print(paste0(name,' <-- ',kbe[i]))
    system(paste0('ln -s ',kbe[i],' ',to,name))
  }
}

dircreations <- function(dir_path,overwriting){
  if(dir.exists(dir_path)){message(paste0(dir_path,' Folder already exist'))
    if(overwriting){message('overwriting...')}
  }else dir.create(dir_path)
}
#========================================================#
########---- get counts DFs by taxonomy ----#############
#======================================================#
get_counts_by_taxonomy <- function(infile,cpmth){
  packages(c('edgeR','tidyverse'),quiet = T)
  #infile <- paste0(opt$lib,j)
  #cpmth <- 5
  if(!file.exists(infile)) stop('in file do not exist. Please check!')
  taxdf <- suppressMessages(read_tsv(infile,col_names=c('Taxa','counts')))
  if(is.na(table(taxdf$counts==0)['TRUE'])) '' else {
    taxdf <- taxdf%>%filter_at(vars(-Taxa), any_vars(. != 0))
    message('ceros removed')
  }
  if(nrow(taxdf)==0) {message('No data in the file');return(NA)}
  
  tt <- suppressMessages(taxdf%>%separate(Taxa,into=c('Domain','Phyllum','Class','Order','Family','Genus','Specie'),sep='\\|'))
  D <- tt%>%filter(is.na(Phyllum))%>%discard(~all(is.na(.x)))%>%map_df(~.x)
  libsize=sum(D$counts)
  D <- D%>%mutate(log2_cpm=cpm(as.data.frame(D)[,2,drop=F],lib.size=libsize,log=TRUE))%>%filter(log2_cpm >= cpmth)%>%arrange(-counts)
  S <- tt%>%filter(!is.na(Specie))
  S <- S%>%mutate(log2_cpm=cpm(as.data.frame(S)[,'counts',drop=F],lib.size=libsize,log=TRUE))%>%filter(log2_cpm >= cpmth)
  G <- tt%>%filter(!is.na(Genus)&is.na(Specie))%>%discard(~all(is.na(.x)))%>%map_df(~.x)
  G <- G%>%mutate(log2_cpm=cpm(as.data.frame(G)[,'counts',drop=F],lib.size=libsize,log=TRUE))%>%filter(log2_cpm >= cpmth)
  f <- tt%>%filter(!is.na(Family)&is.na(Genus))%>%discard(~all(is.na(.x)))%>%map_df(~.x)
  f <- f%>%mutate(log2_cpm=cpm(as.data.frame(f)[,'counts',drop=F],lib.size=libsize,log=TRUE))%>%filter(log2_cpm >= cpmth)
  O <- tt%>%filter(!is.na(Order)&is.na(Family))%>%discard(~all(is.na(.x)))%>%map_df(~.x)
  O <- O%>%mutate(log2_cpm=cpm(as.data.frame(O)[,'counts',drop=F],lib.size=libsize,log=TRUE))%>%filter(log2_cpm >= cpmth)
  C <- tt%>%filter(!is.na(Class)&is.na(Order))%>%discard(~all(is.na(.x)))%>%map_df(~.x)
  C <- C%>%mutate(log2_cpm=cpm(as.data.frame(C)[,'counts',drop=F],lib.size=libsize,log=TRUE))%>%filter(log2_cpm >= cpmth)
  P <- tt%>%filter(!is.na(Phyllum)&is.na(Class))%>%discard(~all(is.na(.x)))%>%map_df(~.x)
  P <- P%>%mutate(log2_cpm=cpm(as.data.frame(P)[,'counts',drop=F],lib.size=libsize,log=TRUE))%>%filter(log2_cpm >= cpmth)
  
  dfs <- list('libsize'=libsize,'Domain'=D,'Phyllum'=P,'Class'=C,'Order'=O,'Family'=f,'Genus'=G,'Specie'=S)
  return(dfs)
}
#==========================================================#
########---- get sample counts by taxonomy ----#############
#==========================================================#
get_sample_CPM_taxonomy <- function(tempData,cpmth,outf){
  packages(c('edgeR','tidyverse'),quiet = T)
  #infile <- paste0(opt$lib,j)
  #cpmth <- 5
  #outf <- paste0(opt$out,name)
  total_libzise <- tempData$F_U$libsize+tempData$R_U$libsize+tempData$PE$libsize
  Domain <- full_join(tempData$F_U$Domain,tempData$R_U$Domain,by='Domain',suffix=c('_F','_R'))%>%full_join(tempData$PE$Domain,by='Domain')
  all_counts <- Domain%>%select(contains('counts'))%>%rowSums(na.rm=T)
  Domain <- Domain%>%mutate(all_counts=all_counts,all_log2_cpm=cpm(as.data.frame(all_counts),lib.size=total_libzise,log=TRUE))%>%filter(all_log2_cpm >= cpmth)
  write_tsv(Domain%>%mutate_if(is.numeric, round,digits=2),paste0(outf,'_Domain.txt'))
  
  Phyllum <- full_join(tempData$F_U$Phyllum,tempData$R_U$Phyllum%>%dplyr::select(-Domain),by='Phyllum',suffix=c('_F','_R'))%>%full_join(tempData$PE$Phyllum%>%dplyr::select(-Domain),by='Phyllum')
  all_counts <- Phyllum%>%select(contains('counts'))%>%rowSums(na.rm=T)
  Phyllum <- Phyllum%>%mutate(all_counts=all_counts,all_log2_cpm=cpm(as.data.frame(all_counts),lib.size=total_libzise,log=TRUE))%>%filter(all_log2_cpm >= cpmth)
  write_tsv(Phyllum%>%mutate_if(is.numeric, round,digits=2),paste0(outf,'_Phyllum.txt'))
  
  Class <- full_join(tempData$F_U$Class,tempData$R_U$Class%>%dplyr::select(-Domain,-Phyllum),by='Class',suffix=c('_F','_R'))%>%full_join(tempData$PE$Class%>%dplyr::select(-Domain,-Phyllum),by='Class')
  all_counts <- Class%>%select(contains('counts'))%>%rowSums(na.rm=T)
  Class <- Class%>%mutate(all_counts=all_counts,all_log2_cpm=cpm(as.data.frame(all_counts),lib.size=total_libzise,log=TRUE))%>%filter(all_log2_cpm >= cpmth)
  write_tsv(Class%>%mutate_if(is.numeric, round,digits=2),paste0(outf,'_Class.txt'))
  
  Order <- full_join(tempData$F_U$Order,tempData$R_U$Order%>%dplyr::select(-Domain,-Phyllum,-Class),by='Order',suffix=c('_F','_R'))%>%full_join(tempData$PE$Order%>%dplyr::select(-Domain,-Phyllum,-Class),by='Order')
  all_counts <- Order%>%select(contains('counts'))%>%rowSums(na.rm=T)
  Order <- Order%>%mutate(all_counts=all_counts,all_log2_cpm=cpm(as.data.frame(all_counts),lib.size=total_libzise,log=TRUE))%>%filter(all_log2_cpm >= cpmth)
  write_tsv(Order%>%mutate_if(is.numeric, round,digits=2),paste0(outf,'_Order.txt'))
  
  Family <- full_join(tempData$F_U$Family,tempData$R_U$Family%>%dplyr::select(-Domain,-Phyllum,-Class,-Order),by='Family',suffix=c('_F','_R'))%>%full_join(tempData$PE$Family%>%dplyr::select(-Domain,-Phyllum,-Class,-Order),by='Family')
  all_counts <- Family%>%select(contains('counts'))%>%rowSums(na.rm=T)
  Family <- Family%>%mutate(all_counts=all_counts,all_log2_cpm=cpm(as.data.frame(all_counts),lib.size=total_libzise,log=TRUE))%>%filter(all_log2_cpm >= cpmth)
  write_tsv(Family%>%mutate_if(is.numeric, round,digits=2),paste0(outf,'_Family.txt'))
  
  Genus <- full_join(tempData$F_U$Genus,tempData$R_U$Genus%>%dplyr::select(-Domain,-Phyllum,-Class,-Order,-Family),by='Genus',suffix=c('_F','_R'))%>%full_join(tempData$PE$Genus%>%dplyr::select(-Domain,-Phyllum,-Class,-Order,-Family),by='Genus')
  all_counts <- Genus%>%select(contains('counts'))%>%rowSums(na.rm=T)
  Genus <- Genus%>%mutate(all_counts=all_counts,all_log2_cpm=cpm(as.data.frame(all_counts),lib.size=total_libzise,log=TRUE))%>%filter(all_log2_cpm >= cpmth)
  write_tsv(Genus%>%mutate_if(is.numeric, round,digits=2),paste0(outf,'_Genus.txt'))
  
  Specie <- full_join(tempData$F_U$Specie,tempData$R_U$Specie%>%dplyr::select(-Domain,-Phyllum,-Class,-Order,-Family,-Genus),by='Specie',suffix=c('_F','_R'))%>%full_join(tempData$PE$Specie%>%dplyr::select(-Domain,-Phyllum,-Class,-Order,-Family,-Genus),by='Specie')
  all_counts <- Specie%>%select(contains('counts'))%>%rowSums(na.rm=T)
  Specie <- Specie%>%mutate(all_counts=all_counts,all_log2_cpm=cpm(as.data.frame(all_counts),lib.size=total_libzise,log=TRUE))%>%filter(all_log2_cpm >= cpmth)
  write_tsv(Specie%>%mutate_if(is.numeric, round,digits=2),paste0(outf,'_Specie.txt'))
  
}

symlinks <- function(){
  library(tidyverse)
  getfiles <- function(pat,searchdir,kbe){
    if(is.null(kbe)){
      kbe <- dir(searchdir,pattern=pat,full.names=T,ignore.case=T,recursive=T,include.dirs=T)%>%str_replace('//','/')
    }
    print(kbe)
    patok <- askYesNo(msg="Is this the expected result?")
    if(isFALSE(patok)){
      if(isFALSE(askYesNo(msg="Use the same pattern?"))){
        pat <- readline(prompt="Enter the regex pattern: ")
      }
      if(isFALSE(askYesNo(msg="Search in the same path?"))){
        searchdir <- readline(prompt="Enter path where to search: ")
      }
      kbe <- dir(searchdir,pattern=pat,full.names=T,ignore.case=T,recursive=T,include.dirs=T)%>%str_replace('//','/')
      getfiles(pat,searchdir,kbe)
    }
    return(kbe)
  }
  
  rename <- function(name,names,n=2){
    if(!name%in%names) {return(name)
    }else{rename(paste0(name,'_',n),names,n=n+1)}
  }
  
  #########################
  ####----Core init----####
  ########################
  message(interactive())
  if(interactive()){
    pat <- readline(prompt="Enter the regex pattern: ")
    message(paste0('',pat))
    searchdir<- readline(prompt="Enter path where to search: ")
    to<- readline(prompt="Enter path of the output folder: ")
    kbe <- getfiles(pat,searchdir,kbe=NULL)
    if(dir.exists(to)){message(paste0(dir_path,' Folder already exist'))}else dir.create(to)
    names <- rep(0,length(kbe))
    for(i in 1:length(kbe)){
      name <- tail(unlist(strsplit(kbe[i],'/')),1)
      name <- rename(name,names)
      names[i] <- name
      print(paste0(name,' <-- ',kbe[i]))
      system(paste0('ln -s ',kbe[i],' ',to,name))
    }
  }
}

symlinks2 <- function(){
  getfiles <- function(pat,searchdir,kbe){
    if(is.null(kbe)){
      kbe <- dir(searchdir,pattern=pat,full.names=T,ignore.case=T,recursive=T,include.dirs=T)%>%str_replace('//','/')
    }
    print(kbe)
    cat("Is this the expected result? yes/no ")
    patok <- scan("stdin",character(),n=1,quiet=T)
    if(str_detect(patok,"n|N|no|NO")){
      cat("Use the same pattern? yes/no ")
      P <- scan("stdin",character(),n=1,quiet=T)
      if(str_detect(P,"n|N|no|NO")){
        cat("Enter the regex pattern: ")
        pat <- scan("stdin",character(),n=1,quiet=T)
      }
      cat("Search in the same path? yes/no ")
      C <- scan("stdin",character(),n=1,quiet=T)
      if(str_detect(C,"n|N|no|NO")){
        cat("Enter path where to search: ")
        searchdir<- scan("stdin",character(),n=1,quiet=T)
      }
      kbe <- dir(searchdir,pattern=pat,full.names=T,ignore.case=T,recursive=T,include.dirs=T)%>%str_replace('//','/')
      getfiles(pat,searchdir,kbe)
    }
    return(kbe)
  }
  
  rename <- function(name,names,n=2){
    if(!name%in%names) {return(name)
    }else{rename(paste0(name,'_',n),names,n=n+1)}
  }
  
  #########################
  ####----Core init----####
  ########################
  cat("Enter the regex pattern: ")
  pat <- scan("stdin",character(),n=1,quiet=T)
  message(pat)
  cat("Enter path where to search: ")
  searchdir<- scan("stdin",character(),n=1,quiet=T)
  message(searchdir)
  cat("Enter path of the output folder: ")
  to<- searchdir<- scan("stdin",character(),n=1,quiet=T)
  message(searchdir)
  kbe <- getfiles(pat,searchdir,kbe=NULL)
  if(dir.exists(to)){message(paste0(dir_path,' Folder already exist'))}else dir.create(to)
  names <- rep(0,length(kbe))
  for(i in 1:length(kbe)){
    name <- tail(unlist(strsplit(kbe[i],'/')),1)
    name <- rename(name,names)
    names[i] <- name
    print(paste0(name,' <-- ',kbe[i]))
    system(paste0('ln -s ',kbe[i],' ',to,name))
  }
  
  
}

