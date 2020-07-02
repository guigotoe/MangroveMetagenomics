argv <- commandArgs(trailingOnly = FALSE)
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
toolbox <- paste0(mypath,'/toolbox.R')
source(toolbox)

requirements <- c('tidyverse','tictoc','phyloseq','vegan','metagenomeSeq','Rlof',"RColorBrewer","RAM",'curatedMetagenomicData','ape','matrixStats','DESeq2','pheatmap')
packages(requirements)
theme_set(theme_bw())
##


results <- 'results/taxonomy/'
if(!dir.exists(results)) dir.create(results,showWarnings=F)
if(!file.exists(paste0(results,'merged_abundance_table.txt'))){
  ##copy from server
  filepath <- '~/hamburg/mangroveCol/taxbin/'
  mtphlanfiles <- list.files(filepath,pattern='taxonomic_profile.tsv',recursive = T,full.names = T)
  
  mergemtphlan <- paste0('merge_metaphlan_tables.py ',paste0(mtphlanfiles,collapse=' '),' > ',results,'merged_abundance_table.txt')
  system(mergemtphlan)
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

mtax_k<- mtax%>%filter(taxlevel==1);data_k <- MRdata(mtax_k,samples_info); ps_k <- pseqGen(data_k) 
mtax_p<- mtax%>%filter(taxlevel==2);data_p <- MRdata(mtax_p,samples_info); ps_p <- pseqGen(data_p) 
mtax_f <- mtax%>%filter(taxlevel==5);data_f <- MRdata(mtax_f,samples_info); ps_f <- pseqGen(data_f)
mtax_g <- mtax%>%filter(taxlevel==6);data_g <- MRdata(mtax_g,samples_info); ps_g <- pseqGen(data_g) 
mtax_s <- mtax%>%filter(taxlevel==7);data_s <- MRdata(mtax_s,samples_info); ps_s <- pseqGen(data_s,tree=T) 

## Kingdom anaylsis
tl <- 'Kingdom'
ps <- ps_k
ps.ord <- ordinate(ps, "NMDS", "bray")
p1 <- plot_ordination(ps, ps.ord, type="samples", color="location", shape="type")+ geom_point(size=5) + ggtitle("Samples")
p1
ggsave(paste0(results,tl,'_NMDS_all_1.pdf'),p1,width=9, height=7)
otu_table(ps)
pdf(paste0(results,tl,'_BarPlot.pdf'),width=7, height=7,useDingbats=F)
plot_bar(ps, fill = tl)
dev.off()
mps <- mergesamplesps_mean(ps,"libname")
pdf(paste0(results,tl,'_BarPlot_merged.pdf'),width=7, height=7,useDingbats=F)
plot_bar(mps, fill = tl)
dev.off()
mps.ord <- ordinate(mps, "NMDS", "bray")
p3 = plot_ordination(mps, mps.ord, type="samples", color="location", shape="type")+ geom_point(size=5) + ggtitle(paste0("Ordinal plot based on ",tl))
ggsave(paste0(results,tl,'_NMDS_merged.pdf'),p3,width=9, height=7)
p4 = plot_ordination(mps, mps.ord, type="split", color=tl, shape="location", label="libname", title="split")
ggsave(paste0(results,tl,'_NMDS_merged_split.pdf'),p4,width=9, height=7)
##
## phyllum anaylsis
tl <- 'Phylum'
ps <- ps_p
ps.ord <- ordinate(ps, "NMDS", "bray")
p1 <- plot_ordination(ps, ps.ord, type="samples", color="location", shape="type")+ geom_point(size=5) + ggtitle("Samples")
p1
ggsave(paste0(results,tl,'_NMDS_all_1.pdf'),p1,width=9, height=7)
mps <- mergesamplesps_mean(ps,"libname")
pdf(paste0(results,tl,'_BarPlot_merged.pdf'),width=7, height=7,useDingbats=F)
plot_bar(mps, fill = tl)
dev.off()
mps.ord <- ordinate(mps, "NMDS", "bray")
p3 = plot_ordination(mps, mps.ord, type="samples", color="location", shape="salinity")+ geom_point(size=5) + ggtitle(paste0("Ordinal plot based on ",tl))
ggsave(paste0(results,tl,'_NMDS_merged.pdf'),p3,width=9, height=7)
p4 = plot_ordination(mps, mps.ord, type="split", color=tl, shape="location", label="libname", title="split")
ggsave(paste0(results,tl,'_NMDS_merged_split.pdf'),p4,width=9, height=7)
##
## family anaylsis

tl <- 'Family'
otu_table(ps_f)%>%colSums()
ps <- ps_f
ps.ord <- ordinate(ps, "NMDS", "bray")
p1 <- plot_ordination(ps, ps.ord, type="samples", color="location", shape="type")+ geom_point(size=5) + ggtitle("Samples")
p1
ggsave(paste0(results,tl,'_NMDS_all_1.pdf'),p1,width=9, height=7)

mps <- mergesamplesps_mean(ps,"libname")
pdf(paste0(results,tl,'_BarPlot_merged.pdf'),width=7, height=7,useDingbats=F)
plot_bar(mps, fill = tl)
dev.off()
mps.ord <- ordinate(mps, "NMDS", "bray")
p3 = plot_ordination(mps, mps.ord, type="samples", color="location", shape="salinity")+ geom_point(size=5) + ggtitle(paste0("Ordinal plot based on ",tl))
ggsave(paste0(results,tl,'_NMDS_merged.pdf'),p3,width=9, height=7)
p4 = plot_ordination(mps, mps.ord, type="split", color=tl, shape="location", label="libname", title="split")+ geom_point(size = 4)
ggsave(paste0(results,tl,'_NMDS_merged_split.pdf'),p4,width=9, height=7)

tax.sum = tapply(taxa_sums(ps), tax_table(ps)[, tl], sum, na.rm=TRUE)
top10 = names(sort(tax.sum, TRUE))[1:10]
psx = prune_taxa((tax_table(ps_f)[, tl] %in% top10), ps_f)

mps <- merge_samples(psx,"libname",fun="mean")
replicates <- sample_data(ps_p)$libname%>%as.vector()%>%table%>%as.tibble()
otu_table(mps) <- otu_table(mps)/replicates$n
newmetdat <- sample_data(psx)%>%as_tibble()%>%distinct(libname, .keep_all = TRUE)%>%dplyr::select(-sample)%>%as.data.frame()
rownames(newmetdat) <- rownames(otu_table(mps))
sample_data(mps) <- newmetdat
pdf(paste0(results,tl,'_BarPlot_merged.pdf'),width=7, height=7,useDingbats=F)
plot_bar(mps, fill = tl)
dev.off()
mps.ord <- ordinate(mps, "NMDS", "bray")
p3 = plot_ordination(mps, mps.ord, type="samples", color="location", shape="salinity")+ geom_point(size=5) + ggtitle("Ordinal plot based on Phylum")
ggsave(paste0(results,tl,'_NMDS_merged.pdf'),p3,width=9, height=7)
p4 = plot_ordination(mps, mps.ord, type="samples", color="location", shape="type")+ geom_point(size=5) + ggtitle("Ordinal plot based on Phylum")
ggsave(paste0(results,tl,'_NMDS_merged_2.pdf'),p4,width=9, height=7)
##
## genus anaylsis
tl <- 'Genus'
ps_f.ord <- ordinate(ps_f, "NMDS", "bray")
p1 <- plot_ordination(ps_f, ps_f.ord, type="samples", color="material", shape="archaeological_site")+ geom_point(size=5) + ggtitle("Samples")
ggsave(paste0('taxbin/chew/',tl,'_NMDS_all.pdf'),p1,width=9, height=7)

tax.sum = tapply(taxa_sums(ps_f), tax_table(ps_f)[, tl], sum, na.rm=TRUE)
top10 = names(sort(tax.sum, TRUE))[1:10]
psx = prune_taxa((tax_table(ps_f)[, tl] %in% top10), ps_f)

mps <- merge_samples(psx,"libname",fun="mean")
replicates <- sample_data(psx)%>%group_by(libname)%>%count
otu_table(mps) <- otu_table(mps)/replicates$n
newmetdat <- sample_data(psx)%>%as_tibble()%>%distinct(libname, .keep_all = TRUE)%>%dplyr::select(-sample,-short.ID)%>%as.data.frame()
rownames(newmetdat) <- rownames(otu_table(mps))
sample_data(mps) <- newmetdat
pdf(paste0('taxbin/chew/',tl,'_BarPlot_merged.pdf'),width=7, height=7,useDingbats=F)
plot_bar(mps, fill = tl)
dev.off()

tl <- 'Genus'
ps <- ps_g
colSums(otu_table(ps))
ps  = transform_sample_counts(ps, function(x) x+1)
counts <- as.data.frame(lapply(otu_table(ps)%>%as.data.frame,function (y) if(class(y)!="integer" ) as.integer(y) else y),stringsAsFactors=F)
rownames(counts) <- rownames(otu_table(ps))
otu_table(ps) <- otu_table(counts,taxa_are_rows=T)
pShannon <- plot_richness(ps, x = "material", color = "archaeological_site",measures='Shannon') + geom_boxplot()+
  theme(axis.text.x = element_text(angle=0, colour = "black"))#+geom_point(position=position_jitterdodge())
ggsave(paste0('taxbin/chew/',tl,'_shannon_all.pdf'),pShannon,width=9, height=7)

##
## Species anaylsis
tl <- 'Species'
ps <- ps_s
colSums(otu_table(ps))
ps  = transform_sample_counts(ps, function(x) x+1)
counts <- as.data.frame(lapply(otu_table(ps)%>%as.data.frame,function (y) if(class(y)!="integer" ) as.integer(y) else y),stringsAsFactors=F)
rownames(counts) <- rownames(otu_table(ps))
otu_table(ps) <- otu_table(counts,taxa_are_rows=T)
pShannon <- plot_richness(ps, x = "material", color = "archaeological_site",measures='Shannon') + geom_boxplot()#+geom_point(position=position_jitterdodge())
ggsave(paste0('taxbin/chew/',tl,'_shannon_all.pdf'),pShannon,width=9, height=7)

#Remove OTUs that do not show appear more than 5 times in more than half the samples
#wh0 = genefilter_sample(ps, filterfun_sample(function(x) x > 5), A=0.023*nsamples(ps))
#ps = prune_taxa(wh0, ps)

###
## Complex plots
##
tl='Genus'
pc <- 0.3
ps <- ps_g
ps <- transform_sample_counts(ps, function(x) x+1)
#methods(class='phyloseq')
sample_data(ps)%>%head

## dufferentially abundant
## dc vs soil
dcsoil = subset_samples(ps, material != "chew")
dcsoil_deseq = phyloseq_to_deseq2(dcsoil, ~ material)
dcsoil_deseq = DESeq2::DESeq(dcsoil_deseq, test="Wald", fitType="parametric")
res = DESeq2::results(dcsoil_deseq, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
## dc vs chew
dcchew = subset_samples(ps, material != "soil")
dcchew_deseq = phyloseq_to_deseq2(dcchew, ~ material)
dcchew_deseq = DESeq2::DESeq(dcchew_deseq, test="Wald", fitType="parametric")
res2 = DESeq2::results(dcchew_deseq, cooksCutoff = FALSE)
alpha = 0.01
sigtab2 = res2[which(res2$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(ps)[rownames(sigtab2), ], "matrix"))
head(sigtab2)
## dc vs soil
chewsoil = subset_samples(ps, material != "dental_calculus")
chewsoil_deseq = phyloseq_to_deseq2(chewsoil, ~ material)
chewsoil_deseq = DESeq2::DESeq(chewsoil_deseq, test="Wald", fitType="parametric")
res3 = DESeq2::results(chewsoil_deseq, cooksCutoff = FALSE)
alpha = 0.01
sigtab3 = res3[which(res3$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(ps)[rownames(sigtab3), ], "matrix"))
head(sigtab3)
####
difftaxa <- c(rownames(sigtab),rownames(sigtab2),rownames(sigtab3))%>%unique
phy_stats <- tibble(!!tl:=taxa_names(ps)%>%str_replace('[a-z]__',''),tax=taxa_names(ps),N=nsamples(ps),PercentPhy=rowMeans(otu_table(ps)),sd=matrixStats::rowSds(otu_table(ps)))%>%mutate(se= sd/sqrt(N))
abund <- subset(phy_stats,tax%in%difftaxa)  # Only take the phyla that are more than 0.01% abundant
abund_ordered <- arrange(abund, desc(PercentPhy))
abund <- abund%>%mutate(!!tl:=factor(abund%>%pull(1),levels=abund_ordered%>%pull(1)%>%rev  ) )
#Relative abundance plot 
abund_plot <- ggplot(abund, aes(y=PercentPhy , x=get(tl) ))  +
  #geom_boxplot(fill = "magenta4", colour = "black") + 
  geom_bar(stat="identity", position=position_dodge(),  fill = "magenta4", colour = "black") +
  theme_bw() + ggtitle(paste0(tl," differentially abundant")) +
  xlab(paste0(tl)) + ylab("Mean Relative Abundance (%)") +
  geom_errorbar(aes(ymin = PercentPhy -se, ymax = PercentPhy +se), width = 0.25) + coord_flip() +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=0, colour = "black", size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        legend.position="none"); 
abund_plot

relabun_plot <- ggplot(abund, aes(y=PercentPhy , x=get(tl) )) + #coord_cartesian(xlim = c(0, 30)) + 
  geom_bar(stat="identity", position=position_dodge(),  fill = "magenta4", colour = "black") +
  ylab("Mean Percent \n Relative \n Abundance (%)") + coord_flip() + theme_bw() + 
  geom_errorbar(aes(ymin = PercentPhy -se, ymax = PercentPhy +se), width = 0.25, color = "black") + 
  scale_y_continuous(expand= c(0,0), limits = c(0,27)) +
  theme(axis.title.x = element_text(face="bold", size=8),
        axis.text.x = element_text(angle=0, colour = "black", size=8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(colour="black", fill = "black"),
        plot.margin = unit(c(0.5, 0.3, 0.4, -0.), "cm"), #top, right, bottom, left
        #plot.margin = unit(c(1.65, 2, 1.35, -1.35), "cm"), #top, right, bottom, left
        legend.position="none"); #relabun_plot

## Ploting ##
mps <- merge_samples(ps,"material",fun="mean")
replicates <- sample_data(ps)%>%group_by(material)%>%count
otu_table(mps) <- otu_table(mps)/replicates$n
newmetdat <- sample_data(ps)%>%as_tibble()%>%distinct(material, .keep_all = TRUE)%>%dplyr::select(-sample,-libname,-short.ID)%>%as.data.frame()
rownames(newmetdat) <- rownames(otu_table(mps))
sample_data(mps) <- newmetdat

dfotus <- otu_table(mps)%>%t()%>%as.data.frame()%>%mutate(tax=colnames(otu_table(mps)))%>%as_tibble()%>%filter(tax%in%difftaxa)%>%gather(-tax,key='material',value='abundance')%>%mutate(Sqrt.abundance=sqrt(abundance),log.abundance=log(abundance+1),taxO=factor(dfotus%>%pull(tax)%>%str_replace('[a-z]__',''),levels=abund_ordered%>%pull(tl)%>%rev))

heat <- ggplot(dfotus, mapping = aes(x=material,y=taxO,fill=Sqrt.abundance))+geom_tile()+
  scale_fill_gradient(name = "Sqrt(Abundance)",low = "#FFFFFF",high = "#012345") +
  theme_bw() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  ggtitle(label = paste0(tl," Class abundance")) +
  scale_y_discrete(limits = rev(levels(dfotus%>%pull(taxO))))


heat <- ggplot(dfotus, mapping = aes(x=material,y=taxO,fill=Sqrt.abundance))+geom_tile() + 
  scale_fill_gradient(name = "Sqrt(Abundance)",low = "#FFFFFF",high = "#012345")  + 
  theme_bw(base_size = 12) + scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  ylab(paste(tl)) + xlab("Sample source") + 
  #facet_grid(cols=vars(material),scales = "free", space = "free") + 
  theme(axis.text.x = element_text(colour="black", size=8, angle = 0, vjust = 1), 
        axis.text.y = element_text(colour="black", vjust=1, size=8),
        axis.title.x = element_text(face="bold", size=8),
        legend.title = element_text(face="bold", size=6),
        legend.text = element_text(size = 6),
        legend.position = c(0.87, 0.13),#"left",top
        axis.title.y = element_text(face="bold", size=8),
        plot.margin = unit(c(0.5, 0.1, 1, 1), "cm"), #top, right, bottom, left
        #plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"), #top, right, bottom, left
        strip.text.x = element_text(size=8, face = "bold", colour = "black"),
        strip.background = element_blank())

pdf(paste0('taxbin/chew/',tl,'_DiffAbundant.pdf'),width=9, height=7,useDingbats=F)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2,width=c(0.83,0.17))))
print(heat, vp=grid::viewport(layout.pos.row=1,layout.pos.col=1))
print(relabun_plot, vp=grid::viewport(layout.pos.row=1,layout.pos.col=2)) 
dev.off()



####  Make a new dataframe with the percent abudance within the entire dataset!
phy_stats <- tibble(!!tl:=taxa_names(ps)%>%str_replace('[a-z]__',''),N=nsamples(ps),PercentPhy=rowMeans(otu_table(ps)),sd=matrixStats::rowSds(otu_table(ps)))%>%mutate(se= sd/sqrt(N))
abund <- subset(phy_stats,PercentPhy > pc)  # Only take the phyla that are more than 0.01% abundant
abund_ordered <- arrange(abund, desc(PercentPhy))
abund <- abund%>%mutate(!!tl:=factor(abund%>%pull(1),levels=abund_ordered%>%pull(1)%>%rev  ) )
#Relative abundance plot 
abund_plot <- ggplot(abund, aes(y=PercentPhy , x=get(tl) ))  +
  #geom_boxplot(fill = "magenta4", colour = "black") + 
  geom_bar(stat="identity", position=position_dodge(),  fill = "magenta4", colour = "black") +
  theme_bw() + ggtitle(paste0(tl," Above ",pc,"% in All Samples")) +
  xlab(paste0(tl)) + ylab("Mean Relative Abundance (%)") +
  geom_errorbar(aes(ymin = PercentPhy -se, ymax = PercentPhy +se), width = 0.25) + coord_flip() +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=0, colour = "black", size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        legend.position="none"); 
abund_plot




