## ---- start time, echo=F----------------------------------------------------------------------------------------------------------------------------------------
start.time <- Sys.time()


## ---- install.load.pkg------------------------------------------------------------------------------------------------------------------------------------------
list.of.packages <- c("gplots","gridExtra","reshape2","tidyverse")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(gridExtra)
library(reshape)
library(gplots)


## ----setwd.1, eval=F--------------------------------------------------------------------------------------------------------------------------------------------

setwd("./") ## setwd("path/to/working/directory")

## ---- eval=T----------------------------------------------------------------------------------------------------------------------------------------------------
main <- "https://github.com/LukeBotanist/Sagebrush/archive/refs/heads/main.zip"

download.file(main, "main.zip")
unzip(zipfile = "main.zip")


## ---- setwd.2, eval=F-------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(getwd(),"/Sagebrush-main"))


## ---- warning=F-------------------------------------------------------------------------------------------------------------------------------------------------
dirs2create <- c("results/figures", "results/tables")
for (i in 1:length(dirs2create)){
  dir.create(dirs2create[i], recursive = T)
}


## ---- functions, echo=T-----------------------------------------------------------------------------------------------------------------------------------------
source("src/functions.R")


## ---- check_filestructure, echo=T-------------------------------------------------------------------------------------------------------------------------------
twee()


## ---- read_files, echo=T----------------------------------------------------------------------------------------------------------------------------------------
gt2x <- read.table("data/EBG_data/dipl_fin-genos.txt") # read diploid genotype file
gt4x <- read.table("data/EBG_data/tetra_fin-genos.txt") # read tetraploid genotype file
ind2x <- read.table("data/EBG_data/dipl_indnames") # read diploid sample names - those must be in the same order as processed by EBG
ind4x <- read.table("data/EBG_data/tetra_indnames") # read tetraploid sample names - those must be in the same order as processed by EBG
shared.vars <- read.table("data/EBG_data/shared-variants.txt") # read shared variants file


## ---- print_dimensions------------------------------------------------------------------------------------------------------------------------------------------
ebg.out.info <- data.frame(ploidy=c("2x", "4x"),
           nr.ind=c(dim(gt2x)[1],dim(gt4x)[1]),
           nr.vars=c(dim(gt2x)[2],dim(gt4x)[2]))
print(ebg.out.info)


## ---- show diploid file content, echo=T-------------------------------------------------------------------------------------------------------------------------
print(head(gt2x[1:5]))


## ---- show tetraploid file content------------------------------------------------------------------------------------------------------------------------------
print(head(gt4x[1:6]))


## ----calculate maf----------------------------------------------------------------------------------------------------------------------------------------------
gt.comb <- rbind(gt2x,gt4x) #combine the two data sets after ensuring the number of loci matches

gt.comb[gt.comb==-9] <- NA #replace missing data value (-9) with NA

#Form sum across columns (=over all individuals), and do not account for NA; divide by the number of diploid alleles (=ind*2) + tetraploid alleles (=ind*4)

maf.weighted <- data.frame("maf"=colSums(gt.comb, na.rm = T)/((nrow(gt2x)*2)+nrow(gt4x)*4))    


## ---- calculate summaries, echo=T-------------------------------------------------------------------------------------------------------------------------------
maf.summ <- maf.weighted %>% 
  summarize(maf_0.01=sum(maf>0.01),
            maf_0.03=sum(maf>0.03),
            maf_0.04=sum(maf>0.04),
            maf_0.05=sum(maf>0.05),
            maf_0.1=sum(maf>0.1)) %>% 
  t()

colnames(maf.summ) <- "variants"
print(maf.summ)


## ---- plot unfiltered maf, echo=F-------------------------------------------------------------------------------------------------------------------------------
dens.plot <- ggplot(maf.weighted, aes(x=maf))+
  geom_density()+
  geom_vline(aes(xintercept= 0.01, linetype = "maf 0.01"), colour= 'red') +
  geom_vline(aes(xintercept= 0.03, linetype = "maf 0.03"), colour= 'darkblue') +
  geom_vline(aes(xintercept= 0.04, linetype = "maf 0.04"), colour= 'midnightblue') +
  geom_vline(aes(xintercept= 0.05, linetype = "maf 0.05"), colour= 'blue') +
  geom_vline(aes(xintercept= 0.1, linetype = "maf 0.1"), colour= 'lightblue') +
  scale_linetype_manual(name = "limit", values = c(2 ,2, 2, 2, 2),
                        guide = guide_legend(override.aes = list(color = c("red","darkblue", 'midnightblue',"blue","lightblue"))))+
  theme_classic()+
  annotation_custom(tableGrob(maf.summ),xmin = 0.75, xmax=1)

print(dens.plot)
ggsave("results/figures/maf_distr_overall.pdf", dens.plot, device="pdf")


## ----plot different maf filters,echo=F, dev="pdf", message=FALSE, warning=F-------------------------------------------------------------------------------------
maf.thresholds <- c(0.01, 0.025, 0.04, 0.05, 0.075, 0.1)
plot.list = list()
for (i in 1:length(maf.thresholds)){
  maf.loop <- maf.weighted %>% 
    filter(maf > maf.thresholds[i])
  var.count <- maf.weighted %>% 
    filter(maf > maf.thresholds[i]) %>% summarize(n())
  
  p <- ggplot(maf.loop, aes(x=maf)) + 
    geom_density() + 
    theme_classic() + 
    labs(title = paste0("maf > ", maf.thresholds[i]), 
         subtitle = paste(var.count,"remaining variants"))
  
  plot.list[[i]] <- p
}


ggsave("results/figures/maf_distr_filtered.pdf", do.call(grid.arrange,plot.list), device="pdf")


## ----write maf, echo=T------------------------------------------------------------------------------------------------------------------------------------------
vars.maf <- cbind(shared.vars, maf.weighted) # append the allele frequencies to the shared variants for subsequent filtering.

names(vars.maf) <- c("CHROM", "POS", "maf") # rename to avoid incompatibilities later on

write.csv(vars.maf, "results/tables/shared_vars_maf.csv", row.names = F, quote = F)


## ----filter by maf, echo=T--------------------------------------------------------------------------------------------------------------------------------------
maf.filter <- 0.04 # set here the best suitable filter based on the figures before
#we can  re-use the previously calculated "maf_weighted" object. as maf is calculated per locus, we can add this info right to the shared variants file

# custom function located in the source file to filter for a specific maf
gt.comb.filtered <- filter_maf(gt.comb, vars.maf, maf.filter) 

ind.comb <- rbind(ind2x,ind4x) 

names(ind.comb) <- c("Sample_ID","Ploidy")
data.comb <- cbind(ind.comb, gt.comb.filtered)
write.csv(data.comb, row.names = F, paste0("results/tables/SNPs_maf",maf.filter,".csv"))


## ---- download polyrelatedness--------------------------------------------------------------------------------------------------------------------------
## PolyRelatedness <- "https://github.com/huangkang1987/polyrelatedness/raw/master/PolyRelatedness_1.11b.zip"
## 
## #path <- "path/to/software/dir"
## path <- "~/Software/"
## 
## download.file(PolyRelatedness, paste0(path,"PolyRelatedness_1.11b.zip"))
## unzip(zipfile = paste0(path,"PolyRelatedness_1.11b.zip"))
## 


## ---- recode2polyrel--------------------------------------------------------------------------------------------------------------------------------------------
polyrel.in <- recode2polyrel(data.comb)


## ----inspect_conversion-----------------------------------------------------------------------------------------------------------------------------------------
head(polyrel.in[1:5])
tail(polyrel.in[1:5])


## ----filter individuals-----------------------------------------------------------------------------------------------------------------------------------------

rm.inds <- c("C4_08","O1_01", "P3_02_rep", "O3_03_rep", "O2_03", "C3_19")

polyrel_input <- polyrel.in %>% 
  select(-Ploidy) %>% 
  mutate(pop=substr(Sample_ID,1,2)) %>% 
  relocate(pop, .after=Sample_ID) %>% 
  filter(!grepl(paste(rm.inds,collapse="|"),.$Sample_ID))



## ---- write polyrel input genotypes-----------------------------------------------------------------------------------------------------------------------------
filename <- "all_maf004"

write_delim(polyrel_input, paste0("~/Software/PolyRelatedness_1/",filename), delim ="\t", quote = "none", eol = "\n") #write tab separated table



## ----comment=''-------------------------------------------------------------------------------------------------------------------------------------------------
cat(readLines('~/Software/PolyRelatedness_1/header.txt'), sep = '\n')


## ---- combine config and data for polyrelatedness---------------------------------------------------------------------------------------------------------------
outdir <- "~/Software/PolyRelatedness_1/"
outfile_name <- paste0(filename,"_polyrel_in.txt")

combine <- paste0("cat ", outdir,"header.txt ",outdir,filename," ",outdir,"end.txt > ",outdir,outfile_name)
system(combine)



## ---- run polyrel from file-------------------------------------------------------------------------------------------------------------------------------------
polyrel_outfile <- paste0(outdir,filename,"_polyrel_out.txt")
command <- paste("~/Software/PolyRelatedness_1/PolyRelatedness",paste0(outdir,outfile_name),polyrel_outfile,"e 5 0")
system(command,intern = T)


## ---- filter.maf.str, echo=T------------------------------------------------------------------------------------------------------------------------------------
data.comb.str <- data.comb %>% 
  filter(!grepl(paste(rm.inds,collapse="|"),.$Sample_ID)) %>% 
  select(c(-Sample_ID,-Ploidy)) # filter individuals/pops/etc.

S_combined_structure <- filter_maf(gt.comb, vars.maf, maf.filter, remove_columns = F) #use function to filter for maf less than desired value


## ---- select.single.snp, echo=T---------------------------------------------------------------------------------------------------------------------------------
sSNP_S_combined_structure <- select_singleSNP.3(S_combined_structure, window_size = 139)# use select single snp function to filter for single snp per locus/bp window


## ---- confirm.dist, echo=T, warning=F, message=F----------------------------------------------------------------------------------------------------------------
sSNP_S_combined_structure %>% 
  group_by(CHROM) %>% 
  mutate(diff=POS-lag(POS)) %>% 
  select(CHROM, diff) %>% 
  summarize(min_dist=min(diff, na.rm=T)) %>% 
  head(9)


## ---- recode2str, echo=T----------------------------------------------------------------------------------------------------------------------------------------
sSNP_S_combined_structure_t <- sSNP_S_combined_structure %>% 
  ungroup %>% 
  select(c(-CHROM,-POS, -maf,-Locus)) %>% 
  t() #remove unwanted columns and transpose data set so it's compatible with Structure (i.e. ind in row, loci in columns)

###now convert to structure genos, which requires each allele of an individual in a row and each locus in a column. 

sSNP_S_combined_structure_t <- cbind(ind.comb, sSNP_S_combined_structure_t) # combine name of samples with data set

data.comb.str <- sSNP_S_combined_structure_t %>% 
  filter(!grepl(paste(rm.inds,collapse="|"),.$Sample_ID)) #%>% dplyr::rename("ploidy"="Ploidy")

sSNP_fin <- recode2structure(data.comb.str, as.pseudotetraploids = F) # utilize function to recode samples to structure.As we're working with a heteroploid data set the diploids need to be encoded either as having missing alleles (i.e. 1,2,-9,9) or as pseudotetraploids, where the homozygote AA or TT need to be encoded as AAAA or TTTT and the heterozygotes as duplex heterozygotes AATT


## ---- eval.diploids.str-----------------------------------------------------------------------------------------------------------------------------------------
head(sSNP_fin[1:5],5)


## ---- eval.tetraplois.str---------------------------------------------------------------------------------------------------------------------------------------
tail(sSNP_fin[1:5],5)


## ---- write str_in file---------------------------------------------------------------------------------------------------------------------------------
#Let's see what we got for our structure analysis - we can also write these results as table as structure requires those numbers for the input file. 
dataset_dims <- dim(sSNP_fin)
filename <- paste0("structure_",dataset_dims[1]/4,"ind",dataset_dims[2]-1,"SNPs")
print(paste("In the final structure data set are",dataset_dims[1]/4,"Individuals and", dataset_dims[2]-1 , "Variants"))
str.info <- data.frame(nr_ind=dataset_dims[1]/4,
                       nr_vars=dataset_dims[2]-1)
write.table(str.info, paste0("results/tables/",filename,".info.txt"), row.names = F, quote=F)

#filename <- "ebg_structure004_correct_sSNP"
write.table(sSNP_fin, file= paste0("results/tables/",filename,".str"), quote=F, col.names = F, sep = "\t", row.names = F)



## ---- print pca, fig.height = 12, fig.width = 12, fig.align = "center", echo=T, message=FALSE, warning=F--------------------------------------------------------

polyrel_out <- read.table(polyrel_outfile, header=T, skip=7)

ind <- as.data.frame(row.names(polyrel_out))
colnames(ind) <- "Sample_ID"
sample_info <- merge(ind, ind.comb, by="Sample_ID", sort = F)
#extract eigenvalues from matrix
eigen <- eigen(polyrel_out)
#extract eigenvalues
ev <- NULL
for (i in 1:5){
  ex_var <- round(eigen$values[i]*100/sum(eigen$values),2)
  ev <- cbind(ev,ex_var)
}



pcoa <- as.data.frame(eigen$vectors)
pca_combined <- cbind(sample_info, pcoa)

pc12 <- pca_combined %>% mutate(Ploidy=as.factor(Ploidy), Population=as.factor(substr(Sample_ID,1,2))) %>% 
  ggplot(aes(x=V1,y=V2, pch=Ploidy, col=as.factor(Population)))+
  geom_point(stroke=.8, size=3)+
  theme_classic()+
  labs(x=paste0("PC1 (", ev[1],"%)"), y=paste0("PC2 (", ev[2],"%)"))+
  scale_shape_manual(values=c(1,15))+
  labs(col="Population")

pc13 <- pca_combined %>% mutate(Ploidy=as.factor(Ploidy), Population=as.factor(substr(Sample_ID,1,2))) %>% 
  ggplot(aes(x=V1,y=V3, pch=Ploidy, col=Population))+
  geom_point(stroke=.8, size=3)+
  theme_classic()+
  labs(x=paste0("PC1 (", ev[1],"%)"), y=paste0("PC2 (", ev[3],"%)"))+
  scale_shape_manual(values=c(1,15))+
  labs(col="Population")


ggpubr::ggarrange(pc12, pc13, common.legend = T, legend = "right")

ggsave(plot = pc12, path="results/figures/", filename = "PCoA12.pdf", device = "pdf")




## ---- save heatmap, fig.height = 24, fig.width = 24, fig.align = "center", echo=T, message=F, warning=FALSE, results="hide"-------------------------------------
polyrel_out <- read.table(polyrel_outfile, header=T, skip=7)
relatedness_matrix <- as.matrix(polyrel_out) # remove diagonal (i.e. individual compared with itself) for better contrast.
diag(relatedness_matrix) <- NA # remove diagonal values
ind <-  row.names(polyrel_out)
heatmap_name <- substr(filename,1,nchar(filename)-4)
#map with dendrogram


## ---- save heatmap.1, fig.height = 24, fig.width = 24, fig.align = "center", echo=F, message=F, warning=FALSE, results="hide"-----------------------------------
pdf(paste0("results/figures/",heatmap_name,"heatmap_combined.pdf"),width=24,height=18)
heatmap.2(relatedness_matrix, trace="none", cexRow=0.2,
          cexCol = 0.2, labRow = ind,
          labCol = ind,
          col= colorRampPalette(c("white",
                                  "white",
                                  "lemonchiffon",
                                  "yellow",
                                  "red",
                                  "darkorchid",
                                  "midnightblue", "black"))(70))

dev.off()


## ---- print heatmap, fig.height = 24, fig.width = 24, fig.align = "center", echo=T, message=F, warning=FALSE----------------------------------------------------
heatmap.2(relatedness_matrix, trace="none", cexRow=0.2,
          cexCol = 0.2, labRow = ind,
          labCol = ind,
          col= colorRampPalette(c("white",
                                  "white",
                                  "lemonchiffon",
                                  "yellow",
                                  "red",
                                  "darkorchid",
                                  "midnightblue", "black"))(70))


## ---- allele freqs, echo=F, warning=F, message=F, results="hide"------------------------------------------------------------------------------------------------
tetraploids <- data.comb %>% filter(Ploidy==4) %>% as.data.frame()

all <- est_inheritance(tetraploids)
ggplot(all, aes(x =allele_freqs , y = gt_freq, group = type, color = type, linetype=type)) + 
  geom_smooth(method='gam',se = F, fullrange=T) +
  labs(x="Allele Frequency", y="Genotype Frequency")

##individual sites

tetraploids_site <- tetraploids %>% 
  mutate(site=as.factor(substr(Sample_ID,1,1)))

inheritance_plots <- NULL

for (i in levels(tetraploids_site$site)){
  print(i)
  site <- tetraploids_site %>% 
    filter(site==i)
  
  estimates_site <- est_inheritance(tetraploids)
  inheritance_plots[[i]] <- ggplot(estimates_site, aes(x =allele_freqs , y = gt_freq, group = type, color = type, linetype=type)) + 
  geom_smooth(method='gam',se = F, fullrange=T)
}

ggpubr::ggarrange(inheritance_plots$C,inheritance_plots$L, inheritance_plots$O, inheritance_plots$P,inheritance_plots$C, common.legend = T)



## ---- echo=F----------------------------------------------------------------------------------------------------------------------------------------------------
end.time <- Sys.time()

elapsed.time <- round((end.time - start.time), 3)

twee()


## ---- warning=F, message=F--------------------------------------------------------------------------------------------------------------------------------------
ssp.assignments <- read.table("data/STRUCTURE_data/K4_assignments.txt", header=T)

assignments.reformat <- ssp.assignments %>% 
  dplyr::mutate(Sample_ID = ifelse(grepl("NA", Sample_ID), 
                                   paste0(substr(Sample_ID,1,2), "_",
                                          substr(Sample_ID,3,4),"_rep"),Sample_ID))


pca_combined <- cbind(sample_info, pcoa) %>% 
  left_join(., assignments.reformat, by="Sample_ID")

pca_combined <-pca_combined %>%  mutate(ssp_pl=paste0(ssp,"_",Ploidy, "x"))

shapes_pca <- c("vas1_2x"=2,
            "vas1_4x"=17,
            "hybrid"=10,
            "tri_2x"=1,
            "tri_4x"=16,
            "wyo1_2x"=0,
            "wyo1_4x"=22,
            "wyo2_4x"=22)
cols_pca <- c("vas1_2x" = "#5b9bd5", # darkblue
              "vas1_4x" = "#5b9bd5", # lightblue
              "tri_2x" = "#5fa659",
              "tri_4x" = "#5fa659",# green
              "wyo1_2x" = "#c6793a", # brown
              "wyo1_4x" = "#c6793a", # brown
              "wyo2_4x" = "#aea642")


labels <- c("vas1_2x"="2x vaseyana",
            "hybrid"="hybrid",
            "vas1_4x"="4x vaseyana",
            "tri_2x"="2x tridentata",
            "wyo1_4x"="4x wyomingensis 1",
            "wyo2_4x"="4x wyomingensis 2")
name <- "subspecies"


pc12 <- pca_combined %>% mutate(Ploidy=as.factor(Ploidy), Population=as.factor(substr(Sample_ID,1,2))) %>% 
  ggplot(aes(x=V1,y=V2, pch=ssp_pl, col=ssp_pl))+
  geom_point(stroke=.8, size=3)+
  theme_classic()+
  labs(x=paste0("PC1 (", ev[1],"%)"), y=paste0("PC2 (", ev[2],"%)"))+
  scale_color_manual(values=cols_pca, limits=force, name=name, labels=labels)+
  scale_fill_manual(values=cols_pca, limits=force, name=name, labels=labels)+
  scale_shape_manual(values=shapes_pca, limits=force, name=name, labels=labels)

pc13 <- pca_combined %>% mutate(Ploidy=as.factor(Ploidy), Population=as.factor(substr(Sample_ID,1,2))) %>% 
  ggplot(aes(x=V1,y=V3, pch=ssp_pl, col=ssp_pl))+
  geom_point(stroke=.8, size=3)+
  theme_classic()+
  labs(x=paste0("PC1 (", ev[1],"%)"), y=paste0("PC2 (", ev[3],"%)"))+
  scale_color_manual(values=cols_pca, limits=force, name=name, labels=labels)+
  scale_fill_manual(values=cols_pca, limits=force, name=name, labels=labels)+
  scale_shape_manual(values=shapes_pca, limits=force, name=name, labels=labels)


ggpubr::ggarrange(pc12, pc13, common.legend = T)




## ---- warning=F, message=F, results="hide"----------------------------------------------------------------------------------------------------------------------

tetraploids_site <- tetraploids %>% 
  mutate(site=as.factor(substr(Sample_ID,1,1)))


for (i in levels(tetraploids_site$site)){
  site <- tetraploids_site %>% 
    filter(site==i)
  
  estimates_site <- est_inheritance(tetraploids)

  p <- ggplot(estimates_site, aes(x = allele_freqs , y = gt_freq, group = type, color = type, linetype=type)) + 
  geom_smooth(method='gam',se = F, fullrange=T)+
  theme_classic()+
  labs(title=paste("Inheritance Site", i))+
  scale_color_manual(values = c("TTTT" = "blue2",
                                "ATTT" = "blue2",
                                "AATT"="cyan3",
                                "AAAT"="cyan", "AAAA"="cyan"), 
                     labels=rev(c("AAAA","AAAT", "AATT", "ATTT", "TTTT")), name="Genotype")+
  scale_linetype_manual(values=c("TTTT" = "solid", 
                                 "ATTT" = "dashed", 
                                 "AATT"="solid", 
                                 "AAAT"="dashed",
                                 "AAAA"="solid"), 
                        labels=rev(c("AAAA","AAAT", "AATT", "ATTT", "TTTT")), name="Genotype")+
  labs(y="Genotype Frequency", x="Allele Frequency")
  
  
  
  inheritance_plots[[i]] <- p
}


ggpubr::ggarrange(inheritance_plots$C,inheritance_plots$L, inheritance_plots$O, inheritance_plots$P,inheritance_plots$C, common.legend = T)


for (i in 1:length(inheritance_plots)){
  pdf(paste0("results/figures/inheritance_plot",names(inheritance_plots)[i],".pdf"))
  print(inheritance_plots[i])
  dev.off()
}


