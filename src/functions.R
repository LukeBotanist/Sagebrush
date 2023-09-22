library(tidyverse)
library(gridExtra)
library(reshape)


#The twee() was developed by Jennifer Bryan and was derived from https://gist.github.com/jennybc/2bf1dbe6eb1f261dfe60
twee <- function(path = getwd(), level = Inf) {
  
  fad <-
    list.files(path = path, recursive = TRUE,no.. = TRUE, include.dirs = TRUE)
  
  fad_split_up <- strsplit(fad, "/")
  
  too_deep <- lapply(fad_split_up, length) > level
  fad_split_up[too_deep] <- NULL
  
  jfun <- function(x) {
    n <- length(x)
    if(n > 1)
      x[n - 1] <- "|__"
    if(n > 2)
      x[1:(n - 2)] <- "   "
    x <- if(n == 1) c("-- ", x) else c("   ", x)
    x
  }
  fad_subbed_out <- lapply(fad_split_up, jfun)
  
  cat(unlist(lapply(fad_subbed_out, paste, collapse = "")), sep = "\n")
}

filter_maf <- function(gt.comb, vars.maf, maf.filter, remove_columns=T){
  gt.comb.t <- t(gt.comb) # transpose gt.comb to get loci in rows and ind in cols. This will make it much faster
  gt.comb.t.maf <- cbind(vars.maf, gt.comb.t)
  head(gt.comb.t.maf)
  if (remove_columns==T){
    gt.comb.filtered <- gt.comb.t.maf %>% filter(maf > maf.filter) %>% select(-CHROM,-POS,-maf)
  } else {
    gt.comb.filtered <- gt.comb.t.maf %>% filter(maf > maf.filter)
  }
  print(paste(dim(gt.comb.filtered)[1], "Variants left after filtering for maf >", maf.filter))
  S.filtered.t <- t(gt.comb.filtered) # transpose back in original format
  if (remove_columns==T){
    return(S.filtered.t)
  } else {
    return(gt.comb.filtered)
  }
}

recode <- function(genos, output){        # this function will automatically detect ploidy from file as prepared by ebg and convert it to character genotypes
  print(paste("Converting genotypes to", output,"genotypes"))
  nallele <- max(genos) #assess max alleles (ploidy) from genotypes
  print(paste0(dim(genos)[1]," ",nallele,"x individuals and ", dim(genos)[2], " variant sited detected..."))
  gts=NULL
  if (output=="character"){
    for (i in 0:nallele){                   #loop through alleles to convert numeric  genotype to character genotype, i.e. 0=AAAA, 1=AAAT, etc.
      gt <- c(rep("A",nallele-i), rep("T",i))
      gts <- rbind.data.frame(gts, gt, stringsAsFactors = F)
    }
  } else {
    for (i in 0:nallele){                   #loop through alleles to convert numeric  genotype to numeric genotype, i.e. 0=1111, 1=1112, etc.
      gt <- c(rep(1,nallele-i), rep(2,i))
      gts <- rbind.data.frame(gts, gt, stringsAsFactors = F)
    }
  }
  gentypes <- as.data.frame(unite(gts,col=gts, sep=",")) # combine all columns to get genotypes
  genos_recode <- genos
  genos_recode <-t(genos_recode) # transpose to make replacing values A LOT faster
  for (j in 0:nallele){ #loop through dataset and replace numeric genotype with character genotype
    genos_recode[genos_recode==j] <- gentypes[j+1,]
  }
  
  genos_recode[genos_recode==-9] <- paste0(rep(9,nallele), collapse = ",") # replace missing data with 9
  genos_recode <-as.data.frame(t(genos_recode)) #backtranspose
  return(genos_recode)
}

recode2polyrel <- function(genos){        # this function will automatically detect ploidy from file as prepared by ebg and convert it to character genotypes
  ploidy_detected <- as.factor(genos$Ploidy)
  print(paste("Dataset with", length(levels(ploidy_detected)),"different ploidy levels detected. Converting accordingly..."))
  counts <- genos %>% group_by(Ploidy) %>% summarize(nind=n())
  genos_combined <- NULL
  for (i in counts$Ploidy){
    genos_filtered <- genos %>% filter(Ploidy==i)
    nallele <- max(genos_filtered$Ploidy) # assess alleles for conversion
    print(paste0(dim(genos_filtered)[1]," ",nallele,"x individuals and ", dim(genos)[2]-2, " variant sited detected. Converting..."))
    gts=NULL
    for (i in 0:nallele){                   #loop through alleles to convert numeric  genotype to character genotype, i.e. 0=1111, 1=1112, etc.
      gt <- c(rep(5,nallele-i), rep(6,i))
      gts <- rbind.data.frame(gts, gt, stringsAsFactors = F)
    }
    gentypes <- as.data.frame(unite(gts,col=gts, sep="")) # combine all columns to get genotypes
    data_info <- genos_filtered[,1:2]
    genos_recode <-t(genos_filtered[,c(-1,-2)]) # transpose to make replacing values A LOT faster
    for (j in 0:nallele){ #loop through dataset and replace numeric genotype with character genotype
      genos_recode[genos_recode==j] <- gentypes[j+1,]
    }
    #genos_recode[genos_recode==-9] <- paste0(rep(-9,nallele), collapse = "") # replace missing data with 9
    genos_recode[is.na(genos_recode)] <- paste0(rep(0,nallele), collapse = "")
    genos_recode <-as.data.frame(t(genos_recode)) #backtranspose
    genos_recode <- cbind(data_info,genos_recode)
    genos_combined <- rbind(genos_combined,genos_recode)
  }
  return(genos_combined)
}

select_singleSNP.3 <- function(combined, window_size=139){
  print(paste("Extracting a single SNP at least", window_size, "bp apart"))
  ploidies_shared_var <- combined
  SNPs <- NULL
  for (j in levels(as.factor(ploidies_shared_var$CHROM))){
    Chm <- ploidies_shared_var %>% filter(CHROM==j)
    k <- 1
    init <- cbind(Chm[1,], Locus=k)
    #print(init)
    SNPs_chm <- init
    if(nrow(Chm)==0){
      next
    }else{
      if(nrow(Chm)==1){
        SNPs <- rbind(SNPs, SNPs_chm)
        next
      } else {
        for (i in 2:nrow(Chm)){
          #print(Chm[i,])
          SNP <- Chm[i,]
          k <- ifelse((SNP[,2]-SNPs_chm[i-1,2])<139,k,k+1)
          SNP$Locus <- k
          SNPs_chm <- rbind(SNPs_chm, SNP)
          
        }
      }
    }
    SNPs <- rbind(SNPs, SNPs_chm)
    print(paste("Extracted", dim(SNPs_chm)[1],"SNPs for Chromosome", j))
  }
  sSNPs <- SNPs %>% group_by(CHROM, Locus) %>% slice_sample(n=1)
  print(paste("Extracted", dim(sSNPs)[1],"single SNPs"))
  
  return(sSNPs)
}


recode2structure <- function(singleSNP_file, as.pseudotetraploids=F){        # this function will automatically detect ploidy from file as prepared by ebg and convert it to character genotypes
  ploidy_detected <- as.factor(singleSNP_file$Ploidy)
  print(paste("Dataset with", length(levels(ploidy_detected)),"different ploidy levels detected. Converting accordingly..."))
  counts <- singleSNP_file %>% group_by(Ploidy) %>% summarize(nind=n())
  genos_combined <- NULL
  as.pseudotetraploids <- F
  for (i in counts$Ploidy){
    genos_filtered <- singleSNP_file %>% filter(Ploidy==i)
    genos_filtered[1:10,1:10]
    nallele <- max(genos_filtered$Ploidy) # assess alleles for conversion
    print(paste0(dim(genos_filtered)[1]," ",nallele,"x individuals and ", dim(genos_filtered)[2], " variant sited detected. Converting..."))
    as.pseudotetraploids <- F
    if(nallele==2 & as.pseudotetraploids==F){
      print("diploid genotypes detected...\n encoding diploids with missing alleles")
      gts=NULL
      for (i in 0:nallele){                   #loop through alleles to convert numeric  genotype to character genotype, i.e. 0=1111, 1=1112, etc.
        gt <- c(rep(1,nallele-i), rep(2,i), rep(-9,nallele))
        gts <- rbind.data.frame(gts, gt, stringsAsFactors = F)
      } 
    } else if (nallele==2 & as.pseudotetraploids==T) {
      print("diploid genotypes detected...\n encoding diploids as pseudotetraploids")
      gts=NULL
      for (i in 0:nallele){                   #loop through alleles to convert numeric  genotype to character genotype, i.e. 0=1111, 1=1112, etc.
        gt <- c(rep("1,1",nallele-i), rep("1,1",i))
        gts <- rbind.data.frame(gts, gt, stringsAsFactors = F)
      } 
    } else {
      print("Polyploid genotypes detected...\n encoding as polyploids")
      gts=NULL
      for (i in 0:nallele){                   #loop through alleles to convert numeric  genotype to character genotype, i.e. 0=1111, 1=1112, etc.
        gt <- c(rep(1,nallele-i), rep(2,i))
        gts <- rbind.data.frame(gts, gt, stringsAsFactors = F)
      } 
    } 
    gentypes <- as.data.frame(unite(gts,col=gts, sep=",")) # combine all columns to get genotypes
    data_info <- genos_filtered[,1:2]
    genos_recode <-t(genos_filtered[,c(-1,-2)]) # transpose to make replacing values A LOT faster
    for (j in 0:nallele){ #loop through dataset and replace numeric genotype with character genotype
      genos_recode[genos_recode==j] <- gentypes[j+1,]
    }
    #genos_recode[genos_recode==-9] <- paste0(rep(-9,nallele), collapse = "") # replace missing data with 9
    if (nallele==2){
      genos_recode[is.na(genos_recode)] <- paste0(rep(-9,nallele*2), collapse = ",")
    }else{
      genos_recode[is.na(genos_recode)] <- paste0(rep(-9,nallele), collapse = ",")
      
    }
    genos_recode <-as.data.frame(t(genos_recode)) #backtranspose
    genos_recode <- cbind(data_info,genos_recode)
    genos_combined <- rbind(genos_combined,genos_recode)
  }
  genos_combined <- as.data.frame(genos_combined)
  genos_unnested <- genos_combined %>% dplyr::select(-2)
  str_in <- genos_unnested %>% separate_rows(2:ncol(.),sep=",")
  return(str_in)
}

est_inheritance <- function(genotypes){
  tetraploid_genos <- genotypes
  
  nallele <- max(tetraploid_genos$Ploidy)
  
  tetraploid_genos <- tetraploid_genos %>% select(-matches("Sample_ID|Ploidy|pop|site"))
  
  gt_freqs <- NULL
  for (i in 0:nallele){                   #loop through alleles to convert numeric  genotype to character genotype, i.e. 0=AAAA, 1=AAAT, etc.
    gt <- paste0(c(rep("A",nallele-i), rep("T",i)),collapse="")
    freqs <- data.frame(colSums(tetraploid_genos == i, na.rm = T)/(nrow(tetraploid_genos)-colSums(is.na(tetraploid_genos))))
    colnames(freqs) <- gt
    gt_freq <- as.data.frame(t(freqs))
    gt_freqs <- rbind(gt_freqs, gt_freq)
  }
  
  allele_counts <- NULL
  i <- 0
  while (i < nallele){                  
    allele_count <- colSums(tetraploid_genos == i, na.rm = T)*(nallele-i)
    allele_counts <- rbind(allele_counts, allele_count)
    i <- i + 1
  }
  
  allele_freqs <- colSums(allele_counts, na.rm = T) / (nrow(tetraploid_genos)*nallele-colSums(is.na(tetraploid_genos)))

  gt_all_freq <- data.frame(rbind(gt_freqs,allele_freqs=allele_freqs))
  gt_all_freq_t <- data.frame(t(gt_all_freq))
  gt_all_freq_t_long <- gather(gt_all_freq_t, type, gt_freq, matches(rownames(gt_freqs)))
  return(gt_all_freq_t_long)
}


