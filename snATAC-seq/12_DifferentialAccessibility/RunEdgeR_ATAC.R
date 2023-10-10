
library(edgeR)

getFilteredSamples = function(cursamples, cellgroup, mincells=50, filter=TRUE){
  rv = c()
  for(cursample in cursamples){
    cursample = gsub("-", "\\.", cursample)
    cursample_cells = cellnumbers[paste0(cursample,"_cells")][cellgroup,]
    
    if((filter==FALSE || cursample != 'pHC15') && !is.na(cursample_cells) && cursample_cells >= mincells){ #One healthy sample was excluded
      rv = c(rv,  gsub("\\.", "-",  cursample))
    }
  }
  return(rv)
}

filterGenes = function(curdata, genepercentages, samples, fraccells=0.00, minsamplesingroup=3){
  headers=rep("", length(samples))
  index=1
  for(cursample in samples){
    headers[index]=gsub("-", "\\.", cursample)
    index = index+1
  }
  tryCatch({
    percentsel = (rowSums(genepercentages[, headers] > fraccells) >= minsamplesingroup)
    percentgenes = rownames(genepercentages)[percentsel]
    
    #mediansel = apply(curdata[headers], 1, function(x) median(x) > 10)
    #mediangenes = rownames(curdata)[mediansel]
    #return(intersect(percentgenes, mediangenes))
    return(percentgenes)
  }, error = function(e){print(e)}
    )
}



runEdgeR = function(curdata, cellgroup, header, c1, c2, c1_label, c2_label, conditionmap, minsamplesingroup=3){
  label = paste0(c1_label, "_vs_", c2_label)
  np1 = 0
  np05 = 0
  
  #Get samples not used for comparison, but included for covariates
  c0=conditionmap$Sample[!(conditionmap$Sample %in% c(c1,c2))]
  
  c0_filtered=getFilteredSamples(c0, cellgroup)
  c1_filtered=getFilteredSamples(c1, cellgroup)
  c2_filtered=getFilteredSamples(c2, cellgroup)
  
  c0samples = paste(c0_filtered, collapse=',')
  c1samples = paste(c1_filtered, collapse=',')
  c2samples = paste(c2_filtered, collapse=',')
  
  genepercentages = read.csv(paste0("/Users/athib/Desktop/CovidRSV/snATAC2/DA/PercentExpressed/",cellgroup,"_countmatrix_samples.txt_percents.txt"), row.names = 1, sep = '\t')
  
  if(length(c1_filtered) >= minsamplesingroup &&  length(c2_filtered) >= minsamplesingroup){
    
    g1 = filterGenes(curdata, genepercentages, c1_filtered)
    g2 = filterGenes(curdata, genepercentages, c2_filtered)
    selectedgenes = unique(c(g1,g2))
    
    c0_data = curdata[,(header %in% c0_filtered)]
    c1_data = curdata[,(header %in% c1_filtered)]
    c2_data = curdata[,(header %in% c2_filtered)]
    combined_data = cbind(c0_data, c1_data, c2_data)
    
    combined_data = combined_data[rownames(combined_data) %in% selectedgenes,]
    
    groupvars = c(rep("NULL", length(c0_filtered)), rep(c1_label, length(c1_filtered)), rep(c2_label, length(c2_filtered)))
    groupvars = factor(groupvars, c("NULL",c1_label, c2_label))
    batchvars = c(conditionmap[match(c0_filtered, conditionmap$Sample,),]$Batch, conditionmap[match(c1_filtered, conditionmap$Sample,),]$Batch, conditionmap[match(c2_filtered, conditionmap$Sample,),]$Batch)
    sexvars = c(conditionmap[match(c0_filtered, conditionmap$Sample,),]$Sex, conditionmap[match(c1_filtered, conditionmap$Sample,),]$Sex, conditionmap[match(c2_filtered, conditionmap$Sample,),]$Sex)
    agevars = c(conditionmap[match(c0_filtered, conditionmap$Sample,),]$Age, conditionmap[match(c1_filtered, conditionmap$Sample,),]$Age, conditionmap[match(c2_filtered, conditionmap$Sample,),]$Age)
    racevars = c(conditionmap[match(c0_filtered, conditionmap$Sample,),]$Race, conditionmap[match(c1_filtered, conditionmap$Sample,),]$Race, conditionmap[match(c2_filtered, conditionmap$Sample,),]$Race)
    steroidvars = c(conditionmap[match(c0_filtered, conditionmap$Sample,),]$Steroid, conditionmap[match(c1_filtered, conditionmap$Sample,),]$Steroid, conditionmap[match(c2_filtered, conditionmap$Sample,),]$Steroid)
    
    y <- DGEList(counts=combined_data,group=groupvars)
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    design <- model.matrix(~batchvars+ sexvars + agevars + groupvars)

    y <- estimateDisp(y,design)
    fit <- glmQLFit(y,design)
    ncoef=ncol(design)
    coefvar = rep(0,ncoef)
    coefvar[ncoef] = 1
    coefvar[ncoef-1] = -1
    lrt <- glmQLFTest(fit,contrast=coefvar)
    
    curtable = lrt$table
    adjustedpval = p.adjust(curtable$PValue, method = 'fdr')
    curtable['p.adjusted'] = adjustedpval
    curtable = curtable[order(curtable$p.adjusted),]
    
    np1 = sum(curtable['p.adjusted'] < 0.1)
    
    np05=0
    np05 = sum(curtable['p.adjusted'] < 0.05)
    
    #write.csv(c1genes, paste0(label,"_", c1_label,".csv"))
    #write.csv(c2genes, paste0(label,"_", c2_label,".csv"))
    write.csv(curtable, paste0(label,".csv"))
    
  }
  
  return(c(label, np1, np05, c1samples, c2samples))
}

ncomparisons=31
conditionmap = read.csv("/Users/athib/Desktop/CovidRSV/snATAC2/DA/SampleGroup.csv")
cellnumbers = read.csv("/Users/athib/Desktop/CovidRSV/snATAC2/DA/ClusterStats.csv", row.names = 1)



for(cursample in rownames(cellnumbers)){
#for(cursample in c("NaiveCD4", 'CD4_TCM', 'CD4_TEM', 'MemoryT_3', 'NaiveTreg', 'MemoryTreg')){
    
  edgermatrix = paste0("/Users/athib/Desktop/CovidRSV/snATAC2/DA/EdgeRMatrixFiles/", cursample, "_countmatrix_samples.txt")
  dir.create(file.path("/Users/athib/Desktop/CovidRSV/snATAC2/DA/EdgeROutputs/", cursample))
  setwd(paste0("/Users/athib/Desktop/CovidRSV/snATAC2/DA/EdgeROutputs/", cursample, "/"))
  
  
  unformattedheader = read.csv(edgermatrix, header = FALSE, nrows = 1, sep='\t')
  unformattedheader = unformattedheader[!is.na(unformattedheader)]
  cellgroup = cursample#strsplit(basename(edgermatrix), '\\.')[[1]][1]
  curdata = read.csv(edgermatrix, row.names=1, sep='\t')
  
  
  summary_comparison = rep("", ncomparisons)
  summary_np1 = rep(NA, ncomparisons)
  summary_np05 = rep(NA, ncomparisons)
  summary_c1samples = rep("", ncomparisons)
  summary_c2samples = rep("", ncomparisons)
  
  
  samples_healthy = conditionmap$Sample[conditionmap$Group %in% c("Healthy")]
  
  samples_cg1 = conditionmap$Sample[conditionmap$Group %in% c("COVID G1")]
  samples_cg2 = conditionmap$Sample[conditionmap$Group %in% c("COVID G2")]
  samples_cg3 = conditionmap$Sample[conditionmap$Group %in% c("COVID G3") & (conditionmap$Steroid == FALSE)]
  samples_cg3_s = conditionmap$Sample[conditionmap$Group %in% c("COVID G3") & (conditionmap$Steroid == TRUE)]
  samples_covid = conditionmap$Sample[conditionmap$Group %in% c("COVID G2", "COVID G3") & (conditionmap$Steroid == FALSE)]
  
  samples_rg1 = conditionmap$Sample[conditionmap$Group %in% c("RSV G1")]
  samples_rg2 = conditionmap$Sample[conditionmap$Group %in% c("RSV G2")]
  samples_rg3 = conditionmap$Sample[conditionmap$Group %in% c("RSV G3")]
  samples_rsv = conditionmap$Sample[conditionmap$Group %in% c("RSV G1", "RSV G2", "RSV G3")]
  
  allcomparisons = c('Healthy', 'COVID G1', 'COVID G2', 'COVID G3', 'COVID G3_S', 'RSV G1', 'RSV G2', 'RSV G3')
  allcomparisons_samples= list(samples_healthy, samples_cg1, samples_cg2, samples_cg3, samples_cg3_s, samples_rg1, samples_rg2, samples_rg3)
  
  comparisonindex = 1
  for(curi in seq(length(allcomparisons))){
    for(curj in seq(curi, length(allcomparisons))){
      if(curi != curj){
        tryCatch({
          rv = runEdgeR(curdata, cellgroup, unformattedheader, allcomparisons_samples[curi][[1]], allcomparisons_samples[curj][[1]], allcomparisons[curi], allcomparisons[curj], conditionmap)
          summary_comparison[comparisonindex] = rv[1]
          summary_np1[comparisonindex] = rv[2]
          summary_np05[comparisonindex] = rv[3]
          summary_c1samples[comparisonindex] = rv[4]
          summary_c2samples[comparisonindex] = rv[5]
        }, error = function(e){print(e)}
        )
        comparisonindex = comparisonindex+1
      }
    }
  }
  
  allcomparisons = c('Healthy', 'COVID (G2+G3)', 'RSV')
  allcomparisons_samples= list(samples_healthy, samples_covid, samples_rsv)
  
  for(curi in seq(length(allcomparisons))){
    for(curj in seq(curi, length(allcomparisons))){
      if(curi != curj){
        tryCatch({
          rv = runEdgeR(curdata, cellgroup, unformattedheader, allcomparisons_samples[curi][[1]], allcomparisons_samples[curj][[1]], allcomparisons[curi], allcomparisons[curj], conditionmap)
          summary_comparison[comparisonindex] = rv[1]
          summary_np1[comparisonindex] = rv[2]
          summary_np05[comparisonindex] = rv[3]
          summary_c1samples[comparisonindex] = rv[4]
          summary_c2samples[comparisonindex] = rv[5]
        }, error = function(e){print(e)}
        )
        comparisonindex = comparisonindex+1
      }
    }
  }
  
  summary = data.frame("Comparisons"=summary_comparison, "Passes 10% FDR"=summary_np1, "Passes 5% FDR"=summary_np05, "C1 Samples"=summary_c1samples, "C2 Samples"=summary_c2samples)
  
  write.csv(summary, "Summary.csv")
}




for(cursample in rownames(cellnumbers)){
  tryCatch({
    edgermatrix = paste0("/Users/athib/Desktop/CovidRSV/snATAC2/DA/EdgeRMatrixFiles/", cursample, ".txt_edgeRMatrix.txt")
    setwd(paste0("/Users/athib/Desktop/CovidRSV/snATAC2/DA/EdgeRMatrixFiles/", cursample,"/"))
    unformattedheader = read.csv(edgermatrix, header = FALSE, nrows = 1, sep = '\t')
    unformattedheader = unformattedheader[!is.na(unformattedheader)]
    cellgroup = cursample#strsplit(basename(edgermatrix), '\\.')[[1]][1]
    curdata = read.csv(edgermatrix, row.names=1, sep = '\t')
    
    s_h = getFilteredSamples(conditionmap$Sample[conditionmap$Group == 'Healthy'], cursample, mincells = 1, filter = F)
    s_c1 = getFilteredSamples(conditionmap$Sample[conditionmap$Group == 'COVID G1'], cursample, mincells = 1, filter = F)
    s_c2 = getFilteredSamples(conditionmap$Sample[conditionmap$Group == 'COVID G2'], cursample, mincells = 1, filter = F)
    s_c3 = getFilteredSamples(conditionmap$Sample[(conditionmap$Group == 'COVID G3') & (conditionmap$Steroid == FALSE)], cursample, mincells = 1, filter = F)
    s_c3_s = getFilteredSamples(conditionmap$Sample[(conditionmap$Group == 'COVID G3') & (conditionmap$Steroid == TRUE)], cursample, mincells = 1, filter = F)
    
    s_r1 = getFilteredSamples(conditionmap$Sample[conditionmap$Group == 'RSV G1'], cursample, mincells = 1, filter = F)
    s_r2 = getFilteredSamples(conditionmap$Sample[conditionmap$Group == 'RSV G2'], cursample, mincells = 1, filter = F)
    s_r3 = getFilteredSamples(conditionmap$Sample[conditionmap$Group == 'RSV G3'], cursample, mincells = 1, filter = F)
    
    combined = c(s_h, s_c1, s_c2, s_c3, s_c3_s, s_r1, s_r2, s_r3)
    sel_data = curdata[,match(combined, unformattedheader)]
    
    y <- DGEList(counts=sel_data)
    y <- calcNormFactors(y)
    
    logcpm = cpm(y, log=T)
    
    write.table(logcpm, paste0(cursample, "_logcpm.txt"))
  }, error = function(e){print(e)}
  )
  
}




