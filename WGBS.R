#!/usr/bin/env Rscript
library("ggplot2")
library("optparse")
library("GenomicRanges")



###Parsing argumets
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="[REQUIRED] WGBS file (bedGraph)", metavar="character"),
  make_option(c("-g", "--genome (hg19 or mm10)"), type="character", default=NULL,
              help="[REQUIRED] build genome annotation", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="[REQUIRED] annotation file (bed)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="out.txt",
              help="[REQUIRED] output directory path [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);

opt = parse_args(opt_parser);


if (is.null(opt$file)) {
  stop("WGBS file needed.")
}
if (is.null(opt$genome)) {
  stop("Build genome annotation needed.")
}
if (is.null(opt$annotation)) {
  stop("Annotation file needed.")
}
if (is.null(opt$outdir)) {
  stop("Output file needed.")
}


file = opt$file
genome = opt$genome
annotation = opt$annotation
outdir = opt$outdir

if (genome != "hg19" | genome != "mm10"){
  stop("Please, only use hg19 or mm10.")
}


###loading annotations from annotatr package
annotation = read.table(annotation, sep= " ", header = TRUE)
type = genome


feature_list_genes = c( paste(type,"genes_1to5kb", sep="_"),
                        paste(type,"genes_promoters", sep="_"),
                        paste(type,"genes_5UTRs", sep="_"),
                        paste(type,"genes_firstexons", sep="_"),
                        paste(type,"genes_exons", sep="_"),
                        paste(type,"genes_introns", sep="_"),
                        paste(type,"genes_3UTRs", sep="_"))

feature_list_inter = c(paste(type,"cpg_islands", sep="_"),
                       paste(type,"cpg_shores", sep="_"),
                       paste(type,"cpg_shelves", sep="_"),
                       paste(type,"enhancers_fantom", sep="_"))

feature_annotation_list = c(unlist(feature_list_genes), unlist(feature_list_inter))

process_analysis <- function(file, CpG_GR, anno, feature_list, selected_chrm, outfolder){
  
  outfolder = paste(outfolder ,gsub('_CpG.bedGraph','', basename(file)), sep ="/")
  dir.create(outfolder)
  
  
  if (selected_chrm == "autosomal"){
    CpG_GR = CpG_GR[seqnames(CpG_GR) != "chrX" & seqnames(CpG_GR) != "chrY",]
  } else if (selected_chrm == "chrX"){
    CpG_GR = CpG_GR[seqnames(CpG_GR) == "chrX",]
  } else if (selected_chrm == "chrY"){
    CpG_GR = CpG_GR[seqnames(CpG_GR) == "chrY",]
  } 
  
  
  ###Process the analysis
  #Object "GRangesObject" with two attributes : a name and a GRange object
  GRangesObject <- function(feature = NULL, GRange = NULL){
    value <- list(name = feature, GR = GRange)
    attr(value, "class") <- "GRangesObject"
    value
  }
  
  #Return the CpGs overlapping the type
  count_overlap <- function(CpG, feature = NULL, anno){
    
    if (!is.null(feature)){
      type_annotation = anno[anno$type == feature,]
    }
    else {
      type_annotation = anno
    }
    
    type_annotation_GR = GRanges(seqnames = type_annotation$seqnames,
                                 ranges = IRanges(as.numeric(type_annotation$start),
                                                  end = as.numeric(type_annotation$end)))
    
    hits_CpG_type = findOverlaps(CpG, type_annotation_GR)
    
    return(unique(CpG[queryHits(hits_CpG_type),]))
  }
  
  #List of GRangesObject. Each GRangesObject contains the name of the analysed type and the CpG overlapping this region
  #print("processing the analysis")
  list_gr = list(GRangesObject())
  for (i in 1:length(feature_list)){
    gr = count_overlap(CpG_GR, feature_list[i], anno)
    list_gr[[i]] = GRangesObject(feature_list[i], gr)
  }
  
  
  #print("creating the table")
  #Creating a table
  feature_list = c(unlist(feature_list), paste(type, "part_of_intergenic", sep="_"), paste(type, "genome_wide", sep="_"))
  usedtype = unlist(strsplit(feature_list, "_"))[1]
  table_res = data.frame(Region = feature_list, Freq = 0, Percentage = 0, MeanBetaValue = 0)
  total_probes_in_regions = GRanges()
  
  for (i in 1:length(list_gr)){
    name_feature = list_gr[[i]]$name
    GR = list_gr[[i]]$GR
    
    table_res[table_res$Region == name_feature,]$Freq = length(GR)
    table_res[table_res$Region == name_feature,]$MeanBetaValue = mean(GR$score)
    
    list = GRangesList(total_probes_in_regions, list_gr[[i]]$GR)
    list2 = unlist(list)
    total_probes_in_regions = unique(list2)
  }
  table_res$Region = gsub(paste0(usedtype,"_"), "", table_res$Region)
  table_res$Region = gsub("genes_", "", table_res$Region)
  
  table_res[table_res$Region =="part_of_intergenic",]$Freq = length(CpG_GR) - length(total_probes_in_regions)
  table_res[table_res$Region =="part_of_intergenic",]$MeanBetaValue = mean(CpG_GR[!(CpG_GR %in% total_probes_in_regions),]$score)
  
  table_res$Percentage = table_res$Freq/length(CpG_GR)
  
  table_res[table_res$Region =="genome_wide",]$Freq = length(CpG_GR)
  table_res[table_res$Region =="genome_wide",]$MeanBetaValue = mean(CpG_GR$score)
  table_res[table_res$Region =="genome_wide",]$Percentage = NA
  
  table_res$Region = factor(table_res$Region, as.character(table_res$Region))
  
  
  
  #saving results
  prefix = gsub('_CpG.bedGraph','', basename(file))
  
  plot_perc <- ggplot(table_res[1:12,], aes(x=Region, y=100*Percentage, fill=Region)) + ylim(0,100) +
    geom_bar(stat="identity")+theme(axis.text.x=element_text(angle = 45, vjust = 1, 
                                                             size = 12, hjust = 1),
                                    legend.position = "none") + ylab("Percentage") +
    geom_text(aes(label=sprintf("%.2f", 100*Percentage) ), vjust=-0.3, size=3.5)
  print(plot_perc)
  
  
  plot_meanValue <- ggplot(table_res, aes(x=Region, y=MeanBetaValue, fill=Region)) + ylim(0,100) +
    geom_bar(stat="identity")+theme(axis.text.x=element_text(angle = 45, vjust = 1, 
                                                             size = 12, hjust = 1),
                                    legend.position = "none") +
    geom_text(aes(label=sprintf("%.2f", MeanBetaValue) ), vjust=-0.3, size=3.5)
  
  print(plot_meanValue)
  
  file_perc_png = sprintf('%s/%s_CpG_perc_%s.png', outfolder, prefix, selected_chrm)
  ggsave(filename = file_perc_png, plot = plot_perc, width = 8, height = 6)
  
  file_mean_png = sprintf('%s/%s_CpG_meanBetaValue_%s.png', outfolder, prefix, selected_chrm)
  ggsave(filename = file_mean_png, plot = plot_meanValue, width = 8, height = 6)
  
  print(file_perc_png)
  print(file_mean_png)
  print(sprintf('%s/%s_CpG_%s.txt', outfolder, prefix, selected_chrm))
  print(table_res)
  write.table(table_res, file = sprintf('%s/%s_CpG_%s.txt', outfolder, prefix, selected_chr))
  
}

###Creating GenomicRanges object from a bedGraph
data_bedGraph = read.table(file, header = FALSE, sep = "\t", comment.char = "#")
allCpG = GRanges(seqnames = data_bedGraph$V1,
                 ranges = IRanges(as.numeric(data_bedGraph$V2),
                                  end = as.numeric(data_bedGraph$V3)),
                 score = data_bedGraph$V4)
process_analysis(file, allCpG, annotation, feature_annotation_list, "autosomal", outdir) 
process_analysis(file, allCpG, annotation, feature_annotation_list, "chrX", outdir)
process_analysis(file, allCpG, annotation, feature_annotation_list, "chrY", outdir)