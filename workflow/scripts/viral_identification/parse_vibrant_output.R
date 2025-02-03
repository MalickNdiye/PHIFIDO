library(data.table)
library(optparse)
library(tidyverse)

# Aguments parser
option_list <- list ( make_option (c("-d","--dir"), type="character",
                                   help="vibrant directory", metavar="character"),
                      make_option (c("-o","--outtab"), type="character",
                                   help="vibrant output table", metavar="character")
)

parser <-OptionParser(option_list=option_list)
args <- parse_args(parser)

# setup, input directory
dir<- args$dir
print(dir)
sample<- basename(args$dir)


# Set path to relevant files
contig_quality<- file.path(dir, "VIBRANT_final.contigs/VIBRANT_results_final.contigs/VIBRANT_genome_quality_final.contigs.tsv")
if (!file.exists(contig_quality)){
  final_df<- data.frame("sample"=NULL, "scaffold"=NULL,"new_scaffold_name"=NULL,
  "type"=NULL,"vibrant_quality"=NULL, "total_genes"=NULL, "length" =NULL, "multiplicity"=NULL)

  write.table(final_df, args$outtab, quote=F, row.names = F, sep="\t" )
  quit(save = "no", status = 0)
} 
contig_genes<- file.path(dir, "/VIBRANT_final.contigs/VIBRANT_results_final.contigs/VIBRANT_summary_results_final.contigs.tsv")


# Parse quality
score<-c("low quality draft"=1,  "medium quality draft"=2, "high quality draft"=3, "complete circular"=4)
quality_df<- fread(contig_quality)%>%
  group_by(scaffold) %>%          
  slice_max(score[Quality]) %>%      
  ungroup()  %>%
  rename("vibrant_quality"=Quality)

if (nrow(quality_df)<1){
  final_df<- data.frame("sample"=NULL, "scaffold"=NULL,"new_scaffold_name"=NULL,
  "type"=NULL,"vibrant_quality"=NULL, "total_genes"=NULL, "length" =NULL, "multiplicity"=NULL)

  write.table(final_df, args$outtab, quote=F, row.names = F, sep="\t" )
  quit(save = "no", status = 0)
} 


# parse gene count
gene_df<- fread(contig_genes)%>%
  select(scaffold, "total genes")%>%
  rename("total_genes"="total genes")

# parse_scaffold
parse_megahit_names<- function(x){
  scaffold_l<- unlist(strsplit(x, split=" "))
  
  len=as.numeric(gsub("len=", "", scaffold_l[4]))
  cov=as.numeric(gsub("multi=", "", scaffold_l[3]))
  
  output=c("len"=len, "cov"=cov)
  
  return(output)
}


# put everything together
final_df<- quality_df %>% left_join(., gene_df) %>%
  mutate(sample=sample,
         length=unlist(lapply(scaffold, function(x) parse_megahit_names(x)["len"])),
         multiplicity=unlist(lapply(scaffold, function(x) parse_megahit_names(x)["cov"])))%>%
  filter(length>5000)%>%
  mutate(new_scaffold_name=paste(sample, "_viral_contig_", row_number(), sep=""))%>%
  relocate(sample) %>%
  relocate(new_scaffold_name, .after = scaffold)%>%
  arrange(-length)

# write table to file
write.table(final_df, args$outtab, quote=F, row.names = F, sep="\t" )
