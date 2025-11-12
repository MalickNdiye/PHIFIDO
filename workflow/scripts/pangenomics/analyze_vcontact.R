library("optparse")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("data.table")
library("GGally")
library("igraph")
library("ggforce")

# this command:"Rscript scripts/pangenomics/analyze_vcontact.R -i {input.vcontact}/c1.ntw -d {input.drep}/dRep_summary_average_updated.tsv -o {output}
# parse the arguments using optparse
option_list = list(
  make_option(c("-i", "--input_vcontact"), type="character", default=NULL,                       
                help="Path to the vContact2 output directory containing the c1.ntw file", metavar="character"), 
    make_option(c("-d", "--input_comm"), type="character", default=NULL,                       
                help="Path to the filtered community data", metavar="character"), 
    make_option(c("-g", "--genome_info"), type="character", default=NULL,                       
                help="Path to the list of genomes to consider", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,                       
                help="Path to the output directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read the dRep summary_average_updated.tsv file
vircom <- fread(opt$input_comm) 

present_genomes<- unique(vircom$genome)
present_vOTUs<- unique(vircom$vOTU)

# genome info
genome_info <- fread(opt$genome_info) 

genome_info_filt <- genome_info %>%
  filter(vOTU %in% present_vOTUs) %>%
  group_by(vOTU) %>%
  reframe(genomes = list(unique(genome)), 
          host_genus = collapse(unique(host_genus), sep = "; "),
          class=unique(class),
          lifestyle=ifelse(all(lifestyle=="Lytic"), "Lytic", "Temperate"),
          origin="This study")

# read the vContact2 c1.ntw file
vcontact_file <- file.path(opt$input_vcontact)
ntw<- fread(vcontact_file, header=FALSE, sep="\t", col.names = c("source", "target", "edge"))
  filter(!grepl("~", source) | !grepl("~", target)) %>%
  mutate(origin_source = ifelse(source %in% genome_info$genome, "This study", "Ndiaye et al. (2025)"),
         origin_source = ifelse(grepl("MT00", source), "Bonilla-Rosso et al. (2020)", origin_source),
          origin_target = ifelse(target %in% genome_info$genome, "This study", "Ndiaye et al. (2025)"),
          origin_target = ifelse(grepl("MT00", target), "Bonilla-Rosso et al. (2020)", origin_target)
) %>%
  filter(!(origin_source =="This study" & !source %in% present_genomes)) %>%
  filter(!(origin_target =="This study" & !target %in% present_genomes))

# add all genomes of Ndiaye et al. (2025) and Bonilla-Rosso et al. (2020) to genome_info_filt
additional_genomes <- unique(c(
  ntw$source[ntw$origin_source != "This study"],
  ntw$target[ntw$origin_target != "This study"]
]))

# create dataframe for additional genomes, all the info except genome and origin will be NA
additional_genome_info <- data.frame(
  vOTU = NA,
  genomes = as.list(additional_genomes),
  host_genus = NA,
  class = NA,
  lifestyle = NA,
  origin = ifelse(grepl("MT00", additional_genomes), "Bonilla-Rosso et al. (2020)", "Ndiaye et al. (2025)")
)

genome_info_complete <- bind_rows(genome_info_filt, additional_genome_info) 

# create an object that contains the network along with the genome information
ntw_info <- ntw %>%
  left_join(genome_info_complete %>% unnest(cols = c(genomes)) %>%
              rename(source = genomes,
                      vOTU_source = vOTU,
                      host_genus_source = host_genus,
                      class_source = class,
                      lifestyle_source = lifestyle,
                      origin_source_info = origin),
            by = "source") %>%
  left_join(genome_info_complete %>% unnest(cols = c(genomes)) %>
              rename(target = genomes,
                      vOTU_target = vOTU,
                      host_genus_target = host_genus,
                      class_target = class,
                      lifestyle_target = lifestyle,
                      origin_target_info = origin),
            by = "target")

# write the output file
fwrite(ntw_info, file.path(opt$output, "vcontact_network_annotated.tsv"), sep="\t")

# Create a ggnet2 object for plotting
ntwk<- ntw_info %>% select(source, target, edge)
nodes <- GGally::ggnet2(ntwk[,-3], 
                mode = "fruchtermanreingold", 
                layout.par = list(list=(niter=2000))) %>% 
  .$data %>% 
  dplyr::rename("Genome" = "label") %>%
  left_join(genome_info_complete %>% unnest(cols = c(genomes)) %>%
              rename(Genome = genomes),
            by = "Genome")

edges <- ntwk %>% 
  mutate(Pair = paste(source, target, sep = ".")) %>% 
  gather(key = "Member", value = "Genome", -Pair, -edge) %>% 
  inner_join(nodes, by = "Genome") 

vcontact_ntw_plt_origin<- vcont_df %>% 
  ggplot() +
  geom_line(data = edges, aes(x, y, group = Pair), alpha = 0.1, size = 0.1) +
  geom_point(data = nodes, aes(x, y, color = origin, shape=host_genus), alpha = 0.8, size = 4) +
  ggforce::geom_circle(data = circle, aes(x0 = 0.5, y0 = 0.5, r = 0.55), color = "gray25", linetype = 2) +
  guides(color = guide_legend(title.hjust = 1, nrow = 3)) +
  theme(text = element_text(size = 12),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        legend.text=element_text(face="bold", size=20),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.title = element_blank())

# plot only This study genomes colored by lifestyle
vcontact_ntw_plt_lifestyle<- vcont_df %>%
  ggplot() +
  geom_line(data = edges, aes(x, y, group = Pair), alpha = 0.1, size = 0.1) +
  geom_point(data = nodes %>% filter(origin == "This study"), 
             aes(x, y, color = lifestyle), alpha = 0.8, size = 4, shape = 16) +
  ggforce::geom_circle(data = circle, aes(x0 = 0.5, y0 = 0.5, r = 0.55), color = "gray25", linetype = 2) +
  guides(color = guide_legend(title.hjust = 1, nrow = 3)) +
  theme(text = element_text(size = 12),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        legend.text=element_text(face="bold", size=20),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.title = element_blank())
  
# plot only This study genomes colored by host_genus
vcontact_ntw_plt_hostgenus<- vcont_df %>%
  ggplot() +
  geom_line(data = edges, aes(x, y, group = Pair), alpha = 0.1, size = 0.1) +
  geom_point(data = nodes %>% filter(origin == "This study"), 
             aes(x, y, color = host_genus), alpha = 0.8, size = 4, shape = 16) +
  ggforce::geom_circle(data = circle, aes(x0 = 0.5, y0 = 0.5, r = 0.55), color = "gray25", linetype = 2) +
  guides(color = guide_legend(title.hjust = 1, nrow = 3)) +
  theme(text = element_text(size = 12),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        panel.grid = element_blank(),
        legend.text=element_text(face="bold", size=20),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",                     
        legend.title = element_blank())

# save all the plots in the output directory in a folder called plots
ggsave(file.path(opt$output, "plots", "vcontact_network_origin.png"), vcontact_ntw_plt_origin, width = 10, height = 8)
ggsave(file.path(opt$output, "plots", "vcontact_network_lifestyle.png"), vcontact_ntw_plt_lifestyle, width = 10, height = 8)
ggsave(file.path(opt$output, "plots", "vcontact_network_hostgenus.png"), vcontact_ntw_plt_hostgenus, width = 10, height = 8)