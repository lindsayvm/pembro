#waar kan ik de informatie vinden op welke positie het centromeer zit?

# aneuploidy score vs TMB
# acrocentric chromosomes: 13,14,15,21,22 dont have a p

#CENTROMEREN CAN BE FOUND HERE
fread("/home/l.leek/pembro/data/chrom_arm_locs.xlsx")


# Library
library(ggplot2)
library(hrbrthemes)

# Create dummy data
data <- data.frame(
  cond = rep(c("condition_1", "condition_2"), each=10), 
  my_x = 1:100 + rnorm(100,sd=9), 
  my_y = 1:100 + rnorm(100,sd=16) 
)

# Basic scatter plot.
p1 <- ggplot(data, aes(x=my_x, y=my_y)) + 
  geom_point( color="#69b3a2") +
  theme_ipsum()

# with linear trend
p2 <- ggplot(data, aes(x=my_x, y=my_y)) +
  geom_point() +
  geom_smooth(method=lm , color="red", se=FALSE) +
  theme_ipsum()

# linear trend + confidence interval
p3 <- ggplot(data, aes(x=my_x, y=my_y)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()


















##############
exit()
purple.qc = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/purple/DRUP01010001T.purple.qc")
purple.qc[purple.qc$QCStatus == "CopyNumberSegments"] #260 CN segments

#purple.cnv has 260  (start end CN)
purple.cnv = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/purple/DRUP01010001T.purple.cnv.somatic.tsv")


cnv.df = purple.cnv %>% 
  mutate(cnv_len = (end-start)) %>% 
  mutate(cnv_weights = cnv_len * copyNumber)
#average cnv per entire genome
avg_cnv = sum(cnv.df$cnv_weights)/sum(cnv.df$cnv_len)


agg_tbl_perChrom = cnv.df %>% group_by(chromosome) %>% 
  dplyr::summarise(sum_cnv_weights = sum(cnv_weights),
            sum_cnv_len = sum(cnv_len),
            .groups = 'drop') %>% 
  mutate(avg_cnv_perChrom = sum_cnv_weights/sum_cnv_len)
#sum of average cnv per chrom
aneuploidy_score = sum(agg_tbl_perChrom$avg_cnv_perChrom)

#sum of average cnv per chrom arm

#obtain centromeric coordinates for hg38:
#   Go to the Table Browser: http://genome.ucsc.edu/cgi-bin/hgTables
# Choose the Mapping and Sequencing group
# Select the "Chromosome Band (Ideogram)" track
# Select filter, and enter "acen" in the gieStain field
# Press "submit" and then "get output"
# Each chromosome will have two entries which overlap.
# They can be simply merged into a single entry.

#I can use centromere positions But that is only under the assumption that the cancer genome length corresponds to reference genome?
#(segmentEndSupport == CENTROMERE and segmentStartSupport CENTROMERE, only once per genome)


#vis_CN has 260 (start end CN)
vis_CN = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/linx/DRUP01010001T.linx.vis_copy_number.tsv")




#only mutation positions, CN (does specify arm)
#purple.driver = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/purple/DRUP01010001T.driver.catalog.somatic.tsv")
#purple.cnvgene = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/purple/DRUP01010001T.purple.cnv.gene.tsv")
#linx.driver = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/linx/DRUP01010001T.linx.driver.catalog.tsv")   #CN per gene

#svId, chrBand (undisruptedCopyNumber, junctionCopyNumber) but no positions
#breakend = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/linx/DRUP01010001T.linx.breakend.tsv")





#Deprecated
##fusion genes
#fusion = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/linx/DRUP01010001T.linx.fusion.tsv")

# #clusterId, arm
# links = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/linx/DRUP01010001T.linx.links.tsv")

#svId, clusterId, vcfId
#svs = fread("/DATA/share/Voesties/data/DRUP/update_5/somatics/DRUP01010001T/linx/DRUP01010001T.linx.svs.tsv")



((3*10)+(4*50)+(3*10))/70



#When the entire short arm of a chromosome is lacking copy number information (and always on chromosome 13,14,15,21, or 22), the copy number of the long arm is extended to the short arm.