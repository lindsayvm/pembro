library(data.table)

#needed: TML, TMB, SV, aneuploidy, clonal

fread("~/pembro_erasmus/CORELR020020T.purple.purity.tsv") #TML TMB 
fread("~/pembro_erasmus/CORELR020002T.purple.purity.tsv") #-
fread("~/pembro_erasmus/CORELR020047T.purple.purity.tsv") #TML TMB
fread("/shared/data/IID/pembro_EMC/CORELR020046T.purple.purity.tsv") #TML, TMB, SV

#calc TML TMB aneuploidy and clonal from somatic vcf files. But structural variants???
fread("~/pembro_erasmus/CORELR020020T.purple.somatic.vcf.gz") #TML TMB 
fread("~/pembro_erasmus/CORELR020002T.purple.somatic.vcf.gz") #-
fread("~/pembro_erasmus/CORELR020047T.purple.somatic.vcf.gz") #TML TMB
fread("/shared/data/IID/pembro_EMC/CORELR020046T.purple.somatic.vcf.gz") #TML, TMB, SV



