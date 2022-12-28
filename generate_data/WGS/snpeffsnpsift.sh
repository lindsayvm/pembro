

#single
java -Xmx8g -jar ~/workspace/snpEff/snpEff.jar -v -stats ex1.html GRCh37.75 \
    /DATA/share/Voesties/data/HMF/update_10/somatics/CPCT02010359T/purple/CPCT02010359T.purple.somatic.vcf.gz \
    > ~/workspace/snpEff/data/pembro/test_ann.vcf
java -Xmx4g -jar ~/workspace/snpEff/SnpSift.jar \
    filter "(ANN[*].IMPACT has 'HIGH') || (ANN[*].IMPACT has 'MODERATE') || (ANN[*].EFFECT has 'splice_region_variant')" \
    ~/workspace/snpEff/data/pembro/test_ann.vcf \
    > ~/workspace/snpEff/data/pembro/test_ann_filt.vcf
cat ~/workspace/snpEff/data/pembro/test_ann_filt.vcf | ~/workspace/SnpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/workspace/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT FILTER AF AC DP MQ "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK" > ~/workspace/snpEff/data/pembro/test_ann_filt_oneLine.vcf

#All HMF patient (CPCT WIDE)
for i in $(cat /home/l.leek/pembro/data/pembro_wgs_patientSelection_somatic_CPCTWIDE.txt); do
    java -Xmx8g -jar ~/workspace/snpEff/snpEff.jar -v -stats ex1.html GRCh37.75 \
        /DATA/share/Voesties/data/HMF/update_10/somatics/$i/purple/$i.purple.somatic.vcf.gz \
        > /home/l.leek/pembro/data/snpeff_output/${i}_ann.vcf
    
    java -Xmx4g -jar ~/workspace/snpEff/SnpSift.jar filter \
        "(ANN[*].IMPACT has 'HIGH') || (ANN[*].IMPACT has 'MODERATE') || (ANN[*].EFFECT has 'splice_region_variant')" \
        /home/l.leek/pembro/data/snpeff_output/${i}_ann.vcf \
        > /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt.vcf

    cat /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt.vcf | ~/workspace/SnpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/workspace/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT FILTER AF AC DP MQ PURPLE_AF "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK" > /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt_oneLine.vcf
done

#Additional DRUP patients
for i in $(cat /home/l.leek/pembro/data/pembro_wgs_patientSelection_somatic_DRUP.txt); do
    java -Xmx8g -jar ~/workspace/snpEff/snpEff.jar -v -stats ex1.html GRCh37.75 \
        /DATA/share/Voesties/data/DRUP/update_3/somatics/$i/purple/$i.purple.somatic.vcf.gz \
        > /home/l.leek/pembro/data/snpeff_output/${i}_ann.vcf
    
    java -Xmx4g -jar ~/workspace/snpEff/SnpSift.jar filter \
        "(ANN[*].IMPACT has 'HIGH') || (ANN[*].IMPACT has 'MODERATE') || (ANN[*].EFFECT has 'splice_region_variant')" \
        /home/l.leek/pembro/data/snpeff_output/${i}_ann.vcf \
        > /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt.vcf

    cat /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt.vcf | ~/workspace/SnpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/workspace/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT FILTER AF AC DP MQ PURPLE_AF "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK" > /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt_oneLine.vcf
done





#All HMF patient (CPCT WIDE)
for i in $(cat /home/l.leek/pembro/data/pembro_wgs_patientSelection_somatic_CPCTWIDE.txt); do
    cat /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt.vcf | ~/workspace/SnpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/workspace/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT FILTER AF AC DP MQ PURPLE_AF "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK" > /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt_oneLine.vcf
done

for i in $(cat /home/l.leek/pembro/data/pembro_wgs_patientSelection_somatic_DRUP.txt); do
    cat /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt.vcf | ~/workspace/SnpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/workspace/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT FILTER AF AC DP MQ PURPLE_AF "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK" > /home/l.leek/pembro/data/snpeff_output/${i}_ann_filt_oneLine.vcf
done