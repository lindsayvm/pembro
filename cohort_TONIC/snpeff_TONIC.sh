#ls >> tonic_wgs_patientSelection.txt
#cd workspace/snpeff/

for i in $(cat ~/TONIC/VCF_files/tonic_wgs_patientSelection.txt); do
    echo $i
    echo ${i/%.vcf/_ann.vcf}
    echo ${i/%.vcf/_ann_filt.vcf}
    echo ${i/%.vcf/_ann_filt_oneLine.vcf}
done


i=m2_5077_1_CF15185_GCGAGTAA_vs_4878_23_CF16646_CAGCGTTA-combined.annotated.vcf
java -Xmx8g -jar ~/workspace/snpEff/snpEff.jar -v -stats ex1.html GRCh38.99 \
    ~/TONIC/VCF_files/$i \
    > ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann.vcf}
java -Xmx4g -jar ~/workspace/snpEff/SnpSift.jar filter \
    "(ANN[*].IMPACT has 'HIGH') || (ANN[*].IMPACT has 'MODERATE') || (ANN[*].EFFECT has 'splice_region_variant')" \
    ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann.vcf} \
    > ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann_filt.vcf}
cat ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann_filt.vcf} | ~/workspace/SnpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/workspace/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT FILTER AF AC DP MQ "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK"> ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann_filt_oneLine.vcf}

for i in $(cat ~/TONIC/VCF_files/tonic_wgs_patientSelection.txt); do
    java -Xmx8g -jar ~/workspace/snpEff/snpEff.jar -v -stats ex1.html GRCh38.99 \
        ~/TONIC/VCF_files/$i \
        > ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann.vcf}
    
    java -Xmx4g -jar ~/workspace/snpEff/SnpSift.jar filter \
        "(ANN[*].IMPACT has 'HIGH') || (ANN[*].IMPACT has 'MODERATE') || (ANN[*].EFFECT has 'splice_region_variant')" \
        ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann.vcf} \
        > ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann_filt.vcf}

    cat ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann_filt.vcf} | ~/workspace/SnpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/workspace/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT FILTER AF AC DP MQ "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "EFF[*].TRID" "EFF[*].RANK"> ~/TONIC/vcf_snpeffsift/${i/%.vcf/_ann_filt_oneLine.vcf}
done
