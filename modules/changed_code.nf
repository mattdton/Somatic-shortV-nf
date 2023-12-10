
gatk Mutect2 \
             -R $refdir/Homo_sapiens_assembly38.fasta \
             -I ${bam_T} \
             -I ${bam_N} \
             -normal ${bam_id}-N \
             --panel-of-normals ${pon_vcf} \
             --germline-resource $refdir/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz \
             --f1r2-tar-gz ${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz \
             -XL chrM \
             -L "./modules/scatter_files/100M_primary_interval_${splitIntervalNumber}.list" \  
             -O ${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz




gatk  --java-options "-Xmx16g" \
                CalculateContamination \
                -I !{bam_id}-T.pileups.table \
                -tumor-segmentation !{bam_id}-T_segments.table \
                -matched !{bam_id}-N.pileups.table \
                -O !{bam_id}-T_!{bam_id}-N_contamination.table