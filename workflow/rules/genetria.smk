rule genetria_concat_genome:
    """
    Concatenate the whole-genome autosomal VCF from QUILT2 and 
    the ligated Chromosome X BCF from GLIMPSE2 into a unified file.
    """
    input:
        quilt=os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "quilt.down{depth}x.mspbwt.genome.vcf.gz"
        ),
        glimpse=os.path.join(
            OUTDIR_GLIMPSE2,
            "refsize{size}",
            "chrX",
            "down{depth}x.chrX.bcf"
        )
    output:
        vcf=os.path.join(
            OUTDIR_GENETRIA,
            "refsize{size}",
            "genetria.down{depth}x.genome.vcf.gz"
        ),
        tbi=os.path.join(
            OUTDIR_GENETRIA,
            "refsize{size}",
            "genetria.down{depth}x.genome.vcf.gz.tbi"
        )
    log:
        os.path.join(
            OUTDIR_GENETRIA,
            "refsize{size}",
            "genetria.down{depth}x.genome.vcf.gz.log"
        )
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        (
        echo "Concatenating Autosomes and Chromosome X..."
        bcftools concat \
            --threads 8 \
            -a \
            -O z \
            -o {output.vcf} \
            {input.quilt} \
            {input.glimpse}

        echo "Indexing unified genome..."
        bcftools index -t {output.vcf}
        ) &>> {log}
        """


rule genetria_split_by_sample:
    """
    Extract individual sample VCFs from the unified mixed-ploidy genome.
    """
    input:
        vcf=rules.genetria_concat_genome.output.vcf,
        tbi=rules.genetria_concat_genome.output.tbi
    output:
        vcf=os.path.join(
            OUTDIR_GENETRIA,
            "refsize{size}",
            "split_samples",
            "{sample}.vcf.gz"
        ),
        tbi=os.path.join(
            OUTDIR_GENETRIA,
            "refsize{size}",
            "split_samples",
            "{sample}.vcf.gz.tbi"
        )
    log:
        os.path.join(
            OUTDIR_GENETRIA,
            "refsize{size}",
            "split_samples",
            "{sample}.split.log"
        )
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        (
        echo "Extracting sample {wildcards.sample} from hybrid genome..."
        bcftools view \
            -s {wildcards.sample} \
            --threads 4 \
            -O z \
            -o {output.vcf} \
            {input.vcf}

        echo "Indexing the sample VCF..."
        bcftools index -t {output.vcf}
        ) &>> {log}
        """


rule genetria_merge_historic:
    input:
        new=rules.genetria_concat_genome.output.vcf,
        historic=config["vcf_qc"]["historic_vcf"]
    output:
        merged_vcf=config["vcf_qc"]["historic_vcf"].replace(".vcf.gz", "_genetria_{size}_{depth}_updated.vcf.gz"),
        indexed_merge=config["vcf_qc"]["historic_vcf"].replace(".vcf.gz", "_genetria_{size}_{depth}_updated.vcf.gz.tbi")
    log:
        config["vcf_qc"]["historic_vcf"].replace(".vcf.gz", "_genetria_{size}_{depth}_updated.log")
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        (
        echo "Merging hybrid pipeline results with historic samples..."
        bcftools merge --force-samples {input.new} {input.historic} -O z -o {output.merged_vcf}
        bcftools index -t {output.merged_vcf}
        ) &>> {log}
        """


rule genetria_stats_for_het_homalt:
    input:
        vcf=rules.genetria_merge_historic.output.merged_vcf
    output:
        stats=os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x_sample_stats.txt"),
        ratios=os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x_clean_ratios.tsv")
    log:
        os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x_clean_ratios.tsv.log")
    conda:
        "../envs/quilt.yaml"
    shell:
        r"""
        (
        echo "Calculating mixed-run allele metrics..."
        bcftools stats -s - {input.vcf} > {output.stats}

        echo -e "Sample_ID\tnRefHom\tnHomAlt(1/1)\tHeterozygous(0/1)\tRatio(Het/Hom-Alt)" > {output.ratios}
        
        grep "^PSC" {output.stats} | \
        awk -F'\t' '{
            if ($5 == 0) { ratio = 0 } else { ratio = $6/$5 }
            printf "%s\t%s\t%s\t%s\t%.2f\n", $3, $4, $5, $6, ratio
        }' >> {output.ratios}
        ) &>> {log}
        """


rule genetria_pruning_and_pca:
    input:    
        vcf=rules.genetria_merge_historic.output.merged_vcf,
        prune_in=config["vcf_qc"]["prune_in"],
        afreq=config["vcf_qc"]["afreq"]
    output:
        eigenvec=os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x_clean_pca.eigenvec"),
        eigenval=os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x_clean_pca.eigenval"),
        kinship=os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x_clean_kinship.kin0")
    params:
        prefix=os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x")
    log:
        os.path.join(OUTDIR_GENETRIA, "refsize{size}", "qcs", f"{config['run_name']}_down{{depth}}x_clean_pca.log")
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        (
        plink2 --vcf {input.vcf} \
            --set-all-var-ids '@:#:$r:$a' \
            --make-pgen \
            --out {params.prefix}_out_temp \
            --autosome
        
        plink2 --pfile {params.prefix}_out_temp \
            --extract {input.prune_in} \
            --make-pgen \
            --out {params.prefix}_clean_batch_for_pca

        plink2 --pfile {params.prefix}_clean_batch_for_pca --rm-dup exclude-all --make-pgen --out {params.prefix}_clean_batch_dedup
        plink2 --pfile {params.prefix}_clean_batch_dedup --read-freq {input.afreq} --pca 10 --out {params.prefix}_clean_pca
        plink2 --pfile {params.prefix}_clean_batch_dedup --make-king-table --out {params.prefix}_clean_kinship
        ) &>> {log}
        """


rule genetria_render_qc_report:
    input:
        kinship=rules.genetria_pruning_and_pca.output.kinship,
        eigenval=rules.genetria_pruning_and_pca.output.eigenval,
        eigenvec=rules.genetria_pruning_and_pca.output.eigenvec,
        ratios=rules.genetria_stats_for_het_homalt.output.ratios
    output:
        report=os.path.join(OUTDIR_GENETRIA, "refsize{size}", f"{config['run_name']}_down{{depth}}x_QC_Report.html")
    params:
        batch=config['run_name']
    log:
        os.path.join(OUTDIR_GENETRIA, "refsize{size}", f"{config['run_name']}_down{{depth}}x_QC_Report.log")
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/imputation_genotyping_qcs.Rmd"