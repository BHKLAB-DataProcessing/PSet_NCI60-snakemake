from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]

mol_data_source = "https://discover.nci.nih.gov/cellminerdata/normalizedArchives/"
sens_data_source = "https://zenodo.org/record/6678327/files/"
annotation = "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/"

rule get_PSet:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "moldata/RNA__Affy_HG_U133_Plus_2.0_RMA.xls"),
        S3.remote(prefix + "moldata/RNA__Agilent_Human_microRNA_V2_GeneSpringGX.xls"),
        S3.remote(prefix + "moldata/RNA__RNA_seq_composite_expression.xls"),
        S3.remote(prefix + "moldata/RNA__RNA_seq_isoforms.xls"),
        S3.remote(prefix + "celldata/NCI60_CELL_LINE_METADATA.txt"),
        S3.remote(prefix + "sensdata/DOSERESP.csv"),
        S3.remote(prefix + "sensdata/profile.sensitivity_v3.rds"),
        S3.remote(prefix + "sensdata/raw.sensitivity_v3.rds"),
        S3.remote(prefix + "annotation/drugs_with_ids.csv"),
        S3.remote(prefix + "annotation/cell_annotation_all.csv")
    resources:
        mem_mb=5000
    shell:
        "Rscript scripts/getNCI60.R {prefix}"

rule extract_mol_data:
    output:
        S3.remote(prefix + "moldata/RNA__Affy_HG_U133_Plus_2.0_RMA.xls"),
        S3.remote(prefix + "moldata/RNA__Agilent_Human_microRNA_V2_GeneSpringGX.xls"),
        S3.remote(prefix + "moldata/RNA__RNA_seq_composite_expression.xls"),
        S3.remote(prefix + "moldata/RNA__RNA_seq_isoforms.xls")
    input:
        S3.remote(prefix + "download/nci60_RNA__Affy_HG_U133_Plus_2.0_RMA.zip"),
        S3.remote(prefix + "download/nci60_RNA__Agilent_Human_microRNA_V2_GeneSpringGX.zip"),
        S3.remote(prefix + "download/nci60_RNA__RNA_seq_composite_expression.zip"),
        S3.remote(prefix + "download/nci60_RNA__RNA_seq_isoforms.zip")
    resources:
        mem_mb=2000
    shell:
        """
        unzip -d {prefix}download/rna_affy/ {prefix}download/nci60_RNA__Affy_HG_U133_Plus_2.0_RMA.zip && \
            mv {prefix}download/rna_affy/output/RNA__Affy_HG_U133_Plus_2.0_RMA.xls {prefix}moldata/RNA__Affy_HG_U133_Plus_2.0_RMA.xls && \
            rm -rf {prefix}download/rna_affy
        unzip -d {prefix}download/micro_rna/ {prefix}download/nci60_RNA__Agilent_Human_microRNA_V2_GeneSpringGX.zip && \
            mv '{prefix}download/micro_rna/output/RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.xls' {prefix}moldata/RNA__Agilent_Human_microRNA_V2_GeneSpringGX.xls && \
            rm -rf {prefix}download/micro_rna
        unzip -d {prefix}download/rna_composite/ {prefix}download/nci60_RNA__RNA_seq_composite_expression.zip && \
            mv {prefix}download/rna_composite/output/RNA__RNA_seq_composite_expression.xls {prefix}moldata/RNA__RNA_seq_composite_expression.xls && \
            rm -rf {prefix}download/rna_composite
        unzip -d {prefix}download/rna_isoform/ {prefix}download/nci60_RNA__RNA_seq_isoforms.zip && \
            mv {prefix}download/rna_isoform/output/RNA__RNA_seq_isoforms.xls {prefix}moldata/RNA__RNA_seq_isoforms.xls && \
            rm -rf {prefix}download/rna_isoform
        """

rule download_sensdata:
    output:
        S3.remote(prefix + "sensdata/DOSERESP.csv"),
        S3.remote(prefix + "sensdata/profile.sensitivity_v3.rds"),
        S3.remote(prefix + "sensdata/raw.sensitivity_v3.rds")
    resources:
        mem_mb=3000
    shell:
        """
        wget {sens_data_source}DOSERESP.csv?download=1 -O {prefix}sensdata/DOSERESP.csv
        wget {sens_data_source}profile.sensitivity_v3.rds?download=1 -O {prefix}sensdata/profile.sensitivity_v3.rds
        wget {sens_data_source}raw.sensitivity_v3.rds?download=1 -O {prefix}sensdata/raw.sensitivity_v3.rds
        """

rule download_data:
    output:
        S3.remote(prefix + "download/nci60_RNA__Affy_HG_U133_Plus_2.0_RMA.zip"),
        S3.remote(prefix + "download/nci60_RNA__Agilent_Human_microRNA_V2_GeneSpringGX.zip"),
        S3.remote(prefix + "download/nci60_RNA__RNA_seq_composite_expression.zip"),
        S3.remote(prefix + "download/nci60_RNA__RNA_seq_isoforms.zip"),
        S3.remote(prefix + "celldata/NCI60_CELL_LINE_METADATA.txt")
    resources:
        mem_mb=2000
    shell:
        """
        wget {mol_data_source}nci60_RNA__Affy_HG_U133_Plus_2.0_RMA.zip -O {prefix}download/nci60_RNA__Affy_HG_U133_Plus_2.0_RMA.zip
        wget '{mol_data_source}nci60_RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.zip' -O {prefix}download/nci60_RNA__Agilent_Human_microRNA_V2_GeneSpringGX.zip
        wget {mol_data_source}nci60_RNA__RNA_seq_composite_expression.zip -O {prefix}download/nci60_RNA__RNA_seq_composite_expression.zip
        wget {mol_data_source}nci60_RNA__RNA_seq_isoforms.zip -O {prefix}download/nci60_RNA__RNA_seq_isoforms.zip
        wget https://discover.nci.nih.gov/cellminer/samples/NCI60_CELL_LINE_METADATA.txt -O {prefix}/celldata/NCI60_CELL_LINE_METADATA.txt
        """ 

rule get_annotation:
    output:
        S3.remote(prefix + "annotation/drugs_with_ids.csv"),
        S3.remote(prefix + "annotation/cell_annotation_all.csv")
    shell:
        """
        wget -O {prefix}annotation/drugs_with_ids.csv {annotation}drugs_with_ids.csv?raw=true
        wget -O {prefix}annotation/cell_annotation_all.csv {annotation}cell_annotation_all.csv?raw=true
        """