import os
import sys
from subprocess import call
import argparse

if __name__ == "__main__":
    samples = {'T2D_1' : ['SRR10247145', 'SRR10247146', 'SRR10247147', 'SRR10247148'],
    'T2D_2' : ['SRR10247149', 'SRR10247150', 'SRR10247151', 'SRR10247152', 'SRR10247153', 'SRR10247154', 'SRR10247155', 'SRR10247156'],
    'T2D_3' : ['SRR10247157', 'SRR10247158', 'SRR10247159', 'SRR10247160'],
    'T2D_4' : ['SRR10247161', 'SRR10247162', 'SRR10247163', 'SRR10247164', 'SRR10247165', 'SRR10247166', 'SRR10247167', 'SRR10247168'],
    'T2D_5' : ['SRR10247169', 'SRR10247170', 'SRR10247171', 'SRR10247172'],
    'T2D_6' : ['SRR10247173', 'SRR10247174', 'SRR10247175', 'SRR10247176', 'SRR10247177', 'SRR10247178', 'SRR10247179', 'SRR10247180'],
    'T2D_7' : ['SRR10247181', 'SRR10247182', 'SRR10247183', 'SRR10247184', 'SRR10247185', 'SRR10247186'],
    'T2D_8' : ['SRR10247187', 'SRR10247188', 'SRR10247189', 'SRR10247190'],
    'Control_1' : ['SRR10247191', 'SRR10247192', 'SRR10247193', 'SRR10247194', 'SRR10247195', 'SRR10247196'],
    'Control_2' : ['SRR10247197', 'SRR10247198', 'SRR10247199', 'SRR10247200'],
    'Control_3' : ['SRR10247201', 'SRR10247202', 'SRR10247203', 'SRR10247204', 'SRR10247205', 'SRR10247206', 'SRR10247207', 'SRR10247208', 'SRR10247209', 'SRR10247210'],
    'Control_4' : ['SRR10247211', 'SRR10247212', 'SRR10247213', 'SRR10247214'],
    'Control_5' : ['SRR10247215', 'SRR10247216', 'SRR10247217', 'SRR10247218', 'SRR10247219', 'SRR10247220', 'SRR10247221', 'SRR10247222'],
    'Control_6' : ['SRR10247223', 'SRR10247224', 'SRR10247225', 'SRR10247226'],
    'Control_7' : ['SRR10247227', 'SRR10247228', 'SRR10247229', 'SRR10247230', 'SRR10247231', 'SRR10247232', 'SRR10247233', 'SRR10247234'],
    'Control_8' : ['SRR10247235', 'SRR10247236', 'SRR10247237', 'SRR10247238'],
    'Control_9' : ['SRR10247239', 'SRR10247240', 'SRR10247241', 'SRR10247242', 'SRR10247243', 'SRR10247244', 'SRR10247245', 'SRR10247246']}

    parser = argparse.ArgumentParser()
    parser.add_argument('--cores_num', help='Number of cores', default=4)
    parser.add_argument('--reference', help='Path to the reference genome folder')
    parser.add_argument('--output_dir', help='Output directory')
    args = parser.parse_args()
    cores_num = int(args.cores_num)

    os.chdir(args.output_dir)

    # Sample processing
    for sample in samples:
        if not os.path.exists(f'{sample}'):
            os.makedirs(f'{sample}')
        srrs = samples[sample]
        os.chdir(f'{sample}')
        cat_files_cmd_1 = "cat "
        cat_files_cmd_2 = "cat "
        for srr in srrs:
            if call(f'fastq-dump --split-files {srr}', shell=True)!=0:
                sys.exit(f"fastq-dump {srr} failed")

            cat_files_cmd_1 += f"{srr}_1.fastq "
            cat_files_cmd_2 += f"{srr}_2.fastq "

        cat_files_cmd_1 += f"> {sample}_1.fastq"
        cat_files_cmd_2 += f"> {sample}_2.fastq"

        if call(cat_files_cmd_1, shell=True)!=0:
            sys.exit(f"cat 1 {sample} failed")

        if call(cat_files_cmd_2, shell=True)!=0:
            sys.exit(f"cat 2 {sample} failed")

        for srr in srrs:
            os.remove(f'{srr}_1.fastq.gz')
            os.remove(f'{srr}_2.fastq.gz')

        if call(f'fastqc {sample}_1.fastq.gz', shell=True)!=0:
            sys.exit(f"FastQC 1 {sample} failed")

        if call(f'fastqc {sample}_2.fastq.gz', shell=True)!=0:
            sys.exit(f"FastQC 2 {sample} failed")

        if call(f'trim_galore --gzip --cores {cores_num} --clip_R1 3 --clip_R2 3 --illumina --paired {sample}_1.fastq.gz {sample}_2.fastq.gz', shell=True)!=0:
            sys.exit(f"trim_galore {sample} failed")

        os.remove(f'{sample}_1.fastq.gz')
        os.remove(f'{sample}_2.fastq.gz')

        if call(f'fastqc {sample}_1_val_1.fq.gz', shell=True)!=0:
            sys.exit(f"FastQC 1 trimmed {sample} failed")

        if call(f'fastqc {sample}_2_val_2.fq.gz', shell=True)!=0:
            sys.exit(f"FastQC 2 trimmed {sample} failed")

        if call(f'bismark -p {cores_num} --gzip {args.reference} -1 {sample}_1_val_1.fq.gz -2 {sample}_2_val_2.fq.gz', shell=True)!=0:
            sys.exit(f"bismark {sample} failed")

        os.remove(f'{sample}_1_val_1.fq.gz')
        os.remove(f'{sample}_2_val_2.fq.gz')

        if call(f'deduplicate_bismark {sample}_1_val_1_bismark_bt2_pe.bam', shell=True)!=0:
            sys.exit(f"deduplicate_bismark {sample} failed")

        if call(f'bismark2report', shell=True)!=0:
            sys.exit(f"bismark2report failed {sample}")

        if call(f'bismark2summary', shell=True)!=0:
            sys.exit(f"bismark2summary failed {sample}")

        os.remove(f'{sample}_1_val_1_bismark_bt2_pe.bam')

        if call(f'bismark_methylation_extractor --bedGraph --parallel {cores_num}  --gzip {sample}_1_val_1_bismark_bt2_pe.deduplicated.bam', shell=True)!=0:
            sys.exit(f"bismark_methylation_extractor failed {sample}")

        if call(f'rm *.bedGraph.gz', shell=True)!=0:
            sys.exit(f"rm *.bedGraph.gz fail")

        if call(f'rm CHG_*', shell=True)!=0:
            sys.exit(f"rm CHG_* fail")
        if call(f'rm CHH_*', shell=True)!=0:
            sys.exit(f"rm CHH_* fail")
        if call(f'rm CpG_*', shell=True)!=0:
            sys.exit(f"rm CpG_* fail")

        os.remove(f'{sample}_1_val_1_bismark_bt2_pe.deduplicated.bam')

        if call(f'mv -f {sample}_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz ../{sample}.cov.gz', shell=True)!=0:
            sys.exit(f"move {sample}.cov failed")

        os.chdir('..')

        if call(f'gunzip -f {sample}.cov.gz', shell=True)!=0:
            sys.exit(f"gunzip {sample}.cov failed")