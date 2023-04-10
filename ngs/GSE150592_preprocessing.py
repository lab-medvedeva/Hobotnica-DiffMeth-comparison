import os
import sys
from subprocess import call
import argparse

if __name__ == "__main__":
    srrs = ['SRR11790875', 'SRR11790876', 'SRR11790877', 'SRR11790878', 'SRR11790879',
    'SRR11790880', 'SRR11790881', 'SRR11790882', 'SRR11790883', 'SRR11790884',
    'SRR11790885', 'SRR11790886', 'SRR11790887', 'SRR11790888', 'SRR11790889',
    'SRR11790890', 'SRR11790891', 'SRR11790892', 'SRR11790893', 'SRR11790894',
    'SRR11790895', 'SRR11790896', 'SRR11790897', 'SRR11790898', 'SRR11790899',
    'SRR11790900', 'SRR11790901', 'SRR11790902', 'SRR11790903', 'SRR11790904']
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--cores_num', help='Number of cores', default=4)
    parser.add_argument('--reference', help='Path to the reference genome folder')
    parser.add_argument('--output_dir', help='Output directory')
    args = parser.parse_args()
    cores_num = int(args.cores_num)

    os.chdir(args.output_dir)

    # Sample processing
    for i in range(len(srrs)):
        srr = srrs[i]

        if not os.path.exists(f'{srr}'):
            os.makedirs(f'{srr}')

        os.chdir(f'{srr}')

        if call(f'fastq-dump --gzip --split-files {srr}', shell=True)!=0:
            sys.exit(f"fastq-dump {srr} failed")

        if call(f'fastqc {srr}_1.fastq.gz', shell=True)!=0:
            sys.exit(f"FastQC 1 {srr} failed")

        if call(f'fastqc {srr}_2.fastq.gz', shell=True)!=0:
            sys.exit(f"FastQC 2 {srr} failed")

        if call(f'trim_galore --gzip --cores {cores_num} --illumina --paired {srr}_1.fastq.gz {srr}_2.fastq.gz', shell=True)!=0:
            sys.exit(f"trim_galore {srr} failed")

        os.remove(f'{srr}_1.fastq.gz')
        os.remove(f'{srr}_2.fastq.gz')

        if call(f'fastqc {srr}_1_val_1.fq.gz', shell=True)!=0:
            sys.exit(f"FastQC 1 trimmed {srr} failed")

        if call(f'fastqc {srr}_2_val_2.fq.gz', shell=True)!=0:
            sys.exit(f"FastQC 2 trimmed {srr} failed")

        if call(f'bismark -p {cores_num} --gzip {args.reference} -1 {srr}_1_val_1.fq.gz -2 {srr}_2_val_2.fq.gz', shell=True)!=0:
            sys.exit(f"bismark {srr} failed")

        os.remove(f'{srr}_1_val_1.fq.gz')
        os.remove(f'{srr}_2_val_2.fq.gz')

        if call(f'bismark2report', shell=True)!=0:
            sys.exit(f"bismark2report failed {srr}")

        if call(f'bismark2summary', shell=True)!=0:
            sys.exit(f"bismark2summary failed {srr}")

        if call(f'bismark_methylation_extractor --bedGraph --parallel {cores_num}  --gzip {srr}_1_val_1_bismark_bt2_pe.bam', shell=True)!=0:
            sys.exit(f"bismark_methylation_extractor failed {srr}")

        if call(f'rm *.bedGraph.gz', shell=True)!=0:
            sys.exit(f"rm *.bedGraph.gz fail")

        if call(f'rm CHG_*', shell=True)!=0:
            sys.exit(f"rm CHG_* fail")
        if call(f'rm CHH_*', shell=True)!=0:
            sys.exit(f"rm CHH_* fail")
        if call(f'rm CpG_*', shell=True)!=0:
            sys.exit(f"rm CpG_* fail")

        os.remove(f'{srr}_1_val_1_bismark_bt2_pe.bam')

        if call(f'mv -f {srr}_1_val_1_bismark_bt2_pe.bismark.cov.gz ../{srr}.cov.gz', shell=True)!=0:
            sys.exit(f"move {srr}.cov failed")

        os.chdir('..')

        if call(f'gunzip -f {srr}.cov.gz', shell=True)!=0:
            sys.exit(f"gunzip {srr}.cov failed")