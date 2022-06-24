import os
import sys
from subprocess import call
import argparse

if __name__ == "__main__":
    srrs = ['SRR7587054', 'SRR7587055', 'SRR7587056', 'SRR7587057', 'SRR7587058', 'SRR7587059', 'SRR7587060',
        'SRR7587061', 'SRR7587062', 'SRR7587063', 'SRR7587064', 'SRR7587065', 'SRR7587066', 'SRR7587067', 'SRR7587068',
        'SRR7587069', 'SRR7587070', 'SRR7587071', 'SRR7587072', 'SRR7587073', 'SRR7587074', 'SRR7587075', 'SRR7587076',
        'SRR7587077', 'SRR7587078', 'SRR7587079', 'SRR7587080', 'SRR7587081', 'SRR7587082', 'SRR7587083', 'SRR7587084',
        'SRR7587085', 'SRR7587086', 'SRR7587087', 'SRR7587088', 'SRR7587089', 'SRR7587090', 'SRR7587091', 'SRR7587092',
        'SRR7587093', 'SRR7587094', 'SRR7587095', 'SRR7587096', 'SRR7587097', 'SRR7587098', 'SRR7587099', 'SRR7587100',
        'SRR7587101', 'SRR7587102', 'SRR7587103']
        
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

        if call(f'fastq-dump --gzip {srr}', shell=True)!=0:
            sys.exit(f"fastq-dump {srr} failed")

        if call(f'fastqc {srr}.fastq.gz', shell=True)!=0:
            sys.exit(f"FastQC {srr} failed")

        if call(f'trim_galore --gzip --cores {cores_num} --illumina {srr}.fastq.gz', shell=True)!=0:
            sys.exit(f"trim_galore {srr} failed")

        os.remove(f'{srr}.fastq.gz')

        if call(f'fastqc {srr}_trimmed.fq.gz', shell=True)!=0:
            sys.exit(f"FastQC trimmed {srr} failed")

        if call(f'bismark -p {cores_num} --gzip {args.reference} {srr}_trimmed.fq.gz', shell=True)!=0:
            sys.exit(f"bismark {srr} failed")

        os.remove(f'{srr}_trimmed.fq.gz')

        if call(f'deduplicate_bismark {srr}_trimmed_bismark_bt2.bam', shell=True)!=0:
            sys.exit(f"deduplicate_bismark {srr} failed")

        if call(f'bismark2report', shell=True)!=0:
            sys.exit(f"bismark2report failed {srr}")

        if call(f'bismark2summary', shell=True)!=0:
            sys.exit(f"bismark2summary failed {srr}")

        os.remove(f'{srr}_trimmed_bismark_bt2.bam')

        if call(f'bismark_methylation_extractor --bedGraph --parallel {cores_num}  --gzip {srr}_trimmed_bismark_bt2.deduplicated.bam', shell=True)!=0:
            sys.exit(f"bismark_methylation_extractor failed {srr}")

        if call(f'rm *.bedGraph.gz', shell=True)!=0:
            sys.exit(f"rm *.bedGraph.gz fail")

        if call(f'rm CHG_*', shell=True)!=0:
            sys.exit(f"rm CHG_* fail")
        if call(f'rm CHH_*', shell=True)!=0:
            sys.exit(f"rm CHH_* fail")
        if call(f'rm CpG_*', shell=True)!=0:
            sys.exit(f"rm CpG_* fail")

        os.remove(f'{srr}_trimmed_bismark_bt2.deduplicated.bam')

        if call(f'mv -f {srr}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz ../{srr}.cov.gz', shell=True)!=0:
            sys.exit(f"move {srr}.cov failed")

        os.chdir('..')

        if call(f'gunzip -f {srr}.cov.gz', shell=True)!=0:
            sys.exit(f"gunzip {srr}.cov failed")