import os
import sys
from subprocess import call
import argparse

if __name__ == "__main__":
    srrs = ['SRR11477199', 'SRR11477200', 'SRR11477201', 'SRR11477202', 'SRR11477203', 'SRR11477204', 'SRR11477205',
        'SRR11477206', 'SRR11477207', 'SRR11477208', 'SRR11477209', 'SRR11477210', 'SRR11477211', 'SRR11477212', 'SRR11477213',
        'SRR11477214', 'SRR11477215', 'SRR11477216', 'SRR11477217', 'SRR11477218', 'SRR11477219', 'SRR11477220', 'SRR11477221',
        'SRR11477222', 'SRR11477223', 'SRR11477224', 'SRR11477225', 'SRR11477226', 'SRR11477227', 'SRR11477228', 'SRR11477229',
        'SRR11477230', 'SRR11477231', 'SRR11477232', 'SRR11477233', 'SRR11477234', 'SRR11477235', 'SRR11477236', 'SRR11477237',
        'SRR11477238', 'SRR11477239', 'SRR11477240', 'SRR11477241', 'SRR11477242', 'SRR11477243', 'SRR11477244', 'SRR11477245',
        'SRR11477246', 'SRR11477247', 'SRR11477248', 'SRR11477249', 'SRR11477250', 'SRR11477251']
        
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

        if call(f'bismark2report', shell=True)!=0:
            sys.exit(f"bismark2report failed {srr}")

        if call(f'bismark2summary', shell=True)!=0:
            sys.exit(f"bismark2summary failed {srr}")

        if call(f'bismark_methylation_extractor --bedGraph --parallel {cores_num}  --gzip {srr}_trimmed_bismark_bt2.bam', shell=True)!=0:
            sys.exit(f"bismark_methylation_extractor failed {srr}")

        if call(f'rm *.bedGraph.gz', shell=True)!=0:
            sys.exit(f"rm *.bedGraph.gz fail")

        if call(f'rm CHG_*', shell=True)!=0:
            sys.exit(f"rm CHG_* fail")
        if call(f'rm CHH_*', shell=True)!=0:
            sys.exit(f"rm CHH_* fail")
        if call(f'rm CpG_*', shell=True)!=0:
            sys.exit(f"rm CpG_* fail")

        os.remove(f'{srr}_trimmed_bismark_bt2.bam')

        if call(f'mv -f {srr}_trimmed_bismark_bt2.bismark.cov.gz ../{srr}.cov.gz', shell=True)!=0:
            sys.exit(f"move {srr}.cov failed")

        os.chdir('..')

        if call(f'gunzip -f {srr}.cov.gz', shell=True)!=0:
            sys.exit(f"gunzip {srr}.cov failed")