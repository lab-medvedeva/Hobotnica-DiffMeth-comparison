import os
import sys
from subprocess import call
import argparse
import re
import shutil
import random
import os.path
from multiprocessing import Pool
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def call_command(param):
	dir, cmd = param
	print(dir)
	os.chdir(dir)
	if call(cmd, shell=True)!=0:
		sys.exit(f"Failed for {cmd}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--cores', help='Number of cores', default=16)
	parser.add_argument('--reference', help='Path to the reference genome folder')
	parser.add_argument('--output_dir', help='Output directory')
	args = parser.parse_args()
	cores = args.cores
	output_dir = args.output_dir
	number_of_simulations = 10
	number_of_samples = 32
	number_of_case_samples = 16
	mCG_level = 0.65
	mCG_level_plus10 = mCG_level + 0.1
	mCG_level_minus10 = mCG_level - 0.1
	mCG_level_plus15 = mCG_level + 0.15
	mCG_level_minus15 = mCG_level - 0.15
	mCG_level_plus20 = mCG_level + 0.2
	mCG_level_minus20 = mCG_level - 0.2
	mCG_level_plus30 = mCG_level + 0.3
	mCG_level_minus30 = mCG_level - 0.3

	os.chdir(output_dir)
	if call("""wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz \
		| gunzip -c \
		| awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4 }' \
		| sort-bed - \
		> cpgIslandExt.hg38.bed""", shell=True)!=0:
		sys.exit("wget hg38 cpgIslandExt failed")

	################ Generate hypo and hyper methylated regions
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.mkdir(dir)
		os.chdir(dir)
		if call("pirs diploid --output-file chr22_diploid.fasta ../chr22.fasta", shell=True)!=0:
			sys.exit("Failed for pirs")
		if call("shuf -n 300 ../cpgIslandExt.hg38.chr22.bed > cpgIslandExt.hg38.chr22.sample_300.bed", shell=True)!=0:
			sys.exit("Failed for shuf -n 300")
		if call("shuf -n 150  cpgIslandExt.hg38.chr22.sample_300.bed > cpgIslandExt.hg38.chr22.hyper.bed", shell=True)!=0:
			sys.exit("Failed for shuf -n 150")
		if call("bedtools subtract -a cpgIslandExt.hg38.chr22.sample_300.bed -b cpgIslandExt.hg38.chr22.hyper.bed > cpgIslandExt.hg38.chr22.hypo.bed", shell=True)!=0:
			sys.exit("Failed for bedtools subtract")
		if call("bedtools getfasta -fi chr22_diploid.fasta -bed cpgIslandExt.hg38.chr22.hypo.bed > cpgIslandExt.hg38.chr22.hypo.fasta", shell=True)!=0:
			sys.exit(f"Failed for bedtools")
		if call("bedtools getfasta -fi chr22_diploid.fasta -bed cpgIslandExt.hg38.chr22.hyper.bed > cpgIslandExt.hg38.chr22.hyper.fasta", shell=True)!=0:
			sys.exit(f"Failed for bedtools")

	################ Generate base FASTQ
	commands = []
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for srr in range(number_of_samples):
			commands.append((dir, f"python RRBSsim --seed {i*100+srr} -f chr22_diploid.fasta -p --cut_site C-CGG -o {srr}_base_chr22"))
	with Pool(cores) as p:
		p.map(call_command, commands)

	################ Generate alt FASTQ
	commands = []
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for srr in range(16, number_of_samples):
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*200+srr} -f cpgIslandExt.hg38.chr22.hyper.fasta --mCG_level {mCG_level_plus10} -o {srr}_plus_10_chr22"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*300+srr} -f cpgIslandExt.hg38.chr22.hyper.fasta --mCG_level {mCG_level_plus30} -o {srr}_plus_30_chr22"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*400+srr} -f cpgIslandExt.hg38.chr22.hyper.fasta --mCG_level {mCG_level_plus20} -o {srr}_plus_20_chr22"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*500+srr} -f cpgIslandExt.hg38.chr22.hyper.fasta --mCG_level {mCG_level_plus15} -o {srr}_plus_15_chr22"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*600+srr} -f cpgIslandExt.hg38.chr22.hypo.fasta --mCG_level {mCG_level_minus10} -o {srr}_minus_10_chr22"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*700+srr} -f cpgIslandExt.hg38.chr22.hypo.fasta --mCG_level {mCG_level_minus30} -o {srr}_minus_30_chr22"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*800+srr} -f cpgIslandExt.hg38.chr22.hypo.fasta --mCG_level {mCG_level_minus20} -o {srr}_minus_20_chr22"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*900+srr} -f cpgIslandExt.hg38.chr22.hypo.fasta --mCG_level {mCG_level_minus15} -o {srr}_minus_15_chr22"))
	with Pool(cores) as p:
		p.map(call_command, commands)

	############ Preprocessing base
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}/rrbssim'
		os.chdir(dir)
		for srr in range(number_of_samples):
			os.mkdir( f'{dir}/{srr}')
			shutil.move(f"{dir}/{srr}_base_chr22.1.fq", f"{dir}/{srr}/{srr}_base_chr22.1.fq")
			shutil.move(f"{dir}/{srr}_base_chr22.2.fq", f"{dir}/{srr}/{srr}_base_chr22.2.fq")
			os.chdir(f'{dir}/{srr}')
			if call(f"trim_galore --rrbs --gzip --cores 8 --illumina --paired {srr}_base_chr22.1.fq {srr}_base_chr22.2.fq" , shell=True)!=0:
				sys.exit(f"Failed for trim_galore")
			if call(f"bismark --parallel 4 -p 4 --gzip {args.reference} -1 {srr}_base_chr22.1_val_1.fq.gz -2 {srr}_base_chr22.2_val_2.fq.gz" , shell=True)!=0:
				sys.exit(f"Failed for bismark")
			if call(f"bismark_methylation_extractor --bedGraph --parallel {cores} --gzip {srr}_base_chr22.1_val_1_bismark_bt2_pe.bam" , shell=True)!=0:
				sys.exit(f"Failed for bismark_methylation_extractor")

	############ Preprocessing alt
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}/rrbssim'
		os.chdir(dir)
		for srr in range(16, number_of_samples):
			for alt_dir in [f'{srr}_plus_10', f'{srr}_minus_10', f'{srr}_plus_15', f'{srr}_minus_15', f'{srr}_plus_20', f'{srr}_minus_20', f'{srr}_plus_30', f'{srr}_minus_30']:
				os.mkdir( f'{dir}/{srr}/{alt_dir}')
				shutil.move(f"{dir}/{alt_dir}_chr22.1.fq", f"{dir}/{srr}/{alt_dir}/{alt_dir}_chr22.1.fq")
				shutil.move(f"{dir}/{alt_dir}_chr22.2.fq", f"{dir}/{srr}/{alt_dir}/{alt_dir}_chr22.2.fq")
				os.chdir( f'{dir}/{srr}/{alt_dir}')
				if call(f"trim_galore --rrbs --gzip --cores 8 --illumina --paired {alt_dir}_chr22.1.fq {alt_dir}_chr22.2.fq" , shell=True)!=0:
					sys.exit(f"Failed for trim_galore")
				if call(f"bismark -p 4 --gzip {args.reference} -1 {alt_dir}_chr22.1_val_1.fq.gz -2 {alt_dir}_chr22.2_val_2.fq.gz" , shell=True)!=0:
					sys.exit(f"Failed for bismark")
				if call(f"bismark_methylation_extractor --bedGraph --parallel {cores} --gzip {alt_dir}_chr22.1_val_1_bismark_bt2_pe.bam" , shell=True)!=0:
					sys.exit(f"Failed for bismark_methylation_extractor")

	############ unpack base cov.gz
	for i in range(number_of_simulations): 
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for srr in range(number_of_samples):
			os.mkdir( f'result_alt')
			os.mkdir( f'result_perm')
			if os.path.isfile(f'rrbssim/{srr}/{srr}_base_chr22.1_val_1_bismark_bt2_pe.bismark.cov.gz'):
				print('unpack cov.gz')
				if call(f"gunzip rrbssim/{srr}/{srr}_base_chr22.1_val_1_bismark_bt2_pe.bismark.cov.gz", shell=True)!=0:
					sys.exit(f"Failed for gunzip")
			if os.path.isfile(f'rrbssim/{srr}/{srr}_base_chr22.1_val_1_bismark_bt2_pe.bismark.cov'):
				shutil.copy(f"rrbssim/{srr}/{srr}_base_chr22.1_val_1_bismark_bt2_pe.bismark.cov", f"result_alt/{srr}_base.cov")
				shutil.move(f"rrbssim/{srr}/{srr}_base_chr22.1_val_1_bismark_bt2_pe.bismark.cov", f"result_perm/{srr}_base.cov")

	############ unpack alt cov.gz
	for i in range(number_of_simulations): 
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for srr in range(16, number_of_samples):
			for alt_dir in [f'{srr}_plus_10', f'{srr}_minus_10', f'{srr}_plus_15', f'{srr}_minus_15', f'{srr}_plus_20', f'{srr}_minus_20', f'{srr}_plus_30', f'{srr}_minus_30']:
				if os.path.isfile(f'rrbssim/{srr}/{alt_dir}/{alt_dir}_chr22.1_val_1_bismark_bt2_pe.bismark.cov.gz'):
					print('unpack cov.gz')
					if call(f"gunzip rrbssim/{srr}/{alt_dir}/{alt_dir}_chr22.1_val_1_bismark_bt2_pe.bismark.cov.gz" , shell=True)!=0:
						sys.exit(f"Failed for gunzip")
				if os.path.isfile(f'rrbssim/{srr}/{alt_dir}/{alt_dir}_chr22.1_val_1_bismark_bt2_pe.bismark.cov'):
					shutil.move(f"rrbssim/{srr}/{alt_dir}/{alt_dir}_chr22.1_val_1_bismark_bt2_pe.bismark.cov", f"result_alt/{alt_dir}.cov")

	############ generate alt cov
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for srr in range(16, number_of_samples):
			if call(f"bedtools intersect -a result_alt/{srr}_base.cov -b cpgIslandExt.hg38.chr22.sample_300.bed -v > result_alt/{srr}_base.filtered.cov", shell=True)!=0:
				sys.exit(f"Failed for bedtools intersect")
			for differ in [f'10', f'15', f'20', f'30']:
				if call(f"cat result_alt/{srr}_plus_{differ}.cov result_alt/{srr}_minus_{differ}.cov result_alt/{srr}_base.filtered.cov > result_alt/{srr}_alt_{differ}.all.cov", shell=True)!=0:
					sys.exit(f"Failed for cat")
				if call(f"bedtools sort -i result_alt/{srr}_alt_{differ}.all.cov > result_alt/{srr}_alt_{differ}.sorted.cov", shell=True)!=0:
					sys.exit(f"Failed for bedtools sort")
				if os.path.isfile(f'result_alt/{srr}_alt_{differ}.sorted.cov'):
					shutil.copy(f"result_alt/{srr}_alt_{differ}.sorted.cov", f"result_perm/{srr}_alt_{differ}.sorted.cov")