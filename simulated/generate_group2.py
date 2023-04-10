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
	GSE103886_srrs = ['SRR6040066', 'SRR6040067', 'SRR6040068', 'SRR6040069', 'SRR6040070', 'SRR6040071', 'SRR6040072', 'SRR6040073', 'SRR6040074', 'SRR6040075', 'SRR6040076', 'SRR6040077', 'SRR6040078', 'SRR6040079', 'SRR6040080', 'SRR6040081', 'SRR6040082', 'SRR6040083', 'SRR6040084', 'SRR6040085', 'SRR6040086', 'SRR6040087', 'SRR6040088']
	GSE103886_srrs_alt = ['SRR6040077', 'SRR6040078', 'SRR6040079', 'SRR6040080', 'SRR6040081', 'SRR6040082', 'SRR6040083', 'SRR6040084', 'SRR6040085', 'SRR6040086', 'SRR6040087', 'SRR6040088']
	GSE103886_depths = [30.965027046892416,
	31.83177141898794, 
	58.44074494587031, 
	57.28513594726502, 
	43.01340247035768, 
	46.16182904931401, 
	39.33012929879564, 
	41.87589938147192, 
	34.30700224628341, 
	31.343917116582723, 
	34.38032426552555, 
	38.72252482108035, 
	36.06931597777508, 
	37.22091363683115, 
	36.92134474957499, 
	45.36860791623172, 
	38.44448995553599, 
	37.096035575098334, 
	33.629546656554666, 
	33.16289014687882, 
	16.642720323729083, 
	29.951886243978745, 
	29.8356264365599]

	GSE103886_means = np.array([0.3931666,
	 0.4188803,
	 0.4060564,
	 0.4055235,
	 0.3889425,
	 0.3889739,
	 0.3991975,
	 0.3924620,
	 0.3862409,
	 0.3867527,
	 0.4017754,
	 0.3922298,
	 0.3948329,
	 0.3979553,
	 0.3974604,
	 0.4005577,
	 0.3939143,
	 0.3994970,
	 0.3843206,
	 0.3935221,
	 0.3832796,
	 0.3984322,
	 0.3915062])

	GSE103886_means_minus_10 = GSE103886_means - 0.1
	GSE103886_means_minus_15 = GSE103886_means - 0.15
	GSE103886_means_minus_20 = GSE103886_means - 0.20
	GSE103886_means_minus_30 = GSE103886_means - 0.30
	GSE103886_means_plus_10 = GSE103886_means + 0.1
	GSE103886_means_plus_15 = GSE103886_means + 0.15
	GSE103886_means_plus_20 = GSE103886_means + 0.20
	GSE103886_means_plus_30 = GSE103886_means + 0.30
	GSE103886_read_length = 50

	os.chdir(output_dir)

	if call("""wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/mm39/database/cpgIslandExt.txt.gz \
	   | gunzip -c \
	   | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4 }' \
	   | sort-bed - \
	   > cpgIslandExt.mm39.bed""", shell=True)!=0:
	    sys.exit("wget mm39 cpgIslandExt failed")

	for j, srr in enumerate(GSE103886_srrs):
		os.chdir(f"{output_dir}/{srr}")
		if call(f"python base_quality_profile_R1 -i  {srr}.1_1.fastq", shell=True)!=0:
			sys.exit("Failed for base_quality_profile_R1")

		if call(f"python base_quality_profile_R2 -i  {srr}.1_2.fastq", shell=True)!=0:
			sys.exit("Failed for base_quality_profile_R2")

	################ Generate hypo and hyper methylated regions
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.mkdir(dir)
		os.chdir(dir)
		if call("pirs diploid --output-file chr19_diploid.fasta ../chr19.fasta", shell=True)!=0:
			sys.exit("Failed for pirs")
		if call("shuf -n 300 ../cpgIslandExt.mm39.chr19.bed > cpgIslandExt.mm39.chr19.sample_300.bed", shell=True)!=0:
			sys.exit("Failed for shuf -n 300")
		if call("shuf -n 150  cpgIslandExt.mm39.chr19.sample_300.bed > cpgIslandExt.mm39.chr19.hyper.bed", shell=True)!=0:
			sys.exit("Failed for shuf -n 150")
		if call("bedtools subtract -a cpgIslandExt.mm39.chr19.sample_300.bed -b cpgIslandExt.mm39.chr19.hyper.bed > cpgIslandExt.mm39.chr19.hypo.bed", shell=True)!=0:
			sys.exit("Failed for bedtools subtract")
		if call("bedtools getfasta -fi chr19_diploid.fasta -bed cpgIslandExt.mm39.chr19.hypo.bed > cpgIslandExt.mm39.chr19.hypo.fasta", shell=True)!=0:
			sys.exit(f"Failed for bedtools")
		if call("bedtools getfasta -fi chr19_diploid.fasta -bed cpgIslandExt.mm39.chr19.hyper.bed > cpgIslandExt.mm39.chr19.hyper.fasta", shell=True)!=0:
			sys.exit(f"Failed for bedtools")

	################ Generate base FASTQ
	commands = []
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for j, srr in enumerate(GSE103886_srrs):
			commands.append((dir, f"python RRBSsim --seed {i*100+j} -f chr19_diploid.fasta -p --mCG_level {GSE103886_means[j]}  -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_base_chr19"))
	with Pool(cores) as p:
		p.map(call_command, commands)

	################ Generate alt FASTQ
	commands = []
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for srr in range(16, number_of_samples):
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*200+j} -f cpgIslandExt.mm39.chr19.hyper.fasta --mCG_level {GSE103886_means_plus_10[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_plus_10_chr19"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*300+j} -f cpgIslandExt.mm39.chr19.hyper.fasta --mCG_level {GSE103886_means_plus_30[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_plus_30_chr19"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*400+j} -f cpgIslandExt.mm39.chr19.hyper.fasta --mCG_level {GSE103886_means_plus_20[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_plus_20_chr19"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*500+j} -f cpgIslandExt.mm39.chr19.hyper.fasta --mCG_level {GSE103886_means_plus_15[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_plus_15_chr19"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*600+j} -f cpgIslandExt.mm39.chr19.hypo.fasta --mCG_level {GSE103886_means_minus_10[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_minus_10_chr19"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*700+j} -f cpgIslandExt.mm39.chr19.hypo.fasta --mCG_level {GSE103886_means_minus_30[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_minus_30_chr19"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*800+j} -f cpgIslandExt.mm39.chr19.hypo.fasta --mCG_level {GSE103886_means_minus_20[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_minus_20_chr19"))
			commands.append((dir, f"python RRBSsim --CG_rate 0.99 --seed {i*900+j} -f cpgIslandExt.mm39.chr19.hypo.fasta --mCG_level {GSE103886_means_minus_15[j]} -l {GSE103886_read_length} -d {round(GSE103886_depths[j])} --cut_site C-CGG  --matrix_qual_R1 {srr}/Base-Calling_Profiles/base_quality_profile.R1 --matrix_qual_R2 {srr}/Base-Calling_Profiles/base_quality_profile.R2  -o {srr}_minus_15_chr19"))
	with Pool(cores) as p:
		p.map(call_command, commands)

	################ Preprocessing base
	for i in range(number_of_simulations): #
		dir = f'{output_dir}/simulation{i}/rrbssim'
		os.chdir(dir)
		for j, srr in enumerate(GSE103886_srrs):
			os.mkdir( f'{dir}/{srr}')
			shutil.move(f"{dir}/{srr}_base_chr19.1.fq", f"{dir}/{srr}/{srr}_base_chr19.1.fq")
			shutil.move(f"{dir}/{srr}_base_chr19.2.fq", f"{dir}/{srr}/{srr}_base_chr19.2.fq")
			os.chdir(f'{dir}/{srr}')
			if call(f"trim_galore --rrbs --gzip --cores 8 --illumina --paired {srr}_base_chr19.1.fq {srr}_base_chr19.2.fq" , shell=True)!=0:
				sys.exit(f"Failed for trim_galore")
			if call(f"bismark --parallel 4 -p 4 --gzip {args.reference} -1 {srr}_base_chr19.1_val_1.fq.gz -2 {srr}_base_chr19.2_val_2.fq.gz" , shell=True)!=0:
				sys.exit(f"Failed for bismark")
			if call(f"bismark_methylation_extractor --bedGraph --parallel {cores} --gzip {srr}_base_chr19.1_val_1_bismark_bt2_pe.bam" , shell=True)!=0:
				sys.exit(f"Failed for bismark_methylation_extractor")

	################ Preprocessing alt
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}/rrbssim'
		os.chdir(dir)
		for j, srr in enumerate(GSE103886_srrs_alt):
			for alt_dir in [f'{srr}_plus_10', f'{srr}_minus_10', f'{srr}_plus_15', f'{srr}_minus_15', f'{srr}_plus_20', f'{srr}_minus_20', f'{srr}_plus_30', f'{srr}_minus_30']:
				os.mkdir( f'{dir}/{srr}/{alt_dir}')
				shutil.move(f"{dir}/{alt_dir}_chr19.1.fq", f"{dir}/{srr}/{alt_dir}/{alt_dir}_chr19.1.fq")
				shutil.move(f"{dir}/{alt_dir}_chr19.2.fq", f"{dir}/{srr}/{alt_dir}/{alt_dir}_chr19.2.fq")
				os.chdir( f'{dir}/{srr}/{alt_dir}')
				if call(f"/home/abudkina/TrimGalore-0.6.6/trim_galore --rrbs --gzip --cores 4 --illumina --paired {alt_dir}_chr19.1.fq {alt_dir}_chr19.2.fq" , shell=True)!=0:
					sys.exit(f"Failed for trim_galore")
				if call(f"~/Bismark-0.22.3/bismark -p 4 --gzip {args.reference} -1 {alt_dir}_chr19.1_val_1.fq.gz -2 {alt_dir}_chr19.2_val_2.fq.gz" , shell=True)!=0:
					sys.exit(f"Failed for bismark")
				if call(f"~/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --parallel 8 --gzip {alt_dir}_chr19.1_val_1_bismark_bt2_pe.bam" , shell=True)!=0:
					sys.exit(f"Failed for bismark_methylation_extractor")

	################ unpack base cov.gz
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for j, srr in enumerate(GSE103886_srrs):
			if os.path.isfile(f'rrbssim/{srr}/{srr}_base_chr19.1_val_1_bismark_bt2_pe.bismark.cov.gz'):
				print('unpack cov.gz')
				if call(f"gunzip rrbssim/{srr}/{srr}_base_chr19.1_val_1_bismark_bt2_pe.bismark.cov.gz" , shell=True)!=0:
					sys.exit(f"Failed for gunzip")
			if os.path.isfile(f'rrbssim/{srr}/{srr}_base_chr19.1_val_1_bismark_bt2_pe.bismark.cov'):
				shutil.copy(f"rrbssim/{srr}/{srr}_base_chr19.1_val_1_bismark_bt2_pe.bismark.cov", f"result_alt/{srr}_base.cov")
				shutil.move(f"rrbssim/{srr}/{srr}_base_chr19.1_val_1_bismark_bt2_pe.bismark.cov", f"result_perm/{srr}_base.cov")

	################ unpack alt cov.gz
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		for j, srr in enumerate(GSE103886_srrs_alt):
			for alt_dir in [f'{srr}_plus_10', f'{srr}_minus_10', f'{srr}_plus_15', f'{srr}_minus_15', f'{srr}_plus_20', f'{srr}_minus_20', f'{srr}_plus_30', f'{srr}_minus_30']:
				if os.path.isfile(f'rrbssim/{srr}/{alt_dir}/{alt_dir}_chr19.1_val_1_bismark_bt2_pe.bismark.cov.gz'):
					print('unpack cov.gz')
					if call(f"gunzip rrbssim/{srr}/{alt_dir}/{alt_dir}_chr19.1_val_1_bismark_bt2_pe.bismark.cov.gz" , shell=True)!=0:
						sys.exit(f"Failed for gunzip")
				if os.path.isfile(f'rrbssim/{srr}/{alt_dir}/{alt_dir}_chr19.1_val_1_bismark_bt2_pe.bismark.cov'):
					shutil.move(f"rrbssim/{srr}/{alt_dir}/{alt_dir}_chr19.1_val_1_bismark_bt2_pe.bismark.cov", f"result_alt/{alt_dir}.cov")
					
	############ generate alt cov
	for i in range(number_of_simulations):
		dir = f'{output_dir}/simulation{i}'
		os.chdir(dir)
		print(dir)
		for j, srr in enumerate(GSE103886_srrs_alt):
			if call(f"bedtools intersect -a result_alt/{srr}_base.cov -b cpgIslandExt.mm39.chr19.sample_300.bed -v > result_alt/{srr}_base.filtered.cov", shell=True)!=0:
				sys.exit(f"Failed for bedtools intersect")
			for differ in [f'10', f'15', f'20', f'30']:
				if call(f"cat result_alt/{srr}_plus_{differ}.cov result_alt/{srr}_minus_{differ}.cov result_alt/{srr}_base.filtered.cov > result_alt/{srr}_alt_{differ}.all.cov", shell=True)!=0:
					sys.exit(f"Failed for cat")
				if call(f"bedtools sort -i result_alt/{srr}_alt_{differ}.all.cov > result_alt/{srr}_alt_{differ}.sorted.cov", shell=True)!=0:
					sys.exit(f"Failed for bedtools sort")
				if os.path.isfile(f'result_alt/{srr}_alt_{differ}.sorted.cov'):
					shutil.copy(f"result_alt/{srr}_alt_{differ}.sorted.cov", f"result_perm/{srr}_alt_{differ}.sorted.cov")
