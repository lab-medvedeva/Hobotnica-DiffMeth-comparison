import os
import sys
from subprocess import call
import argparse
import pathlib

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hg38_reference', help='Path to the hg38 reference genome folder')
    parser.add_argument('--mm39_reference', help='Path to the mm39 reference genome folder')
    parser.add_argument('--hmm_dm_folder', help='Full path to the HMM-DM code folder')
    parser.add_argument('--preprocessing_cores_num', help='Number of cores for preprocessing', default=4)
    parser.add_argument('--methods_cores_num', help='Number of cores for methods test', default=16)
    parser.add_argument('--output_dir', help='Output directory')
    args = parser.parse_args()
    preprocessing_cores_num = int(args.preprocessing_cores_num)
    methods_cores_num = int(args.methods_cores_num)
    source_dir = pathlib.Path(__file__).parent.resolve()
    os.chdir(f'{args.output_dir}')

    # Run genome preparation
    if call(f'bismark_genome_preparation {args.hg38_reference}', shell=True)!=0:
        sys.exit(f"hg38 bismark_genome_preparation failed")

    if call(f'bismark_genome_preparation {args.mm39_reference}', shell=True)!=0:
        sys.exit(f"mm39 bismark_genome_preparation failed")

    datasets = ['GSE149608', 'GSE138598', 'GSE119980', 'GSE117593', 'GSE148060', 'GSE103886']
    for dataset in datasets:
        if not os.path.exists(f'{dataset}'):
            os.makedirs(f'{dataset}')

        # Run dataset preparation
        dataset_preprocessing_py = os.path.join(source_dir,f'{dataset}_preprocessing.py')
        if dataset == 'GSE103886':
            if call(f'python {dataset_preprocessing_py} --cores_num {preprocessing_cores_num} --reference {args.mm39_reference} --output_dir {dataset}', shell=True)!=0:
                sys.exit(f"{dataset_preprocessing_py} failed")
        else:
            if call(f'python {dataset_preprocessing_py} --cores_num {preprocessing_cores_num} --reference {args.hg38_reference} --output_dir {dataset}', shell=True)!=0:
                sys.exit(f"{dataset_preprocessing_py} failed")

        # Run methods for DMC identification
        WGBS_datasets = ['GSE149608', 'GSE138598', 'GSE119980', 'GSE117593']
        
        for dataset in datasets:

            # Initialization
            init_r = os.path.join(source_dir,'init.R')
            if call(f'Rscript --vanilla {init_r} --dataset {dataset} --working_dir {dataset}', shell=True)!=0:
                sys.exit(f"{dataset} initialization failed")

            # methylsig test
            methylsig_r = os.path.join(source_dir,'methylsig.R')
            if call(f'Rscript --vanilla {methylsig_r} --dataset {dataset} --working_dir {dataset} --cores {methods_cores_num} -s {source_dir}', shell=True)!=0:
                sys.exit(f"{dataset} methylsig failed")

            # DSS test
            dss_r = os.path.join(source_dir,'DSS.R')
            if call(f'Rscript --vanilla {dss_r} --dataset {dataset} --working_dir {dataset} --cores {methods_cores_num} -s {source_dir}', shell=True)!=0:
                sys.exit(f"{dataset} DSS failed")

            # methylkit test
            methylkit_r = os.path.join(source_dir,'methylkit.R')
            if call(f'Rscript --vanilla {methylkit_r} --dataset {dataset} --working_dir {dataset} --cores {methods_cores_num} -s {source_dir}', shell=True)!=0:
                sys.exit(f"{dataset} methylkit failed")

            # HMM-DM test
            hmm_dm_r = os.path.join(source_dir,'hmm_dm.R')
            if call(f'Rscript --vanilla {hmm_dm_r} --dataset {dataset} --working_dir {dataset} --cores {methods_cores_num} -s {source_dir} --hmm_dir {hmm_dm_folder}', shell=True)!=0:
                sys.exit(f"{dataset} hmm-dm failed")

            # RADMeth test
            radmeth_r = os.path.join(source_dir,'radmeth.R')
            if call(f'Rscript --vanilla {radmeth_r} --dataset {dataset} --working_dir {dataset} --cores {methods_cores_num} -s {source_dir}', shell=True)!=0:
                sys.exit(f"{dataset} radmeth failed")

            # BSmooth test
            if dataset in WGBS_datasets:
                bsmooth_r = os.path.join(source_dir,'bsmooth.R')
                if call(f'Rscript --vanilla {bsmooth_r} --dataset {dataset} --working_dir {dataset} --cores {methods_cores_num} -s {source_dir}', shell=True)!=0:
                    sys.exit(f"{dataset} bsmooth failed")