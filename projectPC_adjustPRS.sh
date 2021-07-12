#!/bin/bash

## This script takes a sample bed/bim/fam file and projection to the PC (priciple components) space that
## has been established from MGB Biobank dataset. The resulting PC file is saved next to the input file
## Then, an R script takes the above projection PC file (containing 4 PCs), and predict a PRS score
## based on the PCs and the provided model and parameters. The resulting file is saved next to
## the raw_score file

## Note: This script has to be saved with all other required files in the same folder. The folder can
##       be moved to anywhere

## Example of two input files
## $1=/data/pcpgm/projects/prs_pipeline/genome_file/Genomefile_100820/Clin11-00888_final/Clin11-00888_final/Clin11-00888_final/Clin11-00888_final (.bed, .bim, .fam)
## $2=/data/pcpgm/projects/prs_pipeline/genome_file/Genomefile_100820/Clin11-00888_final/Clin11-00888_final/Clin11-00888_final/Clin11-00888_final_Clin11-00888_final/Clin11-00888_final_Afib/Clin11-00888_final/Clin11-00888_final_Afib.sscore
## The format of the $2 file is as below:
## #FID    IID     NMISS_ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_AVG      SCORE1_SUM
## 0       Clin11-00888      13459006        6388466 2.4127e-06      32.4725

### The relative directory of this current script
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Prompt for arguments
usage(){
  echo "You have to provide one sample bed file and 1-6 optional PRS raw_score file(s) to run this script!"
  echo "Usage: $0 sample_bed_file [--af] [--breastC] [--cad] [--colorectalC] [--prostateC] [--t2d]"
  exit 1
}

if [ $# -eq 0 ]; then
  usage
fi

#load two modules
module load plink/2.0a
module load R-mkl/3.0.2

# Handle a positional parameter and six optional parameters
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -a|--af)          af_raw_score="$2";  shift 2;;  # shift twice to pass a pair of argument and value
    -b|--breastC)     bc_raw_score="$2";  shift 2;;
    -c|--cad)         cad_raw_score="$2"; shift 2;;
    -o|--colorectalC) cc_raw_score="$2";  shift 2;;
    -p|--prostateC)   pc_raw_score="$2";  shift 2;;
    -t|--t2d)         t2d_raw_score="$2"; shift 2;;
    *)                POSITIONAL+=("$1"); shift;;
  esac
done

### Read in input file and create an output file name
sample_bed_file=${POSITIONAL}
projection_pc_file=${sample_bed_file}_new_projection


### PC projection of a sample bed file. The resulting file name is ${projection_pc_file}.sscore
plink2 --bfile ${sample_bed_file} \
       --read-freq ${script_dir}/ref_pcs.acount \
       --score ${script_dir}/ref_pcs.eigenvec.var 2 4 header-read \
               no-mean-imputation variance-standardize list-variants \
       --score-col-nums 5-8 \
       --extract ${script_dir}/race_inform_MGBbiobank_id.txt \
       --out ${projection_pc_file}

### Adjust PRS score for six diseases based on the trained model. The restulting
### file is ${sample_raw_score}.adjusted, next to the raw_score file
# Adjust for AF
if [ ! -z "$af_raw_score" ]; then
  Rscript --vanilla ${script_dir}/prs_adjustment_single_sample.R \
                      ${af_raw_score} \
                      ${projection_pc_file}.sscore \
                      ${script_dir}/ref_pcs.eigenval \
                      ${script_dir}/af_pc_model.RDS
  echo "AF score is adjusted!"
fi

# Adjust for cad
if [ ! -z "$cad_raw_score" ]; then
  Rscript --vanilla ${script_dir}/prs_adjustment_single_sample.R \
                      ${cad_raw_score} \
                      ${projection_pc_file}.sscore \
                      ${script_dir}/ref_pcs.eigenval \
                      ${script_dir}/cad_pc_model.RDS
  echo "CAD score is adjusted!"
fi

# Adjust for t2d
if [ ! -z "$t2d_raw_score" ]; then
  Rscript --vanilla ${script_dir}/prs_adjustment_single_sample.R \
                      ${t2d_raw_score} \
                      ${projection_pc_file}.sscore \
                      ${script_dir}/ref_pcs.eigenval \
                      ${script_dir}/t2d_pc_model.RDS
  echo "T2D score is adjusted!"
fi

# Adjust for BC
if [ ! -z "$bc_raw_score" ]; then
  Rscript --vanilla ${script_dir}/prs_adjustment_single_sample.R \
                      ${bc_raw_score} \
                      ${projection_pc_file}.sscore \
                      ${script_dir}/ref_pcs.eigenval \
                      ${script_dir}/breastC_pc_model.RDS
  echo "BreastCA score is adjusted!"
fi

# Adjust for CC
if [ ! -z "$cc_raw_score" ]; then
  Rscript --vanilla ${script_dir}/prs_adjustment_single_sample.R \
                      ${cc_raw_score} \
                      ${projection_pc_file}.sscore \
                      ${script_dir}/ref_pcs.eigenval \
                      ${script_dir}/colorectalC_pc_model.RDS
  echo "ColorectalCA score is adjusted!"
fi

# Adjust for prostateC
if [ ! -z "$pc_raw_score" ]; then
  Rscript --vanilla ${script_dir}/prs_adjustment_single_sample.R \
                      ${pc_raw_score} \
                      ${projection_pc_file}.sscore \
                      ${script_dir}/ref_pcs.eigenval \
                      ${script_dir}/prostateC_pc_model.RDS
  echo "ProstateCA score is adjusted!"
fi
