!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --nodelist=zenith

source /home/opt/Miniconda3/miniconda3/bin/activate singlem_0.18.3

# Define the directory containing your reads
READS_DIR="/home/projects/Agribiome/IN_RICHES/raw_reads/"
OUTPUT_DIR="/home/projects/Agribiome/IN_RICHES/singleM_tax_profiles"

# List of sample prefixes
SAMPLES=(
"CDA_0524_HighDesertCD_23410"
"CDA_0524_HighDesertCD_23472"
"CDA_0524_HighDesertCD_23903"
"CDA_0524_LarimerCD_23145"
"CDA_0524_LarimerCD_23239"
"CDA_0524_LarimerCD_23407"
"CDA_0524_LarimerCD_23729"
"CDA_0524_LarimerCD_23865"
"CDA_0524_LongmontCD_23342"
"CDA_0524_LongmontCD_23841"
"CDA_0524_MancosCD_23551"
"CDA_0524_MancosCD_23618"
"CDA_0524_MancosCD_23716"
"CDA_0524_MancosCD_23752"
"CDA_0524_MancosCD_23786"
"CDA_0524_MesaCD_23135")

# Loop through each sample and run singleM
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Running singleM for $SAMPLE"

    singlem pipe -1 "$READS_DIR/${SAMPLE}_R1.fastq.gz" \
                 -2 "$READS_DIR/${SAMPLE}_R2.fastq.gz" \
                 -p "$OUTPUT_DIR/${SAMPLE}_taxonomic.profile.tsv"

    echo "singleM completed for $SAMPLE"
done


