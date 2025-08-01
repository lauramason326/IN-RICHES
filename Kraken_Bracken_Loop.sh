#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

# --- SETUP ---

# Directories
RAW_INPUT_DIR="/home/projects/Agribiome/IN_RICHES/bbduk_trimmed_reads"
KRAKEN_OUTPUT_DIR="/home/projects/Agribiome/IN_RICHES/kraken_out"
BRACKEN_OUTPUT_DIR="${KRAKEN_OUTPUT_DIR}/bracken_out"
KRONA_OUTPUT_DIR="${BRACKEN_OUTPUT_DIR}/krona_tables"
DB_DIR="/home/projects/Agribiome/Kraken2_fungi_db/kraken_fungi_db"
KRONA_SCRIPT="/home/projects/Agribiome/Kraken2_fungi_db/KrakenTools/kreport2krona.py"

# Make output directories
mkdir -p "$KRAKEN_OUTPUT_DIR" "$BRACKEN_OUTPUT_DIR" "$KRONA_OUTPUT_DIR"

# Sample list
SAMPLES=(
"IN_RICHES_106_ID_Aberdeen_SPM_p1_r1_0"
"IN_RICHES_110_ID_Aberdeen_WPM_p1_r1_0"
"IN_RICHES_114_ID_Aberdeen_BM_p1_r1_0"
"IN_RICHES_118_ID_Aberdeen_SPB_p1_r1_0"
"IN_RICHES_122_ID_Aberdeen_WPB_p1_r1_0"
"IN_RICHES_126_ID_Aberdeen_SPWPB_p1_r1_0"
#"IN_RICHES_4-Jun-2024_CO_AVRC_RM_2_2"#
"IN_RICHES_CO_23-May-2024_Josh_NT10_1_10"
"IN_RICHES_CO_23-May-2024_Josh_NT30_1_11"
"IN_RICHES_CO_23-May-2024_Josh_NT30_1_2"
"IN_RICHES_CO_23-May-2024_Josh_NT30_1_5"
#"IN_RICHES_CO_23-May-2024_Josh_NT30_1_8"#
"IN_RICHES_4-Jun-2024_CO_AVRC_RM_1_1"
"IN_RICHES_4-Jun-2024_CO_AVRC_RT_1_1"
#"CDA_0524_COAllianceCenterEE_23715"#
"CDA_0524_SangreDeCristoAcequiaAssocEE_23519"
"CDA_0524_COFirstCD_23061"
"CDA_0524_COFirstCD_23969"
"CDA_0524_HaxtunCD_23770"
"IN_RICHES_216_WY_Lerwick-2_R2_p1_r1_0")


# --- STEP 1: Kraken2 ---
echo "Running Kraken2 classification..."
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Running Kraken2 for $SAMPLE..."

  R1="${RAW_INPUT_DIR}/${SAMPLE}_R1_trimmed_bbduk.fastq"
  R2="${RAW_INPUT_DIR}/${SAMPLE}_R2_trimmed_bbduk.fastq"
  REPORT="${KRAKEN_OUTPUT_DIR}/${SAMPLE}_fungal_report.txt"
  OUTPUT="${KRAKEN_OUTPUT_DIR}/${SAMPLE}_fungal_output.txt"
  CLASSIFIED="${KRAKEN_OUTPUT_DIR}/${SAMPLE}_classified#.fq"
  UNCLASSIFIED="${KRAKEN_OUTPUT_DIR}/${SAMPLE}_unclassified#.fq"

  # Skip if input files are missing
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "Missing input files for $SAMPLE. Skipping."
    continue
  fi

kraken2 \
    --db "$DB_DIR" \
    --threads 20 \
    --paired "$R1" "$R2" \
    --report "$REPORT" \
    --output "$OUTPUT" \
    --classified-out "$CLASSIFIED" \
    --unclassified-out "$UNCLASSIFIED"
done

# --- STEP 2: Bracken ---
echo "Running Bracken refinement..."
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Running Bracken for $SAMPLE..."

  REPORT="${KRAKEN_OUTPUT_DIR}/${SAMPLE}_fungal_report.txt"
  B_REPORT="${BRACKEN_OUTPUT_DIR}/${SAMPLE}_bracken_report.txt"
  B_OUTPUT="${BRACKEN_OUTPUT_DIR}/${SAMPLE}_bracken_output.txt"

  if [[ ! -f "$REPORT" ]]; then
    echo "Missing Kraken report for $SAMPLE. Skipping Bracken."
    continue
  fi

  bracken \
    -d "$DB_DIR" \
    -i "$REPORT" \
    -r 151 \
    -t 10 \
    -o "$B_OUTPUT" \
    -w "$B_REPORT"
done

# --- STEP 3: kreport2krona + Labeling ---
echo "Generating Krona tables with sample labels..."
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Processing $SAMPLE for Krona..."

  BRACKEN_REPORT="${BRACKEN_OUTPUT_DIR}/${SAMPLE}_bracken_report.txt"
  KRONA_OUT="${KRONA_OUTPUT_DIR}/${SAMPLE}_krona_table.txt"
  KRONA_LABELED="${KRONA_OUTPUT_DIR}/${SAMPLE}_krona_table_with_sample.txt"

  if [[ ! -f "$BRACKEN_REPORT" ]]; then
    echo "Missing Bracken report for $SAMPLE. Skipping Krona."
    continue
  fi

  python "$KRONA_SCRIPT" -r "$BRACKEN_REPORT" -o "$KRONA_OUT"
  awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"} {print sample, $0}' "$KRONA_OUT" > "$KRONA_LABELED"
done

