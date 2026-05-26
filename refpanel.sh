#!/bin/bash
set -e

# 1. Backup the old file
cp config/refpanel.tsv config/refpanel.tsv.bak

# 2. Create a new header
echo -e "chr\tvcf\tstart\tend\tquilt_map" > config/refpanel.tsv

# 3. Loop and populate with REAL start positions from the map files
for CHR in {1..22} X; do
    # Define file paths
    VCF="data/refpanel/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    MAP="data/maps/beagle-map-converted/chr${CHR}.map"

    # Get START from the first line of the map (Skipping header)
    START=$(awk 'NR==2 {print $1}' "$MAP")

    # Get END from the last line of the map
    END=$(tail -n 1 "$MAP" | awk '{print $1}')

    # Append to the new refpanel.tsv
    echo -e "chr${CHR}\t${VCF}\t${START}\t${END}\t${MAP}" >> config/refpanel.tsv

    echo "Updated chr${CHR}: Start=$START, End=$END"
done
