#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------
# Minimal mitogenome CDS collector
# -----------------------------------

INPUT_CSV="Octopus_mitogenomes.csv"
OUTDIR="cds_sequences"
mkdir -p "$OUTDIR"

# Optional: wipe previous outputs (uncomment if you want a clean rebuild every run)
# rm -f "$OUTDIR"/*.fa

# Standardize gene name (GenBank headers vary: cob vs cytb, nd1 vs nad1, etc.)
standardize_gene() {
  local g="${1,,}"          # lowercase
  g="${g//[^a-z0-9]/}"      # strip odd characters
  case "$g" in
    cob)  echo "cytb" ;;
    cytb) echo "cytb" ;;
    nd1)  echo "nad1" ;;
    nd2)  echo "nad2" ;;
    nd3)  echo "nad3" ;;
    nd4)  echo "nad4" ;;
    nd4l) echo "nad4l" ;;
    nd5)  echo "nad5" ;;
    nd6)  echo "nad6" ;;
    cox1|coi) echo "cox1" ;;     # sometimes COI appears
    cox2) echo "cox2" ;;
    cox3) echo "cox3" ;;
    atp6) echo "atp6" ;;
    atp8) echo "atp8" ;;
    *)    echo "$g" ;;
  esac
}

# Parse CSV robustly (handles quoted commas in citation column)
mapfile -t PAIRS < <(python3 - <<'PY' "$INPUT_CSV"
import csv, sys
with open(sys.argv[1], newline="") as f:
    r = csv.reader(f)
    next(r, None)  # header
    for row in r:
        if not row:
            continue
        species = row[0].strip().strip('"')
        acc = row[1].strip().strip('"')
        if species and acc:
            print(f"{species}\t{acc}")
PY
)

echo "Found ${#PAIRS[@]} species/accessions in $INPUT_CSV"
echo "Writing gene FASTAs to: $OUTDIR"
echo "--------------------------------------"

for pair in "${PAIRS[@]}"; do
  IFS=$'\t' read -r species accession <<< "$pair"

  # Clean header label
  species_label="${species// /_}"

  echo "Fetching CDS for: $species_label ($accession)"

  # Fetch CDS FASTA for this accession (fastest: no esearch)
  tmpfile="$(mktemp)"
  if ! efetch -db nuccore -id "$accession" -format fasta_cds_na > "$tmpfile"; then
    echo "  ERROR: efetch failed for $accession"
    rm -f "$tmpfile"
    continue
  fi

  if [[ ! -s "$tmpfile" ]]; then
    echo "  WARNING: empty CDS result for $accession"
    rm -f "$tmpfile"
    continue
  fi

  # Split into gene-specific files:
  # For each CDS record, identify gene from header and append:
  #   >Species_label
  #   ATG...
  #
  # We do this in bash+awk for speed and simplicity.
  awk -v sp="$species_label" -v outdir="$OUTDIR" '
    function write_record() {
      if (gene != "" && seq != "") {
        print ">" sp >> (outdir "/" gene ".fa")
        print seq     >> (outdir "/" gene ".fa")
      }
    }
    function norm_gene(g) {
      # mimic the bash standardize_gene logic for common cases (lowercase + map)
      g = tolower(g)
      gsub(/[^a-z0-9]/, "", g)

      if (g == "cob")  return "cytb"
      if (g == "cytb") return "cytb"
      if (g == "nd1")  return "nad1"
      if (g == "nd2")  return "nad2"
      if (g == "nd3")  return "nad3"
      if (g == "nd4")  return "nad4"
      if (g == "nd4l") return "nad4l"
      if (g == "nd5")  return "nad5"
      if (g == "nd6")  return "nad6"
      if (g == "coi")  return "cox1"
      if (g == "cox1") return "cox1"
      if (g == "cox2") return "cox2"
      if (g == "cox3") return "cox3"
      if (g == "atp6") return "atp6"
      if (g == "atp8") return "atp8"
      return g
    }

    BEGIN { gene=""; seq="" }

    /^>/ {
      # flush previous record
      write_record()

      gene=""
      seq=""

      # Try to grab [gene=XXX]
      if (match($0, /\[gene=([^]]+)\]/, m)) gene = norm_gene(m[1])
      # Fallback: gene=XXX (space-delimited)
      else if (match($0, /gene=([^ ]+)/, m2)) gene = norm_gene(m2[1])

      next
    }

    {
      # sequence line
      if (gene != "") seq = seq $0
    }

    END { write_record() }
  ' "$tmpfile"

  rm -f "$tmpfile"
done

echo "--------------------------------------"
echo "Done."
echo "Per-gene FASTA files are in: $OUTDIR"