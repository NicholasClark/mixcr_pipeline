#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <target_directory>"
    exit 1
fi

# Check that GNU parallel is available
if ! command -v parallel &>/dev/null; then
    echo "Error: GNU parallel is not available"
    echo "Please load the module before running this script:"
    echo ""
    echo "  ml parallel"
    echo ""
    exit 1
fi

TARGET_DIR="$1"

if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: '$TARGET_DIR' is not a valid directory"
    exit 1
fi

CLONES_BASE="$(dirname "$TARGET_DIR")/clones_mouse"

echo "Target directory: $TARGET_DIR"
num_files=$(ls "$TARGET_DIR" | wc -l)
echo "Target directory contains $num_files total files"
num_tsv=$(ls "$TARGET_DIR"/*.tsv.gz 2>/dev/null | wc -l)
echo "Target directory contains $num_tsv .tsv.gz files"
echo "First 10 files in target directory:"
ls "$TARGET_DIR"/*.tsv.gz 2>/dev/null | head | xargs -n1 basename

# Detect available processors and use 80%
total_cores=$(nproc)
jobs=$(echo "$total_cores * 0.8 / 1" | bc)
jobs=$(( jobs < 1 ? 1 : jobs ))
echo ""
echo "Using $jobs of $total_cores available cores"
echo "Moving .tsv.gz files to $CLONES_BASE/mUD##/:"

# Step 1: Parse all filenames with awk into a temp file
tmpfile=$(mktemp)
ls "$TARGET_DIR"/*.tsv.gz 2>/dev/null | awk '{
    filepath = $0
    filename = filepath
    sub(/.*\//, "", filename)

    ud = filename
    if (match(ud, /mUD[0-9]+/))
        ud = substr(ud, RSTART, RLENGTH)
    else
        ud = ""

    tcr_chain = filename
    if (match(tcr_chain, /TCR[ab]/))
        tcr_chain = substr(tcr_chain, RSTART, RLENGTH)
    else
        tcr_chain = ""

    tr_suffix = filename
    if (match(tr_suffix, /TR[AB]/))
        tr_suffix = substr(tr_suffix, RSTART, RLENGTH)
    else
        tr_suffix = ""

    if ((tcr_chain == "TCRa" && tr_suffix == "TRA") || \
        (tcr_chain == "TCRb" && tr_suffix == "TRB"))
        print filepath "\t" ud
}' > "$tmpfile"

echo "Parsed $(wc -l < "$tmpfile") valid files"

# Step 2: Pre-create all unique UD## directories and verify they exist
echo "Creating destination directories..."
awk -F'\t' '{print $2}' "$tmpfile" | sort -u | while read -r ud; do
    dest="$CLONES_BASE/$ud"
    mkdir -p "$dest"
    if [ ! -d "$dest" ]; then
        echo "ERROR: Failed to create directory: $dest"
        rm "$tmpfile"
        exit 1
    fi
    #echo "  Verified: $dest"
done

# Double-check all directories exist before proceeding
echo "Verifying all directories exist..."
all_ok=true
while IFS=$'\t' read -r filepath ud; do
    dest="$CLONES_BASE/$ud"
    if [ ! -d "$dest" ]; then
        echo "ERROR: Directory missing: $dest"
        all_ok=false
    fi
done < "$tmpfile"

if [ "$all_ok" = false ]; then
    echo "Aborting: some destination directories are missing"
    rm "$tmpfile"
    exit 1
fi

echo "All directories verified, starting file moves..."

# Step 3: Move files in parallel using exported function
export CLONES_BASE
move_file() {
    mv "$1" "$CLONES_BASE/$2/"
}
export -f move_file
parallel --colsep '\t' --jobs "$jobs" move_file {1} {2} < "$tmpfile"

# Clean up temp file
rm "$tmpfile"

echo ""
echo "Done. Contents of $CLONES_BASE/:"
for ud_dir in "$CLONES_BASE"/mUD*/; do
    if [ -d "$ud_dir" ]; then
        files=("$ud_dir"/*)
        echo "  $(basename "$ud_dir"): ${#files[@]} file(s)"
    fi
done