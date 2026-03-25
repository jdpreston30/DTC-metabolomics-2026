#!/usr/bin/env bash
# 03_copy.sh  (Step 3)
# Copies the .7z archives and zip_contents.csv produced by 02_compress.sh
# into the raw_zips_dir specified in 00_config.yaml.
# Original files on the external drive are not modified.
#
# Usage: bash R/MWB/Scripts/03_copy.sh
# Config: R/MWB/00_config.yaml

set -uo pipefail

CONFIG="$(dirname "$0")/../00_config.yaml"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: config not found at $CONFIG" >&2
  exit 1
fi

_get() { grep "^$1:" "$CONFIG" | sed "s/^$1:[[:space:]]*//"; }

mzxml_dir=$(_get mzxml_dir)
raw_zips_dir=$(_get raw_zips_dir)
raw_zips_dir="${raw_zips_dir/#\~/$HOME}"

# Source: sibling of mzxml_dir (where 02_compress.sh writes archives)
parent_dir="$(dirname "$mzxml_dir")"

# --- Validate ---
if [[ ! -d "$parent_dir" ]]; then
  echo "ERROR: parent of mzxml_dir not found: $parent_dir" >&2
  exit 1
fi

if [[ ! -d "$(dirname "$raw_zips_dir")" ]]; then
  echo "ERROR: parent of raw_zips_dir not found: $(dirname "$raw_zips_dir")" >&2
  echo "  Check the raw_zips_dir key in 00_config.yaml" >&2
  exit 1
fi

mkdir -p "$raw_zips_dir"
mkdir -p "${raw_zips_dir}/sld"

echo ""
echo "==> Copying archives to $raw_zips_dir"

copied=0
errors=0

# Copy .7z archives
for f in "${parent_dir}"/*.7z; do
  [[ -f "$f" ]] || continue
  fname="$(basename "$f")"
  echo "    Copying: $fname"
  if cp "$f" "${raw_zips_dir}/${fname}"; then
    copied=$((copied + 1))
  else
    echo "    ERROR: failed to copy $fname" >&2
    errors=$((errors + 1))
  fi
done

# Copy zip_contents.csv
csv="${parent_dir}/zip_contents.csv"
if [[ -f "$csv" ]]; then
  echo "    Copying: zip_contents.csv"
  if cp "$csv" "${raw_zips_dir}/zip_contents.csv"; then
    copied=$((copied + 1))
  else
    echo "    ERROR: failed to copy zip_contents.csv" >&2
    errors=$((errors + 1))
  fi
else
  echo "    WARNING: zip_contents.csv not found at $csv — skipping" >&2
fi

# Copy .sld files from parent_dir and its subfolders
echo ""
echo "==> Copying .sld files to ${raw_zips_dir}/sld/"
sld_count=0
while IFS= read -r -d '' sld_file; do
  fname="$(basename "$sld_file")"
  echo "    Copying: $fname"
  if cp "$sld_file" "${raw_zips_dir}/sld/${fname}"; then
    sld_count=$((sld_count + 1))
    copied=$((copied + 1))
  else
    echo "    ERROR: failed to copy $fname" >&2
    errors=$((errors + 1))
  fi
done < <(find "$parent_dir" -iname "*.sld" -print0 2>/dev/null)

if [[ "$sld_count" -eq 0 ]]; then
  echo "    No .sld files found under $parent_dir"
fi

echo ""
echo "Complete: $copied files copied, $errors errors"
echo "Destination: $raw_zips_dir"
