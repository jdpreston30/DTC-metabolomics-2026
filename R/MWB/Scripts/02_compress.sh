#!/usr/bin/env bash
# 02_compress.sh  (Step 2)
# Compresses mzXML files from hilicpos/ and c18neg/ subfolders into .7z archives
# and generates zip_contents.csv. Archives are written as siblings of mzxml_dir.
# REQUIRES: mzxml_dir must contain hilicpos/ and c18neg/ — run 01_organize.sh first.
#
# Usage: bash R/MWB/Scripts/02_compress.sh
# Config: R/MWB/00_config.yaml

set -uo pipefail

CONFIG="$(dirname "$0")/../00_config.yaml"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: config not found at $CONFIG" >&2; exit 1
fi

_get() { grep "^$1:" "$CONFIG" | sed "s/^$1:[[:space:]]*//"; }

mzxml_dir=$(_get mzxml_dir)
[[ -z "$mzxml_dir" ]] && { echo "ERROR: mzxml_dir not set in $CONFIG" >&2; exit 1; }
[[ -d "$mzxml_dir" ]] || { echo "ERROR: mzxml_dir not found: $mzxml_dir" >&2; exit 1; }

if ! command -v 7z &>/dev/null; then
  echo "ERROR: 7z not found — run 00_setup.sh to install p7zip" >&2; exit 1
fi

parent_dir="$(dirname "$mzxml_dir")"

# --- Require subfolders ---
missing=false
[[ -d "${mzxml_dir}/hilicpos" ]] || { echo "ERROR: ${mzxml_dir}/hilicpos not found." >&2; missing=true; }
[[ -d "${mzxml_dir}/c18neg"   ]] || { echo "ERROR: ${mzxml_dir}/c18neg not found."   >&2; missing=true; }
if [[ "$missing" == true ]]; then
  echo "  Run 01_organize.sh first to sort mzXML files into hilicpos/ and c18neg/." >&2
  exit 1
fi

hilicpos_n=$(find "${mzxml_dir}/hilicpos" -maxdepth 1 -type f -iname "*.mzxml" 2>/dev/null | wc -l | tr -d ' ')
c18neg_n=$(find   "${mzxml_dir}/c18neg"   -maxdepth 1 -type f -iname "*.mzxml" 2>/dev/null | wc -l | tr -d ' ')

echo ""
echo "==> mzXML files ready to compress:"
echo "    hilicpos/: $hilicpos_n files"
echo "    c18neg/:   $c18neg_n files"
echo "    Output:    $parent_dir"
echo ""

read -r -p "==> Compress into hilicpos.7z and c18neg.7z? [y/n]: " resp_compress
echo ""
[[ "$resp_compress" =~ ^[Yy]$ ]] || { echo "Aborted."; exit 0; }

echo "==> Archiving hilicpos..."
7z a "${parent_dir}/hilicpos.7z" "${mzxml_dir}/hilicpos/"
echo ""
echo "==> Archiving c18neg..."
7z a "${parent_dir}/c18neg.7z" "${mzxml_dir}/c18neg/"

# --- Generate zip_contents.csv ---
echo ""
echo "==> Writing zip_contents.csv..."
csv_path="${parent_dir}/zip_contents.csv"
printf 'folder\tsubfolder\tfile_name\n' > "$csv_path"

for subdir in hilicpos c18neg; do
  find "${mzxml_dir}/${subdir}" -maxdepth 1 -type f -iname "*.mzxml" 2>/dev/null | sort | while IFS= read -r f; do
    printf '%s\t\t%s\n' "$subdir" "$(basename "$f")"
  done >> "$csv_path"
done

file_count=$(( $(wc -l < "$csv_path" | tr -d ' ') - 1 ))
echo "    hilicpos.7z → ${parent_dir}/hilicpos.7z"
echo "    c18neg.7z   → ${parent_dir}/c18neg.7z"
echo "    Manifest    → $csv_path ($file_count files)"
echo ""
echo "==> Done. Run 03_copy.sh next to move archives to OneDrive."
