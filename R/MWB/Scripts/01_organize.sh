#!/usr/bin/env bash
# 01_organize.sh  (Step 1)
# Organizes .raw and .mzXML files into hilicpos/ and c18neg/ subfolders.
# For flat (unsorted) directories, sorts files by sequential acquisition order
# using odd/even position from sorted filename order.
#
# Usage: bash R/MWB/Scripts/01_organize.sh
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

parent_dir="$(dirname "$mzxml_dir")"

# _odd_even_split <src_dir> <ext>
# Moves flat files with <ext> in <src_dir> into hilicpos/ and c18neg/ subfolders.
# Files are sorted by name (reflects acquisition order), then split by odd/even
# sequential position (1,3,5,... vs 2,4,6,...). User confirms which group is hilicpos.
_odd_even_split() {
  local src_dir="$1"
  local ext="$2"

  local tmp
  tmp=$(mktemp)
  find "$src_dir" -maxdepth 1 -type f -iname "*.$ext" 2>/dev/null | sort > "$tmp"

  local total
  total=$(wc -l < "$tmp" | tr -d ' ')

  if [[ "$total" -eq 0 ]]; then
    echo "    (No .$ext files to sort.)"
    rm -f "$tmp"; return
  fi

  local odd_n even_n
  odd_n=$(( (total + 1) / 2 ))
  even_n=$(( total / 2 ))

  echo "    $total file(s) found."
  echo "      Odd  positions (1, 3, 5, ...): $odd_n files"
  echo "      Even positions (2, 4, 6, ...): $even_n files"
  echo ""
  read -r -p "    Assign odd-position files to hilicpos/ ? (n = even to hilicpos) [y/n]: " resp_odd
  echo ""

  if ! [[ "$resp_odd" =~ ^[YyNn]$ ]]; then
    echo "    Invalid input — skipping sort."
    rm -f "$tmp"; return
  fi

  local hilicpos_gets_odd=true
  [[ "$resp_odd" =~ ^[Nn]$ ]] && hilicpos_gets_odd=false

  mkdir -p "${src_dir}/hilicpos"
  mkdir -p "${src_dir}/c18neg"

  local pos=0 moved=0 dest fname
  while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    pos=$((pos + 1))
    fname="$(basename "$f")"
    if [[ $((pos % 2)) -eq 1 ]]; then
      [[ "$hilicpos_gets_odd" == true ]] && dest="${src_dir}/hilicpos" || dest="${src_dir}/c18neg"
    else
      [[ "$hilicpos_gets_odd" == true ]] && dest="${src_dir}/c18neg" || dest="${src_dir}/hilicpos"
    fi
    mv "$f" "${dest}/${fname}"
    moved=$((moved + 1))
  done < "$tmp"
  rm -f "$tmp"

  local hilicpos_n c18neg_n
  hilicpos_n=$(find "${src_dir}/hilicpos" -maxdepth 1 -type f -iname "*.$ext" 2>/dev/null | wc -l | tr -d ' ')
  c18neg_n=$(find   "${src_dir}/c18neg"   -maxdepth 1 -type f -iname "*.$ext" 2>/dev/null | wc -l | tr -d ' ')
  echo "    Done: hilicpos/=$hilicpos_n, c18neg/=$c18neg_n (total moved: $moved)"
}

# ============================================================
# SECTION A: .raw files  (parent of mzxml_dir)
# ============================================================
echo ""
echo "========================================"
echo " Section A: .raw files"
echo "========================================"
echo "    Location: $parent_dir"
echo ""

raw_flat=0
raw_flat=$(find "$parent_dir" -maxdepth 1 -type f -iname "*.raw" 2>/dev/null | wc -l | tr -d ' ') || true

raw_hilicpos=0
raw_c18neg=0
[[ -d "${parent_dir}/hilicpos" ]] && \
  raw_hilicpos=$(find "${parent_dir}/hilicpos" -maxdepth 1 -type f -iname "*.raw" 2>/dev/null | wc -l | tr -d ' ') || true
[[ -d "${parent_dir}/c18neg" ]] && \
  raw_c18neg=$(find   "${parent_dir}/c18neg"   -maxdepth 1 -type f -iname "*.raw" 2>/dev/null | wc -l | tr -d ' ') || true

echo "    Flat (unsorted) in parent:  $raw_flat .raw"
echo "    Already in hilicpos/:       $raw_hilicpos .raw"
echo "    Already in c18neg/:         $raw_c18neg .raw"

if [[ "$raw_flat" -gt 0 ]]; then
  echo ""
  read -r -p "==> Sort $raw_flat .raw file(s) into hilicpos/ and c18neg/ by acquisition order? [y/n]: " resp_raw
  echo ""
  if [[ "$resp_raw" =~ ^[Yy]$ ]]; then
    _odd_even_split "$parent_dir" "raw"
  else
    echo "    Skipped."
  fi
elif [[ "$raw_hilicpos" -gt 0 || "$raw_c18neg" -gt 0 ]]; then
  echo "    Raw files already organized into subfolders."
else
  echo "    No .raw files found — skipping."
fi

# ============================================================
# SECTION B: .mzXML files  (mzxml_dir)
# ============================================================
echo ""
echo "========================================"
echo " Section B: .mzXML files"
echo "========================================"
echo "    Location: $mzxml_dir"
echo ""

mzxml_flat=0
mzxml_flat=$(find "$mzxml_dir" -maxdepth 1 -type f -iname "*.mzxml" 2>/dev/null | wc -l | tr -d ' ') || true

mzxml_hilicpos=0
mzxml_c18neg=0
[[ -d "${mzxml_dir}/hilicpos" ]] && \
  mzxml_hilicpos=$(find "${mzxml_dir}/hilicpos" -maxdepth 1 -type f -iname "*.mzxml" 2>/dev/null | wc -l | tr -d ' ') || true
[[ -d "${mzxml_dir}/c18neg" ]] && \
  mzxml_c18neg=$(find   "${mzxml_dir}/c18neg"   -maxdepth 1 -type f -iname "*.mzxml" 2>/dev/null | wc -l | tr -d ' ') || true

echo "    Flat (unsorted) in mzxml_dir: $mzxml_flat .mzXML"
echo "    Already in hilicpos/:         $mzxml_hilicpos .mzXML"
echo "    Already in c18neg/:           $mzxml_c18neg .mzXML"

if [[ "$mzxml_flat" -gt 0 ]]; then
  echo ""
  read -r -p "==> Sort $mzxml_flat .mzXML file(s) into hilicpos/ and c18neg/ by acquisition order? [y/n]: " resp_mzxml
  echo ""
  if [[ "$resp_mzxml" =~ ^[Yy]$ ]]; then
    _odd_even_split "$mzxml_dir" "mzxml"
  else
    echo "    Skipped."
  fi
elif [[ "$mzxml_hilicpos" -gt 0 || "$mzxml_c18neg" -gt 0 ]]; then
  echo "    mzXML files already organized into subfolders."
else
  echo "    No .mzXML files found — skipping."
fi

echo ""
echo "==> Organization complete. Run 02_compress.sh next."
