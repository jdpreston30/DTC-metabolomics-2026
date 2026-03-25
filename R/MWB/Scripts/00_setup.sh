#!/usr/bin/env bash
# 00_setup.sh  (Step 0)
# One-time setup for Mac: checks p7zip for compression.
#
# Usage: bash R/MWB/Scripts/00_setup.sh

set -uo pipefail

# --- Check p7zip ---
echo "==> Checking p7zip..."
if command -v 7z &>/dev/null; then
  echo "    OK: already installed"
else
  echo "    Not found. Installing via Homebrew..."
  brew install p7zip
  echo "    OK: p7zip installed"
fi

echo ""
echo "Setup complete — you are ready to run 01_compress.sh"
echo ""
echo "==> R dependencies: (not yet implemented)"
