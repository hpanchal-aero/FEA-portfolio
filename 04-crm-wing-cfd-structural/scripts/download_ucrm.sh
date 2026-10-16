#!/usr/bin/env bash
# Download, verify and extract the uCRM-9 dataset (Brooks, Kenway & Martins, AIAA J. 56(7), 2018).
# Source: https://data.mendeley.com/datasets/gpk4zn73xn  (DOI 10.17632/gpk4zn73xn.1)
# No licence is stated on the dataset page; the files are therefore NOT committed to Git.
set -euo pipefail

URL="https://data.mendeley.com/public-files/datasets/gpk4zn73xn/files/bc64a6be-82a7-442c-978e-442f76a87c72/file_downloaded"
SHA256="231ef07e23d25062b8b0f8595cb4fe1d294d008b84e9b4003a555311c1c82b50"
ZIP="uCRM9.zip"

cd "$(dirname "$0")/../data"

# Re-use a valid existing download; otherwise (re)download.
if [ -f "$ZIP" ] && echo "${SHA256}  ${ZIP}" | sha256sum -c --quiet - 2>/dev/null; then
    echo "ZIP_PRESENT_AND_VERIFIED"
else
    avail_gb=$(df -BG --output=avail . | tail -1 | tr -dc '0-9')
    if [ "$avail_gb" -lt 3 ]; then
        echo "ABORT: only ${avail_gb} GB free in data/ (need at least 3 GB)." >&2
        exit 1
    fi
    wget -c -nv -O "$ZIP" "$URL"
    echo "${SHA256}  ${ZIP}" | sha256sum -c -
fi

mkdir -p ucrm9
python3 - "$ZIP" ucrm9 <<'PY'
import sys, zipfile
zf = zipfile.ZipFile(sys.argv[1])
bad = zf.testzip()
if bad is not None:
    sys.exit(f"ABORT: corrupt member {bad}")
zf.extractall(sys.argv[2])
print(f"EXTRACTED_OK ({len(zf.namelist())} members)")
PY

du -sh ucrm9
find ucrm9 -type f -printf "%10s  %p\n" | sort -k2 | head -60
