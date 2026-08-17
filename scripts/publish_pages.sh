#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

quarto render reporting/EBS_wp_sporc.qmd --to html
mkdir -p docs
cp reporting/EBS_wp_sporc.html docs/index.html
if [ -f reporting/EBS_wp_sporc.pdf ]; then
  cp reporting/EBS_wp_sporc.pdf docs/EBS_wp_sporc.pdf
fi
touch docs/.nojekyll

printf 'Updated docs/index.html for GitHub Pages.\n'
