#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

quarto render reporting/EBS_wp_sporc.qmd --to html
mkdir -p docs
cp reporting/EBS_wp_sporc.html docs/index.html
touch docs/.nojekyll

printf 'Updated docs/index.html for GitHub Pages.\n'
