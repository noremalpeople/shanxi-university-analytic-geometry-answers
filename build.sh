#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for n in 01 02 03; do
  cd "$ROOT/chapters/$n"
  mkdir -p build
  if [[ "$n" == "01" ]]; then
    tex=chapters/chapter-01.tex
  elif [[ "$n" == "02" ]]; then
    tex=chapters/chapter-02.tex
    python3 scripts/build_figures.py
  else
    tex=chapters/chapter-03.tex
  fi
  xelatex -interaction=nonstopmode -halt-on-error -output-directory=build "$tex" >/dev/null
  xelatex -interaction=nonstopmode -halt-on-error -output-directory=build "$tex" >/dev/null
  if [[ "$n" == "01" ]]; then
    cp build/chapter-01.pdf "$ROOT/dist/chapter-01-vector-algebra-illustrated.pdf"
  elif [[ "$n" == "02" ]]; then
    cp build/chapter-02.pdf "$ROOT/dist/chapter-02-spatial-geometry-complete.pdf"
  else
    cp build/chapter-03.pdf "$ROOT/dist/chapter-03-coordinate-transform-conics-illustrated.pdf"
  fi
done
echo "Built PDFs in $ROOT/dist"
