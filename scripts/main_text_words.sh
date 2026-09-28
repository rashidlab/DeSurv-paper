#!/bin/bash
# Main-text word count (Introduction + Results + Discussion prose), one integer.
# The same count scripts/guards.sh reports; used as manuscript.yml word_count_command.
cd "$(dirname "$0")/.." >/dev/null
t=0
for f in 02_introduction_REVISED 04_results_REVISED 05_discussion_REVISED; do
  n=$(python3 code/util_sentence_extract.py paper/$f.Rmd --list | grep -oP '^\S+\s+\S+\s+\K\d+(?=w)' | paste -sd+ | bc)
  t=$((t+n))
done
echo $t
