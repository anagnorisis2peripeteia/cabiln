#!/usr/bin/env python3
"""Check current state of synonyms column in CSV."""
import csv

CSV_PATH = 'src/pyPept/data/monomers.csv'
rows = list(csv.DictReader(open(CSV_PATH, newline='')))
has_syn = [r for r in rows if r.get('synonyms', '').strip()]
print(f"Rows with synonyms: {len(has_syn)} / {len(rows)}")
for r in has_syn[:30]:
    print(f"  {r['token']:20s} synonyms={r['synonyms']}")
