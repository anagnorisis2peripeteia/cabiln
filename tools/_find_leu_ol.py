import pathlib, re
sdf = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
text = sdf.read_bytes().decode('utf-8')
records = [r for r in text.split('$$$$') if r.strip()]

def gp(rec, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    return m.group(1).strip() if m else ''

# Find Leu_ol by abbr or symbol or name
for rec in records:
    abbr   = gp(rec, 'm_abbr')
    symbol = gp(rec, 'symbol')
    name   = gp(rec, 'm_name')
    if 'leu_ol' in abbr.lower() or 'leu_ol' in symbol.lower() or 'leucinol' in name.lower():
        print(f'abbr={abbr!r}  symbol={symbol!r}  name={name!r}')
        print(f'  m_type={gp(rec,"m_type")}  chem_types={gp(rec,"m_chem_types")}')
        print(f'  rgroups={gp(rec,"m_Rgroups")}')
        print()

# Also check for duplicates or near-duplicates
abbrs = [gp(r, 'm_abbr') for r in records]
from collections import Counter
dupes = {k:v for k,v in Counter(abbrs).items() if v > 1}
if dupes:
    print(f'Duplicate abbrs: {dupes}')
else:
    print(f'No duplicate abbrs. Total: {len(abbrs)} records.')
