import pathlib, re
sdf = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
text = sdf.read_bytes().decode('latin-1')
records = [r.strip() for r in text.split('$$$$') if r.strip()]
print(f'Total records: {len(records)}')

def get_prop(rec, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    return m.group(1).strip() if m else ''

for rec in records:
    if get_prop(rec, 'm_abbr') == 'Hsl':
        print(f'Hsl chem_types: {get_prop(rec, "m_chem_types")}')
        print(f'Hsl rgroups:    {get_prop(rec, "m_Rgroups")}')
        break
