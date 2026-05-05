import subprocess, re
result = subprocess.run(
    ['git', 'show', 'HEAD:src/pyPept/data/monomers.sdf'],
    cwd=r'C:\Users\Cameron\repos\pyPept',
    capture_output=True
)
text = result.stdout.decode('latin-1')
records = [r.strip() for r in text.split('$$$$') if r.strip()]
print(f'Git HEAD records: {len(records)}')

def get_prop(rec, name):
    m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    return m.group(1).strip() if m else ''

has_ct = sum(1 for r in records if get_prop(r, 'm_chem_types'))
print(f'With chem_types: {has_ct}')
print(f'Without: {len(records) - has_ct}')
dyes = [get_prop(r, 'm_abbr') for r in records if get_prop(r, 'm_subtype') == 'dye']
print(f'Dye entries: {dyes}')
