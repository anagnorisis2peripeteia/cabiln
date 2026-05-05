#!/usr/bin/env python3
"""Find and fix non-UTF-8 bytes in SDF property values, and report m_Rgroups drift."""
import subprocess, re, pathlib

SDF = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'

raw = SDF.read_bytes()

# Find all non-UTF-8 byte positions
bad_positions = []
i = 0
while i < len(raw):
    try:
        raw[i:i+1].decode('utf-8')
        i += 1
    except UnicodeDecodeError:
        bad_positions.append(i)
        i += 1

print(f'Non-UTF-8 bytes: {len(bad_positions)}')
for pos in bad_positions[:20]:
    ctx = raw[max(0,pos-20):pos+20]
    print(f'  pos {pos}: byte=0x{raw[pos]:02x}  context={ctx!r}')

# Fix: replace soft hyphens (0xad) and other cp1252 chars with UTF-8 equivalents or plain ASCII
fixed = raw.decode('latin-1').encode('utf-8')
SDF.write_bytes(fixed)

# Verify
try:
    fixed.decode('utf-8')
    print(f'\nSDF re-encoded to UTF-8 successfully ({len(fixed)} bytes, was {len(raw)} bytes).')
except UnicodeDecodeError as e:
    print(f'ERROR: still not valid UTF-8: {e}')

# Also check m_Rgroups drift vs git HEAD
print('\n--- m_Rgroups drift check ---')
git_bytes = subprocess.run(
    ['git', 'show', 'HEAD:src/pyPept/data/monomers.sdf'],
    cwd=str(SDF.parent.parent.parent),
    capture_output=True
).stdout
git_text  = git_bytes.decode('latin-1')
curr_text = fixed.decode('utf-8')

def parse(text):
    result = {}
    for rec in text.split('$$$$'):
        def gp(name):
            m = re.search(r'> +<' + re.escape(name) + r'>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
            return m.group(1).strip() if m else ''
        abbr = gp('m_abbr')
        if abbr:
            result[abbr] = {'rg': gp('m_Rgroups'), 'ct': gp('m_chem_types')}
    return result

git_data  = parse(git_text)
curr_data = parse(curr_text)

changed = []
for abbr, curr in curr_data.items():
    git = git_data.get(abbr)
    if git is None or git['ct']:
        continue  # new monomer or already had chem_types — not reprocessed
    if git['rg'] != curr['rg']:
        changed.append((abbr, git['rg'], curr['rg']))

print(f'Reprocessed monomers with changed m_Rgroups: {len(changed)}')
for abbr, old, new in changed[:30]:
    print(f'  {abbr}:')
    print(f'    was: {old}')
    print(f'    now: {new}')
