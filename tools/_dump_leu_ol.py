import pathlib, re
sdf = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
text = sdf.read_bytes().decode('utf-8')
records = text.split('$$$$')
for rec in records:
    m = re.search(r'> +<m_abbr>.*?\r?\n(.+?)(?:\r?\n|$)', rec)
    if m and m.group(1).strip() == 'Leu_ol':
        print(repr(rec[:2000]))
        break
