const { test, expect, render, tile, site, selectChip, connect, capture } = require('./fixtures');

const tabCases = [
  { name: 'later segment', source: 'A%K.{G(4,2)}-A', occurrences: 4, markers: 0, groups: [[3]] },
  { name: 'separate bracket policies', source: 'ac-K.[G(4,2)]-K.{A(4,2)}-am', occurrences: 6, markers: 0, groups: [[4], [5]] },
  { name: 'nested protected arm', source: 'ac-K.{G(4,2)[.A(1,2)]}-D-am', occurrences: 6, markers: 0, groups: [[4, 5], [5]] },
  { name: 'nested sibling', source: 'ac-K.[[K(4,2).A(1,2)].ac(4,2)]-am', occurrences: 6, markers: 0 },
  { name: 'terminal cycle links', source: '!r-C-A-C-!r', occurrences: 3, markers: 2 },
  { name: 'nested arm crosslink', source: 'ac-K.{G(4,2)[.A(1,2).!1(1,4)]}-D.!1(4,1)-am', occurrences: 6, markers: 2 },
  { name: 'marker-only nested group', source: 'ac-K.[G(4,2)[.!1(1,4)]]-D.!1(4,1)-am', occurrences: 5, markers: 2 },
  { name: 'protected link-only group', source: 'C.!r(1,2).!s(4,4)-A-C.{!s(4,4).!r(2,1)}', occurrences: 3, markers: 4 },
  { name: 'three-port scaffold', source: 'ac-C.[TBMB(4,4).!2(5,4).!3(6,4)]-C.!2-C.!3-am', occurrences: 6, markers: 4 },
];

for (const example of tabCases) {
  test(`tabs highlight their atoms: ${example.name}`, async ({ page }, testInfo) => {
    const data = await render(page, example.source);
    await expect(page.locator('#residue-chips [data-residue]')).toHaveCount(example.occurrences);
    await expect(page.locator('#residue-chips .xlink-chip')).toHaveCount(example.markers);
    const actualIds = await page.locator('#residue-chips [data-residue]').evaluateAll(elements =>
      elements.map(el => Number(el.dataset.residue)).sort((a, b) => a - b));
    expect(actualIds).toEqual(Array.from({ length: example.occurrences }, (_, idx) => idx));
    if (example.groups) expect(data.layout.groups.map(group => group.members)).toEqual(example.groups);
    const tabs = page.locator('#residue-chips .branch-chip, #residue-chips .xlink-chip');
    expect(await tabs.count()).toBeGreaterThan(0);
    for (let idx = 0; idx < await tabs.count(); idx++) {
      const tab = tabs.nth(idx);
      const members = JSON.parse(await tab.getAttribute('data-members'));
      if (!members.length) continue;
      await tab.hover();
      const highlighted = await page.locator('#residue-chips [data-residue].hover').evaluateAll(elements =>
        elements.map(el => Number(el.dataset.residue)).sort((a, b) => a - b));
      expect(highlighted).toEqual([...new Set(members)].sort((a, b) => a - b));
      const actual = await page.locator('#render-inner svg').evaluate((svg, { members, ownership }) => {
        const atoms = new Set(members.flatMap(id => ownership[String(id)]));
        const paths = [...svg.querySelectorAll('[class]')].filter(el => /(?:^| )atom-\d+(?: |$)/.test(el.getAttribute('class')));
        return {
          count: svg.querySelectorAll('.res-hl').length,
          correct: paths.every(el => {
            const ids = [...el.getAttribute('class').matchAll(/(?:^| )atom-(\d+)(?= |$)/g)].map(match => Number(match[1]));
            return el.classList.contains('res-hl') === ids.some(id => atoms.has(id));
          }),
        };
      }, { members, ownership: data.residue_map });
      expect(actual.count).toBeGreaterThan(0);
      expect(actual.correct).toBe(true);
    }
    await capture(page, testInfo, 'tabs-and-atom-ownership');
  });
}

test('Examples populate the editable renderer without closing the panel', async ({ page }, testInfo) => {
  await page.locator('#btn-examples').click();
  const example = page.locator('.example-row').filter({ has: page.getByText('Head-to-tail cyclic', { exact: true }) });
  await example.click();
  await expect(page.locator('#cabiln-input')).toHaveValue('!1-A-K-G-E-L-F-!1');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#residue-chips [data-residue]')).toHaveCount(6);
  await expect(page.locator('#examples-panel')).toHaveClass('open');
  await capture(page, testInfo, 'example-selected');
});

test('all main input selectors convert through both output controls', async ({ page }, testInfo) => {
  for (const [format, source] of [
    ['smiles', 'NCC(=O)N[C@@H](C)C(=O)O'],
    ['biln', 'G-A'],
    ['helm', 'PEPTIDE1{G.A}$$$$'],
  ]) {
    for (const output of ['pct', 'bracket']) {
      await page.locator('#notation-select').selectOption(format);
      await expect(page.locator('#cabiln-input')).toHaveValue('');
      await page.locator('#cabiln-input').fill(source);
      await expect(page.locator('#render-inner svg')).toBeVisible();
      await expect(page.locator('#cabiln-input')).toHaveClass('ok');
      const converted = page.waitForResponse(response => new URL(response.url()).pathname === '/to_cabiln');
      await page.locator(`#btn-to-cabiln-${output}`).click();
      const response = await converted;
      expect(response.status(), await response.text()).toBe(200);
      await expect(page.locator('#notation-select')).toHaveValue('cabiln');
      await expect(page.locator('#cabiln-input')).toHaveValue('G-A');
      await expect(page.locator('#cabiln-input')).toHaveClass('ok');
      await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['G', 'A']);
      await expect(page.locator('#btn-mol')).toBeEnabled();
      await capture(page, testInfo, `${format}-to-${output}`);
    }
  }
});

test('Verify, MOL/SDF upload and both notation formatters preserve chemistry', async ({ page }, testInfo) => {
  const data = await render(page, 'ac-K.[G(4,2)]-am');
  await page.locator('#btn-verify').click();
  for (const extension of ['mol', 'sdf']) {
    await page.locator('#mol-upload').setInputFiles({
      name: `peptide.${extension}`, mimeType: 'chemical/x-mdl-molfile',
      buffer: Buffer.from(data.mol_block + (extension === 'sdf' ? '\n$$$$\n' : '')),
    });
    await expect(page.locator('#smiles-inner svg')).toBeVisible();
    await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
    await expect(page.locator('#smiles-input')).not.toHaveValue('');
  }
  for (const target of ['branch', 'bracket']) {
    const converted = page.waitForResponse(response => new URL(response.url()).pathname === '/convert_notation');
    await page.locator(`#btn-to-${target}`).click();
    const response = await converted;
    expect(response.status(), await response.text()).toBe(200);
    await expect(page.locator('#cabiln-input')).toHaveValue((await response.json()).result);
    await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
    await capture(page, testInfo, `verified-${target}`);
  }
  // Reference conversion has its own controls and response lifetime.
  await page.locator('#smiles-input').fill('NCC(=O)O');
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  for (const button of ['#btn-s2c', '#btn-s2c-bracket']) {
    const converted = page.waitForResponse(response => new URL(response.url()).pathname === '/smiles_to_cabiln');
    await page.locator(button).click();
    expect((await converted).status()).toBe(200);
    await expect(page.locator('#cabiln-input')).toHaveValue('G');
    await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['G']);
    await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
  }
});

test('registration refreshes the open library and supplies detected sites to building and recognition', async ({ page, context, app }, testInfo) => {
  await page.locator('#btn-build').click();
  await page.locator('#lib-search').fill('BrowserThiol');
  await expect(page.locator('#lib-list')).toHaveText('No matches');
  const registration = await context.newPage();
  await registration.goto(`${app.url}/register`);
  await registration.locator('#smiles-in').fill('N[C@@H](CCCS)C(=O)O');
  await registration.locator('#btn-preview').click();
  await expect(registration.locator('#preview-canvas svg')).toBeVisible();
  await expect(registration.locator('#detected-display')).toContainText('thiol');
  await expect(registration.locator('#btn-register')).toBeDisabled();
  await registration.locator('#abbr-in').fill('BrowserThiol');
  await registration.locator('#name-in').fill('Browser acceptance thiol');
  await expect(registration.locator('#btn-register')).toBeEnabled();
  await registration.locator('#btn-register').click();
  await expect(registration.locator('#status-msg')).toContainText('registered successfully');
  await expect(registration.locator('#btn-register')).toBeDisabled();
  await registration.screenshot({ path: testInfo.outputPath('registered-monomer.png'), animations: 'disabled' });
  await registration.close();
  // This isolated headless context has no OS window activation; dispatch the
  // same focus event that returning to the renderer tab delivers.
  await page.evaluate(() => window.dispatchEvent(new Event('focus')));
  await expect(page.locator('.lib-row[data-abbr="BrowserThiol"]')).toBeVisible();
  await tile(page, 'BrowserThiol');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await selectChip(page, 0, 'left', 'BrowserThiol');
  await site(page, 'left', 4);
  await page.locator('#btn-rxn-filter').click();
  await tile(page, 'C', 'right');
  await site(page, 'right', 4);
  await connect(page);
  await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['BrowserThiol', 'C']);
  const downloadPromise = page.waitForEvent('download');
  await page.locator('#btn-mol').click();
  const download = await downloadPromise;
  const productPath = testInfo.outputPath('custom-product.mol');
  await download.saveAs(productPath);
  await page.locator('#btn-verify').click();
  await page.locator('#mol-upload').setInputFiles(productPath);
  await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
  const converted = page.waitForResponse(response => new URL(response.url()).pathname === '/smiles_to_cabiln');
  await page.locator('#btn-s2c-bracket').click();
  const response = await converted;
  expect(response.status(), await response.text()).toBe(200);
  expect((await response.json()).cabiln).toContain('BrowserThiol');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#residue-chips [data-residue]')).toContainText(['BrowserThiol']);
  await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
  await capture(page, testInfo, 'registered-built-reimported');
});

test('a read-only instance exposes rendering while registration stays unavailable', async ({ page, readonlyApp }) => {
  await page.goto(readonlyApp.url);
  await expect(page.locator('#register-link')).toBeHidden();
  await render(page, 'A-G');
  const response = await page.goto(`${readonlyApp.url}/register`);
  expect(response.status()).toBe(403);
  await expect(page.locator('body')).toHaveText('This monomer library is read-only.');
});

test('attachment-form choices load the selected family member before its sites', async ({ page }, testInfo) => {
  await render(page, 'G');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'G');
  for (const [label, symbol] of [['N-terminal', 'Bn_'], ['C-terminal', '_Bn']]) {
    await tile(page, 'Bn', 'right');
    await expect(page.locator('#build-right-rgroups button')).toHaveText(['N-terminal: Bn_', 'C-terminal: _Bn']);
    await expect(page.locator('#build-connect')).toBeDisabled();
    const sites = page.waitForResponse(response => {
      const url = new URL(response.url());
      return url.pathname === '/monomer_rgroups' && url.searchParams.get('abbr') === symbol;
    });
    await page.locator('#build-right-rgroups').getByRole('button', { name: `${label}: ${symbol}`, exact: true }).click();
    const response = await sites;
    expect(response.status(), await response.text()).toBe(200);
    const data = await response.json();
    expect(data.rgroups.length).toBeGreaterThan(0);
    await expect(page.locator('#build-right-abbr')).toHaveText(symbol);
    await expect(page.locator('#build-right-rgroups button')).toHaveText(data.rgroups.map(site => `R${site.slot} ${site.chem_type || ''}`));
    await capture(page, testInfo, `attachment-form-${label}`);
  }
});

test('partially recognized input retains warnings, selectable unknowns and backbone editing', async ({ page }, testInfo) => {
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('CC(=O)N[C@@H](CCC(F)(F)F)C(=O)N');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  const converted = page.waitForResponse(response => new URL(response.url()).pathname === '/to_cabiln');
  await page.locator('#btn-to-cabiln-bracket').click();
  const response = await converted;
  expect(response.status(), await response.text()).toBe(200);
  const data = await response.json();
  expect(data.recognition_status).toBe('partial');
  expect(data.assignments.map(item => item.recognized)).toEqual([true, false, true]);
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#conversion-status')).toBeVisible();
  await expect(page.locator('#conversion-status')).toHaveText(data.warning);
  const syntheticLabel = await page.locator('#residue-chips [data-residue="1"]').innerText();
  expect(syntheticLabel).toMatch(/^__SYN_[0-9a-f]+$/);
  await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['ac', syntheticLabel, 'am']);
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', syntheticLabel);
  await expect(page.locator('#build-left-rgroups button').filter({ hasText: /^R1 / })).toHaveClass(/used/);
  await expect(page.locator('#build-left-rgroups button').filter({ hasText: /^R2 / })).toHaveClass(/used/);
  await capture(page, testInfo, 'partial-conversion-selection');
  await selectChip(page, 2, 'right', 'am');
  await page.locator('#build-insert-btn').click();
  await tile(page, 'G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['ac', syntheticLabel, 'G', 'am']);
  await expect(page.locator('#conversion-status')).toBeHidden();
  await page.locator('#btn-verify').click();
  await page.locator('#smiles-input').fill('CC(=O)N[C@@H](CCC(F)(F)F)C(=O)NCC(N)=O');
  await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
  await capture(page, testInfo, 'partial-conversion-edited');
});
