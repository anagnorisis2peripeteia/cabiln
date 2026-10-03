const { test, expect, render, tile, site, selectChip, connect, capture, isCompletedResponse, interceptOnce } = require('./fixtures');
const fs = require('node:fs/promises');

test('Build omits false attachment sites and disables unsupported chemistry', async ({ page }) => {
  const leftSite = slot => page.locator('#build-left-rgroups button')
    .filter({ hasText: new RegExp(`^R${slot} `) });
  await render(page, 'HArg');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'HArg');
  await expect(leftSite(4)).toBeDisabled();
  await expect(leftSite(4)).toContainText('unavailable');
  await expect(leftSite(2)).toBeEnabled();
  for (const symbol of ['Lys_Mtt', 'D_Lys_Mtt', 'Nspe', 'Nrpe']) {
    await render(page, symbol);
    await selectChip(page, 0, 'left', symbol);
    await expect(page.locator('#build-left-abbr')).toHaveText(symbol);
    await expect(leftSite(2)).toBeEnabled();
    await expect(leftSite(4)).toHaveCount(0);
  }
  await render(page, 'AZDye488');
  await selectChip(page, 0, 'left', 'AZDye488');
  await expect(leftSite(7)).toBeEnabled();
  await expect(leftSite(9)).toHaveCount(0);
  await expect(leftSite(10)).toHaveCount(0);
});

for (const width of [1440, 390]) test(`registration requires an explicit attachment choice at ${width}px`, async ({ page, app }, testInfo) => {
  await page.setViewportSize({ width, height: 900 });
  await page.goto(`${app.url}/register`);
  await page.locator('#smiles-in').fill('CNCC(CN)C(=O)O');
  await page.locator('#btn-preview').click();
  const choices = page.getByRole('radio');
  await expect(choices).toHaveCount(2);
  await page.locator('#abbr-in').fill(`BrowserChoice${width}`);
  await page.locator('#name-in').fill('Selected backbone');
  await expect(page.locator('#btn-register')).toBeDisabled();
  await choices.first().focus();
  await choices.first().press('Space');
  await expect(page.locator('#btn-register')).toBeEnabled();
  const first = await page.locator('#chuckles-out').inputValue();
  await page.locator('.attachment-choice').nth(1).locator('.choice-drawing').click();
  await expect(choices.nth(1)).toBeChecked();
  await expect(page.locator('#chuckles-out')).not.toHaveValue(first);
  await page.screenshot({ path: testInfo.outputPath('selected-backbone.png'), fullPage: true, animations: 'disabled' });
  await page.locator('#smiles-in').fill('NCC(=O)O');
  await expect(page.locator('#attachment-choices')).toBeHidden();
  await expect(page.locator('#btn-register')).toBeDisabled();
  await page.locator('#smiles-in').fill('CNCC(CN)C(=O)O');
  await page.locator('#btn-preview').click();
  await expect(choices).toHaveCount(2);
  await expect(page.locator('input[name="attachment-choice"]:checked')).toHaveCount(0);
  await choices.first().check();
  const selected = await page.locator('#chuckles-out').inputValue();
  const submitted = page.waitForRequest(request => new URL(request.url()).pathname === '/register_monomer');
  await page.locator('#btn-register').click();
  const payload = (await submitted).postDataJSON();
  expect(payload.chuckles).toBe(selected);
  expect(payload.activation_policy).toBe('canonical-sites-v1');
  await expect(page.locator('#status-msg')).toContainText('installed');
});

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
      await expect(tab).toHaveRole('button');
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
      await tab.focus();
      await page.keyboard.press('Space');
      expect(await page.locator('#render-inner .res-hl').count()).toBe(actual.count);
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
      await expect(page.locator('#cabiln-input')).toHaveValue(output === 'pct' ? '' : source);
      await page.locator('#cabiln-input').fill(source);
      await expect(page.locator('#render-inner svg')).toBeVisible();
      await expect(page.locator('#cabiln-input')).toHaveClass('ok');
      const converted = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/to_cabiln');
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
    const converted = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/convert_notation');
    await page.locator(`#btn-to-${target}`).click();
    const response = await converted;
    expect(response.status(), await response.text()).toBe(200);
    await expect(page.locator('#cabiln-input')).toHaveValue((await response.json()).result);
    await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    await expect(page.locator('#render-pane')).toHaveAttribute('aria-busy', 'false');
    await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
    await capture(page, testInfo, `verified-${target}`);
  }
  // Reference conversion has its own controls and response lifetime.
  await page.locator('#smiles-input').fill('NCC(=O)O');
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  for (const button of ['#btn-s2c', '#btn-s2c-bracket']) {
    const converted = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/smiles_to_cabiln');
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
  let previews = 0;
  await registration.route('**/preview_monomer', route => {
    if (++previews > 1) return route.continue();
    return route.fulfill({ status: 503, contentType: 'application/json',
      headers: { 'Retry-After': '1' }, body: JSON.stringify({ error: 'Chemistry capacity is busy' }) });
  });
  await registration.goto(`${app.url}/register`);
  await registration.locator('#smiles-in').fill('N[C@@H](CCCS)C(=O)O');
  await registration.locator('#btn-preview').click();
  await expect(registration.locator('#preview-canvas svg')).toBeVisible();
  expect(previews).toBeGreaterThanOrEqual(2);
  await expect(registration.locator('#detected-display')).toContainText('thiol');
  await expect(registration.locator('#btn-register')).toBeDisabled();
  await registration.locator('#abbr-in').fill('BrowserThiol');
  await registration.locator('#name-in').fill('Browser acceptance thiol');
  await expect(registration.locator('#btn-register')).toBeEnabled();
  await registration.locator('#btn-register').click();
  await expect(registration.locator('#status-msg')).toContainText('installed');
  await expect(registration.locator('#status-msg')).toBeVisible();
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
  const converted = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/smiles_to_cabiln');
  await page.locator('#btn-s2c-bracket').click();
  const response = await converted;
  expect(response.status(), await response.text()).toBe(200);
  expect((await response.json()).cabiln).toContain('BrowserThiol');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#residue-chips [data-residue]')).toContainText(['BrowserThiol']);
  await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
  await capture(page, testInfo, 'registered-built-reimported');
});

for (const width of [1440, 390]) test(`public monomer onboarding stays private and travels with projects at ${width}px`, async ({ page, context, readonlyApp }, testInfo) => {
  test.setTimeout(90_000);
  await page.setViewportSize({ width, height: 1000 });
  await page.goto(readonlyApp.url);
  await render(page, 'A-G');
  await page.locator('#btn-lib').click();
  await page.locator('#register-link').click();
  const form = page.frameLocator('#registration-frame');
  await expect(form.locator('#registration-destination')).toBeHidden();
  await expect(form.locator('#registration-scope')).toContainText('this tab');
  await form.locator('#smiles-in').fill('N[C@@H](CCCS)C(=O)O');
  await form.locator('#btn-preview').click();
  await expect(form.locator('#detected-display')).toContainText('thiol');
  const symbol = `TabThiol${width}`;
  await form.locator('#abbr-in').fill(symbol);
  await form.locator('#name-in').fill('Temporary thiol');
  await expect(form.locator('#btn-register')).toHaveText('Add to this tab');
  await form.locator('#btn-register').click();
  await expect(form.locator('#status-msg.ok')).toContainText(`${symbol} added to this tab`);
  await page.screenshot({ path: testInfo.outputPath('temporary-monomer-added.png'), animations: 'disabled' });
  await expect(form.locator('#registration-done')).toBeInViewport();
  await form.locator('#registration-done').click();
  await expect(page.locator('#registration-dialog')).not.toBeVisible();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator(`.lib-row[data-abbr="${symbol}"]`)).toBeVisible();
  await render(page, symbol);
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', symbol);
  await site(page, 'left', 4);
  await tile(page, 'C', 'right');
  await site(page, 'right', 4);
  await connect(page);
  const source = await page.locator('#cabiln-input').inputValue();
  const download = page.waitForEvent('download');
  await page.locator('#btn-project-save').click();
  const project = JSON.parse(await fs.readFile(await (await download).path(), 'utf8'));
  expect(project.monomers.map(item => item.abbr)).toEqual([symbol]);
  expect(await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'))).not.toContain(symbol);
  await page.evaluate(() => {
    const saved = JSON.parse(sessionStorage.getItem('cabiln.monomers.v1'));
    saved.token = 'expired-server-token';
    sessionStorage.setItem('cabiln.monomers.v1', JSON.stringify(saved));
  });
  await page.reload();
  await page.locator('#btn-restore-draft').click();
  await expect(page.locator('#cabiln-input')).toHaveValue(source);
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');

  const other = await context.newPage();
  await other.goto(readonlyApp.url);
  await other.locator('#btn-lib').click();
  await other.locator('#lib-search').fill(symbol);
  await expect(other.locator('#lib-list')).toHaveText('No matches');
  await other.locator('#project-upload').setInputFiles({
    name: 'custom.cabiln.json', mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(project)),
  });
  await expect(other.locator('#project-status')).toContainText('Project opened');
  await expect(other.locator('#cabiln-input')).toHaveValue(source);
  await expect(other.locator('#cabiln-input')).toHaveClass('ok');
  await expect(other.locator(`.lib-row[data-abbr="${symbol}"]`)).toBeVisible();
  const shared = await page.request.get(`${readonlyApp.url}/monomers`);
  expect((await shared.json()).some(item => item.abbr === symbol)).toBe(false);
  await other.close();
});

test('registration explains empty input, failed previews and duplicate names with keyboard recovery', async ({ page, app }) => {
  await page.goto(`${app.url}/register`);
  await page.locator('#btn-preview').click();
  await expect(page.locator('#smiles-in')).toBeFocused();
  await expect(page.locator('#smiles-in')).toHaveAttribute('aria-invalid', 'true');
  await expect(page.locator('#preview-error')).toContainText('Enter a SMILES string');
  await page.locator('#smiles-in').fill('NCC(=O)O');
  await interceptOnce(page, '**/preview_monomer', route => route.abort('failed'));
  await page.locator('#btn-preview').click();
  await expect(page.locator('#preview-error')).toContainText('try Preview again');
  await expect(page.locator('#smiles-in')).toHaveValue('NCC(=O)O');
  await expect(page.locator('#btn-register')).toBeHidden();
  await page.locator('#btn-preview').click();
  await expect(page.locator('#preview-canvas svg')).toBeVisible();
  await expect(page.locator('#preview-error')).toBeEmpty();
  await expect(page.locator('#btn-register')).toBeDisabled();
  await expect(page.locator('#registration-guide')).toContainText('abbreviation and full name');
  await page.locator('#abbr-in').fill('G');
  await page.locator('#name-in').fill('Keyboard glycine');
  await page.locator('#name-in').press('Enter');
  await expect(page.locator('#status-msg.err')).toContainText(/already|exists/i);
  await expect(page.locator('#preview-canvas svg')).toBeVisible();
  await page.locator('#abbr-in').fill('KeyboardGly');
  await page.locator('#abbr-in').press('Enter');
  await expect(page.locator('#status-msg.ok')).toContainText('KeyboardGly installed');
  await expect(page.locator('#btn-register')).toBeDisabled();
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
      return isCompletedResponse(response) && url.pathname === '/monomer_rgroups' && url.searchParams.get('abbr') === symbol;
    });
    await page.locator('#build-right-rgroups').getByRole('button', { name: `${label}: ${symbol}`, exact: true }).click();
    const response = await sites;
    expect(response.status(), await response.text()).toBe(200);
    const data = await response.json();
    expect(data.rgroups.length).toBeGreaterThan(0);
    await expect(page.locator('#build-right-abbr')).toHaveText(symbol);
    await expect(page.locator('#build-right-rgroups button')).toHaveText(data.rgroups.map(site => `R${site.slot} ${(site.chem_type || '').replaceAll('_', ' ')}${site.used ? ' · used' : ''}`));
    await capture(page, testInfo, `attachment-form-${label}`);
  }
});

test('partially recognized input retains warnings, selectable unknowns and backbone editing', async ({ page }, testInfo) => {
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('CC(=O)N[C@@H](CCC(F)(F)F)C(=O)N');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  const converted = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/to_cabiln');
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
