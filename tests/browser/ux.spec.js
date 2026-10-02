const { test, expect, render, selectChip, site, connect, capture, isCompletedResponse, residuePoint } = require('./fixtures');

test('an unknown monomer points to its exact text and explains recovery', async ({ page }) => {
  await render(page, 'A-G');
  const input = page.locator('#cabiln-input');
  await input.fill('A-G-NotAMonomer-G');
  await expect(input).toHaveClass('err');
  await expect(page.locator('#cabiln-status')).toContainText('Library');
  await page.getByRole('button', { name: 'Select problem in sequence', exact: true }).click();
  expect(await input.evaluate(el => el.value.slice(el.selectionStart, el.selectionEnd))).toBe('NotAMonomer');
  await expect(input).toBeFocused();
  await input.fill('  A-🧬-G  ');
  await expect(page.locator('#cabiln-status')).toContainText('🧬');
  await page.getByRole('button', { name: 'Select problem in sequence', exact: true }).click();
  expect(await input.evaluate(el => el.value.slice(el.selectionStart, el.selectionEnd))).toBe('🧬');
  await input.fill('A-G-A-G');
  await expect(input).toHaveClass('ok');
  await expect(page.getByRole('button', { name: 'Select problem in sequence', exact: true })).toHaveCount(0);
});

test('a failed connection check explains the selected sites and retries without losing them', async ({ page }) => {
  await render(page, 'A-G');
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', 'G');
  await page.locator('#lib-search').fill('am');
  await page.locator('.lib-row[data-abbr="am"] .lib-use').click();
  await site(page, 'left', 2);
  await page.route('**/validate_bond', route => route.abort('failed'));
  await site(page, 'right', 1);
  const status = page.locator('#build-status');
  await expect(status).toContainText('G R2');
  await expect(status).toContainText('am R1');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await page.unroute('**/validate_bond');
  await status.getByRole('button', { name: 'Retry', exact: true }).click();
  await expect(page.locator('#build-connect')).toBeEnabled();
  await expect(page.locator('#build-left-rgroups button.selected')).toHaveText(/^R2 /);
  await expect(page.locator('#build-right-rgroups button.selected')).toHaveText(/^R1 /);
  await page.route('**/insert_bond', route => route.abort('failed'));
  await page.locator('#build-connect').click();
  await expect(status).toContainText('Could not connect G R2 with am R1');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await page.unroute('**/insert_bond');
  await status.getByRole('button', { name: 'Retry', exact: true }).click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G-am');
});

test('invalid bracket syntax keeps the last valid drawing without accepting changed chemistry', async ({ page }) => {
  for (const source of [
    'ac-K.[[K(4,2).A(1,2)garbage].ac(4,2)]-am',
    'ac-K.[!x(4,1)]-D.!x(4,4)-am',
  ]) {
    await render(page, 'A-G');
    const drawing = await page.locator('#render-inner').innerHTML();
    await page.locator('#cabiln-input').fill(source);
    await expect(page.locator('#cabiln-input')).toHaveClass('err');
    await expect(page.locator('#cabiln-input')).toHaveValue(source);
    expect(await page.locator('#render-inner').innerHTML()).toBe(drawing);
    await expect(page.locator('#render-progress-label')).toContainText('Previous drawing');
    await expect(page.locator('#btn-mol')).toBeDisabled();
    await expect(page.locator('#residue-chips [data-residue]').first()).toBeDisabled();
  }
});

test('editing retains the drawing and zoom while stale atom controls stay disabled', async ({ page }, testInfo) => {
  const data = await render(page, 'A-G');
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', 'G');
  await page.locator('#render-canvas').hover();
  await page.mouse.wheel(0, -100);
  const before = await page.locator('#render-inner').evaluate(el => ({ svg: el.innerHTML, transform: el.style.transform }));
  let release;
  const held = new Promise(resolve => { release = resolve; });
  await page.route('**/render', async route => {
    if (route.request().postDataJSON().cabiln === 'A-G-A') await held;
    await route.continue();
  });
  try {
    await page.locator('#cabiln-input').fill('A-G-A');
    await expect(page.locator('#render-progress')).toBeVisible();
    await expect(page.locator('#render-pane')).toHaveAttribute('aria-busy', 'true');
    expect(await page.locator('#render-inner').evaluate(el => el.innerHTML)).toBe(before.svg);
    await expect(page.locator('#btn-mol')).toBeDisabled();
    await expect(page.locator('.res-chip[data-residue="1"]')).toBeDisabled();
    await page.locator('.res-chip[data-residue="1"]').dispatchEvent('click');
    const point = await residuePoint(page, data, 1);
    await page.mouse.click(point.x, point.y);
    await expect(page.locator('#build-left-abbr')).toHaveText('—');
    await capture(page, testInfo, 'previous-drawing-updating');
  } finally { release(); }
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  expect(await page.locator('#render-inner').evaluate(el => el.style.transform)).toBe(before.transform);
  await expect(page.locator('#render-progress')).toBeHidden();
  const valid = await page.locator('#render-inner').innerHTML();
  await page.locator('#cabiln-input').fill('NotARegisteredMonomer-?');
  await expect(page.locator('#cabiln-input')).toHaveClass('err');
  expect(await page.locator('#render-inner').innerHTML()).toBe(valid);
  await expect(page.locator('#render-progress-label')).toContainText('Previous drawing');
  await expect(page.locator('#btn-png')).toBeDisabled();
  await render(page, 'A-G-A-G');
  await expect(page.locator('#btn-png')).toBeEnabled();
  await expect(page.locator('#render-progress')).toBeHidden();
});

test('visible Use action and keyboard sites build a peptide with undo and redo', async ({ page }, testInfo) => {
  await render(page, 'A-G');
  await page.locator('#btn-lib').click();
  await page.locator('#lib-search').fill('Alanine');
  const use = page.locator('.lib-row[data-abbr="A"] .lib-use');
  await use.focus();
  await page.keyboard.press('Enter');
  await expect(page.locator('#build-panel')).toBeVisible();
  await expect(page.locator('#build-right-abbr')).toHaveText('A');
  await page.locator('.res-chip[data-residue="1"]').focus();
  await page.keyboard.press('Enter');
  await expect(page.locator('#build-left-abbr')).toHaveText('G');
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await expect(page.locator('#build-left-svg .site-selected').first()).toBeVisible();
  await expect(page.locator('#build-right-svg .site-selected').first()).toBeVisible();
  await capture(page, testInfo, 'keyboard-connection-ready');
  await connect(page);
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G-A');
  await page.locator('#btn-undo').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#build-connect')).toBeDisabled();
  await page.locator('#btn-redo').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G-A');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#cabiln-input').focus();
  await page.keyboard.press('ControlOrMeta+z');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await page.keyboard.press('ControlOrMeta+Shift+z');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G-A');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await capture(page, testInfo, 'keyboard-build-history');
});

test('format drafts and conversion history survive a reload through explicit recovery', async ({ page }) => {
  await render(page, 'K-G');
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('NCC(=O)O');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#notation-select').selectOption('cabiln');
  await expect(page.locator('#cabiln-input')).toHaveValue('K-G');
  await page.locator('#notation-select').selectOption('smiles');
  await expect(page.locator('#cabiln-input')).toHaveValue('NCC(=O)O');
  await page.locator('#btn-to-cabiln-pct').click();
  await expect(page.locator('#notation-select')).toHaveValue('cabiln');
  await expect(page.locator('#cabiln-input')).toHaveValue('G');
  await page.locator('#btn-undo').click();
  await expect(page.locator('#notation-select')).toHaveValue('smiles');
  await expect(page.locator('#cabiln-input')).toHaveValue('NCC(=O)O');
  await page.locator('#btn-verify').click();
  await page.locator('#smiles-input').fill('CC(=O)O');
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  await page.reload();
  await expect(page.locator('#btn-restore-draft')).toBeVisible();
  await page.locator('#btn-restore-draft').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('NCC(=O)O');
  await expect(page.locator('#notation-select')).toHaveValue('smiles');
  await expect(page.locator('#smiles-input')).toHaveValue('CC(=O)O');
  await page.locator('#notation-select').selectOption('cabiln');
  await expect(page.locator('#cabiln-input')).toHaveValue('G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
});

test('conversion feedback disappears when an edit supersedes a delayed result', async ({ page }) => {
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('NCC(=O)O');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  let release;
  const held = new Promise(resolve => { release = resolve; });
  await page.route('**/to_cabiln', async route => {
    const response = await route.fetch();
    await held;
    await route.fulfill({ response });
  });
  try {
    await page.locator('#btn-to-cabiln-pct').click();
    await expect(page.locator('#conversion-progress')).toBeVisible();
    await page.locator('#cabiln-input').fill('CC(=O)O');
    await expect(page.locator('#conversion-progress')).toBeHidden();
  } finally { release(); }
  await expect(page.locator('#cabiln-input')).toHaveValue('CC(=O)O');
  await expect(page.locator('#notation-select')).toHaveValue('smiles');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
});

test('unavailable draft storage leaves editing and history usable', async ({ page }) => {
  await page.addInitScript(() => {
    Object.defineProperty(window, 'localStorage', { get() { throw new Error('Storage blocked'); } });
  });
  await page.reload();
  await render(page, 'A-G');
  await expect(page.locator('#draft-status')).toContainText('unavailable');
  await page.locator('#btn-undo').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('');
  await page.locator('#btn-redo').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
});

test('starting new work while recovery is offered saves the new work', async ({ page }) => {
  await render(page, 'A');
  await page.reload();
  await expect(page.locator('#btn-restore-draft')).toBeVisible();
  await render(page, 'A-G');
  await expect(page.locator('#btn-restore-draft')).toBeHidden();
  await page.reload();
  await page.locator('#btn-restore-draft').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
});

test('inferred stereo warnings travel with format drafts and recovery', async ({ page }) => {
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('NC(C)C(=O)O');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#btn-to-cabiln-pct').click();
  await expect(page.locator('#conversion-status')).toBeVisible();
  const warning = await page.locator('#conversion-status').textContent();
  expect(warning).toMatch(/stereo/i);
  const sequence = await page.locator('#cabiln-input').inputValue();
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#notation-select').selectOption('cabiln');
  await expect(page.locator('#cabiln-input')).toHaveValue(sequence);
  await expect(page.locator('#conversion-status')).toHaveText(warning);
  await page.reload();
  await page.locator('#btn-restore-draft').click();
  await expect(page.locator('#cabiln-input')).toHaveValue(sequence);
  await expect(page.locator('#conversion-status')).toHaveText(warning);
  await expect(page.locator('#conversion-status')).toBeVisible();
});

test('new reference work resolves recovery and Verify uses the recovered reference', async ({ page }) => {
  await render(page, 'G');
  await page.locator('#btn-verify').click();
  await page.locator('#smiles-input').fill('NCC(=O)O');
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  await page.reload();
  await expect(page.locator('#btn-restore-draft')).toBeVisible();
  await page.locator('#btn-verify').click();
  await page.locator('#smiles-input').fill('CC(=O)O');
  await expect(page.locator('#btn-restore-draft')).toBeHidden();
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  await render(page, 'G');
  await page.reload();
  await page.locator('#btn-restore-draft').click();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  const verified = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/verify');
  await page.locator('#btn-verify').click();
  const response = await verified;
  expect(response.request().postDataJSON().smiles).toBe('CC(=O)O');
  await expect(page.locator('#smiles-input')).toHaveValue('CC(=O)O');
});

test('palette refresh preserves unchanged rows and invalidates changed-definition previews', async ({ page }) => {
  await page.locator('#btn-lib').click();
  const row = page.locator('.lib-row[data-abbr="G"]');
  await expect(row).toBeVisible();
  await row.evaluate(el => { el.dataset.identity = 'retained'; });
  await row.hover();
  await expect(page.locator('#lib-preview svg').first()).toBeVisible();
  const refresh = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/monomers');
  await page.evaluate(() => window.dispatchEvent(new Event('focus')));
  await refresh;
  await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  await expect(row).toHaveAttribute('data-identity', 'retained');
  await page.route('**/monomers', async route => {
    const response = await route.fetch();
    await route.fulfill({ response, headers: { ...response.headers(), 'x-library-version': 'changed-definition' } });
  });
  // Replacing the row under the pointer can start a fresh preview immediately.
  const preview = page.waitForRequest(request => {
    const url = new URL(request.url());
    return url.pathname === '/monomer_svg' && url.searchParams.get('abbr') === 'G';
  });
  await page.evaluate(() => window.dispatchEvent(new Event('focus')));
  await expect(row).not.toHaveAttribute('data-identity');
  await expect(page.locator('#lib-preview')).toBeHidden();
  await page.locator('#lib-search').hover();
  await row.hover();
  await preview;
  await expect(page.locator('#lib-preview svg').first()).toBeVisible();
});

test('narrow screens retain reachable panels, inputs, and builder controls', async ({ page }, testInfo) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await render(page, 'A-G');
  await page.locator('#btn-build').click();
  await page.locator('#lib-search').fill('Alanine');
  await page.locator('.lib-row[data-abbr="A"] .lib-use').click();
  await selectChip(page, 1, 'left', 'G');
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await connect(page);
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G-A');
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  await capture(page, testInfo, 'narrow-build');
  await page.locator('#lib-close').click();
  await page.locator('#build-close').click();
  await page.locator('#btn-verify').click();
  await page.locator('#smiles-input').fill('NCC(=O)O');
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  await capture(page, testInfo, 'narrow-verify');
});
