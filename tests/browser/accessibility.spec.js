const { test, expect, render, selectChip, site } = require('./fixtures');

test('workspace, editing, tutorial and registration states pass automated accessibility checks', async ({ page, app }, testInfo) => {
  test.setTimeout(120_000);
  await page.emulateMedia({ reducedMotion: 'reduce' });
  const scans = [];
  async function scan(name) {
    if (!await page.evaluate(() => !!window.axe)) await page.addScriptTag({ path: require.resolve('axe-core/axe.min.js') });
    const result = await page.evaluate(() => axe.run(document, {
      runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21aa', 'wcag22aa', 'best-practice'] },
    }));
    scans.push({ name, violations: result.violations, incomplete: result.incomplete, passes: result.passes.length });
    expect.soft(result.violations.map(v => ({ id: v.id, nodes: v.nodes.map(n => ({ target: n.target, reason: n.failureSummary })) })), name).toEqual([]);
  }
  await scan('empty');
  await render(page, 'ac-K.[G(4,2).ac(1,2)]-am');
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', 'K');
  await expect(page.locator('.lib-row').first()).toBeVisible();
  await scan('library and selected residue');
  await page.locator('#build-action').selectOption('swap');
  await expect(page.locator('#lib-filter-status')).toContainText('Swap');
  await scan('swap');
  await page.locator('#btn-help').click();
  await page.getByText('Start a peptide and add a monomer', { exact: true }).click();
  await scan('help');
  await page.locator('#help-close').click();
  await page.locator('#build-close').click();
  await page.locator('#lib-close').click();
  await page.locator('#btn-examples').click();
  await expect(page.locator('.example-row').first()).toBeVisible();
  await scan('examples');
  await page.locator('#examples-close').click();
  await render(page, 'G');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'G');
  await page.locator('#build-action').selectOption('connect');
  await page.locator('#lib-search').fill('Alanine');
  await page.locator('.lib-row[data-abbr="A"] .lib-use').click();
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-inner svg')).toBeVisible();
  await scan('connection preview');
  await page.locator('#btn-verify').click();
  await page.locator('#smiles-input').fill('NCC(=O)O');
  await expect(page.locator('#compare-bar .match')).toBeVisible();
  await scan('verification');
  await page.locator('#btn-dark').click();
  await scan('light drawings');
  await page.goto(`${app.url}/?tutorial=swap`);
  await expect(page.locator('#tutorial-panel')).toBeVisible();
  await scan('swap tutorial');
  await page.setViewportSize({ width: 320, height: 640 });
  await scan('narrow tutorial');
  await page.goto(`${app.url}/register`);
  await page.locator('#btn-preview').click();
  await scan('registration error');
  await page.locator('#smiles-in').fill('NCC(=O)O');
  await page.locator('#btn-preview').click();
  await expect(page.locator('#preview-canvas svg')).toBeVisible();
  await scan('registration preview');
  await page.locator('#smiles-in').fill('CNCC(CN)C(=O)O');
  await page.locator('#btn-preview').click();
  await expect(page.getByRole('radio')).toHaveCount(2);
  await scan('registration attachment choice');
  await testInfo.attach('accessibility-results', { body: JSON.stringify(scans, null, 2), contentType: 'application/json' });
});
