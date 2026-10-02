const { test, expect, render, selectChip, tile, isCompletedResponse } = require('./fixtures');

async function selectSwap(page, index, abbr) {
  await page.locator('#btn-build').click();
  await page.locator('#build-action').selectOption('swap');
  const options = page.waitForResponse(r => isCompletedResponse(r) && new URL(r.url()).pathname === '/replacement_options');
  await selectChip(page, index, 'left', abbr);
  const response = await options;
  expect(response.status(), await response.text()).toBe(200);
  return response.json();
}

async function previewSwap(page) {
  await expect(page.locator('#build-preview-button')).toBeEnabled();
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-inner svg')).toBeVisible();
  await expect(page.locator('#build-connect')).toBeEnabled();
  return page.locator('#build-preview-source').textContent();
}

for (const sample of [
  { source: 'ac-K.{G(4,2)[.A(1,2)]}-am', selected: 'K', to: 'Orn', expected: 'ac-Orn.{G(4,2)[.A(1,2)]}-am', excluded: 'A' },
  { source: 'ac-C.!r(4,4)-A-C.!r-am', selected: 'C', to: 'dC', expected: 'ac-dC.!r(4,4)-A-C.!r-am', excluded: 'S' },
]) {
  test(`Swap preserves ${sample.selected === 'K' ? 'nested branches' : 'a disulfide'}, preview, exports and history`, async ({ page }) => {
    await render(page, sample.source);
    const options = await selectSwap(page, 1, sample.selected);
    expect(options.candidates.map(c => c.abbr)).not.toContain(sample.excluded);
    await tile(page, sample.to);
    await expect(page.locator('#build-right-svg svg')).toBeVisible();
    await expect(page.locator('#build-swap-sites select')).toHaveCount(3);
    await expect(page.locator('#build-connect')).toBeDisabled();
    expect(await previewSwap(page)).toBe(sample.expected);
    await expect(page.locator('#build-preview-changes')).toContainText(`Replace ${sample.selected} with ${sample.to}`);
    const bonds = await page.locator('#build-preview-inner .edit-connection').evaluateAll(elements =>
      [...new Set(elements.flatMap(el => [...el.classList].filter(name => /^bond-\d+$/.test(name))))]);
    expect(bonds).toHaveLength(3);
    const marked = page.locator('#build-preview-inner .edit-residue').first();
    const before = await marked.boundingBox();
    await page.getByRole('button', { name: 'Focus change', exact: true }).click();
    const after = await marked.boundingBox();
    expect(Math.hypot(after.width, after.height)).toBeGreaterThan(Math.hypot(before.width, before.height));
    await expect(marked).toBeInViewport({ ratio: 1 });
    const focused = await page.locator('#build-preview-inner').evaluate(el => el.style.transform);
    await page.getByRole('button', { name: 'Focus change', exact: true }).click();
    expect(await page.locator('#build-preview-inner').evaluate(el => el.style.transform)).toBe(focused);
    await page.locator('#build-preview [data-zoom="reset"]').click();
    expect(await marked.boundingBox()).toEqual(before);
    await expect(page.locator('#cabiln-input')).toHaveValue(sample.source);
    if (sample.selected === 'K') {
      await page.locator('#lib-close').click();
      await page.getByRole('button', { name: 'Change replacement monomer', exact: true }).click();
      await expect(page.locator('#lib-panel')).toBeVisible();
      await expect(page.locator('#build-left-abbr')).toHaveText('K');
      await expect(page.locator('#residue-chips [data-residue="1"]')).toHaveCSS('outline-color', 'rgb(90, 154, 224)');
      await expect(page.locator('#build-preview')).toBeHidden();
      await expect(page.locator('#build-connect')).toBeDisabled();
      await expect(page.locator('#build-right-abbr')).toHaveText('—');
      await tile(page, sample.to);
      await previewSwap(page);
      // Two individually compatible choices cannot consume the same site.
      await page.getByLabel('Replacement site for R4', { exact: true }).selectOption('1');
      await expect(page.locator('#build-preview')).toBeHidden();
      await expect(page.locator('#build-preview-button')).toBeDisabled();
      await expect(page.locator('#build-connect')).toBeDisabled();
      await page.getByLabel('Replacement site for R4', { exact: true }).selectOption('4');
      await previewSwap(page);
    }
    await page.locator('#build-connect').click();
    await expect(page.locator('#cabiln-input')).toHaveValue(sample.expected);
    await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    await expect(page.locator('#btn-mol')).toBeEnabled();
    await expect(page.locator('#btn-png')).toBeEnabled();
    await page.locator('#residue-chips [data-residue="1"]').hover();
    await expect(page.locator('#render-inner .res-hl').first()).toHaveCount(1);
    await page.locator('#btn-undo').click();
    await expect(page.locator('#cabiln-input')).toHaveValue(sample.source);
    await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    await page.locator('#btn-redo').click();
    await expect(page.locator('#cabiln-input')).toHaveValue(sample.expected);
  });
}

test('new monomer with a different site number is discovered, mapped and previewed before apply', async ({ page, app }) => {
  const source = 'ac-E.[G(4,1)]-am';
  await render(page, source);
  await selectSwap(page, 1, 'E');
  const registration = await page.request.post(`${app.url}/register_monomer`, { data: {
    chuckles: '[1*]N[C@@H](CC(=O)[7*])C(=O)[2*]',
    abbr: 'SwapD', name: 'Aspartate with sidechain R7', type: 'aa', subtype: 'modified',
    chem_types: {1:'backbone_n',2:'backbone_c',7:'carboxyl'},
    leaving: {1:'[H]',2:'[OH]',7:'[OH]'},
  } });
  expect(registration.status(), await registration.text()).toBe(200);
  const options = page.waitForResponse(r => isCompletedResponse(r) && new URL(r.url()).pathname === '/replacement_options');
  await selectChip(page, 1, 'left', 'E');
  expect((await (await options).json()).candidates.map(c=>c.abbr)).toContain('SwapD');
  await tile(page, 'SwapD');
  await expect(page.getByLabel('Replacement site for R4', { exact: true })).toHaveValue('7');
  const proposed = await previewSwap(page);
  expect(proposed).toContain('SwapD');
  await expect(page.locator('#build-preview-reaction')).toContainText('R4 → R7');
  await expect(page.locator('#cabiln-input')).toHaveValue(source);
  for (const width of [1280, 390]) {
    await page.setViewportSize({ width, height: 900 });
    expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(width);
    await page.locator('#build-connect').scrollIntoViewIfNeeded();
    await expect(page.locator('#build-connect')).toBeInViewport();
  }
  await page.locator('#build-connect').click();
  await expect(page.locator('#cabiln-input')).toHaveValue(proposed);
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  const proof = await page.request.post(`${app.url}/verify`, { data: {
    cabiln: proposed, smiles: 'CC(=O)N[C@@H](CC(=O)NCC(=O)O)C(N)=O',
  } });
  expect(proof.status(), await proof.text()).toBe(200);
  expect((await proof.json()).match).toBe(true);
});
