const { test, expect, render, selectChip, site, capture, residuePoint } = require('./fixtures');

async function choosePracticeMonomer(page, abbr) {
  await page.goto(`${page.url()}?tutorial=1`);
  await page.locator('#tutorial-load').click();
  await page.locator('#tutorial-next').click();
  await page.locator('#tutorial-show').click();
  await page.locator('#tutorial-next').click();
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', 'G');
  await page.locator('#tutorial-next').click();
  await page.locator(`.lib-row[data-abbr="${abbr}"] .lib-use`).click();
  await expect(page.locator('#build-right-abbr')).toHaveText(abbr);
}

test('Show control preserves the browsed monomer and does not open an offscreen preview', async ({ page }) => {
  await choosePracticeMonomer(page, 'G');
  await page.locator('.lib-row[data-abbr="G"]').hover();
  await expect(page.locator('#lib-preview')).toBeVisible();
  await page.locator('.lib-row[data-abbr="DAla"] .lib-use').click();
  await expect(page.locator('#build-right-abbr')).toHaveText('DAla');
  const scroll = await page.locator('#lib-list').evaluate(el => el.scrollTop);
  expect(scroll).toBeGreaterThan(0);
  await page.locator('#tutorial-show').click();
  // The library preview opens after a 200 ms hover/focus debounce.
  await page.waitForTimeout(300);
  await expect(page.locator('#lib-preview')).toBeHidden();
  await expect(page.locator('#lib-list')).toBeFocused();
  expect(await page.locator('#lib-list').evaluate(el => el.scrollTop)).toBe(scroll);
  await expect(page.locator('#build-right-abbr')).toHaveText('DAla');
});

for (const preview of [false, true]) test(`completed connection supports Back and Undo review (${preview ? 'preview, guide closed during apply' : 'direct Connect'})`, async ({ page }) => {
  await choosePracticeMonomer(page, 'DAla');
  const next = page.locator('#tutorial-next');
  await next.click();
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  if (preview) {
    await page.locator('#build-preview-button').click();
    await expect(page.locator('#build-preview-source')).toHaveText('A-G-DAla');
    let release, arrived;
    const held = new Promise(resolve => { arrived = resolve; });
    const gate = new Promise(resolve => { release = resolve; });
    await page.route('**/insert_bond', async route => {
      const response = await route.fetch();
      arrived();
      await gate;
      await route.fulfill({ response });
    }, { times: 1 });
    await page.locator('#build-connect').click();
    await held;
    await page.locator('#tutorial-close').click();
    release();
    await expect(page.locator('#cabiln-input')).toHaveValue('A-G-DAla');
    await page.locator('#btn-tutorial').click();
  } else await page.locator('#build-connect').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G-DAla');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(next).toBeEnabled();
  for (const undone of [false, true]) {
    while (await page.locator('#tutorial-back').isEnabled()) {
      await page.locator('#tutorial-back').click();
      await expect(next).toBeEnabled();
    }
    for (let step = 1; step < 8; step++) {
      await page.locator('#tutorial-show').click();
      await next.click();
      await expect(page.locator('#tutorial-progress')).toHaveText(`Step ${step + 1} of 9`);
    }
    await expect(page.locator('#cabiln-input')).toHaveValue(undone ? 'A-G' : 'A-G-DAla');
    if (!undone) {
      await expect(next).toBeDisabled();
      await page.locator('#btn-undo').click();
      await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    }
    await expect(next).toBeEnabled();
  }
  await page.locator('#btn-redo').click();
  await expect(next).toBeDisabled();
  await page.locator('#btn-undo').click();
  await next.click();
  await next.click();
  await expect(page.locator('#tutorial-panel')).toBeHidden();
});

for (const width of [1440, 390]) test(`Retatrutide Swap lesson preserves its branch, mapping review, Back and saved work at ${width}px`, async ({ page, context }, testInfo) => {
  await page.setViewportSize({ width, height: 900 });
  await render(page, 'ac-K-am');
  await expect.poll(() => page.evaluate(() => JSON.parse(localStorage.getItem('cabiln.draft.v1'))?.document.text)).toBe('ac-K-am');
  const saved = await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'));
  await page.locator('#btn-build').click();
  await page.locator('#build-action').selectOption('swap');
  const opened = context.waitForEvent('page');
  await page.locator('#build-swap-tutorial').click();
  const practice = await opened;
  await practice.setViewportSize({ width, height: 900 });
  await expect(practice.locator('#tutorial-lesson')).toHaveValue('swap');
  await expect(practice.locator('#draft-status')).toContainText('Practice tab');
  const rendered = practice.waitForResponse(response => response.url().endsWith('/render'));
  await practice.getByRole('button', { name: 'Load Retatrutide', exact: true }).click();
  const original = await (await rendered).json();
  const source = await practice.locator('#cabiln-input').inputValue();
  const next = practice.locator('#tutorial-next');
  await next.click();
  await practice.locator('#tutorial-show').click();
  await practice.locator('#build-action').selectOption('swap');
  await practice.locator('#tutorial-show').click();
  await expect(practice.locator('.res-chip[data-residue="16"]')).toBeFocused();
  await selectChip(practice, 16, 'left', 'K');
  await next.click();
  await practice.locator('#lib-close').click();
  await practice.locator('#tutorial-show').click();
  await practice.locator('#lib-search').fill('Orn');
  await practice.getByRole('button', { name: 'Select Orn as replacement', exact: true }).click();
  await next.click();
  await expect(practice.locator('#build-swap-sites select')).toHaveCount(3);
  await expect(practice.locator('#build-swap-sites')).toContainText('AEEA');
  await practice.locator('#tutorial-show').click();
  await expect(practice.locator('#build-swap-mapping')).toBeFocused();
  if (width === 1440) {
    const mapping = await practice.locator('#build-swap-sites').boundingBox();
    const panel = await practice.locator('#build-panel').boundingBox();
    expect(mapping.y + mapping.height).toBeLessThanOrEqual(panel.y + panel.height);
  }
  await practice.getByLabel('Replacement site for R4', { exact: true }).selectOption('1');
  await expect(next).toBeDisabled();
  await expect(practice.locator('#build-preview-button')).toBeDisabled();
  await practice.getByLabel('Replacement site for R4', { exact: true }).selectOption('4');
  await next.click();
  await expect(practice.locator('#build-connect')).toBeDisabled();
  await practice.getByRole('button', { name: 'Preview swap', exact: true }).click();
  const expected = source.replace('K.[AEEA', 'Orn.[AEEA');
  expect(expected).not.toBe(source);
  await expect(practice.locator('#build-preview-source')).toHaveText(expected);
  await expect(practice.locator('#build-preview-reaction')).toHaveText('R1 → R1 · R2 → R2 · R4 → R4');
  await capture(practice, testInfo, 'retatrutide-swap-preview');
  await next.click();
  const product = practice.waitForResponse(response => response.url().endsWith('/render'));
  await practice.locator('#build-connect').click();
  const replaced = await (await product).json();
  expect(replaced.layout.groups).toEqual(original.layout.groups);
  expect(replaced.residues.find(residue => residue.idx === 16).abbr).toBe('Orn');
  await expect(practice.locator('#cabiln-input')).toHaveValue(expected);
  await expect(next).toBeEnabled();
  while (await practice.locator('#tutorial-back').isEnabled()) {
    await practice.locator('#tutorial-back').click();
    await expect(next).toBeEnabled();
  }
  for (let step = 1; step < 7; step++) await next.click();
  await expect(practice.locator('#tutorial-title')).toHaveText('Restore Retatrutide');
  await practice.locator('#btn-undo').click();
  await next.click();
  await expect(practice.locator('#cabiln-input')).toHaveValue(source);
  await expect(practice.locator('#cabiln-input')).toHaveClass('ok');
  expect(await practice.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await practice.locator('#tutorial-lesson').selectOption('connect');
  await expect(practice.locator('#tutorial-load')).toHaveText('Load A–G');
  await expect(practice.locator('#cabiln-input')).toHaveValue(source);
  await practice.close();
  await expect(page.locator('#cabiln-input')).toHaveValue('ac-K-am');
  expect(await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'))).toBe(saved);
});

test('example loading can fail, retry and finish after switching lessons without replacing the document', async ({ page }) => {
  await page.goto(`${page.url()}?tutorial=swap`);
  await page.route('**/examples', route => route.fulfill({ status: 500,
    contentType: 'application/json', body: '{"error":"Temporary failure"}' }), { times: 1 });
  await page.locator('#tutorial-load').click();
  await expect(page.locator('#tutorial-status')).toContainText('Could not load Retatrutide');
  await expect(page.locator('#tutorial-next')).toBeDisabled();
  let release, arrived;
  const held = new Promise(resolve => { arrived = resolve; });
  const gate = new Promise(resolve => { release = resolve; });
  await page.route('**/examples', async route => {
    const response = await route.fetch();
    arrived();
    await gate;
    await route.fulfill({ response });
  }, { times: 1 });
  await page.locator('#tutorial-load').click();
  await held;
  await page.locator('#tutorial-lesson').selectOption('connect');
  await page.locator('#tutorial-load').click();
  const completed = page.waitForResponse(response => response.url().endsWith('/examples'));
  release();
  await completed;
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#tutorial-next')).toBeEnabled();
  await page.locator('#tutorial-lesson').selectOption('swap');
  await page.locator('#tutorial-load').click();
  await expect(page.locator('#tutorial-next')).toBeEnabled();
  await expect(page.locator('#cabiln-input')).toHaveValue(/^Y-Aib-Q-G-T/);
});

for (const width of [1440, 390]) {
  test(`guided practice builds, previews and undoes without changing saved work at ${width}px`, async ({ page, context }, testInfo) => {
    await page.setViewportSize({ width, height: 900 });
    await expect(page.locator('#tutorial-panel')).toBeHidden();
    await expect(page.locator('#btn-tutorial')).toBeHidden();
    await render(page, 'ac-K-am');
    await page.locator('#btn-verify').click();
    await page.locator('#smiles-input').fill('NCC(=O)O');
    await expect(page.locator('#smiles-input')).toHaveClass('ok');
    await expect.poll(() => page.evaluate(() => JSON.parse(localStorage.getItem('cabiln.draft.v1'))?.reference)).toBe('NCC(=O)O');
    const saved = await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'));
    await expect(page.locator('#tutorial-launch')).toBeVisible();
    const opened = context.waitForEvent('page');
    await page.locator('#tutorial-launch').click();
    const practice = await opened;
    await practice.setViewportSize({ width, height: 900 });
    const next = practice.locator('#tutorial-next');
    await expect(next).toBeDisabled();
    await expect(practice.locator('#btn-restore-draft')).toBeHidden();
    await expect(practice.locator('#draft-status')).toContainText('Practice tab');
    await practice.locator('#tutorial-load').click();
    await next.click();
    await expect(practice.locator('#tutorial-title')).toHaveText('Find a residue in the drawing');
    await practice.locator('#tutorial-show').click();
    await expect(practice.locator('#render-inner .res-hl').first()).toBeVisible();
    await next.click();
    await expect(next).toBeDisabled();
    await practice.locator('#btn-build').click();
    await selectChip(practice, 1, 'left', 'G');
    await next.click();
    if (width === 1440) await practice.locator('#lib-search').fill('Alanine');
    await practice.locator('.lib-row[data-abbr="A"] .lib-use').click();
    await next.click();
    await expect(next).toBeDisabled();
    await site(practice, 'left', 2);
    await site(practice, 'right', 1);
    await next.click();
    await practice.locator('#build-preview-button').click();
    await expect(practice.locator('#build-preview-source')).toHaveText('A-G-A');
    await capture(practice, testInfo, 'tutorial-preview');
    await next.click();
    await practice.locator('#build-connect').click();
    await next.click();
    await practice.locator('#btn-undo').click();
    await next.click();
    await expect(practice.locator('#cabiln-input')).toHaveValue('A-G');
    await expect(practice.locator('#cabiln-input')).toHaveClass('ok');
    await expect(next).toHaveText('Finish');
    await next.click();
    await expect(practice.locator('#tutorial-panel')).toBeHidden();
    await expect(practice.locator('#btn-tutorial')).toBeFocused();
    expect(await practice.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
    await practice.locator('#btn-help').click();
    await expect(practice.locator('#btn-clear-draft')).toBeHidden();
    await practice.reload();
    await expect(practice.locator('#btn-restore-draft')).toBeHidden();
    await practice.close();
    await expect(page.locator('#cabiln-input')).toHaveValue('ac-K-am');
    await expect(page.locator('#smiles-input')).toHaveValue('NCC(=O)O');
    expect(await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'))).toBe(saved);
  });
}

test('empty-canvas entry, failed drawing, dismissal and restart work without browser storage', async ({ page, context }) => {
  await context.addInitScript(() => {
    Object.defineProperty(window, 'localStorage', { get() { throw new Error('Storage blocked'); } });
  });
  const opened = context.waitForEvent('page');
  await page.locator('#render-inner .tutorial-start').click();
  const practice = await opened;
  await practice.route('**/render', route => route.fulfill({ status: 500,
    contentType: 'application/json', body: '{"error":"Temporary failure"}' }), { times: 1 });
  await practice.locator('#tutorial-load').click();
  const retry = practice.locator('#cabiln-status').getByRole('button', { name: 'Retry' });
  await expect(retry).toBeVisible();
  await expect(practice.locator('#tutorial-next')).toBeDisabled();
  await retry.click();
  await expect(practice.locator('#tutorial-next')).toBeEnabled();
  await practice.locator('#tutorial-title').focus();
  await practice.keyboard.press('Escape');
  await expect(practice.locator('#tutorial-panel')).toBeHidden();
  await expect(practice.locator('#btn-tutorial')).toBeFocused();
  await practice.keyboard.press('Enter');
  await practice.locator('#tutorial-next').click();
  await render(practice, 'G');
  await practice.locator('#tutorial-restart').click();
  await expect(practice.locator('#tutorial-progress')).toHaveText('Step 1 of 9');
  await expect(practice.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(practice.locator('#tutorial-next')).toBeEnabled();
  await expect(practice.locator('#draft-status')).toContainText('not saved automatically');
  await practice.close();
});

for (const width of [1440, 390]) test(`tutorial recovers closed panels and follows a browsed monomer and changed selection at ${width}px`, async ({ page }) => {
  await page.setViewportSize({ width, height: 900 });
  await page.goto(`${page.url()}?tutorial=1`);
  const next = page.locator('#tutorial-next');
  const rendered = page.waitForResponse(response => response.url().endsWith('/render') && response.request().postDataJSON()?.cabiln === 'A-G');
  await page.locator('#tutorial-load').click();
  const data = await (await rendered).json();
  await next.click();
  await page.locator('#tutorial-show').click();
  await next.click();
  await page.locator('#btn-build').click();
  if (width === 1440) {
    await page.locator('.lib-row').first().waitFor();
    const point = await residuePoint(page, data, 1);
    await page.mouse.click(point.x, point.y);
    await expect(page.locator('#build-left-abbr')).toHaveText('G');
  } else await selectChip(page, 1, 'left', 'G');
  await next.click();
  await page.locator('#lib-close').click();
  await page.locator('#tutorial-show').click();
  await expect(page.locator('#lib-panel')).toBeVisible();
  await page.locator('.lib-row[data-abbr="L"] .lib-use').click();
  await expect(page.locator('#lib-search')).toHaveValue('');
  await next.click();
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await expect(next).toBeEnabled();
  await page.locator('#build-close').click();
  await expect(next).toBeDisabled();
  await page.locator('#tutorial-show').click();
  await expect(page.locator('#build-panel')).toBeVisible();
  await expect(page.locator('#tutorial-status')).toContainText('glycine');
  await selectChip(page, 1, 'left', 'G');
  await page.locator('.lib-row[data-abbr="L"] .lib-use').click();
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  if (width === 390) await next.click();
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-source')).toHaveText('A-G-L');
  await expect(next).toBeEnabled();
  await page.locator('.lib-row[data-abbr="V"] .lib-use').click();
  await expect(next).toBeDisabled();
  await site(page, 'right', 1);
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-source')).toHaveText('A-G-V');
  await page.locator('#build-preview-close').click();
  if (width === 390) await expect(next).toBeDisabled();
  await page.locator('#tutorial-show').click();
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-source')).toHaveText('A-G-V');
  if (width === 390) await next.click();
  await page.locator('#build-connect').click();
  if (width === 1440) {
    // Applying before advancing the guide already satisfies the connection and preview steps.
    await next.click();
    await next.click();
  }
  await next.click();
  await page.locator('#btn-undo').click();
  await next.click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await next.click();
  await expect(page.locator('#tutorial-panel')).toBeHidden();
});
