const { test, expect, render, tile, site, selectChip, capture, residuePoint } = require('./fixtures');

test('panels return focus on dismissal and Escape closes the active panel independently', async ({ page, browserName }) => {
  // Safari on macOS uses Option+Tab to include links in keyboard navigation.
  await page.keyboard.press(browserName === 'webkit' && process.platform === 'darwin' ? 'Alt+Tab' : 'Tab');
  await expect(page.getByRole('link', { name: 'Skip to peptide drawing' })).toBeFocused();
  await page.keyboard.press('Enter');
  await expect(page.locator('#render-canvas')).toBeFocused();
  for (const [button, panel, focused] of [
    ['btn-lib', 'lib-panel', 'lib-search'],
    ['btn-examples', 'examples-panel', 'examples-close'],
    ['btn-verify', 'verify-pane', 'smiles-input'],
    ['btn-help', 'help-panel', 'help-close'],
    ['btn-build', 'build-panel', 'build-action'],
  ]) {
    await page.locator(`#${button}`).click();
    await expect(page.locator(`#${focused}`)).toBeFocused();
    await page.keyboard.press('Escape');
    await expect(page.locator(`#${panel}`)).toBeHidden();
    await expect(page.locator(`#${button}`)).toBeFocused();
  }
  // Build opened Library; dismissing Build did not dismiss its independent panel.
  await expect(page.locator('#lib-panel')).toBeVisible();
  await page.locator('#lib-search').focus();
  await page.keyboard.press('Escape');
  await expect(page.locator('#btn-lib')).toBeFocused();
  await page.locator('#btn-examples').click();
  await expect(page.locator('.example-row').first()).toBeVisible();
  await page.locator('.example-row').first().focus();
  await page.keyboard.press('Enter');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#examples-panel')).toBeVisible();
});

test('all three drawings share visible zoom controls and keyboard pan and reset', async ({ page }, testInfo) => {
  await render(page, 'G');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'G');
  await tile(page, 'A', 'right');
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-inner svg')).toBeVisible();
  await page.locator('#btn-verify').click();
  await page.locator('#smiles-input').fill('NCC(=O)O');
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  for (const [canvas, inner] of [
    ['render-canvas', 'render-inner'], ['smiles-canvas', 'smiles-inner'],
    ['build-preview-canvas', 'build-preview-inner'],
  ]) {
    const controls = page.locator(`[data-viewport="${canvas}"]`);
    const view = page.locator(`#${inner}`);
    await controls.getByRole('button', { name: 'Zoom in', exact: true }).click();
    await expect(view).toHaveCSS('transform', 'matrix(1.2, 0, 0, 1.2, 0, 0)');
    await controls.getByRole('button', { name: 'Zoom out', exact: true }).click();
    await page.locator(`#${canvas}`).focus();
    await page.keyboard.press('ArrowRight');
    await page.keyboard.press('ArrowDown');
    await expect(view).toHaveCSS('transform', 'matrix(1, 0, 0, 1, 30, 30)');
    await page.keyboard.press('+');
    await expect(view).toHaveCSS('transform', 'matrix(1.2, 0, 0, 1.2, 30, 30)');
    await page.keyboard.press('0');
    await expect(view).toHaveCSS('transform', 'matrix(1, 0, 0, 1, 0, 0)');
    await controls.getByRole('button', { name: 'Reset view' }).click();
  }
  await page.locator('#build-preview-close').click();
  await expect(page.locator('#build-preview-button')).toBeFocused();
  await expect(page.locator('#build-connect')).toBeEnabled();
  await page.locator('#verify-close').click();
  await expect(page.locator('#btn-verify')).toBeFocused();
  await capture(page, testInfo, 'shared-drawing-controls');
});

test('touch pan and pinch preserve selection and work at a narrow viewport', async ({ browser, browserName, app }, testInfo) => {
  const context = await browser.newContext({ hasTouch: true, viewport: { width: 390, height: 844 } });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  try {
    await page.goto(app.url);
    const data = await render(page, 'G-A');
    const canvas = page.locator('#render-canvas');
    await canvas.scrollIntoViewIfNeeded();
    const box = await canvas.boundingBox();
    const x = Math.round(box.x + box.width / 2), y = Math.round(box.y + box.height / 2);
    const view = page.locator('#render-inner');
    if (browserName === 'chromium') {
      const session = await context.newCDPSession(page);
      const touch = (type, points) => session.send('Input.dispatchTouchEvent', {
        type, touchPoints: points.map(([id, x, y]) => ({ id, x, y })),
      });
      await touch('touchStart', [[0, x, y]]);
      await touch('touchMove', [[0, x + 35, y + 25]]);
      await touch('touchEnd', []);
      await expect(view).toHaveCSS('transform', 'matrix(1, 0, 0, 1, 35, 25)');
      await touch('touchStart', [[0, x - 30, y], [1, x + 30, y]]);
      await touch('touchMove', [[0, x - 60, y], [1, x + 60, y]]);
      await touch('touchEnd', []);
      await expect.poll(() => view.evaluate(el => new DOMMatrix(getComputedStyle(el).transform).a)).toBeCloseTo(2);
      await page.locator('[data-viewport="render-canvas"]').getByRole('button', { name: 'Reset view' }).tap();
      await expect(view).toHaveCSS('transform', 'matrix(1, 0, 0, 1, 0, 0)');
      await session.detach();
    }
    await page.locator('#btn-build').tap();
    await canvas.scrollIntoViewIfNeeded();
    const point = await residuePoint(page, data, 0);
    await page.touchscreen.tap(point.x, point.y);
    await expect(page.locator('#build-left-abbr')).toHaveText('G');
    await page.locator('#build-left-change').tap();
    await page.locator('#residue-chips [data-residue="0"]').tap();
    await expect(page.locator('#build-left-abbr')).toHaveText('G');
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
    await capture(page, testInfo, 'touch-builder');
    expect(errors).toEqual([]);
  } finally { await context.close(); }
});

test('failed library, example and drawing requests have working Retry actions in place', async ({ page }) => {
  for (const [url, button, container, success] of [
    ['monomers', 'btn-lib', 'lib-list', '.lib-row'],
    ['examples', 'btn-examples', 'examples-list', '.example-row'],
  ]) {
    await page.route(`**/${url}`, route => route.fulfill({ status: 500, contentType: 'application/json', body: '{"error":"Temporary failure"}' }), { times: 1 });
    await page.locator(`#${button}`).click();
    await expect(page.locator(`#${container}`).getByRole('button', { name: 'Retry' })).toBeVisible();
    if (url === 'monomers') {
      await page.locator('#lib-search').fill('Glycine');
      await page.locator('#btn-rxn-filter').click();
      await expect(page.locator(`#${container}`).getByRole('button', { name: 'Retry' })).toBeVisible();
    }
    await page.locator(`#${container}`).getByRole('button', { name: 'Retry' }).click();
    await expect(page.locator(success).first()).toBeVisible();
    if (url === 'monomers') {
      await page.locator('#build-close').click();
      await page.locator('#lib-search').fill('');
    }
    await page.locator(`#${button}`).click();
  }
  await page.route('**/render', route => route.fulfill({ status: 500, contentType: 'application/json', body: '{"error":"Temporary renderer failure"}' }), { times: 1 });
  await page.locator('#cabiln-input').fill('G-A');
  await expect(page.locator('#cabiln-input')).toHaveValue('G-A');
  await page.locator('#cabiln-status').getByRole('button', { name: 'Retry' }).click();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#btn-mol')).toBeEnabled();
  await page.locator('#btn-verify').click();
  await page.route('**/render_reference', route => route.fulfill({ status: 500, contentType: 'application/json', body: '{"error":"Temporary reference failure"}' }), { times: 1 });
  await page.locator('#smiles-input').fill('NCC(=O)O');
  await page.locator('#smiles-status').getByRole('button', { name: 'Retry' }).click();
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  await page.route('**/verify', route => route.abort('failed'), { times: 1 });
  await page.locator('#smiles-input').fill('NCC(=O)N[C@@H](C)C(=O)O');
  await page.locator('#compare-bar').getByRole('button', { name: 'Retry' }).click();
  await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
});

test('opening all panels at laptop width keeps both drawings usable', async ({ page }) => {
  await page.setViewportSize({ width: 1024, height: 768 });
  await render(page, 'A-G');
  for (const id of ['btn-lib', 'btn-examples', 'btn-verify', 'btn-build']) await page.locator('#'+id).click();
  await page.locator('#smiles-input').fill('N[C@@H](C)C(=O)NCC(=O)O');
  await expect(page.locator('#compare-bar .match')).toBeVisible();
  for (const id of ['render-canvas', 'smiles-canvas']) {
    expect((await page.locator('#'+id).boundingBox()).width).toBeGreaterThan(300);
  }
  expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBe(1024);
});

test('losing window focus ends a drag instead of leaving the drawing attached to the pointer', async ({ page }) => {
  await render(page, 'G-A');
  const box = await page.locator('#render-canvas').boundingBox();
  const x = box.x + 30, y = box.y + 30;
  await page.mouse.move(x, y);
  await page.mouse.down();
  await page.mouse.move(x + 20, y + 20);
  const view = page.locator('#render-inner');
  const before = await view.evaluate(el => el.style.transform);
  await page.evaluate(() => window.dispatchEvent(new Event('blur')));
  await page.mouse.move(x + 60, y + 60);
  expect(await view.evaluate(el => el.style.transform)).toBe(before);
  await page.mouse.up();
  await page.locator('[data-viewport="render-canvas"]').getByRole('button', { name: 'Reset view' }).click();
  await expect(view).toHaveCSS('transform', 'matrix(1, 0, 0, 1, 0, 0)');
});

test('properties, visible warnings and expandable normalization details remain distinct', async ({ page }, testInfo) => {
  await page.route('**/render', async route => {
    const response = await route.fetch();
    const data = await response.json();
    if (response.ok()) {
      data.normalization_note = 'Detailed normalization explanation <img src=x onerror=alert(1)>';
      data.warnings = ['Check the inferred stereochemistry before using this structure.'];
    }
    await route.fulfill({ response, json: data });
  });
  await render(page, 'G-A');
  const status = page.locator('#cabiln-status');
  await expect(status.locator('.structure-info')).toContainText('MW');
  await expect(status.locator('.structure-warning')).toBeVisible();
  await expect(status.locator('.structure-details p')).toBeHidden();
  await status.getByText('Notation normalized · Details', { exact: true }).click();
  await expect(status.locator('.structure-details p')).toHaveText('Detailed normalization explanation <img src=x onerror=alert(1)>');
  await expect(status.locator('img')).toHaveCount(0);
  await capture(page, testInfo, 'structured-drawing-status');
});
