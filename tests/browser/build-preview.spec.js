const { test, expect, render, selectChip, site, tile, connect, capture, isCompletedResponse } = require('./fixtures');

for (const example of [
  { name: 'backbone addition', source: 'A-G', host: 1, left: 'G', right: 'A', slots: [2, 1], groups: ['[OH]', '[H]'], reaction: /amide/i },
  { name: 'N-terminal cap', source: 'K-G', host: 0, left: 'K', right: 'ac', slots: [1, 2], groups: ['[H]', '[OH]'], reaction: /amide/i },
  { name: 'existing disulfide', source: 'C-A-C', host: 0, left: 'C', right: 'C', target: 2, slots: [4, 4], groups: ['[H]', '[H]'], reaction: /disulfide/ },
  { name: 'maleimide addition', source: 'C', host: 0, left: 'C', right: 'Mal', slots: [4, 4], groups: ['[H]', '[H]'], reaction: /thiol maleimide/ },
  { name: 'alkyne–azide cycloaddition', source: 'Pra', host: 0, left: 'Pra', right: 'AzK', slots: [4, 4], groups: ['[H]', '[H]'], reaction: /cuaac 1 4 triazole/ },
]) {
  test(`connection preview shows ${example.name} without applying it`, async ({ page }, testInfo) => {
    await render(page, example.source);
    await page.locator('#btn-build').click();
    await selectChip(page, example.host, 'left', example.left);
    if (example.target === undefined) await tile(page, example.right, 'right');
    else await selectChip(page, example.target, 'right', example.right);
    await site(page, 'left', example.slots[0]);
    await site(page, 'right', example.slots[1]);
    for (const [index, side] of ['left', 'right'].entries()) {
      await expect(page.locator(`#build-${side}-site`)).toHaveText(`R${example.slots[index]} · Free-site group: ${example.groups[index]}`);
      // A vertical bond has zero bounding-box width even when its stroke is
      // painted. Check the attachment label and its glow instead.
      const marker = page.locator(`#build-${side}-svg .site-selected:not([class*="bond-"])`).first();
      await expect(marker).toBeVisible();
      await expect(marker).toHaveCSS('filter', /drop-shadow/);
    }
    await expect(page.locator('#build-preview-button')).toBeEnabled();
    const before = await page.locator('#render-inner').innerHTML();
    const history = await page.evaluate(() => JSON.stringify(editor));
    const drawing = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/render');
    await page.locator('#build-preview-button').click();
    const product = await (await drawing).json();
    await expect(page.locator('#build-preview-inner svg')).toBeVisible();
    await expect(page.locator('#build-preview-reaction')).toHaveText(example.reaction);
    const left = `${example.left} (residue ${example.host + 1}) R${example.slots[0]}`;
    const right = `${example.right} (${example.target === undefined ? 'new' : `residue ${example.target + 1}`}) R${example.slots[1]}`;
    await expect(page.locator('#build-preview-status')).toContainText(`${left} ↔ ${right}`);
    await expect(page.locator('#cabiln-input')).toHaveValue(example.source);
    expect(await page.locator('#render-inner').innerHTML()).toBe(before);
    expect(await page.evaluate(() => JSON.stringify(editor))).toBe(history);
    const proposed = await page.locator('#build-preview-source').textContent();
    expect(product.cabiln_echo).toBe(proposed);
    if (example.name === 'backbone addition') {
      const palette = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/monomers');
      await page.evaluate(() => window.dispatchEvent(new Event('focus')));
      await palette;
      await expect(page.locator('#build-preview-inner svg')).toBeVisible();
      for (const viewport of [{ width: 1280, height: 800 }, { width: 1024, height: 768 }, { width: 1440, height: 1000 }]) {
        await page.setViewportSize(viewport);
        for (const id of ['build-left', 'build-right', 'build-preview-button', 'build-connect', 'build-preview-canvas']) {
          await expect(page.locator(`#${id}`)).toBeInViewport({ ratio: 1 });
        }
        expect((await page.locator('#render-pane').boundingBox()).height).toBeGreaterThan(0);
      }
      await capture(page, testInfo, 'connection-preview');
      await page.setViewportSize({ width: 390, height: 844 });
      await page.locator('#build-preview').scrollIntoViewIfNeeded();
      expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);
      await expect(page.locator('#build-preview-inner svg')).toBeVisible();
      await expect(page.locator('#build-preview-canvas')).toBeInViewport({ ratio: 1 });
      await capture(page, testInfo, 'narrow-connection-preview');
    }
    const applied = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/render');
    const { data } = await connect(page);
    expect(data.result).toBe(proposed);
    expect((await (await applied).json()).residues).toEqual(product.residues);
    await expect(page.locator('#build-preview')).toBeHidden();
    await page.locator('#btn-undo').click();
    await expect(page.locator('#cabiln-input')).toHaveValue(example.source);
  });
}

test('a failed preview keeps the peptide and lets the user retry or connect', async ({ page }) => {
  await render(page, 'G');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'G');
  await tile(page, 'A', 'right');
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await page.route('**/render', route => route.fulfill({ status: 500, contentType: 'application/json', body: JSON.stringify({ error: 'Drawing failed' }) }));
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-status')).toContainText('Drawing failed');
  await expect(page.locator('#cabiln-input')).toHaveValue('G');
  await expect(page.locator('#build-connect')).toBeEnabled();
  await expect(page.locator('#build-preview-button')).toBeEnabled();
  await page.unroute('**/render');
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-inner svg')).toBeVisible();
  await page.locator('#build-left-rgroups button.selected').click();
  await expect(page.locator('#build-preview')).toBeHidden();
  await expect(page.locator('#build-preview-button')).toBeDisabled();
});
