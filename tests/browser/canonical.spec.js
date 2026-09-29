const { test, expect, render, isCompletedResponse } = require('./fixtures');

test('canonical formatting preserves chemistry, recursive tabs and undo history', async ({ page }) => {
  const original = 'K.[G(4,2).ac(1,2)]-A';
  const initial = await render(page, original);
  await expect(page.locator('#notation-policy')).toHaveValue('layout');
  await page.locator('#notation-policy').selectOption('canonical');
  for (const [target, expected] of [
    ['bracket', 'K.[G(4,2).[ac(1,2)]]-A'],
    ['branch', 'K.!1(4,2)-A%ac-G.!1(2,4)'],
  ]) {
    const conversion = page.waitForResponse(response =>
      isCompletedResponse(response) && new URL(response.url()).pathname === '/convert_notation');
    await page.locator(`#btn-to-${target}`).click();
    const response = await conversion;
    expect(response.request().postDataJSON().canonical).toBe(true);
    expect(response.status(), await response.text()).toBe(200);
    await expect(page.locator('#cabiln-input')).toHaveValue(expected);
    await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    await expect(page.locator('#residue-chips [data-residue]')).toHaveCount(4);
    if (target === 'bracket') {
      await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['K', 'G', 'ac', 'A']);
      const inner = page.locator('#residue-chips .branch-chip').filter({ hasText: /^\[$/ }).nth(1);
      await expect(inner).toHaveAttribute('data-members', '[3]');
      await inner.hover();
      await expect(page.locator('#render-inner svg')).toHaveClass(/has-highlight/);
      await expect(page.locator('#residue-chips [data-residue="3"]')).toHaveClass(/hover/);
      await expect(page.locator('#residue-chips [data-residue="2"]')).not.toHaveClass(/hover/);
      await page.locator('#btn-undo').click();
      await expect(page.locator('#cabiln-input')).toHaveValue(original);
      await page.locator('#btn-redo').click();
      await expect(page.locator('#cabiln-input')).toHaveValue(expected);
      await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    }
  }
  await page.locator('#btn-verify').click();
  await page.locator('#mol-upload').setInputFiles({
    name: 'original.mol', mimeType: 'chemical/x-mdl-molfile',
    buffer: Buffer.from(initial.mol_block),
  });
  await expect(page.locator('#compare-bar .match')).toHaveText('✓ EXACT MATCH');
});
