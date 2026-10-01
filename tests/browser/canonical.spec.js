const fs = require('node:fs');
const { test, expect, render, isCompletedResponse } = require('./fixtures');

test('canonical formatting reuses the drawing with correct tabs, exports and undo history', async ({ page }) => {
  const original = 'K.[G(4,2).ac(1,2)]-A';
  const initial = await render(page, original);
  let current = initial;
  await expect(page.locator('#notation-policy')).toHaveValue('layout');
  await page.locator('#notation-policy').selectOption('canonical');
  for (const [target, expected] of [
    ['bracket', 'K.[G(4,2).[ac(1,2)]]-A'],
    ['branch', 'K.!1(4,2)-A%ac-G.!1(2,4)'],
  ]) {
    await page.locator('#cabiln-input').hover();
    const geometry = await page.locator('#render-inner svg').innerHTML();
    const requests = [];
    const listen = request => {
      const path = new URL(request.url()).pathname;
      if (['/render', '/convert_notation'].includes(path)) requests.push(path);
    };
    page.on('request', listen);
    const conversion = page.waitForResponse(response =>
      isCompletedResponse(response) && new URL(response.url()).pathname === '/convert_notation');
    await page.locator(`#btn-to-${target}`).click();
    const response = await conversion;
    expect(response.request().postDataJSON().canonical).toBe(true);
    expect(response.status(), await response.text()).toBe(200);
    await expect(page.locator('#cabiln-input')).toHaveValue(expected);
    await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    await expect(page.locator('#render-pane')).toHaveAttribute('aria-busy', 'false');
    await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
    page.off('request', listen);
    expect(requests).toEqual(['/convert_notation']);
    expect(await page.locator('#render-inner svg').innerHTML()).toBe(geometry);
    await expect(page.locator('#residue-chips [data-residue]')).toHaveCount(4);
    const downloaded = page.waitForEvent('download');
    await page.locator('#btn-mol').click();
    expect(fs.readFileSync(await (await downloaded).path(), 'utf8')).toBe(current.mol_block);
    // The cap's known source atoms remain its atoms after canonical reordering.
    const cap = current.residues.find(residue => residue.abbr === 'ac');
    await page.locator('#residue-chips [data-residue]').filter({ hasText: /^ac$/ }).hover();
    const correctlyPainted = await page.locator('#render-inner svg').evaluate(async (svg, atoms) => {
      await Promise.all(svg.getAnimations({ subtree: true }).map(animation => animation.finished.catch(() => {})));
      return [...svg.querySelectorAll('path,text')].every(element => {
        const indices = [...(element.getAttribute('class') || '').matchAll(/\batom-(\d+)\b/g)].map(match => Number(match[1]));
        if (!indices.length) return true;
        const selected = indices.some(index => atoms.includes(index));
        return element.classList.contains('res-hl') === selected &&
          Math.abs(Number(getComputedStyle(element).opacity) - (selected ? 1 : 0.12)) < 0.005;
      });
    }, current.residue_map[String(cap.idx)]);
    expect(correctlyPainted).toBe(true);
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
      await expect(page.locator('#cabiln-input')).toHaveClass('ok');
      const redraw = page.waitForResponse(response => isCompletedResponse(response) &&
        new URL(response.url()).pathname === '/render' && response.request().postDataJSON().cabiln === expected);
      await page.locator('#btn-redo').click();
      current = await (await redraw).json();
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
