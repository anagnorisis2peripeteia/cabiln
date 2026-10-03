const fs = require('node:fs');
const { test, expect, render, tile, site, selectChip, connect, capture, isCompletedResponse, residuePoint } = require('./fixtures');

test('panels persist independently and the library starts a peptide', async ({ page }, testInfo) => {
  await expect(page.locator('#btn-mol')).toBeDisabled();
  await expect(page.locator('#btn-png')).toBeDisabled();
  await page.locator('#btn-lib').click();
  await expect(page.locator('#lib-panel')).toHaveClass('open');
  await tile(page, 'K');
  await expect(page.locator('#cabiln-input')).toHaveValue('K');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#btn-build').click();
  await expect(page.locator('#build-panel')).toHaveClass('open');
  await page.locator('#btn-examples').click();
  await expect(page.locator('.example-row').first()).toBeVisible();
  await page.locator('#btn-verify').click();
  await expect(page.locator('#verify-pane')).toBeVisible();
  await expect(page.locator('#lib-panel')).toHaveClass('open');
  await expect(page.locator('#build-panel')).toHaveClass('open');
  await capture(page, testInfo, 'panels-open');
  await page.locator('#examples-close').click();
  await page.locator('#build-close').click();
  await page.locator('#lib-close').click();
  await page.locator('#btn-verify').click();
  await expect(page.locator('#examples-panel')).not.toHaveClass('open');
  await expect(page.locator('#build-panel')).not.toHaveClass('open');
  await expect(page.locator('#lib-panel')).not.toHaveClass('open');
  await expect(page.locator('#verify-pane')).toBeHidden();
  await expect(page.locator('#cabiln-input')).toHaveValue('K');
  await page.locator('#btn-build').click();
  await expect(page.locator('#lib-panel')).toHaveClass('open');
  await expect(page.locator('#build-left-abbr')).toHaveText('—');
  await expect(page.locator('#build-connect')).toBeDisabled();
  await capture(page, testInfo, 'panels-reopened');
});

test('right-click tiles and detected sites build backbone, branch and cap', async ({ page }, testInfo) => {
  await page.locator('#btn-build').click();
  await tile(page, 'K');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await selectChip(page, 0, 'left', 'K');
  await expect(page.locator('#build-left-rgroups button')).toHaveText([
    'R1 primary amine', 'R2 carboxyl', 'R3 primary amine',
    'R4 primary amine', 'R5 primary amine',
  ]);
  await tile(page, 'G', 'right');
  await expect(page.locator('#build-right-abbr')).toHaveText('G');
  await expect(page.locator('#build-connect')).toBeDisabled();
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await capture(page, testInfo, 'backbone-ready');
  const backbone = await connect(page);
  expect(backbone.submitted).toMatchObject({ host_residue_idx: 0, new_abbr: 'G', r_host: 2, r_new: 1, target_residue_idx: -1 });
  await expect(page.locator('#cabiln-input')).toHaveValue('K-G');

  await selectChip(page, 0, 'left', 'K');
  await expect(page.locator('#build-left-rgroups button').filter({ hasText: /^R2 / })).toHaveClass(/used/);
  await site(page, 'left', 4);
  await page.locator('#btn-rxn-filter').click();
  await expect(page.locator('#btn-rxn-filter')).toHaveClass(/active/);
  await tile(page, 'A', 'right');
  await site(page, 'right', 2);
  await connect(page);
  await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['K', 'A', 'G']);

  await selectChip(page, 0, 'left', 'K');
  await site(page, 'left', 1);
  await tile(page, 'ac', 'right');
  await site(page, 'right', 2);
  await connect(page);
  await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['ac', 'K', 'A', 'G']);
  await expect(page.locator('#build-connect')).toBeDisabled();
  await capture(page, testInfo, 'branched-capped-product');
});

test('repeat occurrences, clear controls and insertion keep the selected position', async ({ page }, testInfo) => {
  await render(page, 'A-A');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'A');
  await selectChip(page, 1, 'right', 'A');
  await expect(page.locator('#build-insert-row')).toBeVisible();
  await page.locator('#build-insert-btn').click();
  await expect(page.locator('#build-insert-btn')).toHaveText('✕ Cancel');
  const inserted = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/insert_backbone');
  await tile(page, 'G');
  const response = await inserted;
  expect(response.request().postDataJSON()).toEqual({ cabiln: 'A-A', after_idx: 0, new_abbr: 'G' });
  expect(response.status(), await response.text()).toBe(200);
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G-A');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#build-insert-row')).toBeHidden();
  await selectChip(page, 2, 'left', 'A');
  await expect(page.locator('#build-left-rgroups button').filter({ hasText: /^R1 / })).toHaveClass(/used/);
  await expect(page.locator('#build-left-rgroups button').filter({ hasText: /^R2 / })).not.toHaveClass(/used/);
  await capture(page, testInfo, 'repeated-last-occurrence');
  for (const side of ['left', 'right']) {
    await tile(page, 'K', 'right');
    await page.locator(`#build-${side}-change`).click();
    await expect(page.locator('#build-left-abbr')).toHaveText('—');
    await expect(page.locator('#build-right-abbr')).toHaveText('—');
    await expect(page.locator('#build-connect')).toBeDisabled();
    await selectChip(page, 2, 'left', 'A');
  }
});

test('successive substitutions update chemistry without renumbering and Undo restores it', async ({ page }) => {
  await render(page, 'K');
  await page.locator('#btn-build').click();
  const left = slot => page.locator('#build-left-rgroups button').filter({ hasText: new RegExp(`^R${slot} `) });
  await selectChip(page, 0, 'left', 'K');
  await expect(left(5)).toHaveText('R5 primary amine');
  for (const slot of [4, 5]) {
    if (slot === 5) {
      await selectChip(page, 0, 'left', 'K');
      await expect(left(4)).toBeDisabled();
      await expect(left(5)).toHaveText('R5 secondary amine');
      await expect(left(5)).toBeEnabled();
    }
    await site(page, 'left', slot);
    if (slot === 4) await page.locator('#btn-rxn-filter').click();
    await tile(page, 'TBMB', 'right');
    const checked = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/validate_bond');
    await site(page, 'right', 4);
    const response = await checked;
    expect(response.request().postDataJSON()).toMatchObject({
      cabiln: await page.locator('#cabiln-input').inputValue(), residue_idx_a: 0,
      residue_idx_b: -1, slot_a: slot, slot_b: 4,
    });
    expect((await response.json()).valid).toBe(true);
    await page.locator('#build-preview-button').click();
    await expect(page.locator('#build-preview-inner svg')).toBeVisible();
    await connect(page);
  }
  await selectChip(page, 0, 'left', 'K');
  await expect(left(4)).toBeDisabled();
  await expect(left(5)).toBeDisabled();
  await page.locator('#btn-undo').click();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await selectChip(page, 0, 'left', 'K');
  await expect(left(5)).toHaveText('R5 secondary amine');
  await expect(left(5)).toBeEnabled();
  await render(page, 'K.ac(4,2)');
  await selectChip(page, 0, 'left', 'K');
  await expect(left(5)).toHaveText('R5 amide N');
  await expect(left(5)).toBeEnabled();
});

test('sulfonamide controls reject amine coupling and retain alkylation', async ({ page }) => {
  await render(page, 'K.pbf(4,2)');
  await page.locator('#btn-build').click();
  const left = slot => page.locator('#build-left-rgroups button').filter({ hasText: new RegExp(`^R${slot} `) });
  await selectChip(page, 0, 'left', 'K');
  await expect(left(5)).toHaveText('R5 sulfonamide N');
  await site(page, 'left', 5);
  await tile(page, 'ac', 'right');
  await site(page, 'right', 2);
  await expect(page.locator('#build-status')).toContainText('No reaction');
  await expect(page.locator('#build-connect')).toBeDisabled();
  await tile(page, 'TBMB', 'right');
  await site(page, 'right', 4);
  await expect(page.locator('#build-connect')).toBeEnabled();
  await connect(page);
  await selectChip(page, 0, 'left', 'K');
  await expect(left(5)).toHaveText('R5 sulfonamide N · used');
  await expect(left(5)).toBeDisabled();
  await page.locator('#btn-undo').click();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await selectChip(page, 0, 'left', 'K');
  await expect(left(5)).toHaveText('R5 sulfonamide N');
  await expect(left(5)).toBeEnabled();
  await render(page, 'C._Me(4,2)');
  await selectChip(page, 0, 'left', 'C');
  await expect(left(4)).toHaveText('R4 thioether · used');
  await expect(left(4)).toBeDisabled();
});

async function clickSvgOccurrence(page, rendered, occurrence) {
  // Find a visible point whose topmost painted SVG element belongs to the
  // requested occurrence. This exercises the browser's SVG hit testing.
  const point = await page.evaluate(({ ownership, occurrence }) => {
    const owned = new Set(ownership[String(occurrence)]);
    for (const element of document.querySelectorAll('#render-inner svg [class]')) {
      const atom = (element.getAttribute('class') || '').match(/atom-(\d+)/);
      if (!atom || !owned.has(Number(atom[1]))) continue;
      const box = element.getBoundingClientRect();
      for (let x = 0; x < 5; x++) for (let y = 0; y < 5; y++) {
        const px = box.x + box.width * (x + 0.5) / 5;
        const py = box.y + box.height * (y + 0.5) / 5;
        const top = document.elementFromPoint(px, py);
        const hit = top?.getAttribute('class')?.match(/atom-(\d+)/);
        if (top?.closest('#render-inner') && hit && owned.has(Number(hit[1]))) return { x: px, y: py };
      }
    }
    return null;
  }, { ownership: rendered.residue_map, occurrence });
  expect(point, 'A painted atom can be selected through the molecular drawing').not.toBeNull();
  await page.mouse.click(point.x, point.y);
}

test('SVG selection connects the intended existing cysteines', async ({ page }, testInfo) => {
  await page.locator('#btn-build').click();
  const data = await render(page, 'C-A-C');
  await clickSvgOccurrence(page, data, 0);
  await expect(page.locator('#build-left-abbr')).toHaveText('C');
  await expect(page.locator('#build-left-rgroups button')).toHaveCount(4);
  const point = await residuePoint(page, data, 2);
  await page.mouse.click(point.x, point.y);
  await expect(page.locator('#build-right-abbr')).toHaveText('C');
  await site(page, 'left', 4);
  await site(page, 'right', 4);
  const { submitted } = await connect(page);
  expect(submitted).toMatchObject({ host_residue_idx: 0, target_residue_idx: 2, r_host: 4, r_new: 4 });
  await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['C', 'A', 'C']);
  await expect(page.locator('#residue-chips .xlink-chip')).toHaveCount(2);
  await capture(page, testInfo, 'disulfide-through-svg-and-chip');
});

test('drawing selection accepts bond margins and label interiors with highlight off and after zoom', async ({ page }) => {
  const data = await render(page, 'A-G');
  await page.locator('#btn-build').click();
  await page.locator('.lib-row').first().waitFor();
  await page.locator('#btn-hl').click();
  for (const label of [false, true]) {
    const point = await residuePoint(page, data, 1, label);
    await page.mouse.click(point.x, point.y);
    await expect(page.locator('#build-left-abbr')).toHaveText('G');
    await expect(page.locator('#build-left-rgroups button')).toHaveCount(3);
    await page.locator('#build-left-change').click();
  }
  const canvas = page.locator('#render-canvas');
  await canvas.click({ position: { x: 12, y: 12 } });
  await expect(page.locator('#build-left-abbr')).toHaveText('—');
  const point = await residuePoint(page, data, 1);
  await page.mouse.move(point.x, point.y);
  await page.mouse.down();
  await page.mouse.move(point.x + 45, point.y + 25, { steps: 5 });
  await page.mouse.up();
  await expect(page.locator('#build-left-abbr')).toHaveText('—');
  await canvas.focus();
  await page.keyboard.press('+');
  const zoomed = await residuePoint(page, data, 1);
  await page.mouse.click(zoomed.x, zoomed.y);
  await expect(page.locator('#build-left-abbr')).toHaveText('G');
});

test('backbone end sites close a cycle through the builder', async ({ page }, testInfo) => {
  await render(page, 'A-G-A');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'A');
  await selectChip(page, 2, 'right', 'A');
  await site(page, 'left', 1);
  await site(page, 'right', 2);
  const { submitted } = await connect(page);
  expect(submitted).toMatchObject({ host_residue_idx: 0, target_residue_idx: 2, r_host: 1, r_new: 2 });
  await expect(page.locator('#residue-chips .xlink-chip')).toHaveCount(2);
  await capture(page, testInfo, 'backbone-cycle');
});

test('three scaffold ports connect through the same tile and site controls', async ({ page }, testInfo) => {
  await render(page, 'ac-C-C-C-am');
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', 'C');
  await tile(page, 'TBMB', 'right');
  await site(page, 'left', 4);
  await site(page, 'right', 4);
  await connect(page);
  await expect(page.locator('#residue-chips [data-residue="5"]')).toHaveText('TBMB');
  for (const [scaffoldSlot, cysteine] of [[5, 2], [6, 3]]) {
    await selectChip(page, 5, 'left', 'TBMB');
    await selectChip(page, cysteine, 'right', 'C');
    await site(page, 'left', scaffoldSlot);
    await site(page, 'right', 4);
    const { submitted } = await connect(page);
    expect(submitted).toMatchObject({ host_residue_idx: 5, target_residue_idx: cysteine, r_host: scaffoldSlot, r_new: 4 });
  }
  await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['ac', 'C', 'TBMB', 'C', 'C', 'am']);
  await expect(page.locator('#residue-chips .xlink-chip')).toHaveCount(4);
  await capture(page, testInfo, 'three-port-scaffold-built');
});

test('library hover preview appears and closes when the pointer leaves', async ({ page }, testInfo) => {
  await page.locator('#btn-lib').click();
  await page.locator('#lib-search').fill('G');
  await page.locator('.lib-row[data-abbr="G"]').hover();
  await expect(page.locator('#lib-preview svg').first()).toBeVisible();
  await capture(page, testInfo, 'library-preview');
  await page.locator('#lib-preview').hover();
  await page.waitForTimeout(200);
  await expect(page.locator('#lib-preview')).toBeVisible();
  await page.locator('#lib-preview').focus();
  await page.keyboard.press('Escape');
  await expect(page.locator('#lib-preview')).toBeHidden();
  await expect(page.locator('#lib-search')).toBeFocused();
  await expect(page.locator('#lib-panel')).toBeVisible();
  await page.locator('.lib-row[data-abbr="G"]').hover();
  await expect(page.locator('#lib-preview')).toBeVisible();
  await page.locator('#cabiln-input').hover();
  await expect(page.locator('#lib-preview')).toBeHidden();
});

test('library insertion replaces selected text and restores the input caret', async ({ page }, testInfo) => {
  await render(page, 'A-G');
  await page.locator('#btn-lib').click();
  await page.locator('#cabiln-input').focus();
  await page.locator('#cabiln-input').press('End');
  await page.locator('#cabiln-input').press('Shift+ArrowLeft');
  // Search remains outside the text area, as in the user's normal palette flow.
  await tile(page, 'A');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-A');
  await expect(page.locator('#cabiln-input')).toBeFocused();
  expect(await page.locator('#cabiln-input').evaluate(el => [el.selectionStart, el.selectionEnd])).toEqual([3, 3]);
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await capture(page, testInfo, 'tile-insert-caret');
});

test('atom and stereo glyphs and carbon joins retain highlighting under the pointer', async ({ page }, testInfo) => {
  const rendered = await render(page, 'A');
  const nitrogen = rendered.mol_block.split('\n').slice(4).findIndex(line => line.slice(31, 34).trim() === 'N');
  expect(nitrogen).toBeGreaterThanOrEqual(0);
  const svg = page.locator('#render-inner svg');
  const chip = page.locator('#residue-chips [data-residue="0"]');
  const targets = await svg.evaluate((drawing, nitrogen) => {
    return [...drawing.querySelectorAll('path')].flatMap((path, index) => {
      const kind = path.classList.contains('CIP_Code') ? 'stereo glyph'
        : path.classList.contains(`atom-${nitrogen}`)
        && !/\bbond-\d+\b/.test(path.getAttribute('class') || '') ? 'NH2 glyph'
        : (path.getAttribute('d').match(/[a-z]/gi) || []).join('') === 'MLL' ? 'carbon join' : null;
      if (!kind) return [];
      const box = path.getBoundingClientRect();
      for (let x = 0; x < 15; x++) for (let y = 0; y < 15; y++) {
        const px = box.x + box.width * (x + 0.5) / 15;
        const py = box.y + box.height * (y + 0.5) / 15;
        if (document.elementFromPoint(px, py) === path) return [{ index, kind, x: px, y: py }];
      }
      throw new Error(`No painted point found for ${kind} path ${index}`);
    });
  }, nitrogen);
  expect(targets.filter(target => target.kind === 'NH2 glyph').length).toBeGreaterThanOrEqual(3);
  expect(targets.filter(target => target.kind === 'stereo glyph').length).toBe(3);
  expect(targets.filter(target => target.kind === 'carbon join').length).toBeGreaterThan(0);

  await chip.hover();
  await expect(svg).toHaveClass(/has-highlight/);
  expect.soft(await svg.locator('path:not(.res-hl)').count(), 'Every alanine path belongs to its selected residue').toBe(0);
  const observed = [];
  for (const target of targets) {
    await page.mouse.move(target.x, target.y);
    await expect.soft(chip, `${target.kind} hover selects alanine`).toHaveClass(/hover/, { timeout: 500 });
    const path = svg.locator('path').nth(target.index);
    await expect.soft(path, `${target.kind} remains visible`).toHaveCSS('opacity', '1', { timeout: 500 });
    observed.push({ ...target, ...await path.evaluate(element => ({
      classes: element.getAttribute('class'), opacity: getComputedStyle(element).opacity,
      highlighted: element.closest('svg').classList.contains('has-highlight'),
    })) });
  }
  await testInfo.attach('painted-highlight-targets', { body: JSON.stringify(observed, null, 2), contentType: 'application/json' });
});

test('canvas toggles, zoom, pan, reset, reroll and downloads remain usable', async ({ page }, testInfo) => {
  await render(page, 'ac-K-G-am');
  const canvas = page.locator('#render-canvas');
  await expect(canvas).toHaveClass(/dark/);
  await page.locator('#btn-dark').click();
  await expect(canvas).not.toHaveClass(/dark/);
  await page.locator('#btn-dark').click();
  await page.locator('#residue-chips [data-residue="1"]').hover();
  await expect(page.locator('#render-inner svg')).toHaveClass(/has-highlight/);
  await page.locator('#btn-hl').click();
  await page.locator('#residue-chips [data-residue="1"]').hover();
  await expect(page.locator('#render-inner svg')).not.toHaveClass(/has-highlight/);
  await page.locator('#btn-hl').click();
  await canvas.hover();
  await page.mouse.wheel(0, -100);
  await expect(page.locator('#render-inner')).toHaveCSS('transform', /1\.12/);
  const box = await canvas.boundingBox();
  await page.mouse.move(box.x + 15, box.y + 15);
  await page.mouse.down();
  await page.mouse.move(box.x + 55, box.y + 45);
  await page.mouse.up();
  await expect(page.locator('#render-inner')).toHaveCSS('transform', /40, 30/);
  await canvas.dblclick({ position: { x: 15, y: 15 } });
  await expect(page.locator('#render-inner')).toHaveCSS('transform', 'matrix(1, 0, 0, 1, 0, 0)');
  for (const [seed, label] of [[2, 'CoordGen'], [3, 'Indigo']]) {
    const rendered = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/render');
    await page.locator('#btn-reroll').click();
    const response = await rendered;
    expect(response.request().postDataJSON().seed).toBe(seed);
    expect(response.status()).toBe(200);
    await expect(page.locator('#btn-reroll')).toHaveText(`⟳ ${label}`);
    await expect(page.locator('#render-inner svg')).toBeVisible();
    await expect(page.locator('#btn-reroll')).toBeEnabled();
  }
  for (const extension of ['mol', 'png']) {
    const downloaded = page.waitForEvent('download');
    await page.locator(`#btn-${extension}`).click();
    const download = await downloaded;
    expect(download.suggestedFilename()).toBe(`structure.${extension}`);
    const destination = testInfo.outputPath(download.suggestedFilename());
    await download.saveAs(destination);
    const bytes = fs.readFileSync(destination);
    if (extension === 'mol') expect(bytes.toString()).toContain('M  END');
    else expect(bytes.subarray(0, 8)).toEqual(Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]));
  }
  await capture(page, testInfo, 'canvas-controls');
});
