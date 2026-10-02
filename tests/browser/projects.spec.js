const fs = require('node:fs/promises');
const { test, expect, render, selectChip, isCompletedResponse } = require('./fixtures');

async function saveProject(page) {
  const download = page.waitForEvent('download');
  await page.locator('#btn-project-save').click();
  const saved = await download;
  expect(saved.suggestedFilename()).toBe('peptide.cabiln.json');
  const text = await fs.readFile(await saved.path(), 'utf8');
  await expect(page.locator('#project-status')).toContainText('Project saved');
  return JSON.parse(text);
}

async function openProject(page, project) {
  await page.locator('#project-upload').setInputFiles({
    name: 'saved.cabiln.json', mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(project)),
  });
}

test('portable projects preserve import quality, notation drafts and the exact uploaded reference', async ({ page }) => {
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('NC(C)C(=O)O');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  const rendered = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/render');
  await page.locator('#btn-to-cabiln-pct').click();
  const structure = await (await rendered).json();
  await expect(page.locator('#quality-summary')).toContainText('stereochemistry inferred');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  const source = await page.locator('#cabiln-input').inputValue();
  await page.locator('#btn-verify').click();
  const original = structure.mol_block.replace(/\r?\n/g, '\r\n');
  await page.locator('#mol-upload').setInputFiles({
    name: 'original.mol', mimeType: 'chemical/x-mdl-molfile', buffer: Buffer.from(original),
  });
  await expect(page.locator('#smiles-input')).toHaveClass('ok');
  const project = await saveProject(page);
  expect(project.document.text).toBe(source);
  expect(project.document.quality.inferred_stereo).toBe(true);
  expect(project.document.quality.assignments.length).toBeGreaterThan(0);
  expect(project.drafts.smiles.text).toBe('NC(C)C(=O)O');
  expect(project.reference.original).toEqual({ kind: 'mol', name: 'original.mol', content: original });
  expect(project.context.canonical.format).toBe('cabiln-graph-v1');
  expect(project.context.library_binding.monomers).toBeTruthy();
  expect(project.reference.context.library_binding).toEqual(project.context.library_binding);
  await render(page, 'A-G');
  await openProject(page, project);
  await expect(page.locator('#project-status')).toContainText('Project opened');
  await expect(page.locator('#cabiln-input')).toHaveValue(source);
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#quality-summary')).toContainText('stereochemistry inferred');
  await page.locator('#notation-select').selectOption('smiles');
  await expect(page.locator('#cabiln-input')).toHaveValue('NC(C)C(=O)O');
  await page.locator('#notation-select').selectOption('cabiln');
  await expect(page.locator('#quality-summary')).toContainText('stereochemistry inferred');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  expect((await saveProject(page)).reference.original.content).toBe(original);
});

test('binding mismatches explain recovery without replacing current work', async ({ page }) => {
  await render(page, 'G');
  const project = await saveProject(page);
  project.context.library_binding.reactions = 'different-reaction-rules';
  await render(page, 'A-G');
  await openProject(page, project);
  await expect(page.locator('#project-status')).toContainText('current work is unchanged');
  await expect(page.locator('#project-status')).toContainText(/rules|binding|library/i);
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
});

test('a delayed project validation cannot overwrite a newer sequence edit', async ({ page }) => {
  await render(page, 'G');
  const project = await saveProject(page);
  let release;
  const held = new Promise(resolve => { release = resolve; });
  await page.route('**/validate_project', async route => {
    const response = await route.fetch();
    await held;
    await route.fulfill({ response });
  });
  try {
    await openProject(page, project);
    await expect(page.locator('#project-status')).toContainText('Checking the project');
    await render(page, 'A-K');
  } finally { release(); }
  await expect(page.locator('#project-status')).toContainText('changed during the project check');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-K');
});

test('quality clears on edits and canonical reordering, and Undo restores the original evidence', async ({ page }) => {
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('NC(C)C(=O)O');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#btn-to-cabiln-pct').click();
  await expect(page.locator('#import-quality')).toBeVisible();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await render(page, 'G');
  await expect(page.locator('#import-quality')).toBeHidden();
  await page.locator('#btn-undo').click();
  await expect(page.locator('#import-quality')).toBeVisible();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#notation-policy').selectOption('canonical');
  await page.locator('#btn-to-bracket').click();
  await expect(page.locator('#canonical-status')).toContainText('cabiln-graph-v1');
  await expect(page.locator('#import-quality')).toBeHidden();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  const project = await saveProject(page);
  expect(project.document.canonical.format).toBe('cabiln-graph-v1');
  expect(project.document.quality).toBeNull();
  await page.locator('#btn-undo').click();
  await expect(page.locator('#import-quality')).toBeVisible();
});

test('synthetic regions stay editable and library notes are distinct from chemistry warnings', async ({ page }) => {
  const data = await render(page, 'A-<[1*]NCC([2*])=O>');
  const synthetic = page.locator('#residue-chips [data-quality="synthetic"]');
  await expect(synthetic).toHaveCount(1);
  await expect(synthetic).toBeEnabled();
  await expect(synthetic).toHaveAccessibleName(/Synthetic preserved region.*editable/);
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', data.residues[1].abbr);
  await expect(page.locator('#build-left-rgroups button').filter({ hasText: /^R2 / })).toBeEnabled();
  await expect(page.locator('#lib-panel')).toHaveClass('open');
  await page.locator('#lib-search').fill('meC');
  const row = page.locator('.lib-row[data-abbr="meC"]');
  await expect(row.locator('.lib-quality')).toHaveCount(0);
  await row.hover();
  await expect(page.locator('#lib-preview')).toContainText('Library notes:');
  await page.locator('#lib-search').fill('aMeLeu');
  const uncertain = page.locator('.lib-row[data-abbr="aMeLeu"]');
  await expect(uncertain.locator('.lib-quality')).toContainText('stereo unspecified');
  await uncertain.hover();
  await expect(page.locator('#lib-preview .prev-warn')).toContainText('stereo unspecified');
});

test('help explains data handling and clearing a browser draft preserves the current session only', async ({ page }) => {
  await render(page, 'A-G');
  await expect(page.locator('#draft-status')).toContainText('Draft saved');
  await page.locator('#btn-help').click();
  await expect(page.locator('#help-panel')).toContainText('server hosting this page');
  await page.getByText('Privacy and reporting a problem', { exact: true }).click();
  await expect(page.getByRole('link', { name: 'Report a problem' })).toHaveAttribute('href', 'https://github.com/anagnorisis2peripeteia/pyPept/issues');
  await page.getByText('Save, restore and export your work', { exact: true }).click();
  await page.locator('#btn-clear-draft').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#btn-undo')).toBeEnabled();
  expect(await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'))).toBeNull();
  await page.keyboard.press('Escape');
  await expect(page.locator('#help-panel')).toBeHidden();
  await page.reload();
  await expect(page.locator('#btn-restore-draft')).toBeHidden();
  await expect(page.locator('#cabiln-input')).toHaveValue('');
});

test('a cleared recovery draft can be recovered without overwriting new work', async ({ page }) => {
  await render(page, 'A-G');
  await expect(page.locator('#draft-status')).toContainText('Draft saved');
  const original = await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'));
  await page.reload();
  await page.locator('#btn-help').click();
  await page.getByText('Save, restore and export your work', { exact: true }).click();
  await page.locator('#btn-clear-draft').click();
  await page.locator('#btn-undo-clear-draft').click();
  expect(await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'))).toBe(original);
  await expect(page.locator('#cabiln-input')).toHaveValue('');
  await page.locator('#btn-restore-draft').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#btn-clear-draft').click();
  await render(page, 'A-G-L');
  await expect(page.locator('#btn-undo-clear-draft')).toBeHidden();
  await expect.poll(() => page.evaluate(() => JSON.parse(localStorage.getItem('cabiln.draft.v1')).document.text)).toBe('A-G-L');
});

test('a busy project save offers Retry without losing the document', async ({ page }) => {
  await render(page, 'A-G');
  await page.route('**/prepare_project', route => route.fulfill({
    status: 503, contentType: 'application/json', headers: { 'Retry-After': '1' },
    body: JSON.stringify({ error: 'Chemistry capacity is busy; retry shortly' }),
  }));
  await page.locator('#btn-project-save').click();
  const retry = page.locator('#project-status').getByRole('button', { name: 'Retry', exact: true });
  await expect(retry).toBeVisible();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await expect(page.locator('#btn-project-save')).toBeEnabled();
  await page.unroute('**/prepare_project');
  const download = page.waitForEvent('download');
  await retry.click();
  const saved = await download;
  const project = JSON.parse(await fs.readFile(await saved.path(), 'utf8'));
  expect(project.document.text).toBe('A-G');
  await expect(page.locator('#project-status')).toContainText('Project saved');
});

test('short server admission pressure retries a current render without losing input', async ({ page }) => {
  let attempts = 0;
  await page.route('**/render', async route => {
    attempts++;
    if (attempts === 1) {
      await route.fulfill({ status: 503, headers: { 'Retry-After': '1' },
        contentType: 'application/json', body: JSON.stringify({ error: 'Chemistry workers are busy' }) });
    } else await route.continue();
  });
  await page.locator('#cabiln-input').fill('A-G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  expect(attempts).toBe(2);
});
