const fs = require('node:fs/promises');
const { test, expect, render, interceptOnce } = require('./fixtures');

async function openRegistration(page) {
  await page.locator('#register-link').click();
  const form = page.frameLocator('#registration-frame');
  await expect(form.locator('#smiles-in')).toBeFocused();
  return form;
}

async function prepareMonomer(form, symbol, smiles = 'N[C@@H](CCCS)C(=O)O') {
  await form.locator('#smiles-in').fill(smiles);
  await form.locator('#btn-preview').click();
  await expect(form.locator('#preview-canvas svg')).toBeVisible();
  await form.locator('#abbr-in').fill(symbol);
  await form.locator('#name-in').fill(symbol + ' test monomer');
  await expect(form.locator('#btn-register')).toBeEnabled();
}

async function saveProject(page) {
  const downloaded = page.waitForEvent('download');
  await page.locator('#btn-project-save').click();
  const project = JSON.parse(await fs.readFile(await (await downloaded).path(), 'utf8'));
  await expect(page.locator('#project-status')).toContainText('Project saved');
  return project;
}

for (const width of [1440, 390]) test(`onboarding survives closing pending forms at ${width}px`, async ({ page, readonlyApp }) => {
  await page.setViewportSize({ width, height: 1000 });
  await page.goto(readonlyApp.url);
  await render(page, 'A-G');
  await page.locator('#btn-lib').click();
  let releasePreview, previewStarted;
  const previewGate = new Promise(resolve => { releasePreview = resolve; });
  const pendingPreview = new Promise(resolve => { previewStarted = resolve; });
  await interceptOnce(page, '**/preview_monomer', async route => {
    previewStarted();
    await previewGate;
    await route.continue();
  });
  let form = await openRegistration(page);
  await form.locator('#smiles-in').fill('NCC(=O)O');
  await form.locator('#btn-preview').click();
  await pendingPreview;
  await form.locator('#smiles-in').press('Escape');
  await expect(page.locator('#registration-dialog')).not.toBeVisible();
  form = await openRegistration(page);
  await expect(form.locator('#smiles-in')).toHaveValue('');
  releasePreview();
  await prepareMonomer(form, 'InterruptedThiol');
  await expect(form.locator('#detected-display')).toContainText('thiol');

  let releaseAddition, additionStored;
  const additionGate = new Promise(resolve => { releaseAddition = resolve; });
  const stored = new Promise(resolve => { additionStored = resolve; });
  await interceptOnce(page, '**/session_library', async route => {
    const response = await route.fetch();
    expect(response.status()).toBe(200);
    additionStored();
    await additionGate;
    await route.fulfill({ response });
  });
  await form.locator('#btn-register').click();
  await stored;
  await page.locator('#registration-close').click();
  await expect(page.locator('#registration-dialog')).not.toBeVisible();
  form = await openRegistration(page);
  releaseAddition();
  await expect.poll(() => page.evaluate(() => window.CabilnLibrary.monomers.map(item => item.abbr)))
    .toEqual(['InterruptedThiol']);
  await expect(form.locator('#smiles-in')).toHaveValue('');
  await expect(form.locator('#btn-register')).toBeDisabled();
  await form.locator('#smiles-in').press('Escape');
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
  await page.locator('#lib-search').fill('InterruptedThiol');
  await expect(page.locator('.lib-row[data-abbr="InterruptedThiol"]')).toBeVisible();
  await render(page, 'InterruptedThiol-G');
});

test('onboarding retries failures and keeps several monomers through project conflicts and expiry', async ({ page, context, readonlyApp }) => {
  test.setTimeout(90_000);
  await page.goto(readonlyApp.url);
  await render(page, 'A-G');
  await page.locator('#btn-lib').click();
  const form = await openRegistration(page);
  await prepareMonomer(form, 'RetryThiol');
  await interceptOnce(page, '**/session_library', route => route.fulfill({
    status: 500, contentType: 'application/json', body: JSON.stringify({ error: 'Temporary service failure. Try again.' }),
  }));
  await form.locator('#btn-register').click();
  await expect(form.locator('#status-msg')).toContainText('Temporary service failure');
  await expect(form.locator('#abbr-in')).toHaveValue('RetryThiol');
  expect(await page.evaluate(() => window.CabilnLibrary.monomers)).toEqual([]);
  await form.locator('#btn-register').click();
  await expect(form.locator('#status-msg.ok')).toContainText('RetryThiol added');
  await prepareMonomer(form, 'RetryAlcohol', 'N[C@@H](CCCO)C(=O)O');
  await form.locator('#btn-register').click();
  await expect(form.locator('#status-msg.ok')).toContainText('RetryAlcohol added');
  await form.locator('#registration-done').click();
  await render(page, 'ac-RetryThiol-RetryAlcohol-am');
  const project = await saveProject(page);
  expect(project.monomers.map(item => item.abbr)).toEqual(['RetryThiol', 'RetryAlcohol']);
  const current = project.document.text;
  const conflict = structuredClone(project);
  conflict.monomers[0].name = 'Another definition';
  await page.locator('#project-upload').setInputFiles({ name: 'conflict.cabiln.json', mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(conflict)) });
  await expect(page.locator('#project-status')).toContainText('another definition');
  await expect(page.locator('#cabiln-input')).toHaveValue(current);
  expect(await page.evaluate(() => window.CabilnLibrary.monomers)).toEqual(project.monomers);

  const other = await context.newPage();
  await other.goto(readonlyApp.url);
  await other.locator('#project-upload').setInputFiles({ name: 'shared.cabiln.json', mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(project)) });
  await expect(other.locator('#project-status')).toContainText('Project opened');
  await expect(other.locator('#cabiln-input')).toHaveClass('ok');
  for (let attempt = 0; attempt < 3; attempt++) {
    await other.evaluate(() => {
      const library = JSON.parse(sessionStorage.getItem('cabiln.monomers.v1'));
      library.token = 'expired';
      sessionStorage.setItem('cabiln.monomers.v1', JSON.stringify(library));
    });
    await other.reload();
    await other.locator('#btn-restore-draft').click();
    await expect(other.locator('#cabiln-input')).toHaveValue(current);
    await expect(other.locator('#cabiln-input')).toHaveClass('ok');
    expect(await other.evaluate(() => window.CabilnLibrary.monomers)).toEqual(project.monomers);
  }
  await other.close();
});

test('onboarding without tab storage still saves a portable project', async ({ page, readonlyApp }) => {
  await page.addInitScript(() => {
    const setItem = Storage.prototype.setItem;
    Storage.prototype.setItem = function (...args) {
      if (this === sessionStorage) throw new DOMException('Storage unavailable', 'QuotaExceededError');
      return setItem.apply(this, args);
    };
  });
  await page.goto(readonlyApp.url);
  await render(page, 'A-G');
  await page.locator('#btn-lib').click();
  const form = await openRegistration(page);
  await prepareMonomer(form, 'MemoryThiol');
  await form.locator('#btn-register').click();
  await expect(form.locator('#status-msg.ok')).toContainText('Tab storage is unavailable');
  await form.locator('#registration-done').click();
  await render(page, 'MemoryThiol-G');
  const project = await saveProject(page);
  expect(project.monomers.map(item => item.abbr)).toEqual(['MemoryThiol']);
  expect(project.document.text).toBe('MemoryThiol-G');
});
