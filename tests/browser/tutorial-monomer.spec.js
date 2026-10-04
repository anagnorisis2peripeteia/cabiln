const fs = require('node:fs/promises');
const { test, expect, render, site, selectChip, interceptOnce, capture } = require('./fixtures');

async function nextStep(page, title) {
  await page.locator('#tutorial-next').click();
  await expect(page.locator('#tutorial-title')).toHaveText(title);
}

async function startMonomerLesson(page, url) {
  await page.goto(`${url}/?tutorial=monomer`);
  await page.locator('#tutorial-load').click();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await nextStep(page, 'Enter a practice molecule');
  await page.locator('#tutorial-show').click();
  expect(await page.evaluate(() => {
    const text = document.createRange();
    text.selectNodeContents(document.getElementById('lib-title'));
    const title = text.getBoundingClientRect(), add = document.getElementById('register-link').getBoundingClientRect();
    return title.right <= add.left || title.bottom <= add.top || add.bottom <= title.top;
  }), 'Library title and Add monomer do not overlap').toBe(true);
  await page.locator('#register-link').click();
  const form = page.frameLocator('#registration-frame');
  await expect(form.locator('#smiles-in')).toBeFocused();
  await expect(page.locator('#registration-dialog #tutorial-panel')).toBeVisible();
  await expect(form.locator('#registration-destination')).toBeHidden();
  await page.locator('#tutorial-example').click();
  await nextStep(page, 'See where the block can connect');
  return form;
}

async function addPracticeBlock(page, form, abbr) {
  await form.locator('#btn-preview').click();
  await expect(form.locator('#preview-canvas svg')).toBeVisible();
  await expect(form.locator('#detected-display')).toContainText('R1');
  await expect(form.locator('#detected-display')).toContainText('R2');
  await expect(form.locator('#detected-display')).toBeInViewport({ ratio: .9 });
  await expect(page.locator('#tutorial-cue-label')).toHaveText('Your block’s connection points');
  await nextStep(page, 'Add the block to this tab');
  await form.locator('#abbr-in').fill(abbr);
  await form.locator('#btn-register').click();
  await expect(form.locator('#status-msg.ok')).toContainText(`${abbr} added to this tab`);
  await nextStep(page, 'Find your block in Library');
  await form.locator('#registration-done').click();
  await expect(page.locator('#registration-dialog')).not.toBeVisible();
  await expect(page.locator('#cabiln-input')).toHaveValue('A-G');
}

for (const width of [1440, 390]) test(`custom monomer lesson builds, saves and restores in an isolated tab at ${width}px`, async ({ page, context, readonlyApp }, testInfo) => {
  test.setTimeout(120_000);
  await page.setViewportSize({ width, height: 1000 });
  await page.goto(readonlyApp.url);
  await render(page, 'ac-K-am');
  await expect.poll(() => page.evaluate(() => JSON.parse(localStorage.getItem('cabiln.draft.v1'))?.document.text)).toBe('ac-K-am');
  const draft = await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'));
  const form = await startMonomerLesson(page, readonlyApp.url);
  const abbr = `Lesson${width}`;
  await page.locator('#tutorial-show').click();
  await expect(page.locator('#tutorial-cue-label')).toHaveText('Click Preview & detect R-groups');
  const ring = await page.locator('#tutorial-cue-ring').boundingBox();
  const button = await form.locator('#btn-preview').boundingBox();
  expect(Math.abs(ring.x - button.x + 4)).toBeLessThan(2);
  expect(Math.abs(ring.y - button.y + 4)).toBeLessThan(2);
  await capture(page, testInfo, 'monomer-form-guide');
  await addPracticeBlock(page, form, abbr);
  await nextStep(page, 'Choose the marked G tile');
  await page.locator('#btn-build').click();
  await selectChip(page, 1, 'left', 'G');
  await nextStep(page, 'Choose a block to add');
  await expect(page.locator('#tutorial-body')).toContainText(abbr);
  await page.locator(`.lib-row[data-abbr="${abbr}"] .lib-use`).click();
  await nextStep(page, 'Choose the connection');
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await nextStep(page, 'Check the chain before adding');
  await page.locator('#build-preview-button').click();
  await expect(page.locator('#build-preview-source')).toHaveText(`A-G-${abbr}`);
  await nextStep(page, 'Add the new block');
  await page.locator('#build-connect').click();
  await expect(page.locator('#cabiln-input')).toHaveValue(`A-G-${abbr}`);
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');

  // Reviewing completed instructions never asks for the consumed Build selections.
  for (let i = 0; i < 6; i++) {
    await page.locator('#tutorial-back').click();
    await expect(page.locator('#tutorial-next')).toBeEnabled();
  }
  for (let i = 0; i < 7; i++) await page.locator('#tutorial-next').click();
  await expect(page.locator('#tutorial-title')).toHaveText('Keep the block with your project');
  await expect(page.locator('#tutorial-next')).toBeDisabled();
  await interceptOnce(page, '**/prepare_project', route => route.fulfill({ status: 500,
    contentType: 'application/json', body: JSON.stringify({ error: 'Practice save interrupted' }) }));
  await page.locator('#btn-project-save').click();
  await expect(page.locator('#project-status')).toContainText('Practice save interrupted');
  await expect(page.locator('#tutorial-next')).toBeDisabled();
  const download = page.waitForEvent('download');
  await page.locator('#btn-project-save').click();
  const file = await (await download).path();
  const project = JSON.parse(await fs.readFile(file, 'utf8'));
  expect(project.document.text).toBe(`A-G-${abbr}`);
  expect(project.monomers.map(item => item.abbr)).toEqual([abbr]);
  await nextStep(page, 'Restore your custom block');
  const opened = context.waitForEvent('page');
  await page.locator('#tutorial-open-tab').click();
  const fresh = await opened;
  await fresh.setViewportSize({ width, height: 1000 });
  await expect(fresh.locator('#tutorial-title')).toHaveText('Restore your custom block');
  expect(await fresh.evaluate(() => window.opener)).toBe(null);
  expect(await fresh.evaluate(() => window.CabilnLibrary.monomers)).toEqual([]);
  await expect(fresh.locator('#tutorial-next')).toBeDisabled();
  await fresh.locator('#project-upload').setInputFiles({ name: 'broken.json', mimeType: 'application/json', buffer: Buffer.from('{') });
  await expect(fresh.locator('#project-status')).toContainText('not valid project JSON');
  await expect(fresh.locator('#tutorial-next')).toBeDisabled();
  await fresh.locator('#project-upload').setInputFiles({ name: 'peptide.cabiln.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(project)) });
  await expect(fresh.locator('#project-status')).toContainText('Project opened');
  await expect(fresh.locator('#cabiln-input')).toHaveClass('ok');
  expect(await fresh.evaluate(() => window.CabilnLibrary.monomers)).toEqual(project.monomers);
  await expect(fresh.locator('#tutorial-next')).toBeEnabled();
  await capture(fresh, testInfo, 'monomer-project-restored');
  await nextStep(fresh, 'Your project carries its building blocks');
  await expect(fresh.locator('#tutorial-body')).toContainText('closing it can lose them');
  await fresh.locator('#tutorial-next').click();
  await expect(fresh.locator('#tutorial-panel')).toBeHidden();
  await fresh.locator('#btn-tutorial').click();
  await fresh.locator('#tutorial-restart').click();
  await expect(fresh.locator('#tutorial-progress')).toHaveText('Step 1 of 13');
  await expect(fresh).toHaveURL(/\?tutorial=monomer$/);
  await expect(fresh.locator('#cabiln-input')).toHaveClass('ok');
  expect(await page.evaluate(() => localStorage.getItem('cabiln.draft.v1'))).toBe(draft);
  await fresh.close();
});

test('custom monomer lesson recovers closed forms, changed inputs, reloads and lesson changes on a writable local app', async ({ page, app }) => {
  let form = await startMonomerLesson(page, app.url);
  await form.locator('#smiles-in').click();
  await form.locator('#smiles-in').fill('invalid');
  await expect(form.locator('#smiles-in')).toHaveValue('invalid');
  await form.locator('#btn-preview').click();
  await expect(form.locator('#preview-error')).not.toBeEmpty();
  await expect(page.locator('#tutorial-next')).toBeDisabled();
  await page.locator('#registration-close').click();
  await expect(page.locator('body > #tutorial-panel')).toBeVisible();
  await expect(page.locator('#tutorial-instruction')).toHaveText('Click Add monomer');
  await page.locator('#register-link').click();
  form = page.frameLocator('#registration-frame');
  await page.locator('#tutorial-example').click();
  await addPracticeBlock(page, form, 'ReusedBlock');
  // WebKit reports navigation-cancelled fetches as page errors. This check
  // reloads a completed addition; the separate pending-form journey holds its response.
  await expect.poll(() => page.evaluate(() => requests.has('main-render', 'library', 'library-focus', 'reactions'))).toBe(false);
  await page.reload();
  await page.locator('#tutorial-load').click();
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  for (let i = 0; i < 4; i++) {
    await page.locator('#tutorial-next').click();
    await expect(page.locator('#tutorial-next')).toBeEnabled();
  }
  await expect(page.locator('#tutorial-status')).toContainText('ReusedBlock');
  await page.locator('#tutorial-lesson').selectOption('connect');
  await expect(page.locator('#tutorial-progress')).toHaveText('Step 1 of 9');
  await expect(page).toHaveURL(/\?tutorial=1$/);
  await page.locator('#tutorial-lesson').selectOption('monomer');
  await expect(page.locator('#tutorial-progress')).toHaveText('Step 1 of 13');
  await expect(page).toHaveURL(/\?tutorial=monomer$/);
  expect(await page.evaluate(() => window.CabilnLibrary.monomers.map(item => item.abbr))).toEqual(['ReusedBlock']);
  const shared = await page.request.get(`${app.url}/monomers`);
  expect((await shared.json()).some(item => item.abbr === 'ReusedBlock')).toBe(false);
});

test('custom monomer lesson records an addition completed after closing the guide and form', async ({ page, readonlyApp }) => {
  await page.emulateMedia({ reducedMotion: 'reduce' });
  const form = await startMonomerLesson(page, readonlyApp.url);
  await form.locator('#btn-preview').click();
  await expect(page.locator('#tutorial-next')).toBeEnabled();
  await nextStep(page, 'Add the block to this tab');
  await form.locator('#abbr-in').fill('LateBlock');
  await page.locator('#tutorial-show').click();
  expect(await page.locator('#tutorial-cue').evaluate(el => el.getAnimations({ subtree: true }).length)).toBe(0);
  // The form can finish while both its guide and iframe are dismissed.
  await page.locator('#tutorial-close').click();
  let arrived, release;
  const pending = new Promise(resolve => { arrived = resolve; });
  const gate = new Promise(resolve => { release = resolve; });
  await interceptOnce(page, '**/session_library', async route => {
    const response = await route.fetch();
    arrived();
    await gate;
    await route.fulfill({ response });
  });
  await form.locator('#btn-register').click();
  await pending;
  await form.locator('#abbr-in').click();
  await form.locator('#abbr-in').fill('UnusedEdit');
  await page.locator('#registration-close').click();
  release();
  await expect.poll(() => page.evaluate(() => window.CabilnLibrary.monomers.map(item => item.abbr))).toEqual(['LateBlock']);
  await page.locator('#btn-tutorial').click();
  await expect(page.locator('#tutorial-status')).toContainText('LateBlock is in this tab');
  await expect(page.locator('#tutorial-next')).toBeEnabled();
  await nextStep(page, 'Find your block in Library');
  await expect(page.locator('#tutorial-next')).toBeEnabled();
});
