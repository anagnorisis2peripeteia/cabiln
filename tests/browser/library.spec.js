const { test, expect, render, selectChip, site } = require('./fixtures');

test('Filter guides selection and changes the library for the chosen attachment site', async ({ page }) => {
  await render(page, 'K');
  await page.locator('#btn-lib').click();
  const filter = page.locator('#btn-rxn-filter');
  const scope = page.locator('#lib-filter-status');
  await expect(filter).toBeEnabled();
  await filter.click();
  await expect(page.locator('#build-panel')).toBeVisible();
  await expect(scope).toContainText('Select a residue');
  await selectChip(page, 0, 'left', 'K');
  await expect(scope).toContainText('any free site');
  await site(page, 'left', 4);
  await expect(scope).toContainText('K · R4');
  await expect(page.locator('.lib-row[data-abbr="ac"]')).toHaveCount(1);
  await expect(page.locator('.lib-row[data-abbr="am"]')).toHaveCount(0);
  await site(page, 'left', 2);
  await expect(scope).toContainText('K · R2');
  await expect(page.locator('.lib-row[data-abbr="am"]')).toHaveCount(1);
  await expect(page.locator('.lib-row[data-abbr="ac"]')).toHaveCount(0);
  await filter.click();
  await expect(page.locator('.lib-row[data-abbr="ac"]')).toHaveCount(1);
  await expect(page.locator('.lib-row[data-abbr="am"]')).toHaveCount(1);
  await filter.click();
  await page.locator('#lib-search').fill('unlikely-monomer-name');
  await expect(page.locator('#lib-list')).toHaveText('No matches');
  await page.locator('#lib-search').fill('');
  await page.locator('#build-close').click();
  await expect(filter).toHaveAttribute('aria-pressed', 'false');
  await expect(page.locator('.lib-row[data-abbr="ac"]')).toHaveCount(1);
  await expect(scope).toContainText('Filter');
});

test('a failed reaction filter stays explained across refresh and Retry restores matching', async ({ page }) => {
  await render(page, 'K');
  await page.route('**/reactions', route => route.fulfill({ status: 500,
    contentType: 'application/json', body: '{"error":"Temporary failure"}' }));
  await page.locator('#btn-build').click();
  const scope = page.locator('#lib-filter-status');
  await expect(scope).toContainText('unavailable');
  await selectChip(page, 0, 'left', 'K');
  await site(page, 'left', 4);
  await page.locator('#btn-rxn-filter').click();
  await expect(scope).toContainText('unavailable');
  await expect(page.locator('.lib-row')).toHaveCount(0);
  const refreshed = page.waitForResponse(response => response.url().endsWith('/monomers'));
  await page.evaluate(() => window.dispatchEvent(new Event('focus')));
  await refreshed;
  await expect(scope.getByRole('button', { name: 'Retry' })).toBeVisible();
  await page.unroute('**/reactions');
  await scope.getByRole('button', { name: 'Retry' }).click();
  await expect(scope).toContainText('K · R4');
  await expect(page.locator('.lib-row[data-abbr="ac"]')).toHaveCount(1);
  await expect(page.locator('.lib-row[data-abbr="am"]')).toHaveCount(0);
});

test('closing Build while reaction rules load cancels the filter intent', async ({ page }) => {
  await render(page, 'K');
  let release;
  const held = new Promise(resolve => { release = resolve; });
  await page.route('**/reactions', async route => {
    const response = await route.fetch();
    await held;
    await route.fulfill({ response });
  });
  const rules = page.waitForResponse(response => response.url().endsWith('/reactions'));
  try {
    await page.locator('#btn-build').click();
    await selectChip(page, 0, 'left', 'K');
    await site(page, 'left', 4);
    await page.locator('#btn-rxn-filter').click();
    await expect(page.locator('#lib-count')).toHaveText('Loading matches…');
    await expect(page.locator('.lib-row')).toHaveCount(0);
    await page.locator('#build-close').click();
  } finally { release(); }
  await rules;
  await expect(page.locator('#btn-rxn-filter')).toHaveAttribute('aria-pressed', 'false');
  await expect(page.locator('.lib-row[data-abbr="ac"]')).toHaveCount(1);
  await expect(page.locator('.lib-row[data-abbr="am"]')).toHaveCount(1);
});

test('Filter explains conversion when the input has no editable residue tiles', async ({ page }) => {
  await page.locator('#notation-select').selectOption('smiles');
  await page.locator('#cabiln-input').fill('NCC(=O)O');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await page.locator('#btn-lib').click();
  await page.locator('#btn-rxn-filter').click();
  await expect(page.locator('#lib-filter-status')).toContainText('Convert this input to CABILN');
  await expect(page.locator('#build-hint')).toContainText('Convert this input to CABILN');
  await page.locator('#btn-to-cabiln-pct').click();
  await expect(page.locator('#cabiln-input')).toHaveValue('G');
  await expect(page.locator('#cabiln-input')).toHaveClass('ok');
  await selectChip(page, 0, 'left', 'G');
  await page.locator('#btn-rxn-filter').click();
  await expect(page.locator('#lib-filter-status')).toContainText('Compatible with G');
  await expect(page.locator('.lib-row[data-abbr="ac"]')).toHaveCount(1);
});
