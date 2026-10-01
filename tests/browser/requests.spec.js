const { test, expect, render, tile, site, selectChip, capture, isCompletedResponse } = require('./fixtures');

function deferred() {
  let resolve;
  const promise = new Promise(done => { resolve = done; });
  return { promise, resolve };
}

for (const previous of ['', 'G']) {
  test(`immediate notation conversion draws its new input from ${previous || 'an empty document'}`, async ({ page }) => {
    if (previous) await render(page, previous);
    const source = 'K.[G(4,2).ac(1,2)]-A';
    const requests = [];
    page.on('request', request => {
      const path = new URL(request.url()).pathname;
      if (['/render', '/convert_notation'].includes(path)) requests.push(path);
    });
    const converted = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/convert_notation');
    const rendered = page.waitForResponse(response => isCompletedResponse(response) && new URL(response.url()).pathname === '/render');
    // Both events occur before the browser can start the pending drawing timer.
    await page.evaluate(source => {
      const input = document.getElementById('cabiln-input');
      input.value = source;
      input.dispatchEvent(new Event('input', { bubbles: true }));
      document.getElementById('btn-to-branch').click();
    }, source);
    const conversion = await converted;
    expect(conversion.status(), await conversion.text()).toBe(200);
    const result = (await conversion.json()).result;
    const response = await rendered;
    expect(response.status(), await response.text()).toBe(200);
    const drawing = await response.json();
    expect(response.request().postDataJSON().cabiln).toBe(result);
    await expect(page.locator('#cabiln-input')).toHaveValue(result);
    await expect(page.locator('#cabiln-input')).toHaveClass('ok');
    await expect(page.locator('#render-pane')).toHaveAttribute('aria-busy', 'false');
    await expect(page.locator('#residue-chips [data-residue]')).toHaveCount(4);
    await expect(page.locator('#btn-mol')).toBeEnabled();
    expect(drawing.residues.map(residue => residue.abbr).sort()).toEqual(['A', 'G', 'K', 'ac']);
    expect(requests).toEqual(['/convert_notation', '/render']);
  });
}

test('an older render cannot replace a newer peptide', async ({ page }, testInfo) => {
  const arrived = deferred();
  const release = deferred();
  const finished = deferred();
  await page.route('**/render', async route => {
    if (route.request().postDataJSON().cabiln !== 'A') return route.continue();
    const response = await route.fetch();
    arrived.resolve();
    await release.promise;
    try { await route.fulfill({ response }); } finally { finished.resolve(); }
  });
  try {
    await page.locator('#cabiln-input').fill('A');
    await arrived.promise;
    await render(page, 'G');
    release.resolve();
    await finished.promise;
    await expect(page.locator('#cabiln-input')).toHaveValue('G');
    await expect(page.locator('#residue-chips [data-residue]')).toHaveText(['G']);
    await expect(page.locator('#btn-mol')).toBeEnabled();
    await capture(page, testInfo, 'newer-render-retained');
  } finally { release.resolve(); }
});

test('closing Build discards a delayed residue selection', async ({ page }, testInfo) => {
  await render(page, 'K-G');
  await page.locator('#btn-build').click();
  const arrived = deferred();
  const release = deferred();
  const finished = deferred();
  await page.route('**/monomer_rgroups?*', async route => {
    const response = await route.fetch();
    arrived.resolve();
    await release.promise;
    try { await route.fulfill({ response }); } finally { finished.resolve(); }
  });
  try {
    await page.locator('#residue-chips [data-residue="0"]').click();
    await arrived.promise;
    await page.locator('#build-close').click();
    release.resolve();
    await finished.promise;
    await page.locator('#btn-build').click();
    await expect(page.locator('#build-left-abbr')).toHaveText('—');
    await expect(page.locator('#build-left-rgroups button')).toHaveCount(0);
    await expect(page.locator('#build-connect')).toBeDisabled();
    await capture(page, testInfo, 'closed-build-discards-selection');
  } finally { release.resolve(); }
});

test('failed bond validation leaves a usable retry through the same site controls', async ({ page }, testInfo) => {
  await render(page, 'A');
  await page.locator('#btn-build').click();
  await selectChip(page, 0, 'left', 'A');
  await tile(page, 'G', 'right');
  await page.route('**/validate_bond', route => route.fulfill({
    status: 503, contentType: 'application/json', body: JSON.stringify({ error: 'Temporary test outage' }),
  }), { times: 1 });
  await site(page, 'left', 2);
  await site(page, 'right', 1);
  await expect(page.locator('#build-status')).toHaveText('Temporary test outage');
  await expect(page.locator('#build-connect')).toBeDisabled();
  const right = page.locator('#build-right-rgroups button').filter({ hasText: /^R1 / });
  await right.click();
  await right.click();
  await expect(page.locator('#build-connect')).toBeEnabled();
  await expect(page.locator('#build-status')).toHaveClass(/valid/);
  await capture(page, testInfo, 'bond-validation-recovered');
});

test('editing SMILES invalidates a delayed registration preview', async ({ page, app }, testInfo) => {
  await page.goto(`${app.url}/register`);
  const arrived = deferred();
  const release = deferred();
  const finished = deferred();
  await page.route('**/preview_monomer', async route => {
    const response = await route.fetch();
    arrived.resolve();
    await release.promise;
    try { await route.fulfill({ response }); } finally { finished.resolve(); }
  });
  try {
    await page.locator('#smiles-in').fill('NCC(=O)O');
    await page.locator('#btn-preview').click();
    await arrived.promise;
    await page.locator('#smiles-in').fill('N[C@@H](C)C(=O)O');
    release.resolve();
    await finished.promise;
    await expect(page.locator('#preview-section')).toBeHidden();
    await expect(page.locator('#btn-register')).toBeDisabled();
    await expect(page.locator('#btn-preview')).toBeEnabled();
    await expect(page.locator('#chuckles-out')).toHaveValue('');
    await page.screenshot({ path: testInfo.outputPath('stale-registration-preview-discarded.png'), animations: 'disabled' });
  } finally { release.resolve(); }
});

test('a delayed library preview stays hidden after the library closes', async ({ page }) => {
  const arrived = deferred();
  const release = deferred();
  const finished = deferred();
  await page.route('**/monomer_svg?abbr=G', async route => {
    const response = await route.fetch();
    arrived.resolve();
    await release.promise;
    try { await route.fulfill({ response }); } finally { finished.resolve(); }
  });
  try {
    await page.locator('#btn-lib').click();
    await page.locator('#lib-search').fill('G');
    await page.locator('.lib-row[data-abbr="G"]').hover();
    await arrived.promise;
    await page.locator('#lib-close').click();
    release.resolve();
    await finished.promise;
    // A fixed client aborts this request, so a response event is not guaranteed.
    // The Node test also resolves a cancelled fetch to cover already-queued data.
    await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
    await expect(page.locator('#lib-preview')).toBeHidden();
  } finally { release.resolve(); }
});
