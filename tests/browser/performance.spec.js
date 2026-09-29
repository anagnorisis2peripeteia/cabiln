const fs = require('node:fs');
const { test, expect, tile, site, selectChip } = require('./fixtures');

// Run this file alone for a fresh server. CABILN_APP_ROOT can select a frozen
// revision while the browser, Python environment, and timing code stay fixed.
test('measure edit and connection latency through the next painted result', async ({ page }, testInfo) => {
  async function measure(selector, eventName, source, action) {
    await page.evaluate(({ selector, eventName, source }) => {
      window.paintMeasurement = new Promise((resolve, reject) => {
        let started;
        const observer = new MutationObserver(() => {
          const input = document.querySelector('#cabiln-input');
          if (started === undefined || input.value !== source || !input.classList.contains('ok') ||
              document.querySelectorAll('#residue-chips [data-residue]').length !== source.split('-').length) return;
          observer.disconnect();
          clearTimeout(timeout);
          requestAnimationFrame(() => requestAnimationFrame(() => resolve(performance.now() - started)));
        });
        const timeout = setTimeout(() => { observer.disconnect(); reject(new Error('Render timed out')); }, 15000);
        observer.observe(document.body, { subtree: true, childList: true, attributes: true });
        document.querySelector(selector).addEventListener(eventName, () => { started = performance.now(); }, { once: true, capture: true });
      });
    }, { selector, eventName, source });
    await action();
    return page.evaluate(() => window.paintMeasurement);
  }

  const edits = [];
  let source;
  for (source of ['A-G', 'A-G-A', 'A-G-A-G', 'A-G-A-G-A', 'A-G-A-G-A-G']) {
    edits.push({ source, milliseconds: await measure('#cabiln-input', 'input', source,
      () => page.locator('#cabiln-input').fill(source)) });
  }
  await page.locator('#btn-build').click();
  const connections = [];
  for (const abbr of ['A', 'G', 'A', 'G', 'A']) {
    const residues = source.split('-');
    await selectChip(page, residues.length - 1, 'left', residues.at(-1));
    await tile(page, abbr, 'right');
    await site(page, 'left', 2);
    await site(page, 'right', 1);
    await expect(page.locator('#build-connect')).toBeEnabled();
    source += `-${abbr}`;
    connections.push({ source, milliseconds: await measure('#build-connect', 'click', source,
      () => page.locator('#build-connect').click()) });
  }
  const result = {
    application: process.env.CABILN_APP_ROOT || 'current checkout',
    browser: process.env.CABILN_BROWSER_CHANNEL || 'chromium',
    viewport: testInfo.project.use.viewport,
    edits, connections,
    note: 'Browser event through completed drawing/chips and two animation frames. Excludes automation polling. Run this file alone for a cold first edit. Small local samples, not production percentiles.',
  };
  const filename = testInfo.outputPath('browser-timings.json');
  fs.writeFileSync(filename, JSON.stringify(result, null, 2));
  await testInfo.attach('browser-timings', { path: filename, contentType: 'application/json' });
});
