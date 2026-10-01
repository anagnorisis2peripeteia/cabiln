const fs = require('node:fs');
const { test, expect, tile, site, selectChip } = require('./fixtures');

// Run this file alone for a fresh server. CABILN_APP_ROOT can select a frozen
// revision while the browser, Python environment, and timing code stay fixed.
test('measure feedback and completed chemistry for edits, connections and larger peptides', async ({ page }, testInfo) => {
  test.setTimeout(120_000);
  async function measure(selector, eventName, source, action, residues = source.split('-').length) {
    await page.evaluate(({ selector, eventName, source, residues }) => {
      window.paintMeasurement = new Promise((resolve, reject) => {
        let started, feedbackMs, completing = false;
        const afterPaint = callback => requestAnimationFrame(() => requestAnimationFrame(callback));
        const observer = new MutationObserver(() => {
          const input = document.querySelector('#cabiln-input');
          if (started === undefined || completing || input.value !== source || !input.classList.contains('ok') ||
              document.querySelectorAll('#residue-chips [data-residue]').length !== residues) return;
          completing = true;
          observer.disconnect();
          clearTimeout(timeout);
          afterPaint(() => resolve({ feedbackMs, completedMs: performance.now() - started }));
        });
        const timeout = setTimeout(() => { observer.disconnect(); reject(new Error('Render timed out')); }, 30000);
        observer.observe(document.body, { subtree: true, childList: true, attributes: true });
        document.querySelector(selector).addEventListener(eventName, () => {
          started = performance.now();
          afterPaint(() => {
            const busy = selector === '#cabiln-input'
              ? document.querySelector('#render-pane').getAttribute('aria-busy') === 'true'
              : document.querySelector('#build-connect').disabled && document.querySelector('#build-status').textContent;
            if (busy || completing) feedbackMs = performance.now() - started;
          });
        }, { once: true, capture: true });
      });
    }, { selector, eventName, source, residues });
    await action();
    const timing = await page.evaluate(() => window.paintMeasurement);
    expect(timing.feedbackMs, 'The action paints immediate feedback or its completed result').toBeGreaterThanOrEqual(0);
    return timing;
  }

  const edits = [];
  let source;
  for (source of ['A-G', 'A-G-A', 'A-G-A-G', 'A-G-A-G-A', 'A-G-A-G-A-G']) {
    edits.push({ source, ...await measure('#cabiln-input', 'input', source,
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
    connections.push({ source, ...await measure('#build-connect', 'click', source,
      () => page.locator('#build-connect').click()) });
  }
  await page.locator('#build-close').click();
  await page.locator('#lib-close').click();
  const largerPeptides = [];
  for (const sample of [
    { name: '60 residues with terminal caps', residues: 62,
      source: 'ac-' + Array(12).fill('A-G-K-S-D').join('-') + '-am' },
    { name: 'six nested branches', residues: 20,
      source: 'ac-' + Array(6).fill('K.{G(4,2)[.A(1,2)]}').join('-') + '-am' },
  ]) {
    largerPeptides.push({ ...sample, ...await measure('#cabiln-input', 'input', sample.source,
      () => page.locator('#cabiln-input').fill(sample.source), sample.residues) });
  }
  const result = {
    application: process.env.CABILN_APP_ROOT || 'current checkout',
    browser: process.env.CABILN_BROWSER_CHANNEL || 'chromium',
    viewport: testInfo.project.use.viewport,
    edits, connections, largerPeptides,
    note: 'Event to visible feedback and to completed drawing/chips; both sampled after two animation frames. Includes network and chemistry in completedMs. Excludes automation polling. First edit is cold only when run alone. Local samples, not field INP or production percentiles.',
  };
  const filename = testInfo.outputPath('browser-timings.json');
  fs.writeFileSync(filename, JSON.stringify(result, null, 2));
  await testInfo.attach('browser-timings', { path: filename, contentType: 'application/json' });
});
