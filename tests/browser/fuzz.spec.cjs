const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const { execFileSync } = require('node:child_process');
const fc = require('fast-check');
const { test, expect, render, tile, site, selectChip, connect, isCompletedResponse } = require('./fixtures');

const profile = process.env.CABILN_FUZZ_PROFILE || 'ci';
assert.ok(['ci', 'deep'].includes(profile), 'CABILN_FUZZ_PROFILE must be ci or deep');
const runs = Number(process.env.CABILN_FUZZ_RUNS || (profile === 'deep' ? 40 : 8));
assert.ok(Number.isSafeInteger(runs) && runs > 0, 'CABILN_FUZZ_RUNS must be a positive integer');
const seed = !process.env.CABILN_FUZZ_SEED ? undefined : Number(process.env.CABILN_FUZZ_SEED);
assert.ok(seed === undefined || (/^-?\d+$/.test(process.env.CABILN_FUZZ_SEED) && Number.isSafeInteger(seed) && seed >= -2147483648 && seed <= 2147483647), 'CABILN_FUZZ_SEED must be a signed 32-bit integer');
assert.ok(!process.env.CABILN_FUZZ_PATH || /^\d+(?::\d+)*$/.test(process.env.CABILN_FUZZ_PATH), 'CABILN_FUZZ_PATH must contain colon-separated non-negative integers');
const artifacts = process.env.CABILN_FUZZ_ARTIFACTS || path.join(os.tmpdir(), 'cabiln-browser-fuzz');
fs.mkdirSync(artifacts, { recursive: true });
assert.ok(!process.env.CABILN_FUZZ_CASE || ['browser-branch-graphs', 'browser-editable-history', 'browser-delayed-ownership'].includes(process.env.CABILN_FUZZ_CASE), 'Unknown CABILN_FUZZ_CASE');
const observedApplications = new Map();
const root = path.resolve(__dirname, '../..');
const sourceFiles = ['builder.js', 'document.js', 'project.js', 'requests.js', 'register.js'].map(name => `src/pyPept/web/static/${name}`).concat(['tests/browser/fixtures.js', 'tests/browser/fuzz.spec.cjs', 'tests/browser/fuzz.config.js', 'tests/browser/package-lock.json']);
const source = { commit: execFileSync('git', ['rev-parse', 'HEAD'], { cwd: root, encoding: 'utf8' }).trim(),
  sha256: Object.fromEntries(sourceFiles.map(name => [name, crypto.createHash('sha256').update(fs.readFileSync(path.join(root, name))).digest('hex')])) };

async function freshPage(browser, app, name, body) {
  const context = await browser.newContext({ viewport: { width: 1440, height: 1000 } });
  await context.tracing.start({ screenshots: true, snapshots: true, sources: true });
  const page = await context.newPage();
  const errors = [];
  let failed = false;
  const assetHashes = {};
  const assets = [];
  page.on('response', response => {
    const asset = new URL(response.url()).pathname;
    if (/\/(builder|document|project|requests|register)\.js$/.test(asset)) {
      assets.push(response.body().then(body => { assetHashes[asset] = crypto.createHash('sha256').update(body).digest('hex'); }));
    }
  });
  page.on('pageerror', error => errors.push(error.message));
  try {
    await page.goto(app.url);
    await body(page);
    expect(errors).toEqual([]);
  } catch (error) {
    failed = true;
    // Every failing shrink replaces this capture; the final file is the last
    // failing concrete browser state, alongside fast-check's minimized case.
    await page.screenshot({ path: path.join(artifacts, `${name}-failure.png`) }).catch(() => {});
    fs.writeFileSync(path.join(artifacts, `${name}-failure.html`), await page.content().catch(() => ''));
    throw error;
  } finally {
    await Promise.all(assets);
    const binding = await page.evaluate(() => typeof projectSnapshot === 'function' ? projectSnapshot().context : null).catch(() => null);
    if (!observedApplications.has(name)) observedApplications.set(name, { appRoot: app.appRoot,
      installed: process.env.CABILN_BROWSER_INSTALLED === '1', python: process.env.CABILN_PYTHON || 'python3',
      execution: process.env.CABILN_EXECUTION || 'default', workerMemoryMb: process.env.CABILN_WORKER_MEMORY_MB || 'default',
      librarySha256: Object.fromEntries(['monomers.sdf', 'monomers.csv'].map(file => [file,
        crypto.createHash('sha256').update(fs.readFileSync(path.join(app.temporary, file))).digest('hex')])),
      servedAssetsSha256: {}, libraryBindings: [] });
    const observed = observedApplications.get(name);
    Object.assign(observed.servedAssetsSha256, assetHashes);
    if (binding && !observed.libraryBindings.some(item => JSON.stringify(item) === JSON.stringify(binding))) observed.libraryBindings.push(binding);
    await context.tracing.stop(failed ? { path: path.join(artifacts, `${name}-failure.trace.zip`) } : {});
    await context.close();
  }
}

async function campaign(name, arbitrary, body, browser) {
  const counters = {};
  const uniqueHistories = new Set();
  const count = key => { counters[key] = (counters[key] || 0) + 1; };
  const started = Date.now();
  const result = await fc.check(fc.asyncProperty(arbitrary, async value => {
    count('historiesIncludingShrinks');
    try { await body(value, count); }
    finally { uniqueHistories.add(fc.stringify(value)); }
  }), { numRuns: runs, ...(seed === undefined ? {} : { seed }),
    ...(process.env.CABILN_FUZZ_PATH ? { path: process.env.CABILN_FUZZ_PATH } : {}), verbose: 2,
    endOnFailure: !!process.env.CABILN_FUZZ_PATH,
    plugins: [fc.interruptAfterTimeLimit(profile === 'deep' ? 840_000 : 150_000, { failOnInterrupt: true })] });
  const report = { property: name, profile, seed: result.seed, path: result.counterexamplePath,
    replayPath: result.counterexample?.map(value => fc.stringify(value)).join('\n').match(/replayPath="([^"]*)"/)?.[1] || null,
    counterexample: result.counterexample?.map(value => fc.stringify(value)), failed: result.failed,
    runs: result.numRuns, shrinks: result.numShrinks, durationMs: Date.now() - started, counters, uniqueHistoriesIncludingShrinks: uniqueHistories.size,
    source, application: observedApplications.get(name), versions: { node: process.version, fastCheck: fc.__version, playwright: require('@playwright/test/package.json').version,
      browser: browser.version(), platform: process.platform },
    error: result.errorInstance?.stack || null };
  fs.writeFileSync(path.join(artifacts, `${name}.json`), JSON.stringify(report, null, 2) + '\n');
  expect(result.failed, JSON.stringify(report, null, 2)).toBe(false);
}

const amino = fc.constantFrom('A', 'G');
const graph = fc.record({ prefix: fc.array(amino, { maxLength: 3 }), suffix: fc.array(amino, { maxLength: 3 }),
  arm: fc.array(amino, { minLength: 1, maxLength: 4 }), nested: fc.boolean(), protected: fc.boolean(), capped: fc.boolean() });

test('generated branch graphs keep exact chip and highlighted atom ownership', async ({ browser, app }) => {
  const name = 'browser-branch-graphs';
  test.skip(process.env.CABILN_FUZZ_CASE && process.env.CABILN_FUZZ_CASE !== name);
  await campaign(name, graph, async (example, count) => freshPage(browser, app, name, async page => {
    const prefix = [...(example.capped ? ['ac'] : []), ...example.prefix];
    const suffix = [...example.suffix, ...(example.capped ? ['am'] : [])];
    const tail = example.arm.slice(1).map(symbol => `.${symbol}(1,2)`).join('');
    const opening = example.protected ? '{' : '[';
    const closing = example.protected ? '}' : ']';
    const arm = `${example.arm[0]}(4,2)` + (example.nested && tail ? `[${tail}]` : tail);
    const source = [...prefix, `K.${opening}${arm}${closing}`, ...suffix].join('-');
    const data = await render(page, source);
    count('renderGraph'); if (example.nested && tail) count('nestedBranch');
    if (example.protected) count('protectedBranch');
    const mainLength = prefix.length + 1 + suffix.length;
    const armIds = example.arm.map((_, index) => mainLength + index);
    const expectedGroups = [armIds, ...(example.nested && tail ? [armIds.slice(1)] : [])];
    expect(data.layout.groups.map(group => group.members)).toEqual(expectedGroups);
    await expect(page.locator('#residue-chips [data-residue]')).toHaveText([...prefix, 'K', ...example.arm, ...suffix]);
    const actualIds = await page.locator('#residue-chips [data-residue]').evaluateAll(elements => elements.map(el => Number(el.dataset.residue)).sort((a, b) => a - b));
    expect(actualIds).toEqual(Array.from({ length: mainLength + example.arm.length }, (_, index) => index));
    const tabs = page.locator('#residue-chips .branch-chip');
    for (let index = 0; index < await tabs.count(); index++) {
      const tab = tabs.nth(index);
      const members = JSON.parse(await tab.getAttribute('data-members'));
      await tab.hover(); count('highlightGroup');
      const highlighted = await page.locator('#residue-chips [data-residue].hover').evaluateAll(elements => elements.map(el => Number(el.dataset.residue)).sort((a, b) => a - b));
      expect(highlighted).toEqual([...members].sort((a, b) => a - b));
      const actual = await page.locator('#render-inner svg').evaluate((svg, { members, ownership }) => {
        const atoms = new Set(members.flatMap(id => ownership[String(id)]));
        const paths = [...svg.querySelectorAll('[class]')].filter(el => /(?:^| )atom-\d+(?: |$)/.test(el.getAttribute('class')));
        return { count: svg.querySelectorAll('.res-hl').length, correct: paths.every(el => {
          const ids = [...el.getAttribute('class').matchAll(/(?:^| )atom-(\d+)(?= |$)/g)].map(match => Number(match[1]));
          return el.classList.contains('res-hl') === ids.some(id => atoms.has(id));
        }) };
      }, { members, ownership: data.residue_map });
      expect(actual.count).toBeGreaterThan(0); expect(actual.correct).toBe(true);
    }
    await page.locator('#cabiln-input').hover();
    await expect(page.locator('#render-inner svg')).not.toHaveClass(/has-highlight/);
    await expect(page.locator('#btn-mol')).toBeEnabled();
  }), browser);
});

function remember(model, source) {
  if (source === model.source) return;
  model.undo.push(model.source); model.redo = []; model.source = source;
}
async function checkHistory(model, real) {
  await expect(real.page.locator('#cabiln-input')).toHaveValue(model.source);
  if (model.undo.length) await expect(real.page.locator('#btn-undo')).toBeEnabled();
  else await expect(real.page.locator('#btn-undo')).toBeDisabled();
  if (model.redo.length) await expect(real.page.locator('#btn-redo')).toBeEnabled();
  else await expect(real.page.locator('#btn-redo')).toBeDisabled();
  if (model.source) {
    await expect(real.page.locator('#cabiln-input')).toHaveClass('ok');
    await expect(real.page.locator('#btn-mol')).toBeEnabled();
  } else await expect(real.page.locator('#btn-mol')).toBeDisabled();
}
class BrowserEdit {
  constructor(symbols) { this.source = symbols.join('-'); }
  check() { return true; }
  async run(model, real) {
    real.count('edit');
    await real.page.clock.fastForward(1000);
    const response = real.page.waitForResponse(r => isCompletedResponse(r) && new URL(r.url()).pathname === '/render');
    await real.page.locator('#cabiln-input').fill(this.source);
    await real.page.clock.runFor(200);
    expect((await response).status()).toBe(200);
    remember(model, this.source); await checkHistory(model, real);
  }
  toString() { return `Edit(${JSON.stringify(this.source)})`; }
}
class BrowserTravel {
  constructor(direction) { this.direction = direction; }
  check(model) { return model[this.direction].length > 0; }
  async run(model, real) {
    real.count(this.direction);
    model[this.direction === 'undo' ? 'redo' : 'undo'].push(model.source);
    model.source = model[this.direction].pop();
    await real.page.locator(`#btn-${this.direction}`).click();
    await checkHistory(model, real);
  }
  toString() { return this.direction; }
}
class BrowserClear {
  check(model) { return !!model.source; }
  async run(model, real) {
    real.count('clear');
    await real.page.clock.fastForward(1000);
    await real.page.locator('#cabiln-input').fill('');
    remember(model, ''); await checkHistory(model, real);
  }
  toString() { return 'Clear'; }
}
class BrowserConnect {
  constructor(symbol) { this.symbol = symbol; }
  check(model) { return !!model.source; }
  async run(model, real) {
    const page = real.page;
    real.count('connect');
    if (!await page.locator('#build-panel').evaluate(el => el.classList.contains('open'))) await page.locator('#btn-build').click();
    await page.locator('#build-left-change').click();
    const symbols = model.source.split('-');
    await selectChip(page, symbols.length - 1, 'left', symbols.at(-1));
    await tile(page, this.symbol, 'right');
    await site(page, 'left', 2); await site(page, 'right', 1);
    const expected = `${model.source}-${this.symbol}`;
    const connected = await connect(page);
    expect(connected.data.result).toBe(expected);
    remember(model, expected); await checkHistory(model, real);
  }
  toString() { return `Connect(${this.symbol})`; }
}
class BrowserPanels {
  check() { return true; }
  async run(model, real) {
    real.count('panelLifetime');
    await real.page.locator('#btn-build').click();
    await real.page.locator('#btn-examples').click();
    await checkHistory(model, real);
  }
  toString() { return 'TogglePanels'; }
}
class BrowserSave {
  check(model) { return !!model.source; }
  async run(model, real) {
    real.count('saveProject');
    const downloaded = real.page.waitForEvent('download');
    await real.page.locator('#btn-project-save').click();
    const download = await downloaded;
    real.saved = fs.readFileSync(await download.path());
    expect(JSON.parse(real.saved.toString()).document.text).toBe(model.source);
    model.saved = model.source;
    await checkHistory(model, real);
  }
  toString() { return 'SaveProject'; }
}
class BrowserOpen {
  check(model) { return !!model.saved && model.saved !== model.source; }
  async run(model, real) {
    real.count('openProject');
    await real.page.locator('#project-upload').setInputFiles({ name: 'generated.cabiln.json', mimeType: 'application/json', buffer: real.saved });
    remember(model, model.saved);
    await checkHistory(model, real);
  }
  toString() { return 'OpenSavedProject'; }
}
class BrowserReload {
  check(model) { return !!model.source; }
  async run(model, real) {
    real.count('reload');
    await real.page.reload();
    await real.page.locator('#btn-restore-draft').click();
    model.undo = ['']; model.redo = [];
    await checkHistory(model, real);
  }
  toString() { return 'ReloadAndRestore'; }
}
const browserCommands = fc.commands([
  fc.array(amino, { minLength: 1, maxLength: 4 }).map(symbols => new BrowserEdit(symbols)),
  fc.constant(new BrowserTravel('undo')), fc.constant(new BrowserTravel('redo')),
  fc.constant(new BrowserClear()), amino.map(symbol => new BrowserConnect(symbol)), fc.constant(new BrowserPanels()),
  fc.constant(new BrowserSave()), fc.constant(new BrowserOpen()), fc.constant(new BrowserReload()),
], { maxCommands: profile === 'deep' ? 18 : 10, size: 'large',
  ...(process.env.CABILN_FUZZ_REPLAY_PATH ? { replayPath: process.env.CABILN_FUZZ_REPLAY_PATH } : {}) });

test('generated browser edit and build histories preserve the independent timeline', async ({ browser, app }) => {
  const name = 'browser-editable-history';
  test.skip(process.env.CABILN_FUZZ_CASE && process.env.CABILN_FUZZ_CASE !== name);
  await campaign(name, browserCommands, async (history, count) => freshPage(browser, app, name, async page => {
    await page.clock.install({ time: new Date('2026-09-30T12:00:00Z') });
    await fc.asyncModelRun(() => ({ model: { source: '', undo: [], redo: [] }, real: { page, count } }), history);
  }), browser);
});

test('generated delayed completions preserve immediate conversion and current source', async ({ browser, app }) => {
  const name = 'browser-delayed-ownership';
  test.skip(process.env.CABILN_FUZZ_CASE && process.env.CABILN_FUZZ_CASE !== name);
  await campaign(name, fc.record({ source: fc.constantFrom('NCC(=O)O', 'N[C@@H](C)C(=O)O'),
    editDuring: fc.boolean(), cancelPhase: fc.constantFrom('drawing', 'conversion'), useMol: fc.boolean() }), async (example, count) => freshPage(browser, app, name, async page => {
    let release;
    let observed;
    const arrived = new Promise(resolve => { observed = resolve; });
    const gate = new Promise(resolve => { release = resolve; });
    const endpoint = example.cancelPhase === 'drawing' ? '/render_reference' : '/to_cabiln';
    let held = false;
    let originalMol;
    await page.route(`**${endpoint}`, async route => {
      if (held) { await route.continue(); return; }
      held = true;
      // Obtain the real server result first, then hold delivery. This exercises
      // an already queued response even when the client cancels its request.
      const response = await route.fetch();
      observed(); await gate;
      await route.fulfill({ response }).catch(() => {});
    });
    if (example.useMol) {
      const molecule = await render(page, 'A');
      originalMol = molecule.mol_block;
      await page.locator('#btn-verify').click();
      await page.locator('#mol-upload').setInputFiles({ name: 'original.mol', mimeType: 'chemical/x-mdl-molfile', buffer: Buffer.from(molecule.mol_block) });
      await expect(page.locator('#smiles-input')).toHaveClass('ok');
      count('originalMol');
    }
    await page.locator('#notation-select').selectOption('smiles');
    await page.locator('#cabiln-input').fill(example.source);
    if (example.cancelPhase === 'drawing') await arrived;
    await page.locator('#btn-to-cabiln-bracket').click();
    if (example.cancelPhase === 'conversion') await arrived;
    if (example.editDuring) {
      await page.locator('#cabiln-input').fill('CCO');
      count('cancelByEdit');
    } else count('immediateConversion');
    release();
    if (example.editDuring) {
      await expect(page.locator('#notation-select')).toHaveValue('smiles');
      await expect(page.locator('#cabiln-input')).toHaveValue('CCO');
      await expect(page.locator('#cabiln-input')).toHaveClass('ok');
      await expect(page.locator('#conversion-progress')).toBeHidden();
    } else {
      await expect(page.locator('#notation-select')).toHaveValue('cabiln');
      await expect(page.locator('#cabiln-input')).toHaveValue(example.source === 'NCC(=O)O' ? 'G' : 'A');
      await expect(page.locator('#btn-mol')).toBeEnabled();
      await expect(page.locator('#conversion-status')).toBeHidden();
    }
    if (example.useMol) {
      const reference = await page.evaluate(() => JSON.parse(window.localStorage.getItem('cabiln.draft.v1')).reference_original);
      expect(reference.kind).toBe('mol'); expect(reference.name).toBe('original.mol');
      expect(reference.content).toBe(originalMol);
    }
  }), browser);
});
