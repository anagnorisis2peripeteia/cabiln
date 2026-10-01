// npm --prefix tests/browser run fuzz:frontend
// Replay one property with CABILN_FUZZ_CASE, CABILN_FUZZ_SEED and CABILN_FUZZ_PATH.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const crypto = require('node:crypto');
const { execFileSync } = require('node:child_process');
const test = require('node:test');
const fc = require('./browser/node_modules/fast-check');
const { page } = require('./_frontend_harness');

const profile = process.env.CABILN_FUZZ_PROFILE || 'ci';
assert.ok(['ci', 'deep'].includes(profile), 'CABILN_FUZZ_PROFILE must be ci or deep');
const runs = Number(process.env.CABILN_FUZZ_RUNS || (profile === 'deep' ? 1000 : 100));
assert.ok(Number.isSafeInteger(runs) && runs > 0, 'CABILN_FUZZ_RUNS must be a positive integer');
const seed = !process.env.CABILN_FUZZ_SEED ? undefined : Number(process.env.CABILN_FUZZ_SEED);
assert.ok(seed === undefined || (/^-?\d+$/.test(process.env.CABILN_FUZZ_SEED) && Number.isSafeInteger(seed) && seed >= -2147483648 && seed <= 2147483647), 'CABILN_FUZZ_SEED must be a signed 32-bit integer');
assert.ok(!process.env.CABILN_FUZZ_PATH || /^\d+(?::\d+)*$/.test(process.env.CABILN_FUZZ_PATH), 'CABILN_FUZZ_PATH must contain colon-separated non-negative integers');
const completionOrder = process.env.CABILN_FUZZ_SCHEDULE ? JSON.parse(process.env.CABILN_FUZZ_SCHEDULE) : null;
assert.ok(completionOrder === null || (Array.isArray(completionOrder) && completionOrder.length > 0 &&
  completionOrder.every(id => Number.isSafeInteger(id) && id > 0) && new Set(completionOrder).size === completionOrder.length), 'CABILN_FUZZ_SCHEDULE must be a JSON array of distinct positive task IDs');
const artifacts = process.env.CABILN_FUZZ_ARTIFACTS || path.join(os.tmpdir(), 'cabiln-frontend-fuzz');
const root = path.resolve(__dirname, '..');
assert.ok(!process.env.CABILN_FUZZ_CASE || ['editable-history', 'render-completion-order', 'foreground-admission-and-errors', 'registration-preview-and-write-ownership', 'normalization-keeps-redo-and-reference', 'selection-and-library-revisions'].includes(process.env.CABILN_FUZZ_CASE), 'Unknown CABILN_FUZZ_CASE');
const sourceFiles = ['builder.js', 'document.js', 'project.js', 'requests.js', 'register.js'].map(name => `src/pyPept/web/static/${name}`).concat(['tests/_frontend_harness.js', 'tests/test_frontend_fuzz.js', 'tests/browser/package-lock.json']);
const source = { commit: execFileSync('git', ['rev-parse', 'HEAD'], { cwd: root, encoding: 'utf8' }).trim(),
  sha256: Object.fromEntries(sourceFiles.map(name => [name, crypto.createHash('sha256').update(fs.readFileSync(path.join(root, name))).digest('hex')])) };
const tick = () => new Promise(setImmediate);
const binding = { project_version: 1, library_binding: { monomers: 'test-library', aliases: 'a', reactions: 'r', caps: 'c' }, canonical: { format: 'test-v1', rdkit: 'test' } };
const sequences = fc.constantFrom('A', 'G', 'K', 'A-G', 'K-G-A', 'K.[G(4,2)]', 'K.!1(4,2)%G.!1(2,4)', 'D.(4,1)-G%G-A', '');
const referenceText = fc.constantFrom('NCC(=O)O', 'CCO', 'O', '');
const finished = new WeakSet();

function drawing(source) {
  return { svg: `<svg>${source}</svg>`, mol_block: `MOL:${source}`, info: source,
    context: binding, warnings: [], residues: [], residue_map: {}, layout: { segments: [], groups: [], markers: [] } };
}
function pending(ui, url) { return ui.requests.filter(r => r.url === url && !finished.has(r)); }
function resolve(request, data, ok = true, status, headers) {
  assert.ok(request, 'Expected an actual handler request');
  finished.add(request);
  request.resolve(data, ok, status, headers);
}
function reject(request) { finished.add(request); request.reject(new Error('controlled network failure')); }
function snapshot(ui) { return JSON.parse(ui.run('JSON.stringify(projectSnapshot())')); }
async function settle(ui) {
  // Complete only background calculations. Foreground operations retain explicit
  // responses in their command so a missing or unexpected endpoint cannot pass.
  await ui.timers();
  for (let turn = 0; turn < 12; turn++) {
    await tick();
    const requests = ui.requests.filter(r => !finished.has(r) && ['/render', '/render_reference', '/render_mol', '/verify', '/monomers', '/reactions'].includes(r.url));
    if (!requests.length) return;
    for (const request of requests) {
      const body = request.options.body ? JSON.parse(request.options.body) : {};
      if (request.url === '/monomers' || request.url === '/reactions') resolve(request, []);
      else if (request.url === '/verify') resolve(request, { match: true });
      else resolve(request, { ...drawing(body.cabiln || body.input || 'MOL'), format: 'SMILES', smiles: body.input || 'NCC(=O)O' });
    }
  }
  assert.fail('Background calculations did not settle');
}

function property(name, arbitraries, body) {
  test(name, { timeout: profile === 'deep' ? 125_000 : 35_000,
    skip: !!process.env.CABILN_FUZZ_CASE && process.env.CABILN_FUZZ_CASE !== name }, async () => {
    const counters = {};
    const count = key => { counters[key] = (counters[key] || 0) + 1; };
    const started = Date.now();
    let failureSchedule = null;
    const result = await fc.check(fc.asyncProperty(...arbitraries, async (...args) => {
      count('historiesIncludingShrinks');
      if (completionOrder && typeof args[0]?.schedule === 'function') args[0] = fc.schedulerFor(completionOrder);
      try { await body(...args, count); }
      catch (error) {
        if (typeof args[0]?.report === 'function') failureSchedule = {
          order: args[0].report().map(task => task.taskId), replay: fc.stringify(args[0]),
        };
        throw error;
      }
    }), { numRuns: runs, ...(seed === undefined ? {} : { seed }),
      ...(process.env.CABILN_FUZZ_PATH ? { path: process.env.CABILN_FUZZ_PATH } : {}), verbose: 2,
      endOnFailure: !!process.env.CABILN_FUZZ_PATH,
      plugins: [fc.interruptAfterTimeLimit(profile === 'deep' ? 120_000 : 30_000, { failOnInterrupt: true })] });
    fs.mkdirSync(artifacts, { recursive: true });
    const report = { property: name, profile, seed: result.seed, path: result.counterexamplePath,
      replayPath: result.counterexample?.map(value => fc.stringify(value)).join('\n').match(/replayPath="([^"]*)"/)?.[1] || null,
      counterexample: result.counterexample?.map((value, index) => index === 0 && failureSchedule ? failureSchedule.replay : fc.stringify(value)),
      completionOrder: failureSchedule?.order || null,
      failed: result.failed, runs: result.numRuns, shrinks: result.numShrinks,
      durationMs: Date.now() - started, counters,
      source, versions: { node: process.version, fastCheck: fc.__version,
        playwright: require('./browser/node_modules/@playwright/test/package.json').version, platform: process.platform },
      error: result.errorInstance?.stack || null };
    fs.writeFileSync(path.join(artifacts, `${name}.json`), JSON.stringify(report, null, 2) + '\n');
    assert.equal(result.failed, false, JSON.stringify(report, null, 2));
  });
}

// A sequence timeline is the independent oracle. It knows no request IDs,
// document-class internals, rendering implementation, or private undo objects.
function commit(model, source) {
  if (source === model.source) return;
  model.undo.push(model.source);
  model.redo = [];
  model.source = source;
  model.storageCleared = false;
}
function checkDocument(model, real) {
  const ui = real.ui;
  assert.equal(ui.element('cabiln-input').value, model.source);
  assert.equal(ui.element('notation-select').value, 'cabiln');
  assert.equal(ui.element('btn-undo').disabled, model.undo.length === 0);
  assert.equal(ui.element('btn-redo').disabled, model.redo.length === 0);
  assert.equal(ui.element('btn-mol').disabled, model.source === '');
  assert.equal(ui.element('btn-png').disabled, model.source === '');
  if (model.source) assert.equal(ui.run('lastMolBlock'), `MOL:${model.source}`);
  const current = snapshot(ui);
  assert.equal(current.document.text, model.source);
  assert.equal(current.document.warning, '');
  assert.equal(current.reference.text, model.reference);
  assert.deepEqual(current.reference.original, model.original);
}

class Edit {
  constructor(values) { this.values = values; }
  check() { return true; }
  async run(model, real) {
    real.count('edit'); real.clock += 1000;
    let changed = false;
    const previous = model.source;
    for (const value of this.values) {
      real.clock += 50;
      if (value !== model.source) changed = true;
      model.source = value;
      await real.ui.input('cabiln-input', value);
    }
    if (changed) { model.undo.push(previous); model.redo = []; model.storageCleared = false; }
    await settle(real.ui); checkDocument(model, real);
  }
  toString() { return `Edit(${JSON.stringify(this.values)})`; }
}
class Travel {
  constructor(direction) { this.direction = direction; }
  check(model) { return model[this.direction].length > 0; }
  async run(model, real) {
    real.count(this.direction);
    model[this.direction === 'undo' ? 'redo' : 'undo'].push(model.source);
    model.source = model[this.direction].pop();
    model.storageCleared = false;
    await real.ui.element(`btn-${this.direction}`).click();
    await settle(real.ui); checkDocument(model, real);
  }
  toString() { return this.direction; }
}
class Format {
  check(model) { return !!model.source; }
  async run(model, real) {
    real.count('convert');
    const sequence = model.source === 'K.[G(4,2)]' ? 'K.!1(4,2)%G.!1(2,4)' : model.source === 'K.!1(4,2)%G.!1(2,4)' ? 'K.[G(4,2)]' : model.source;
    const task = real.ui.element('btn-to-bracket').click();
    await tick();
    resolve(pending(real.ui, '/convert_notation')[0], { result: sequence });
    await task;
    commit(model, sequence);
    await settle(real.ui); checkDocument(model, real);
  }
  toString() { return 'Format'; }
}
class Normalize {
  check(model) { return model.source === 'D.(4,1)-G%G-A'; }
  async run(model, real) {
    real.count('normalize');
    const normalized = 'D.[G(4,1).A(2,1)]-G';
    await real.ui.element('btn-reroll').click();
    resolve(pending(real.ui, '/render')[0], { ...drawing(normalized), normalized_cabiln: normalized });
    await tick();
    model.source = normalized;
    model.storageCleared = false;
    checkDocument(model, real);
  }
  toString() { return 'Normalize'; }
}
class Connect {
  check(model) { return /^[AGK](?:-[AGK])*$/.test(model.source); }
  async run(model, real) {
    real.count('connect');
    // Selection geometry is generated against real DOM in fuzz.spec.cjs. Here
    // the actual Connect handler owns the asynchronous document commit.
    real.ui.run(`buildLeftRIdx=${model.source.split('-').length - 1}; buildRightRIdx=null; buildLeft={abbr:'K',selectedSlot:2}; buildRight={abbr:'G',selectedSlot:1}`);
    const task = real.ui.element('build-connect').click();
    const sequence = `${model.source}-G`;
    resolve(pending(real.ui, '/insert_bond')[0], { result: sequence });
    await task;
    commit(model, sequence);
    await settle(real.ui); checkDocument(model, real);
  }
  toString() { return 'ConnectBackboneG'; }
}
class Reference {
  constructor(text, mol) { this.text = text; this.mol = mol; }
  check() { return true; }
  async run(model, real) {
    real.count(this.mol ? 'referenceMol' : 'referenceText');
    model.storageCleared = false;
    if (this.mol) {
      const content = `ORIGINAL\r\n${this.text}\nM  END\n`;
      const task = real.ui.element('mol-upload').dispatchEvent({ type: 'change', target: { files: [{ name: 'fuzz.mol', text: async () => content }] } });
      await tick();
      resolve(pending(real.ui, '/render_mol')[0], { ...drawing('reference'), smiles: 'NCC(=O)O' });
      await task;
      model.reference = 'NCC(=O)O';
      model.original = { kind: 'mol', content, name: 'fuzz.mol' };
    } else {
      await real.ui.input('smiles-input', this.text);
      model.reference = this.text;
      model.original = this.text ? { kind: 'text', content: this.text, name: '' } : null;
    }
    await settle(real.ui); checkDocument(model, real);
  }
  toString() { return `Reference(${JSON.stringify(this.text)},mol=${this.mol})`; }
}
class Save {
  check() { return true; }
  async run(model, real) {
    real.count('save');
    const task = real.ui.element('btn-project-save').click();
    const request = pending(real.ui, '/prepare_project')[0];
    const project = JSON.parse(request.options.body).project;
    assert.equal(project.document.text, model.source);
    assert.equal(project.reference.text, model.reference);
    assert.deepEqual(project.reference.original, model.original);
    project.context = binding;
    resolve(request, { project });
    await task;
    real.saved = JSON.parse(await real.ui.downloads.at(-1).text());
    model.saved = { source: model.source, reference: model.reference, original: model.original };
    checkDocument(model, real);
  }
  toString() { return 'SaveProject'; }
}
class Open {
  check(model) { return !!model.saved && model.saved.source !== model.source; }
  async run(model, real) {
    real.count('open');
    const text = JSON.stringify(real.saved);
    const task = real.ui.element('project-upload').dispatchEvent({ type: 'change', target: { files: [{ size: text.length, text: async () => text }] } });
    await tick();
    resolve(pending(real.ui, '/validate_project')[0], { valid: true, context: binding });
    await task;
    commit(model, model.saved.source);
    Object.assign(model, { reference: model.saved.reference, original: model.saved.original });
    await settle(real.ui); checkDocument(model, real);
  }
  toString() { return 'OpenSavedProject'; }
}
class Reload {
  check(model) { return !!(model.source || model.reference) && !model.storageCleared; }
  async run(model, real) {
    real.count('reload');
    await real.ui.run('window.dispatchEvent(new Event("pagehide"))');
    const storedDraft = JSON.parse(real.ui.storage.get('cabiln.draft.v1'));
    assert.equal(storedDraft.document.text, model.source);
    assert.deepEqual(storedDraft.reference_original, model.original);
    real.ui = page('builder.js', { storedDraft, now: () => real.clock });
    const task = real.ui.element('btn-restore-draft').click();
    const request = pending(real.ui, '/prepare_project')[0];
    if (request) resolve(request, { project: JSON.parse(request.options.body).project });
    await task;
    model.undo = storedDraft.document.text || storedDraft.document.context ? [''] : [];
    model.redo = [];
    await settle(real.ui); checkDocument(model, real);
  }
  toString() { return 'ReloadAndRestore'; }
}
class ClearStorage {
  check() { return true; }
  async run(model, real) {
    real.count('clearStorage');
    await real.ui.element('btn-clear-draft').click();
    await real.ui.run('window.dispatchEvent(new Event("pagehide"))');
    assert.equal(real.ui.storage.has('cabiln.draft.v1'), false);
    model.storageCleared = true;
    checkDocument(model, real);
  }
  toString() { return 'ClearBrowserStorage'; }
}

const commands = fc.commands([
  fc.array(sequences, { minLength: 1, maxLength: 3 }).map(values => new Edit(values)),
  fc.constant(new Travel('undo')), fc.constant(new Travel('redo')), fc.constant(new Format()),
  fc.constant(new Normalize()), fc.constant(new Connect()),
  fc.tuple(referenceText, fc.boolean()).map(([text, mol]) => new Reference(text, mol)),
  fc.constant(new Save()), fc.constant(new Open()), fc.constant(new Reload()), fc.constant(new ClearStorage()),
], { maxCommands: profile === 'deep' ? 50 : 25, size: 'large',
  ...(process.env.CABILN_FUZZ_REPLAY_PATH ? { replayPath: process.env.CABILN_FUZZ_REPLAY_PATH } : {}) });

property('editable-history', [commands], async (history, count) => {
  const real = { clock: 10000, count };
  real.ui = page('builder.js', { now: () => real.clock });
  await fc.asyncModelRun(() => ({ model: { source: '', undo: [], redo: [], reference: '', original: null }, real }), history);
});

property('render-completion-order', [fc.scheduler(), fc.array(fc.record({ source: sequences.filter(Boolean), reference: referenceText,
  outcome: fc.constantFrom('success', 'network', 'invalid'), referenceOutcome: fc.constantFrom('success', 'network', 'invalid'),
  releaseOne: fc.boolean() }), { minLength: 2, maxLength: 8 })], async (scheduler, edits, count) => {
  const ui = page('builder.js');
  const tasks = [];
  const scheduled = new WeakSet();
  let revision = -1;
  const reference = { text: '', svg: '', status: '', className: '', lastSmiles: '', context: null };
  function checkReference() {
    assert.equal(ui.element('smiles-input').value, reference.text);
    assert.equal(ui.element('smiles-input').className, reference.className);
    assert.equal(ui.element('smiles-inner').innerHTML, reference.svg);
    assert.equal(ui.element('smiles-status').textContent, reference.status);
    assert.equal(ui.run('lastSmiles'), reference.lastSmiles);
    assert.deepEqual(snapshot(ui).reference.context, reference.context);
  }
  for (const [index, edit] of edits.entries()) {
    revision = index;
    count('editMainAndReference');
    await ui.input('cabiln-input', edit.source);
    await ui.input('smiles-input', edit.reference);
    Object.assign(reference, { text: edit.reference, status: edit.reference ? 'Updating reference…' : '', className: '', lastSmiles: '' });
    if (!edit.reference) {
      count('emptyReference');
      reference.svg = '<div class="placeholder">Paste SMILES, BILN, or HELM…</div>';
      reference.context = null;
    }
    await ui.timers();
    for (const request of ui.requests.filter(r => !scheduled.has(r))) {
      scheduled.add(request);
      const body = JSON.parse(request.options.body);
      const isReference = request.url === '/render_reference';
      const outcome = isReference ? edit.referenceOutcome : edit.outcome;
      const source = body.cabiln || body.input;
      const error = `controlled input error ${index}`;
      const context = { ...binding, library_binding: { ...binding.library_binding, monomers: `reference:${index}` } };
      tasks.push(scheduler.schedule(Promise.resolve(), `${request.url}:${index}:${source}:${outcome}`).then(async () => {
        count(`${isReference ? 'reference' : 'main'}.${outcome}`);
        // This oracle admits only a completion of the latest independently
        // edited reference revision. It does not inspect the app's request map.
        if (isReference && revision === index) {
          if (outcome === 'network') reference.status = 'Could not reach the renderer. Your reference input is preserved.';
          else if (outcome === 'invalid') {
            reference.status = error; reference.className = 'err';
            reference.svg = `<div class="placeholder err">${error}</div>`;
          } else Object.assign(reference, { svg: `<svg>${source}</svg>`, lastSmiles: source,
            status: `SMILES: ${source}`, className: 'ok', context });
        }
        if (outcome === 'network') reject(request);
        else if (outcome === 'invalid') resolve(request, { error }, false, 400);
        else resolve(request, { ...drawing(source), context: isReference ? context : undefined, smiles: body.input, format: 'SMILES' });
        await tick();
      }));
    }
    if (edit.releaseOne && scheduler.count()) {
      await scheduler.waitNext(1); await tick(); checkReference();
    }
  }
  await scheduler.waitAll(); await Promise.all(tasks);
  checkReference();
  const last = edits.at(-1);
  assert.equal(ui.element('cabiln-input').value, last.source);
  assert.equal(ui.element('btn-mol').disabled, last.outcome !== 'success');
  assert.equal(ui.element('cabiln-input').className, last.outcome === 'success' ? 'ok' : last.outcome === 'invalid' ? 'err' : '');
  if (last.outcome === 'success') {
    assert.equal(ui.run('lastMolBlock'), `MOL:${last.source}`);
    assert.equal(ui.element('render-inner').innerHTML, `<svg>${last.source}</svg>`);
  } else assert.equal(snapshot(ui).document.warning, '');
  assert.equal(ui.run('buildLeft'), null); assert.equal(ui.run('buildRight'), null);
});

const foregroundCase = fc.record({
  format: fc.boolean(), drawingStarted: fc.boolean(),
  presentation: fc.constantFrom('missing', 'valid', 'changed-binding'),
  drawingBusy: fc.integer({ min: 0, max: 2 }), conversionBusy: fc.integer({ min: 0, max: 2 }),
  outcome: fc.constantFrom('success', 'network', 'error'),
}).chain(example => fc.constantFrom('none', 'conversion', 'reference',
  ...(!example.format || example.drawingStarted ? ['drawing'] : []),
  ...(example.format && !example.drawingStarted && example.outcome !== 'success' ? ['fallback'] : []),
).map(cancel => ({ ...example, cancel })));

property('foreground-admission-and-errors', [foregroundCase], async (example, count) => {
  const ui = page('builder.js');
  const source = example.format ? 'K.[G(4,2)]' : 'NCC(=O)O';
  if (!example.format) {
    ui.element('notation-select').value = 'smiles';
    await ui.element('notation-select').dispatchEvent({ type: 'change' });
  }
  await ui.input('cabiln-input', source);
  if (example.drawingStarted) await ui.timers();
  const endpoint = example.format ? '/convert_notation' : '/to_cabiln';
  const renderEndpoint = example.format ? '/render' : '/render_reference';
  const button = example.format ? 'btn-to-bracket' : 'btn-to-cabiln-bracket';
  const conversion = ui.element(button).click();
  await tick();
  const waitsForDrawing = !example.format || example.drawingStarted;
  const residues = [{ idx: 0, abbr: 'K' }, { idx: 1, abbr: 'G' }];
  const residue_map = { 0: [0, 1], 1: [2, 3] };
  count(`${example.format ? 'format' : 'foreign'}.${example.drawingStarted ? 'activeDrawing' : 'pendingTimer'}`);
  count(waitsForDrawing ? 'waitForDrawing' : 'skipUnsentFormatDrawing');

  async function finishDrawing(url, text, blocksConversion = false) {
    let rendering = pending(ui, url)[0];
    assert.ok(rendering, 'The current source has an actual drawing request');
    assert.equal(JSON.parse(rendering.options.body).cabiln || JSON.parse(rendering.options.body).input, text);
    for (let i = 0; i < example.drawingBusy; i++) {
      count('drawingRetry');
      resolve(rendering, { error: 'busy' }, false, 503, { 'Retry-After': '1' });
      await tick();
      if (blocksConversion) assert.equal(pending(ui, endpoint).length, 0);
      await ui.timers(); await tick();
      rendering = pending(ui, url)[0];
    }
    resolve(rendering, { ...drawing(text), format: 'SMILES',
      ...(example.format ? { residues, residue_map } : {}) });
    await tick();
  }

  if (waitsForDrawing) {
    assert.equal(pending(ui, endpoint).length, 0, 'Foreground work waits for an active or required drawing');
  } else {
    await ui.timers();
    assert.equal(pending(ui, renderEndpoint).length, 0, 'Formatting consumes the unsent drawing timer');
    assert.equal(pending(ui, endpoint).length, 1, 'Unsent CABILN can be formatted directly');
  }
  if (example.cancel === 'drawing') {
    count('cancelDuringDrawing');
    await ui.input('cabiln-input', 'C');
    resolve(pending(ui, renderEndpoint)[0], drawing('obsolete'));
    await conversion;
    await settle(ui);
    assert.equal(ui.element('cabiln-input').value, 'C');
    assert.equal(pending(ui, endpoint).length, 0);
    assert.equal(ui.element('conversion-progress').hidden, true);
    return;
  }
  if (waitsForDrawing) await finishDrawing(renderEndpoint, source, true);
  let request = pending(ui, endpoint)[0];
  assert.ok(request);
  if (example.cancel === 'conversion' || example.cancel === 'reference') {
    const referenceEdit = example.cancel === 'reference';
    count(referenceEdit ? 'cancelByReferenceEdit' : 'cancelDuringConversion');
    await ui.input(referenceEdit ? 'smiles-input' : 'cabiln-input', referenceEdit ? 'CCO' : 'C');
    resolve(request, { cabiln: 'OLD', result: 'OLD' });
    await conversion; await settle(ui);
    const expected = referenceEdit ? source : 'C';
    assert.equal(ui.element('cabiln-input').value, expected);
    assert.equal(ui.element('cabiln-input').className, 'ok');
    assert.equal(ui.element('render-inner').innerHTML, `<svg>${expected}</svg>`);
    if (referenceEdit && !waitsForDrawing) count('cancelRestoresUnsentDrawing');
    assert.equal(snapshot(ui).document.warning, '');
    return;
  }
  for (let i = 0; i < example.conversionBusy; i++) {
    count('conversionRetry');
    resolve(request, { error: 'busy' }, false, 503, { 'Retry-After': '1' });
    await tick(); await ui.timers(); await tick();
    request = pending(ui, endpoint)[0];
  }
  const result = example.format ? 'K.!1(4,2)%G.!1(2,4)' : 'G';
  const success = { cabiln: result, result, from: 'SMILES', warnings: [], assignments: [] };
  if (example.format && example.presentation !== 'missing') Object.assign(success, {
    source_echo: source, occurrence_order: [0, 1],
    context: example.presentation === 'valid' ? binding : {
      ...binding, library_binding: { ...binding.library_binding, monomers: 'changed' } },
    presentation: { cabiln_echo: result, residues,
      layout: { segments: [{ roots: [0, 1], members: [0, 1] }], groups: [], markers: [] },
      crosslink_groups: [], warnings: [] },
  });
  count(example.outcome);
  if (example.outcome === 'network') reject(request);
  else resolve(request, example.outcome === 'error' ? { error: 'controlled conversion failure' } : success);
  await tick();
  if (!waitsForDrawing) {
    if (example.cancel === 'fallback') {
      count('cancelDuringFallbackDrawing');
      const fallback = pending(ui, '/render')[0];
      assert.equal(JSON.parse(fallback.options.body).cabiln, source);
      await ui.input('cabiln-input', 'C');
      resolve(fallback, drawing('obsolete'));
      await conversion; await settle(ui);
      assert.equal(ui.element('cabiln-input').value, 'C');
      assert.equal(ui.element('render-inner').innerHTML, '<svg>C</svg>');
      return;
    }
    count(example.outcome === 'success' ? 'drawOnlyFormattedResult' : 'drawAfterFormatFailure');
    await finishDrawing('/render', example.outcome === 'success' ? result : source);
  }
  await conversion;
  function checkRetainedDrawing() {
    count('retainedDrawing');
    assert.equal(pending(ui, '/render').length, 0, 'The existing drawing needs no new request');
    assert.equal(ui.element('render-inner').innerHTML, `<svg>${source}</svg>`);
    assert.equal(ui.run('lastMolBlock'), `MOL:${source}`);
    assert.equal(ui.run('lastCabiln'), result);
    assert.deepEqual(JSON.parse(ui.run('JSON.stringify(residueMap)')), residue_map);
  }
  if (example.format && example.presentation === 'valid' &&
      waitsForDrawing && example.outcome === 'success') checkRetainedDrawing();
  assert.equal(snapshot(ui).document.warning, '');
  if (example.outcome !== 'success') {
    assert.equal(ui.element('cabiln-input').value, source);
    assert.equal(ui.element('cabiln-input').className, 'ok');
    assert.match(ui.element('cabiln-status').textContent, /failure|failed|Could not/);
    const retry = ui.element(button).click();
    await tick();
    resolve(pending(ui, endpoint)[0], success);
    await retry;
    if (example.format && example.presentation === 'valid') checkRetainedDrawing();
    count('successfulRetry');
  }
  await settle(ui);
  assert.equal(ui.element('cabiln-input').value, result);
  assert.doesNotMatch(ui.element('cabiln-status').textContent, /failure|failed|busy/);
  await ui.run('window.dispatchEvent(new Event("pagehide"))');
  assert.equal(JSON.parse(ui.storage.get('cabiln.draft.v1')).document.warning, '');
  await ui.element('btn-undo').click(); await settle(ui);
  assert.equal(ui.element('cabiln-input').value, source);
  assert.equal(ui.element('notation-select').value, example.format ? 'cabiln' : 'smiles');
  await ui.element('btn-redo').click(); await settle(ui);
  assert.equal(ui.element('cabiln-input').value, result);
  assert.equal(ui.element('notation-select').value, 'cabiln');
});

property('registration-preview-and-write-ownership', [fc.scheduler(), fc.record({
  sources: fc.array(fc.constantFrom('NCC(=O)O', 'CC(=O)O', 'NCCC(=O)O'), { minLength: 2, maxLength: 6 }),
  cancelLast: fc.boolean(), failure: fc.boolean(), writeFailure: fc.boolean(), editDuringWrite: fc.boolean(),
})], async (scheduler, example, count) => {
  const ui = page('register.js');
  await ui.input('abbr-in', 'FuzzMonomer');
  await ui.input('name-in', 'Fuzz monomer');
  const tasks = [];
  const responses = [];
  for (const [index, source] of example.sources.entries()) {
    await ui.input('smiles-in', source);
    tasks.push(ui.element('btn-preview').click());
    const request = pending(ui, '/preview_monomer').at(-1);
    responses.push(scheduler.schedule(Promise.resolve(), `preview:${index}:${source}`).then(async () => {
      count('previewCompletion');
      if (index === example.sources.length - 1 && example.failure) reject(request);
      else resolve(request, { svg: `<svg>PREVIEW:${index}</svg>`, chuckles: `CHUCKLES:${index}`,
        chem_types: { 1: 'backbone_n' }, leaving: { 1: '[H]' } });
      await tick();
    }));
  }
  if (example.cancelLast) { count('previewCancelledByEdit'); await ui.input('smiles-in', 'O'); }
  await scheduler.waitAll(); await Promise.all(responses); await Promise.all(tasks);
  assert.equal(ui.element('btn-register').disabled, example.cancelLast || example.failure);
  if (example.cancelLast || example.failure) {
    assert.equal(ui.run('detectedData'), null);
    return;
  }
  assert.equal(ui.element('chuckles-out').value, `CHUCKLES:${example.sources.length - 1}`);
  ui.element('chuckles-out').value = 'UNTRUSTED DISPLAY EDIT';
  const writing = ui.element('btn-register').click();
  await ui.element('btn-register').click();
  const write = pending(ui, '/register_monomer');
  assert.equal(write.length, 1, 'A double click cannot duplicate a library write');
  assert.equal(JSON.parse(write[0].options.body).chuckles, `CHUCKLES:${example.sources.length - 1}`);
  if (example.editDuringWrite) await ui.input('abbr-in', 'NextMonomer');
  resolve(write[0], example.writeFailure ? { error: 'busy' } : { total: 1129 }, !example.writeFailure,
    example.writeFailure ? 503 : 200, { 'Retry-After': '1' });
  await writing; await ui.timers(); await tick();
  count(example.writeFailure ? 'failedWriteNoRetry' : 'write');
  assert.equal(ui.requests.filter(r => r.url === '/register_monomer').length, 1, 'Library writes are never automatically retried');
  if (example.editDuringWrite) assert.equal(ui.element('status-msg').style.display, 'none');
  assert.equal(ui.element('btn-register').disabled, !example.writeFailure && !example.editDuringWrite);
});

property('normalization-keeps-redo-and-reference', [fc.array(fc.constantFrom('A', 'K', 'G'), { minLength: 1, maxLength: 5 }), referenceText], async (symbols, reference, count) => {
  let clock = 10000;
  const ui = page('builder.js', { now: () => clock });
  const before = symbols.join('-');
  await ui.input('cabiln-input', before); await settle(ui);
  await ui.input('smiles-input', reference); await settle(ui);
  const original = snapshot(ui).reference;
  clock += 1000;
  await ui.input('cabiln-input', 'D.(4,1)-G%G-A');
  await ui.timers();
  const normalized = 'D.[G(4,1).A(2,1)]-G';
  resolve(pending(ui, '/render')[0], { ...drawing(normalized), normalized_cabiln: normalized });
  await tick(); count('normalization');
  await ui.element('btn-undo').click();
  resolve(pending(ui, '/render')[0], { ...drawing(before), normalized_cabiln: before });
  await tick(); count('normalizationWithRedo');
  assert.equal(ui.element('cabiln-input').value, before);
  assert.equal(ui.element('btn-redo').disabled, false);
  await ui.element('btn-redo').click(); await settle(ui);
  assert.equal(ui.element('cabiln-input').value, normalized);
  assert.deepEqual(snapshot(ui).reference, original);
});

property('selection-and-library-revisions', [fc.scheduler(), fc.array(fc.integer({ min: 0, max: 1 }), { minLength: 2, maxLength: 7 }),
  fc.constantFrom('none', 'close', 'edit'), fc.boolean()], async (scheduler, selections, invalidate, changedLibrary, count) => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'K-A'); await ui.timers();
  resolve(pending(ui, '/render')[0], { ...drawing('K-A'), residue_map: { 0: [0], 1: [1] },
    residues: [{ idx: 0, abbr: 'K' }, { idx: 1, abbr: 'A' }],
    layout: { segments: [{ roots: [0, 1], members: [0, 1] }], groups: [], markers: [] } });
  await tick();
  await ui.element('btn-build').click();
  const initialLibrary = [
    { abbr: 'K', name: 'Lysine', type: 'aa', chem_types: '' },
    { abbr: 'A', name: 'Alanine', type: 'aa', chem_types: '' },
  ];
  resolve(pending(ui, '/monomers')[0], initialLibrary, true, 200, { 'X-Library-Version': 'initial' });
  await settle(ui);
  const chips = ui.element('residue-chips').children;
  const tasks = [];
  const responses = [];
  for (const [revision, id] of selections.entries()) {
    count('select');
    tasks.push(chips[id].click());
    const sites = ui.requests.at(-1);
    assert.match(sites.url, /^\/monomer_rgroups\?/);
    await ui.run('window.dispatchEvent(new Event("focus"))');
    const library = pending(ui, '/monomers').at(-1);
    responses.push(scheduler.schedule(Promise.resolve(), `sites:${revision}:${id}`).then(async () => {
      resolve(sites, { svg: `<svg>SITE:${revision}</svg>`, rgroups: [{ slot: 1, chem_type: 'backbone_n' }] });
      await tick();
    }));
    responses.push(scheduler.schedule(Promise.resolve(), `library:${revision}`).then(async () => {
      const palette = changedLibrary
        ? [{ abbr: `Revision${revision}`, name: 'Generated library revision', type: 'aa', chem_types: '' }]
        : initialLibrary;
      resolve(library, palette, true, 200, { 'X-Library-Version': changedLibrary ? String(revision) : 'initial' });
      await tick();
    }));
  }
  if (invalidate === 'close') { count('closeBuild'); await ui.element('build-close').click(); }
  if (invalidate === 'edit') { count('editCancelsSelection'); await ui.input('cabiln-input', 'G'); }
  await ui.element('lib-close').click();
  await scheduler.waitAll(); await Promise.all(responses); await Promise.all(tasks);
  assert.equal(ui.run('allMonomers[0].abbr'), changedLibrary ? `Revision${selections.length - 1}` : 'K');
  assert.equal(ui.element('lib-panel').classList.contains('open'), false);
  assert.equal(ui.element('lib-preview').style.display, 'none');
  assert.equal(ui.element('build-connect').disabled, true);
  count(changedLibrary ? 'changedLibrary' : 'unchangedLibrary');
  if (invalidate !== 'none' || changedLibrary) assert.equal(ui.run('buildLeft'), null);
  else {
    assert.equal(ui.run('buildLeft.abbr'), selections.at(-1) === 0 ? 'K' : 'A');
    assert.equal(ui.element('build-left-svg').innerHTML, `<svg>SITE:${selections.length - 1}</svg>`);
  }
});
