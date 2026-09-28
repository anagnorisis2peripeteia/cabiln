// Run with: node --test tests/test_frontend.js
// The VM runs the shipped scripts. Deferred fetches deliberately ignore aborts
// to cover responses already queued when an input is edited.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const test = require('node:test');
const vm = require('node:vm');

class Element {
  constructor() {
    this.value = '';
    this.innerHTML = '';
    this.textContent = '';
    this.style = {};
    this.dataset = {};
    this.children = [];
    this.listeners = new Map();
    this.disabled = false;
    this.hidden = true;
    const classes = new Set();
    this.classList = {
      add: name => classes.add(name),
      remove: name => classes.delete(name),
      contains: name => classes.has(name),
      toggle(name, enabled = !classes.has(name)) {
        if (enabled) classes.add(name);
        else classes.delete(name);
      },
    };
  }
  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, []);
    this.listeners.get(type).push(listener);
  }
  dispatchEvent(event) {
    event.target ??= this;
    return Promise.all((this.listeners.get(event.type) || []).map(fn => fn(event)));
  }
  querySelectorAll() { return []; }
  querySelector() { return null; }
  appendChild(element) { this.children.push(element); }
  focus() {}
}

function page(script, { registration = false } = {}) {
  const elements = new Map();
  const element = id => {
    if (!elements.has(id)) elements.set(id, new Element());
    return elements.get(id);
  };
  element('notation-select').value = 'cabiln';
  element('type-in').value = 'aa';
  element('subtype-in').value = 'modified';
  for (const id of ['btn-register', 'btn-png', 'btn-mol']) element(id).disabled = true;
  const requests = [];
  const timers = new Map();
  let timerId = 0;
  const context = vm.createContext({
    console, AbortController,
    Event: class { constructor(type) { this.type = type; } },
    document: {
      getElementById: element,
      querySelectorAll: () => [],
      querySelector: () => null,
      createElement: () => new Element(),
    },
    window: new Element(),
    setTimeout(callback) { timers.set(++timerId, callback); return timerId; },
    clearTimeout(id) { timers.delete(id); },
    setInterval() {},
    fetch(url, options = {}) {
      if (url === '/capabilities') {
        return Promise.resolve({ ok: true, json: async () => ({ registration }) });
      }
      return new Promise((resolve, reject) => {
        requests.push({
          url, options, reject,
          resolve(data, ok = true) { resolve({ ok, json: async () => data }); },
        });
      });
    },
  });
  const filename = path.join(__dirname, '../src/pyPept/web/static', script);
  vm.runInContext(fs.readFileSync(filename, 'utf8'), context, { filename });
  return {
    element, requests,
    run: source => vm.runInContext(source, context),
    async input(id, value) {
      element(id).value = value;
      await element(id).dispatchEvent({ type: 'input' });
    },
    async timers() {
      const pending = [...timers.values()];
      timers.clear();
      for (const callback of pending) callback();
      await Promise.resolve();
    },
  };
}

const rendered = {
  svg: '<svg>OLD STRUCTURE</svg>', mol_block: 'OLD MOL', info: 'old',
  residue_map: {}, residues: [], chains: [], bracket_groups: [], crosslink_groups: [],
};
const preview = {
  svg: '<svg>PREVIEW</svg>', chuckles: '[1*]NCC([2*])=O',
  chem_types: { 1: 'backbone_n', 2: 'backbone_c' },
  leaving: { 1: '[H]', 2: '[OH]' },
};

function tabs(ui, symbols, layout, crosslinks = []) {
  const residues = symbols.map((abbr, idx) => ({ abbr, idx }));
  const atoms = Object.fromEntries(symbols.map((_, idx) => [idx, [idx * 2, idx * 2 + 1]]));
  ui.run(`buildResidueUI(${JSON.stringify(atoms)}, ${JSON.stringify(residues)}, ${JSON.stringify(layout)}, ${JSON.stringify(crosslinks)})`);
  return ui.element('residue-chips').children;
}

function branch(id, host, roots, opening, parent = null, members = roots) {
  return { id, host, roots, members, parent, opening, closing: opening === '{' ? '}' : ']' };
}

test('each branch keeps its own delimiters and occurrence IDs', () => {
  const ui = page('builder.js');
  const chips = tabs(ui, ['ac', 'K', 'K', 'am', 'G', 'A'], {
    segments: [{ roots: [0, 1, 2, 3], members: [0, 1, 2, 3, 4, 5] }],
    groups: [branch(0, 1, [4], '['), branch(1, 2, [5], '{')], markers: [],
  });
  assert.deepEqual(chips.map(chip => chip.textContent), ['ac', 'K', '$', '[', 'G', ']', 'K', '$', '{', 'A', '}', 'am']);
  assert.deepEqual(chips.filter(chip => chip.dataset.residue !== undefined).map(chip => chip.dataset.residue), [0, 1, 4, 2, 5, 3]);
});

test('nested group tabs appear once and highlight their own members', async () => {
  const ui = page('builder.js');
  const chips = tabs(ui, ['ac', 'K', 'D', 'am', 'G', 'A'], {
    segments: [{ roots: [0, 1, 2, 3], members: [0, 1, 2, 3, 4, 5] }],
    groups: [branch(0, 1, [4], '{', null, [4, 5]), branch(1, 4, [5], '[', 0)], markers: [],
  });
  assert.deepEqual(chips.map(chip => chip.textContent), ['ac', 'K', '$', '{', 'G', '$', '[', 'A', ']', '}', 'D', 'am']);
  ui.run('let hovered; highlightGroup = ids => { hovered = ids; }');
  await chips.find(chip => chip.textContent === '{').dispatchEvent({ type: 'mouseenter' });
  assert.equal(ui.run('JSON.stringify(hovered)'), '[4,5]');
  await chips.find(chip => chip.textContent === '[').dispatchEvent({ type: 'mouseenter' });
  assert.equal(ui.run('JSON.stringify(hovered)'), '[5]');
});

test('branches and backbone insertion work on the second explicit segment', () => {
  const ui = page('builder.js');
  const chips = tabs(ui, ['A', 'K', 'A', 'G'], {
    segments: [{ roots: [0], members: [0] }, { roots: [1, 2], members: [1, 2, 3] }],
    groups: [branch(0, 1, [3], '{')], markers: [],
  });
  assert.deepEqual(chips.map(chip => chip.textContent), ['%', 'A', '%', 'K', '$', '{', 'G', '}', 'A']);
  assert.equal(ui.run('buildLeftRIdx = 1; buildRightRIdx = 2; checkAdjacentBackbone()'), true);
  assert.equal(ui.run('buildLeftRIdx = 0; checkAdjacentBackbone()'), false);
  assert.equal(ui.run('buildLeftRIdx = 3; checkAdjacentBackbone()'), false);
});

test('terminal link tabs use recorded placement and highlight both ends', async () => {
  const ui = page('builder.js');
  const members = [0, 2];
  const chips = tabs(ui, ['C', 'A', 'C'], {
    segments: [{ roots: [0, 1, 2], members: [0, 1, 2] }], groups: [],
    markers: [
      { residue: 0, tag: '!ring', members, group: null, before: true },
      { residue: 2, tag: '!ring', members, group: null, before: false },
    ],
  }, [{ tag: '!ring', members }]);
  assert.deepEqual(chips.map(chip => chip.textContent), ['!ring', 'C', 'A', 'C', '!ring']);
  ui.run('let hovered; highlightGroup = ids => { hovered = ids; }');
  for (const chip of chips.filter(chip => chip.textContent === '!ring')) {
    await chip.dispatchEvent({ type: 'mouseenter' });
    assert.equal(ui.run('JSON.stringify(hovered)'), '[0,2]');
  }
});

test('a marker-only protected group keeps its links and their atom owners', () => {
  const ui = page('builder.js');
  const members = [0, 2];
  const chips = tabs(ui, ['C', 'A', 'C'], {
    segments: [{ roots: [0, 1, 2], members: [0, 1, 2] }],
    groups: [branch(0, 2, [], '{')],
    markers: [
      { residue: 0, tag: '!r', members, group: null, before: false },
      { residue: 0, tag: '!s', members, group: null, before: false },
      { residue: 2, tag: '!s', members, group: 0, before: false },
      { residue: 2, tag: '!r', members, group: 0, before: false },
    ],
  }, [{ tag: '!r', members }, { tag: '!s', members }]);
  assert.deepEqual(chips.map(chip => chip.textContent), ['C', '$', '!r', '!s', 'A', 'C', '$', '{', '!s', '!r', '}']);
  assert.deepEqual(JSON.parse(chips.find(chip => chip.textContent === '{').dataset.members), members);
});

test('a crosslink inside a nested arm renders both endpoints', () => {
  const ui = page('builder.js');
  const members = [2, 5];
  const chips = tabs(ui, ['ac', 'K', 'D', 'am', 'G', 'A'], {
    segments: [{ roots: [0, 1, 2, 3], members: [0, 1, 2, 3, 4, 5] }],
    groups: [branch(0, 1, [4], '{', null, [4, 5]), branch(1, 4, [5], '[', 0)],
    markers: [
      { residue: 5, tag: '!1', members, group: 1, before: false },
      { residue: 2, tag: '!1', members, group: null, before: false },
    ],
  }, [{ tag: '!1', members }]);
  assert.deepEqual(chips.map(chip => chip.textContent), ['ac', 'K', '$', '{', 'G', '$', '[', 'A', '!1', ']', '}', 'D', '!1', 'am']);
});

test('an arm containing only a marker keeps its brackets', () => {
  const ui = page('builder.js');
  const members = [2, 4];
  const chips = tabs(ui, ['ac', 'K', 'D', 'am', 'G'], {
    segments: [{ roots: [0, 1, 2, 3], members: [0, 1, 2, 3, 4] }],
    groups: [branch(0, 1, [4], '[', null, [4]), branch(1, 4, [], '[', 0)],
    markers: [
      { residue: 4, tag: '!1', members, group: 1, before: false },
      { residue: 2, tag: '!1', members, group: null, before: false },
    ],
  }, [{ tag: '!1', members }]);
  assert.deepEqual(chips.map(chip => chip.textContent), ['ac', 'K', '$', '[', 'G', '$', '[', '!1', ']', ']', 'D', '!1', 'am']);
});

test('reopening the palette discovers new monomers and preserves its search', async () => {
  const ui = page('builder.js');
  const gly = { abbr: 'G', name: 'Glycine', type: 'aa', subtype: 'natural', chem_types: '' };
  const fresh = { ...gly, abbr: 'NewThiol', name: 'New thiol' };
  const first = ui.run('loadMonomers()');
  ui.requests[0].resolve([gly]);
  await first;
  ui.element('lib-search').value = 'new';
  ui.run('openLib()');
  const refresh = ui.requests.find((r, i) => i > 0 && r.url === '/monomers');
  assert.ok(refresh);
  refresh.resolve([gly, fresh]);
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(ui.element('lib-search').value, 'new');
  assert.match(ui.element('lib-list').innerHTML, /NewThiol/);
  assert.equal(ui.run('allMonomers.length'), 2);
});

test('returning from registration refreshes an open palette', async () => {
  const ui = page('builder.js');
  ui.run("libPanel.classList.add('open')");
  await ui.run("window.dispatchEvent(new Event('focus'))");
  assert.equal(ui.requests[0].url, '/monomers');
});

test('older library refresh cannot replace a newer library', async () => {
  const ui = page('builder.js');
  const old = ui.run('loadMonomers()');
  const current = ui.run('loadMonomers()');
  ui.requests[1].resolve([{ abbr: 'New', name: 'New', type: 'aa', subtype: '', chem_types: '' }]);
  await current;
  ui.requests[0].resolve([]);
  await old;
  assert.equal(ui.run('allMonomers[0].abbr'), 'New');
});

test('a discovered cap family asks for its attachment form before loading slots', async () => {
  const ui = page('builder.js');
  ui.run("allMonomers = [{abbr:'NovelCap', degenerate:true, nterm_abbr:'NovelCap_', cterm_abbr:'_NovelCap'}]");
  await ui.run("loadBuildRight('NovelCap')");
  assert.equal(ui.requests.length, 0);
  const choices = ui.element('build-right-rgroups').children;
  assert.equal(choices.length, 2);
  const pending = choices[1].dispatchEvent({ type: 'click' });
  assert.equal(ui.requests[0].url, '/monomer_rgroups?abbr=_NovelCap');
  ui.requests[0].resolve({ svg: '<svg/>', rgroups: [] });
  await pending;
});

test('conversion warnings survive rendering and clear when the sequence is edited', async () => {
  const ui = page('builder.js');
  ui.run("useConvertedCabiln('A', 'Input stereochemistry was inferred')");
  const pending = ui.run("doRenderCabiln('A')");
  ui.requests[0].resolve(rendered);
  await pending;
  assert.equal(ui.element('conversion-status').hidden, false);
  assert.match(ui.element('conversion-status').textContent, /stereochemistry/);
  await ui.input('cabiln-input', 'G');
  assert.equal(ui.element('conversion-status').hidden, true);
});

test('starting detection makes the preview visible over its CSS default', async () => {
  const ui = page('register.js');
  await ui.input('smiles-in', 'NCC(=O)O');
  const pending = ui.element('btn-preview').dispatchEvent({ type: 'click' });
  assert.equal(ui.element('preview-section').style.display, 'block');
  ui.requests[0].resolve(preview);
  await pending;
});

test('legacy normalization updates the source before exposing selectable residues', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'D.(4,1)-G%G-A';
  const pending = ui.run("doRenderCabiln(cabilnInput.value)");
  const normalized = 'D.[G(4,1).A(2,1)]-G';
  ui.requests[0].resolve({ ...rendered, normalized_cabiln: normalized, cabiln_echo: normalized });
  await pending;
  assert.equal(ui.element('cabiln-input').value, normalized);
  assert.equal(ui.run('lastCabiln'), normalized);
});

test('clearing the main input discards an already pending render and exports', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doRenderCabiln("A-G")');
  await ui.input('cabiln-input', '');
  ui.requests[0].resolve(rendered);
  await pending;
  assert.doesNotMatch(ui.element('render-inner').innerHTML, /OLD STRUCTURE/);
  assert.equal(ui.run('lastCabiln'), '');
  assert.equal(ui.element('btn-mol').disabled, true);
});

test('editing main input invalidates rendering before the debounce fires', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doRenderCabiln("A-G")');
  await ui.input('cabiln-input', 'G-A');
  ui.requests[0].resolve(rendered);
  await pending;
  assert.doesNotMatch(ui.element('render-inner').innerHTML, /OLD STRUCTURE/);
  assert.equal(ui.element('btn-png').disabled, true);
});

test('clearing the reference invalidates its pending render', async () => {
  const ui = page('builder.js');
  ui.element('smiles-input').value = 'NCC(=O)O';
  const pending = ui.run('doRenderRef(smilesInput.value)');
  await ui.input('smiles-input', '');
  ui.requests[0].resolve({ ...rendered, smiles: 'NCC(=O)O', format: 'SMILES' });
  await pending;
  assert.equal(ui.run('lastSmiles'), '');
  assert.doesNotMatch(ui.element('smiles-inner').innerHTML, /OLD STRUCTURE/);
});

test('closing build mode discards pending R-group selections', async () => {
  const ui = page('builder.js');
  const pending = ui.run('loadBuildLeft("K", 0)');
  ui.run('clearBuild()');
  ui.requests[0].resolve({ svg: '<svg>OLD</svg>', rgroups: [{ slot: 4 }] });
  await pending;
  assert.equal(ui.run('buildLeft'), null);
  assert.equal(ui.element('build-left-abbr').textContent, '—');
});

test('notation conversion cannot replace an input edited during the request', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('convertNotation("bracket")');
  await ui.input('cabiln-input', 'K-C');
  ui.requests[0].resolve({ result: 'OLD CONVERSION' });
  await pending;
  assert.equal(ui.element('cabiln-input').value, 'K-C');
});

test('editing monomer SMILES invalidates detected chemistry and registration', async () => {
  const ui = page('register.js');
  ui.element('smiles-in').value = 'NCC(=O)O';
  const pending = ui.element('btn-preview').dispatchEvent({ type: 'click' });
  ui.requests[0].resolve(preview);
  await pending;
  await ui.input('smiles-in', 'CC(=O)O');
  assert.equal(ui.run('detectedData'), null);
  assert.equal(ui.element('btn-register').disabled, true);
  assert.equal(ui.element('chuckles-out').value, '');
});

test('a monomer preview arriving after an edit is discarded', async () => {
  const ui = page('register.js');
  ui.element('smiles-in').value = 'NCC(=O)O';
  const pending = ui.element('btn-preview').dispatchEvent({ type: 'click' });
  await ui.input('smiles-in', 'CC(=O)O');
  ui.requests[0].resolve(preview);
  await pending;
  assert.equal(ui.run('detectedData'), null);
  assert.equal(ui.element('btn-register').disabled, true);
});

test('the latest successful render enables exports; an older response cannot replace it', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const old = ui.run('doRenderCabiln("A-G")');
  await ui.input('cabiln-input', 'G-A');
  const latest = ui.run('doRenderCabiln("G-A")');
  ui.requests[1].resolve({ ...rendered, svg: '<svg>NEW</svg>', cabiln_echo: 'G-A' });
  await latest;
  ui.requests[0].resolve(rendered);
  await old;
  assert.equal(ui.element('render-inner').innerHTML, '<svg>NEW</svg>');
  assert.equal(ui.element('btn-mol').disabled, false);
  assert.equal(ui.run('lastCabiln'), 'G-A');
});

test('changing notation cancels a pending foreign render', async () => {
  const ui = page('builder.js');
  ui.element('notation-select').value = 'smiles';
  ui.element('cabiln-input').value = 'NCC(=O)O';
  const pending = ui.run('doRenderForeign(cabilnInput.value)');
  ui.element('notation-select').value = 'cabiln';
  await ui.element('notation-select').dispatchEvent({ type: 'change' });
  ui.requests[0].resolve({ ...rendered, format: 'SMILES' });
  await pending;
  assert.doesNotMatch(ui.element('render-inner').innerHTML, /OLD STRUCTURE/);
  assert.equal(ui.element('btn-png').disabled, true);
});

test('editing either structure discards a pending verification result', async () => {
  for (const field of ['cabiln-input', 'smiles-input']) {
    const ui = page('builder.js');
    ui.run('verifyMode = true; lastCabiln = "A-G"; lastSmiles = "NCC(=O)O"');
    const pending = ui.run('triggerVerify()');
    await ui.input(field, '');
    ui.requests[0].resolve({
      match: true, smiles_canonical: 'OLD', cabiln_canonical: 'OLD',
    });
    await pending;
    assert.equal(ui.element('compare-bar').innerHTML, '');
  }
});

test('competing build selections keep the latest R-group response on each side', async () => {
  for (const side of ['Left', 'Right']) {
    const ui = page('builder.js');
    const first = ui.run(`loadBuild${side}("K", 0)`);
    const second = ui.run(`loadBuild${side}("C", 1)`);
    ui.requests[1].resolve({ svg: '<svg>C</svg>', rgroups: [{ slot: 4 }] });
    await second;
    ui.requests[0].resolve({ svg: '<svg>K</svg>', rgroups: [{ slot: 5 }] });
    await first;
    assert.equal(ui.run(`build${side}.abbr`), 'C');
    assert.equal(ui.run(`build${side}.rgroups[0].slot`), 4);
  }
});

test('a cleared R-group selection cannot be enabled by old bond validation', async () => {
  const ui = page('builder.js');
  ui.run(`
    buildLeft = {abbr: 'K', selectedSlot: 4, rgroups: [{slot: 4, chem_type: 'amine_primary'}]};
    buildRight = {abbr: 'D', selectedSlot: 4, rgroups: [{slot: 4, chem_type: 'carboxyl'}]};
  `);
  const pending = ui.run('checkBuildValidity()');
  ui.run('selectRgroup("left", 4)');
  ui.requests[0].resolve({ valid: true, reaction: 'amide' });
  await pending;
  assert.equal(ui.element('build-connect').disabled, true);
  assert.equal(ui.element('build-status').textContent, '');
});

test('cancelled SMILES conversion restores the button without replacing main input', async () => {
  const ui = page('builder.js');
  ui.element('btn-s2c').textContent = '→ %';
  ui.element('smiles-input').value = 'NCC(=O)O';
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doS2c("percent")');
  await ui.input('smiles-input', 'CC(=O)O');
  assert.equal(ui.element('btn-s2c').disabled, false);
  assert.equal(ui.element('btn-s2c').textContent, '→ %');
  ui.requests[0].resolve({ cabiln: 'OLD', details: [] });
  await pending;
  assert.equal(ui.element('cabiln-input').value, 'A-G');
});

test('conversion results cannot override a notation switch', async () => {
  const ui = page('builder.js');
  ui.element('notation-select').value = 'smiles';
  ui.element('cabiln-input').value = 'NCC(=O)O';
  const pending = ui.run('doToCabiln("bracket")');
  ui.element('notation-select').value = 'helm';
  await ui.element('notation-select').dispatchEvent({ type: 'change' });
  ui.requests[0].resolve({ cabiln: 'OLD', from: 'SMILES' });
  await pending;
  assert.equal(ui.element('notation-select').value, 'helm');
  assert.equal(ui.element('cabiln-input').value, '');
});

test('a pending connection cannot overwrite manual sequence edits', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-K';
  ui.run(`
    buildLeft = {abbr: 'K', selectedSlot: 4};
    buildRight = {abbr: 'D', selectedSlot: 4};
  `);
  const pending = ui.element('build-connect').dispatchEvent({ type: 'click' });
  await ui.input('cabiln-input', 'A-C');
  ui.requests[0].resolve({ result: 'OLD' });
  await pending;
  assert.equal(ui.element('cabiln-input').value, 'A-C');
});

test('cancelling backbone insertion discards the pending insertion result', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-K';
  ui.run(`
    chainData = [{residues: [0, 1]}];
    buildLeftRIdx = 0; buildRightRIdx = 1; insertBetweenActive = true;
  `);
  const pending = ui.run('doInsertBetween("G")');
  await ui.element('build-insert-btn').dispatchEvent({ type: 'click' });
  ui.requests[0].resolve({ result: 'A-G-K' });
  await pending;
  assert.equal(ui.element('cabiln-input').value, 'A-K');
  assert.equal(ui.run('insertBetweenActive'), false);
});

test('file reading is cancelled by reference input edits before upload rendering', async () => {
  const ui = page('builder.js');
  let finishReading;
  const pending = ui.element('mol-upload').dispatchEvent({
    type: 'change',
    target: { files: [{ name: 'old.mol', text: () => new Promise(resolve => { finishReading = resolve; }) }] },
  });
  await ui.input('smiles-input', 'NEW');
  finishReading('OLD MOL');
  await pending;
  assert.equal(ui.requests.length, 0);
  assert.equal(ui.element('smiles-input').value, 'NEW');
});

test('a failed second preview cannot reuse the earlier chemistry', async () => {
  const ui = page('register.js');
  ui.element('smiles-in').value = 'NCC(=O)O';
  let pending = ui.element('btn-preview').dispatchEvent({ type: 'click' });
  ui.requests[0].resolve(preview);
  await pending;
  pending = ui.element('btn-preview').dispatchEvent({ type: 'click' });
  ui.requests[1].reject(new Error('Network error'));
  await pending;
  assert.equal(ui.run('detectedData'), null);
  assert.equal(ui.element('btn-register').disabled, true);
});

test('registration uses preview chemistry and permits a second record after success', async () => {
  const ui = page('register.js');
  ui.element('smiles-in').value = 'NCC(=O)O';
  await ui.input('abbr-in', 'TestGly');
  await ui.input('name-in', 'Test glycine');
  let pending = ui.element('btn-preview').dispatchEvent({ type: 'click' });
  ui.requests[0].resolve(preview);
  await pending;
  assert.equal(ui.element('btn-register').disabled, false);
  ui.element('chuckles-out').value = 'UNRELATED EDIT';
  pending = ui.element('btn-register').dispatchEvent({ type: 'click' });
  assert.equal(JSON.parse(ui.requests[1].options.body).chuckles, preview.chuckles);
  ui.requests[1].resolve({ total: 1000 });
  await pending;
  assert.equal(ui.element('btn-register').disabled, true);
  assert.equal(ui.element('btn-register').textContent, '✓ Registered');
  await ui.input('smiles-in', 'CC(=O)O');
  await ui.input('abbr-in', 'TestAc');
  pending = ui.element('btn-preview').dispatchEvent({ type: 'click' });
  ui.requests[2].resolve({ ...preview, chuckles: 'CC([2*])=O' });
  await pending;
  assert.equal(ui.element('btn-register').disabled, false);
  assert.equal(ui.element('btn-register').textContent, 'Register monomer');
});

test('Register is shown only when the server explicitly enables registration', async () => {
  for (const registration of [false, true]) {
    const ui = page('builder.js', { registration });
    await new Promise(setImmediate);
    assert.equal(ui.element('register-link').hidden, !registration);
  }
  const html = fs.readFileSync(path.join(__dirname, '../src/pyPept/web/static/index.html'), 'utf8');
  assert.match(html, /<a\b[^>]*id="register-link"[^>]*\bhidden\b/);
});
