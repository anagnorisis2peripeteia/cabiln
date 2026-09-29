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
    this.attributes = new Map();
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
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
  getAttribute(name) { return this.attributes.get(name) ?? null; }
  querySelector() { return null; }
  appendChild(element) { this.children.push(element); }
  click() { return this.dispatchEvent({ type: 'click' }); }
  focus() {}
}

function page(script, { registration = false, storedDraft = null, now = Date.now } = {}) {
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
  const downloads = [];
  const storage = new Map(storedDraft ? [['cabiln.draft.v1', JSON.stringify(storedDraft)]] : []);
  const window = new Element();
  window.localStorage = {
    getItem: key => storage.get(key) || null,
    setItem: (key, value) => storage.set(key, value),
    removeItem: key => storage.delete(key),
  };
  const context = vm.createContext({
    console, AbortController, Blob,
    Date: class extends Date { static now() { return now(); } },
    URL: { createObjectURL(blob) { downloads.push(blob); return 'blob:project'; }, revokeObjectURL() {} },
    Event: class { constructor(type) { this.type = type; } },
    document: {
      getElementById: element,
      querySelectorAll: () => [],
      querySelector: () => null,
      createElement: () => new Element(),
    },
    window,
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
          resolve(data, ok = true, status = ok ? 200 : 400, headers = {}) {
            resolve({ ok, status, headers: { get: key => headers[key] || null }, json: async () => data });
          },
        });
      });
    },
  });
  const filename = path.join(__dirname, '../src/pyPept/web/static', script);
  if (script === 'builder.js') {
    for (const name of ['project.js', 'document.js']) {
      const dependency = path.join(path.dirname(filename), name);
      vm.runInContext(fs.readFileSync(dependency, 'utf8'), context, { filename: dependency });
    }
  }
  vm.runInContext(fs.readFileSync(filename, 'utf8'), context, { filename });
  return {
    element, requests, downloads, storage,
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

test('older library refresh cannot replace a newer library', async () => {
  const ui = page('builder.js');
  const old = ui.run('loadMonomers()');
  ui.run("libPanel.classList.add('open')");
  ui.run("window.dispatchEvent(new Event('focus'))");
  ui.requests[1].resolve([{ abbr: 'New', name: 'New', type: 'aa', subtype: '', chem_types: '' }]);
  await new Promise(setImmediate);
  ui.requests[0].resolve([]);
  await old;
  assert.equal(ui.run('allMonomers[0].abbr'), 'New');
  assert.match(ui.element('lib-list').innerHTML, /New/);
});

function previewRow(ui) {
  ui.run(`
    window.innerHeight = 1000;
    const previewRow = { getBoundingClientRect: () => ({top: 50, right: 300}) };
  `);
}

test('hiding a library preview discards its queued response', async () => {
  const ui = page('builder.js');
  previewRow(ui);
  ui.run('startPreview("G", previewRow)');
  await ui.timers();
  ui.run('hidePreview()');
  ui.requests[0].resolve({svg: '<svg>OLD PREVIEW</svg>'});
  await new Promise(setImmediate);
  assert.equal(ui.element('lib-preview').style.display, 'none');
  assert.doesNotMatch(ui.element('lib-preview').innerHTML, /OLD PREVIEW/);
});

test('an older library hover cannot replace the latest preview', async () => {
  const ui = page('builder.js');
  previewRow(ui);
  ui.run('startPreview("G", previewRow)');
  await ui.timers();
  ui.run('startPreview("A", previewRow)');
  await ui.timers();
  ui.requests[1].resolve({svg: '<svg>CURRENT PREVIEW</svg>'});
  await new Promise(setImmediate);
  ui.requests[0].resolve({svg: '<svg>OLD PREVIEW</svg>'});
  await new Promise(setImmediate);
  assert.equal(ui.element('lib-preview').style.display, 'block');
  assert.match(ui.element('lib-preview').innerHTML, /CURRENT PREVIEW/);
  assert.doesNotMatch(ui.element('lib-preview').innerHTML, /OLD PREVIEW/);
});

test('closing a cached preview before its debounce leaves it hidden', async () => {
  const ui = page('builder.js');
  previewRow(ui);
  ui.run('previewCache.G = {svg: "<svg>CACHED PREVIEW</svg>"}; startPreview("G", previewRow)');
  ui.run('hidePreview()');
  await ui.timers();
  assert.equal(ui.requests.length, 0);
  assert.equal(ui.element('lib-preview').style.display, 'none');
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

test('canonical formatting is explicit and the original spelling remains undoable', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'Ala-Gly');
  ui.element('notation-policy').value = 'canonical';
  const pending = ui.run('convertNotation("bracket")');
  const request = ui.requests.find(request => request.url === '/convert_notation');
  assert.deepEqual(JSON.parse(request.options.body), {
    cabiln: 'Ala-Gly', target: 'bracket', canonical: true,
  });
  request.resolve({ result: 'A-G' });
  await pending;
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  ui.run('travelHistory("undo")');
  assert.equal(ui.element('cabiln-input').value, 'Ala-Gly');
  ui.run('travelHistory("redo")');
  assert.equal(ui.element('cabiln-input').value, 'A-G');
});

test('typing coalesces Undo while server normalization preserves Redo', async () => {
  let clock = 1000;
  const ui = page('builder.js', { now: () => clock });
  await ui.input('cabiln-input', 'Ala');
  clock += 300;
  await ui.input('cabiln-input', 'Ala-Gly');
  await ui.element('btn-undo').click();
  assert.equal(ui.element('cabiln-input').value, '');
  await ui.element('btn-redo').click();
  assert.equal(ui.element('cabiln-input').value, 'Ala-Gly');

  clock += 1000;
  await ui.input('cabiln-input', 'Ala-Gly-Lys');
  await ui.element('btn-undo').click();
  assert.equal(ui.element('cabiln-input').value, 'Ala-Gly');
  latestRequest(ui, '/render').resolve({ ...rendered, normalized_cabiln: 'A-G' });
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  assert.equal(ui.element('btn-redo').disabled, false);
  assert.equal(JSON.parse(ui.storage.get('cabiln.draft.v1')).document.text, 'A-G');
  await ui.element('btn-redo').click();
  assert.equal(ui.element('cabiln-input').value, 'Ala-Gly-Lys');
  await ui.element('btn-undo').click();
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  clock += 1000;
  await ui.input('cabiln-input', 'S');
  assert.equal(ui.element('btn-redo').disabled, true);
});

test('layout formatting stays the default and errors retain the document', async () => {
  const ui = page('builder.js');
  const original = 'K.{G(4,2).ac(1,2)}-A';
  await ui.input('cabiln-input', original);
  const rendering = ui.run('doRenderCabiln(cabilnInput.value)');
  const pending = ui.run('convertNotation("branch")');
  const request = ui.requests.find(request => request.url === '/convert_notation');
  assert.equal(JSON.parse(request.options.body).canonical, false);
  request.resolve({ error: 'Cannot preserve the requested attachment' }, false);
  await pending;
  assert.equal(ui.element('cabiln-input').value, original);
  assert.equal(ui.element('conversion-status').hidden, false);
  assert.equal(ui.element('conversion-status').textContent, 'Cannot preserve the requested attachment');
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await rendering;
  assert.equal(ui.element('conversion-status').hidden, false);
  assert.equal(ui.element('conversion-status').textContent, 'Cannot preserve the requested attachment');
  assert.equal(JSON.parse(ui.storage.get('cabiln.draft.v1')).document.warning,
    'Cannot preserve the requested attachment');
});

test('recursive sibling tabs return to their parent and highlight exact descendants', async () => {
  const ui = page('builder.js');
  const chips = tabs(ui, ['K', 'A', 'K', 'G', 'ac', 'ac'], {
    segments: [{ roots: [0, 1], members: [0, 1, 2, 3, 4, 5] }],
    groups: [
      branch(0, 0, [2], '[', null, [2, 3, 4, 5]),
      branch(1, 2, [3], '[', 0, [3, 4]),
      branch(2, 3, [4], '[', 1),
      branch(3, 2, [5], '[', 0),
    ], markers: [],
  });
  assert.deepEqual(chips.filter(chip => chip.dataset.residue !== undefined)
    .map(chip => chip.dataset.residue), [0, 2, 3, 4, 5, 1]);
  ui.run('var hovered = null; highlightGroup = ids => { hovered = ids; }');
  const brackets = chips.filter(chip => chip.textContent === '[');
  for (const [index, members] of [[0, [2, 3, 4, 5]], [1, [3, 4]], [2, [4]], [3, [5]]]) {
    await brackets[index].dispatchEvent({ type: 'mouseenter' });
    assert.equal(ui.run('JSON.stringify(hovered)'), JSON.stringify(members));
  }
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
  assert.equal(ui.element('build-status').textContent, 'Choose a site on the current residue');
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
  assert.deepEqual(JSON.parse(ui.requests[0].options.body), {
    input: 'NCC(=O)O', input_format: 'smiles', notation: 'bracket',
  });
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
  assert.equal(ui.element('preview-section').style.display, 'block');
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

const binding = {
  project_version: 1,
  library_binding: { monomers: 'library-one', aliases: 'aliases', reactions: 'rules', caps: 'caps' },
  canonical: { format: 'cabiln-graph-v1', labeling: 'rdkit-colored-port-graph-v1', rdkit: '2026.03.1' },
};
const imported = {
  cabiln: 'G-<NCC(=O)O>', context: binding, recognition_status: 'partial',
  search_complete: false, inferred_stereo: true, synthetic_components: [0],
  warnings: ['Input stereochemistry was inferred'],
  assignments: [
    { symbol: 'G', recognized: true, residue_index: 0, source_atoms: [0, 1], attachments: [[2, 1]] },
    { symbol: '_SYN0', recognized: false, residue_index: 1, source_atoms: [2, 3], attachments: [[1, 2]] },
  ],
};
function portableProject() {
  const quality = { ...imported };
  delete quality.cabiln;
  delete quality.context;
  const document = { text: imported.cabiln, notation: 'cabiln',
    warning: imported.warnings[0], quality, canonical: null, context: binding };
  return { format: 'cabiln-project', version: 1, document,
    drafts: { cabiln: document, smiles: { text: 'NCC(=O)O', notation: 'smiles' } },
    reference: { text: 'NCC(=O)O', original: { kind: 'mol', content: 'EXACT ORIGINAL\r\nMOL\n', name: 'source.mol' }, context: binding },
    context: binding, saved_at: '2026-09-29T12:00:00Z' };
}
function latestRequest(ui, url) {
  return [...ui.requests].reverse().find(request => request.url === url);
}
function uploadProject(ui, project) {
  const content = typeof project === 'string' ? project : JSON.stringify(project);
  return ui.element('project-upload').dispatchEvent({ type: 'change', target: { files: [{
    size: Buffer.byteLength(content), text: async () => content,
  }] } });
}
function setImported(ui) {
  ui.run(`useConvertedCabiln(${JSON.stringify(imported.cabiln)}, ${JSON.stringify(imported.warnings[0])}, ${JSON.stringify(imported)})`);
}

test('conversion records structured quality through rendering, notation drafts, Undo and Redo', async () => {
  const ui = page('builder.js');
  ui.element('notation-select').value = 'smiles';
  await ui.input('cabiln-input', 'NCC(=O)O');
  const converting = ui.run('doToCabiln("bracket")');
  latestRequest(ui, '/to_cabiln').resolve({ ...imported, from: 'SMILES', warning: imported.warnings[0] });
  await converting;
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
  assert.equal(ui.element('conversion-status').hidden, false);
  assert.equal(ui.element('conversion-status').textContent, imported.warnings[0]);
  assert.match(ui.element('quality-summary').textContent, /Partial recognition.*1\/2 library matches.*search limited.*stereochemistry inferred/);
  assert.equal(ui.element('import-quality').hidden, false);
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(editor.present.quality.assignments)')), imported.assignments);
  ui.element('notation-select').value = 'smiles';
  await ui.element('notation-select').dispatchEvent({ type: 'change' });
  assert.equal(ui.element('import-quality').hidden, true);
  assert.equal(ui.element('cabiln-input').value, 'NCC(=O)O');
  ui.element('notation-select').value = 'cabiln';
  await ui.element('notation-select').dispatchEvent({ type: 'change' });
  assert.equal(ui.element('import-quality').hidden, false);
  await ui.input('cabiln-input', 'A-G');
  assert.equal(ui.element('conversion-status').hidden, true);
  assert.equal(ui.element('import-quality').hidden, true);
  assert.equal(ui.run('editor.present.quality'), null);
  await ui.element('btn-undo').click();
  assert.equal(ui.element('import-quality').hidden, false);
  assert.equal(ui.element('conversion-status').textContent, imported.warnings[0]);
  await ui.element('btn-redo').click();
  assert.equal(ui.run('editor.present.quality'), null);
});

test('formatting clears occurrence assignments while preserving warning and canonical convention', async () => {
  const ui = page('builder.js');
  setImported(ui);
  ui.element('notation-policy').value = 'canonical';
  const pending = ui.run('convertNotation("bracket")');
  latestRequest(ui, '/convert_notation').resolve({ result: '<NCC(=O)O>-G', context: binding,
    canonical: { ...binding.canonical, binding: binding.library_binding } });
  await pending;
  assert.equal(ui.run('editor.present.quality'), null);
  assert.match(ui.element('canonical-status').textContent, /cabiln-graph-v1.*2026.03.1/);
  assert.equal(ui.element('conversion-status').textContent, imported.warnings[0]);
  await ui.element('btn-undo').click();
  assert.equal(ui.run('editor.present.quality.assignments.length'), 2);
  assert.equal(ui.element('canonical-status').hidden, true);
});

test('new library or convention clears import evidence and cancels an outdated project save', async () => {
  for (const current of [
    { ...binding, library_binding: { ...binding.library_binding, monomers: 'library-two' } },
    { ...binding, canonical: { ...binding.canonical, rdkit: '2026.09.1' } },
  ]) {
    const ui = page('builder.js');
    setImported(ui);
    const saving = ui.element('btn-project-save').click();
    const request = latestRequest(ui, '/prepare_project');
    latestRequest(ui, '/render').resolve({ ...rendered, context: current });
    await new Promise(setImmediate);
    request.resolve({ project: portableProject() });
    await saving;
    assert.equal(ui.downloads.length, 0);
    assert.equal(request.options.signal.aborted, true);
    assert.equal(ui.run('editor.present.quality'), null);
    assert.equal(ui.element('conversion-status').hidden, true);
    assert.match(ui.element('project-status').textContent, /Import evidence was cleared/);
  }
});

test('synthetic and opaque occurrences remain selectable and library symbols are never guessed from names', async () => {
  const ui = page('builder.js');
  const residues = [
    { idx: 0, abbr: '_SYNcustom', kind: 'library' },
    { idx: 1, abbr: '_SYN0', kind: 'synthetic' },
    { idx: 2, abbr: '_SYN1', kind: 'opaque' },
  ];
  ui.run(`buildResidueUI({}, ${JSON.stringify(residues)}, {segments: [{roots: [0, 1, 2], members: [0, 1, 2]}], groups: [], markers: []})`);
  ui.run('let selected = []; selectBuildResidue = residue => selected.push(residue.idx)');
  const chips = ui.element('residue-chips').children;
  assert.equal(chips[0].dataset.quality, undefined);
  assert.equal(chips[1].dataset.quality, 'synthetic');
  assert.equal(chips[2].dataset.quality, 'opaque');
  assert.match(chips[2].getAttribute('aria-label'), /Opaque preserved fragment.*editable/);
  for (const chip of chips) {
    assert.equal(chip.disabled, false);
    await chip.click();
  }
  assert.equal(ui.run('JSON.stringify(selected)'), '[0,1,2]');
});

test('palette and preview show library quality issues as escaped text', async () => {
  const ui = page('builder.js');
  const monomer = { abbr: 'G', name: 'Glycine', type: 'aa', subtype: 'natural', chem_types: '',
    quality: { status: 'review', issues: [{ code: 'curation', message: 'Review <atom> mapping' }] } };
  const pending = ui.run('loadMonomers()');
  latestRequest(ui, '/monomers').resolve([monomer]);
  await pending;
  assert.match(ui.element('lib-list').innerHTML, /lib-quality.*Review &lt;atom&gt; mapping/);
  assert.match(ui.element('lib-list').innerHTML, /aria-label="G: Glycine. Library quality: Review/);
  previewRow(ui);
  ui.run(`showPreview(${JSON.stringify({ ...preview, quality: monomer.quality })}, previewRow)`);
  assert.match(ui.element('lib-preview').innerHTML, /Library quality: Review &lt;atom&gt; mapping/);
  assert.doesNotMatch(ui.element('lib-preview').innerHTML, /<atom>/);
});

test('project save downloads only server-prepared source, drafts, evidence and exact reference', async () => {
  const ui = page('builder.js');
  setImported(ui);
  ui.run(`referenceOriginal = ${JSON.stringify(portableProject().reference.original)}; referenceContext = ${JSON.stringify(binding)}; smilesInput.value = 'NCC(=O)O'; editor.drafts.smiles = {text: 'NCC(=O)O', notation: 'smiles', warning: '', quality: null, canonical: null, context: null}`);
  const pending = ui.element('btn-project-save').click();
  const request = latestRequest(ui, '/prepare_project');
  const submitted = JSON.parse(request.options.body).project;
  assert.equal(ui.downloads.length, 0);
  assert.deepEqual(submitted.reference, portableProject().reference);
  assert.deepEqual(submitted.document.quality.assignments, imported.assignments);
  const context = { ...binding, resolutions: { document: { source: 'source', resolved: 'chemistry' } } };
  const prepared = { ...submitted, context,
    document: { ...submitted.document, context },
    drafts: Object.fromEntries(Object.entries(submitted.drafts).map(([mode, state]) => [mode, { ...state, context }])) };
  request.resolve({ project: prepared });
  await pending;
  assert.equal(ui.element('btn-project-save').disabled, false);
  assert.equal(ui.downloads.length, 1);
  assert.deepEqual(JSON.parse(await ui.downloads[0].text()), prepared);
  const draft = JSON.parse(ui.storage.get('cabiln.draft.v1'));
  assert.deepEqual(draft.document.context, context);
  assert.deepEqual(draft.reference_original, portableProject().reference.original);
});

test('validated project open restores all document metadata and accepts proved unrelated library additions', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'A-G');
  const pending = uploadProject(ui, portableProject());
  await new Promise(setImmediate);
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  const current = { ...binding, library_binding: { ...binding.library_binding, monomers: 'with-unrelated-addition' } };
  latestRequest(ui, '/validate_project').resolve({ valid: true, context: current });
  await pending;
  assert.equal(ui.element('cabiln-input').value, imported.cabiln);
  assert.equal(ui.element('import-quality').hidden, false);
  assert.equal(ui.element('smiles-input').value, 'NCC(=O)O');
  assert.equal(ui.run('referenceOriginal.content'), portableProject().reference.original.content);
  latestRequest(ui, '/render').resolve({ ...rendered, context: current });
  await new Promise(setImmediate);
  assert.equal(ui.element('import-quality').hidden, false);
  assert.equal(ui.run('editor.drafts.smiles.text'), 'NCC(=O)O');
  await ui.element('btn-undo').click();
  assert.equal(ui.element('cabiln-input').value, 'A-G');
});

test('binding mismatch or incompatible conventions leave current work and history unchanged', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'A-G');
  const history = ui.run('editor.past.length');
  const pending = uploadProject(ui, portableProject());
  await new Promise(setImmediate);
  latestRequest(ui, '/validate_project').resolve({ code: 'project_binding_mismatch',
    error: 'Selected definitions changed. Use the original library, or open the JSON to recover the source text.' }, false);
  await pending;
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  assert.equal(ui.run('editor.past.length'), history);
  assert.equal(ui.downloads.length, 0);
  assert.match(ui.element('project-status').textContent, /original library.*current work is unchanged/);
});

test('invalid project JSON, versions and malformed quality cannot mutate current state', async () => {
  const invalid = portableProject();
  invalid.document.quality.assignments = [null];
  for (const project of ['not json', { ...portableProject(), version: 2 }, invalid]) {
    const ui = page('builder.js');
    await ui.input('cabiln-input', 'A-G');
    await uploadProject(ui, project);
    assert.equal(ui.element('cabiln-input').value, 'A-G');
    assert.equal(latestRequest(ui, '/validate_project'), undefined);
    assert.match(ui.element('project-status').textContent, /current work is unchanged/);
  }
});

test('project validation and preparation responses cannot overwrite an intervening edit', async () => {
  for (const action of ['open', 'save']) {
    const ui = page('builder.js');
    await ui.input('cabiln-input', 'A-G');
    const pending = action === 'open' ? uploadProject(ui, portableProject()) : ui.element('btn-project-save').click();
    await new Promise(setImmediate);
    const request = latestRequest(ui, action === 'open' ? '/validate_project' : '/prepare_project');
    await ui.input('cabiln-input', 'A-K');
    request.resolve(action === 'open' ? { valid: true, context: binding } : { project: portableProject() });
    await pending;
    assert.equal(request.options.signal.aborted, true);
    assert.equal(ui.element('cabiln-input').value, 'A-K');
    assert.equal(ui.downloads.length, 0);
    assert.match(ui.element('project-status').textContent, /changed during the project check/);
  }
});

test('a project file still being read cannot override a reference edit', async () => {
  const ui = page('builder.js');
  let finishRead;
  const pending = ui.element('project-upload').dispatchEvent({ type: 'change', target: { files: [{
    size: 100, text: () => new Promise(resolve => { finishRead = resolve; }),
  }] } });
  await ui.input('smiles-input', 'NEW REFERENCE');
  finishRead(JSON.stringify(portableProject()));
  await pending;
  assert.equal(latestRequest(ui, '/validate_project'), undefined);
  assert.equal(ui.element('smiles-input').value, 'NEW REFERENCE');
});

test('MOL originals survive normalized rendering and rejected reference input', async () => {
  for (const response of [{ svg: '<svg/>', smiles: 'NCC(=O)O' }, { error: 'Unsupported enhanced stereochemistry' }]) {
    const ui = page('builder.js');
    const original = 'EXACT ORIGINAL\r\nMOL WITH METADATA\n';
    const pending = ui.element('mol-upload').dispatchEvent({ type: 'change', target: { files: [{
      name: 'original.mol', text: async () => original,
    }] } });
    await new Promise(setImmediate);
    latestRequest(ui, '/render_mol').resolve(response);
    await pending;
    const snapshot = JSON.parse(ui.run('JSON.stringify(projectSnapshot())'));
    assert.equal(snapshot.reference.original.content, original);
    assert.equal(snapshot.reference.original.name, 'original.mol');
    assert.equal(snapshot.reference.text, response.smiles || '');
  }
});

test('normalizing a pending MOL reference cancels a save of the previous reference state', async () => {
  const ui = page('builder.js');
  const rendering = ui.run('renderMolReference("ORIGINAL", "original.mol")');
  const saving = ui.element('btn-project-save').click();
  latestRequest(ui, '/render_mol').resolve({ svg: '<svg/>', smiles: 'NCC(=O)O' });
  await rendering;
  latestRequest(ui, '/prepare_project').resolve({ project: portableProject() });
  await saving;
  assert.equal(ui.downloads.length, 0);
  assert.equal(ui.element('smiles-input').value, 'NCC(=O)O');
});

test('restored original text controls reference rendering and survives normalization', async () => {
  const ui = page('builder.js');
  const project = portableProject();
  project.reference = { text: 'NORMALIZED OR DIFFERENT',
    original: { kind: 'text', content: '  NC(C)C(=O)O  ', name: '' } };
  const pending = uploadProject(ui, project);
  await new Promise(setImmediate);
  latestRequest(ui, '/validate_project').resolve({ valid: true, context: binding });
  await pending;
  assert.equal(ui.element('smiles-input').value, project.reference.original.content);
  ui.run('restoreReferenceDrawing()');
  assert.equal(JSON.parse(latestRequest(ui, '/render_reference').options.body).input, 'NC(C)C(=O)O');
  latestRequest(ui, '/render_reference').resolve({ svg: '<svg/>', smiles: 'CC(N)C(=O)O', format: 'SMILES' });
  await new Promise(setImmediate);
  assert.equal(ui.run('lastSmiles'), 'CC(N)C(=O)O');
  assert.equal(ui.run('referenceOriginal.content'), project.reference.original.content);
});

test('clearing browser storage preserves work and Undo, cancels recovery, and stays cleared until editing', async () => {
  const project = portableProject();
  const ui = page('builder.js', { storedDraft: { version: 1, document: project.document,
    drafts: project.drafts, reference: project.reference.text,
    reference_original: project.reference.original, context: project.context } });
  assert.equal(ui.element('btn-restore-draft').hidden, false);
  const restoring = ui.element('btn-restore-draft').click();
  const validation = latestRequest(ui, '/prepare_project');
  await ui.element('btn-clear-draft').click();
  validation.resolve({ project });
  await restoring;
  assert.equal(ui.element('cabiln-input').value, '');
  assert.equal(ui.storage.has('cabiln.draft.v1'), false);
  ui.run(`acceptDocumentContext(${JSON.stringify(binding)})`);
  await ui.run('window.dispatchEvent(new Event("pagehide"))');
  assert.equal(ui.storage.has('cabiln.draft.v1'), false);
  await ui.input('cabiln-input', 'A-G');
  const history = ui.run('editor.past.length');
  ui.run('saveDraft()');
  assert.equal(ui.storage.has('cabiln.draft.v1'), true);
  await ui.element('btn-clear-draft').click();
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  assert.equal(ui.run('editor.past.length'), history);
});

test('help is keyboard dismissible and production script cannot force reload on server changes', async () => {
  const ui = page('builder.js');
  await ui.element('btn-help').click();
  assert.equal(ui.element('help-panel').hidden, false);
  assert.equal(ui.element('btn-help').getAttribute('aria-expanded'), 'true');
  await ui.run('window.dispatchEvent({type: "keydown", key: "Escape"})');
  assert.equal(ui.element('help-panel').hidden, true);
  const source = fs.readFileSync(path.join(__dirname, '../src/pyPept/web/static/builder.js'), 'utf8');
  assert.doesNotMatch(source, /location\.reload|\/server_id/);
});

test('a busy calculation retries after Retry-After and uses the same current document', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doRenderCabiln("A-G")');
  latestRequest(ui, '/render').resolve({ error: 'Chemistry workers are busy' }, false, 503, { 'Retry-After': '1' });
  await new Promise(setImmediate);
  assert.equal(ui.requests.length, 1);
  await ui.timers();
  assert.equal(ui.requests.length, 2);
  assert.equal(ui.requests[1].options.body, ui.requests[0].options.body);
  ui.requests[1].resolve(rendered);
  await pending;
  assert.equal(ui.element('cabiln-input').className, 'ok');
});

test('editing during a busy retry cancels the wait without submitting obsolete chemistry', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doRenderCabiln("A-G")');
  latestRequest(ui, '/render').resolve({ error: 'Busy' }, false, 503, { 'Retry-After': '1' });
  await new Promise(setImmediate);
  await ui.input('cabiln-input', '');
  await pending;
  await ui.timers();
  assert.equal(ui.requests.filter(request => request.url === '/render').length, 1);
  assert.equal(ui.element('cabiln-input').value, '');
});

test('persistent overload is visible after two retries and other errors are never retried', async () => {
  for (const [status, header, attempts] of [[503, '1', 3], [503, '30', 1], [503, '', 1], [500, '1', 1]]) {
    const ui = page('builder.js');
    ui.element('cabiln-input').value = 'A-G';
    const pending = ui.run('doRenderCabiln("A-G")');
    for (let attempt = 0; attempt < attempts; attempt++) {
      latestRequest(ui, '/render').resolve({ error: 'Server remains busy; try again shortly' }, false, status, { 'Retry-After': header });
      await new Promise(setImmediate);
      if (attempt < attempts - 1) await ui.timers();
    }
    await pending;
    assert.equal(ui.requests.length, attempts);
    assert.match(ui.element('cabiln-status').textContent, /Server remains busy/);
    assert.equal(ui.element('cabiln-input').value, 'A-G');
  }
});

test('reaction loading failure stays visible and does not poison the retry cache', async () => {
  const ui = page('builder.js');
  let pending = ui.run('loadReactions()');
  latestRequest(ui, '/reactions').resolve({ error: 'Server remains busy' }, false, 503);
  await pending;
  assert.equal(ui.run('reactionPairs'), null);
  assert.equal(ui.element('lib-status').hidden, false);
  assert.match(ui.element('lib-status').textContent, /Server remains busy/);
  pending = ui.run('loadReactions()');
  latestRequest(ui, '/reactions').resolve([['backbone_n', 'backbone_c']]);
  await pending;
  assert.equal(ui.run('reactionPairs.length'), 1);
  assert.equal(ui.element('lib-status').hidden, true);
});

test('ordinary browser drafts are prepared against their saved binding before recovery', async () => {
  const project = portableProject();
  for (const mismatch of [false, true]) {
    const ui = page('builder.js', { storedDraft: { version: 1, document: project.document,
      drafts: project.drafts, reference: project.reference.text,
      reference_original: project.reference.original, context: binding } });
    const restoring = ui.element('btn-restore-draft').click();
    assert.equal(latestRequest(ui, '/validate_project'), undefined);
    const prepare = latestRequest(ui, '/prepare_project');
    if (mismatch) prepare.resolve({ error: 'Library definitions changed. Use the original library.' }, false);
    else prepare.resolve({ project });
    await restoring;
    assert.equal(ui.element('cabiln-input').value, mismatch ? '' : imported.cabiln);
    if (mismatch) {
      assert.equal(ui.element('btn-restore-draft').hidden, false);
      assert.match(ui.element('project-status').textContent, /Library definitions changed.*unchanged/);
    } else assert.equal(ui.element('import-quality').hidden, false);
  }
});

test('foreign input renders bind the current document and its recoverable notation draft', async () => {
  const ui = page('builder.js');
  ui.element('notation-select').value = 'biln';
  await ui.input('cabiln-input', 'A-G');
  const pending = ui.run('doRenderForeign("A-G")');
  assert.equal(JSON.parse(latestRequest(ui, '/render_reference').options.body).input_format, 'biln');
  latestRequest(ui, '/render_reference').resolve({ ...rendered, format: 'BILN', context: binding });
  await pending;
  const saved = JSON.parse(ui.storage.get('cabiln.draft.v1'));
  assert.deepEqual(saved.document.context, binding);
  assert.deepEqual(saved.drafts.biln.context, binding);
  const checking = ui.element('btn-project-save').click();
  const submitted = JSON.parse(latestRequest(ui, '/prepare_project').options.body).project;
  assert.deepEqual(submitted.document.context, binding);
  latestRequest(ui, '/prepare_project').resolve({ error: 'Library changed; render the source again' }, false);
  await checking;
  assert.equal(ui.downloads.length, 0);
});

test('a main-only library change preserves the reference binding and saved chemistry proof', async () => {
  const ui = page('builder.js');
  const oldContext = { ...binding, resolutions: { reference: { source: 'A', resolved: 'original-alanine' } } };
  const current = { ...binding, library_binding: { ...binding.library_binding, monomers: 'new-library' } };
  ui.run(`commitDocument('G', 'cabiln', '', {context: ${JSON.stringify(binding)}});
    projectContext = ${JSON.stringify(oldContext)};
    referenceOriginal = {kind: 'text', content: 'A', name: ''}; smilesInput.value = 'A'; lastSmiles = 'ORIGINAL';`);
  latestRequest(ui, '/render').resolve({ ...rendered, context: current });
  await new Promise(setImmediate);
  const snapshot = JSON.parse(ui.run('JSON.stringify(projectSnapshot())'));
  assert.deepEqual(snapshot.document.context, current);
  assert.deepEqual(snapshot.context, current);
  assert.deepEqual(snapshot.reference.context, oldContext);
  assert.equal(ui.run('lastSmiles'), 'ORIGINAL');
  const checking = ui.element('btn-project-save').click();
  const submitted = JSON.parse(latestRequest(ui, '/prepare_project').options.body).project;
  assert.deepEqual(submitted.reference.context, oldContext);
  latestRequest(ui, '/prepare_project').resolve({ error: 'The reference library definition changed' }, false);
  await checking;
  assert.equal(ui.downloads.length, 0);
  assert.equal(ui.element('smiles-input').value, 'A');
});

test('reference rendering records its own context and reference edits invalidate only its proof', async () => {
  const ui = page('builder.js');
  await ui.input('smiles-input', 'A-G');
  const pending = ui.run('doRenderRef("A-G")');
  latestRequest(ui, '/render_reference').resolve({ ...rendered, smiles: 'CHEMISTRY', format: 'BILN', context: binding });
  await pending;
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(referenceContext)')), binding);
  const full = { ...binding, resolutions: { reference: { source: 'A-G', resolved: 'CHEMISTRY' } } };
  ui.run(`referenceContext = ${JSON.stringify(full)}; acceptReferenceContext(${JSON.stringify(binding)})`);
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(referenceContext)')), full);
  await ui.input('smiles-input', 'A-K');
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(referenceContext)')), binding);
  ui.run('saveDraft()');
  assert.deepEqual(JSON.parse(ui.storage.get('cabiln.draft.v1')).reference_context, binding);
  await ui.input('smiles-input', '');
  assert.equal(ui.run('referenceContext'), null);
});

test('validated project proof upgrades reference context and is retained by later matching renders', async () => {
  const ui = page('builder.js');
  const project = portableProject();
  const current = { ...binding, library_binding: { ...binding.library_binding, monomers: 'unrelated-addition' },
    resolutions: { reference: { source: 'original', resolved: 'same-reference' } } };
  const pending = uploadProject(ui, project);
  await new Promise(setImmediate);
  latestRequest(ui, '/validate_project').resolve({ valid: true, context: current });
  await pending;
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(referenceContext)')), current);
  const rendering = ui.run('renderMolReference(referenceOriginal.content, referenceOriginal.name)');
  latestRequest(ui, '/render_mol').resolve({ svg: '<svg/>', smiles: project.reference.text, context: current });
  await rendering;
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(referenceContext)')), current);
  assert.equal(ui.run('referenceOriginal.content'), project.reference.original.content);
});
