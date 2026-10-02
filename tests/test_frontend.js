// Run with: node --test tests/test_frontend.js
// The VM runs the shipped scripts. Deferred fetches deliberately ignore aborts
// to cover responses already queued when an input is edited.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const test = require('node:test');
const { page, flush, openBuilder, loadSelection, chooseSite, prepareConnection } = require('./_frontend_harness');

const rendered = {
  svg: '<svg>OLD STRUCTURE</svg>', mol_block: 'OLD MOL', info: 'old',
  residue_map: {}, residues: [], chains: [], bracket_groups: [], crosslink_groups: [],
};
const glycineDrawing = {
  svg: '<svg>GLYCINE</svg>', format: 'SMILES', smiles: 'NCC(=O)O', info: '5 atoms',
};
const preview = {
  svg: '<svg>PREVIEW</svg>', chuckles: '[1*]NCC([2*])=O',
  chem_types: { 1: 'backbone_n', 2: 'backbone_c' },
  leaving: { 1: '[H]', 2: '[OH]' },
};

function tabs(ui, symbols, layout, crosslinks = []) {
  const residues = symbols.map((abbr, idx) => ({ abbr, idx }));
  const atoms = Object.fromEntries(symbols.map((_, idx) => [idx, [idx * 2, idx * 2 + 1]]));
  ui.run(`residueView.render({residue_map: ${JSON.stringify(atoms)}, residues: ${JSON.stringify(residues)}, layout: ${JSON.stringify(layout)}, crosslink_groups: ${JSON.stringify(crosslinks)}}, editor.present.quality)`);
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
  ui.run('let hovered; residueView.highlightGroup = ids => { hovered = ids; }');
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
  assert.equal(ui.run('residueView.insertionAnchor(1, 2)'), 1);
  assert.equal(ui.run('residueView.insertionAnchor(0, 2)'), null);
  assert.equal(ui.run('residueView.insertionAnchor(3, 2)'), null);
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
  ui.run('let hovered; residueView.highlightGroup = ids => { hovered = ids; }');
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

test('library metadata waits for the current palette before using the shared worker', async () => {
  const ui = page('builder.js');
  ui.run('library.open()');
  assert.deepEqual(ui.requests.map(request => request.url), ['/monomers']);

  ui.run("window.dispatchEvent(new Event('focus'))");
  const latest = ui.requests[1];
  ui.requests[0].resolve([]);
  await new Promise(setImmediate);
  assert.deepEqual(ui.requests.map(request => request.url), ['/monomers', '/monomers']);

  latest.resolve([{ abbr: 'G', name: 'Glycine', type: 'aa', subtype: 'natural', chem_types: '' }]);
  await new Promise(setImmediate);
  assert.deepEqual(ui.requests.map(request => request.url), ['/monomers', '/monomers', '/reactions']);
  ui.requests[2].resolve([['backbone_c', 'backbone_n']]);
  await new Promise(setImmediate);
  assert.equal(ui.run('library.monomers[0].abbr'), 'G');
  assert.equal(ui.run('library.reactionPairs[0].join(",")'), 'backbone_c,backbone_n');
});

test('reopening the palette discovers new monomers and preserves its search', async () => {
  const ui = page('builder.js');
  const gly = { abbr: 'G', name: 'Glycine', type: 'aa', subtype: 'natural', chem_types: '' };
  const fresh = { ...gly, abbr: 'NewThiol', name: 'New thiol' };
  const first = ui.run('library.load()');
  ui.requests[0].resolve([gly]);
  await first;
  ui.element('lib-search').value = 'new';
  ui.run('library.open()');
  const refresh = ui.requests.find((r, i) => i > 0 && r.url === '/monomers');
  assert.ok(refresh);
  refresh.resolve([gly, fresh]);
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(ui.element('lib-search').value, 'new');
  assert.match(ui.element('lib-list').innerHTML, /NewThiol/);
  assert.equal(ui.run('library.monomers.length'), 2);
});

test('older library refresh cannot replace a newer library', async () => {
  const ui = page('builder.js');
  const old = ui.run('library.load()');
  ui.run("library.panel.classList.add('open')");
  ui.run("window.dispatchEvent(new Event('focus'))");
  ui.requests[1].resolve([{ abbr: 'New', name: 'New', type: 'aa', subtype: '', chem_types: '' }]);
  await new Promise(setImmediate);
  ui.requests[0].resolve([]);
  await old;
  assert.equal(ui.run('library.monomers[0].abbr'), 'New');
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
  ui.run('library.startPreview("G", previewRow)');
  await ui.timers();
  ui.run('library.hidePreview()');
  ui.requests[0].resolve({svg: '<svg>OLD PREVIEW</svg>'});
  await new Promise(setImmediate);
  assert.equal(ui.element('lib-preview').style.display, 'none');
  assert.doesNotMatch(ui.element('lib-preview').innerHTML, /OLD PREVIEW/);
});

test('an older library hover cannot replace the latest preview', async () => {
  const ui = page('builder.js');
  previewRow(ui);
  ui.run('library.startPreview("G", previewRow)');
  await ui.timers();
  ui.run('library.startPreview("A", previewRow)');
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
  ui.run('library.previewCache.G = {svg: "<svg>CACHED PREVIEW</svg>"}; library.startPreview("G", previewRow)');
  ui.run('library.hidePreview()');
  await ui.timers();
  assert.equal(ui.requests.length, 0);
  assert.equal(ui.element('lib-preview').style.display, 'none');
});

test('a discovered cap family asks for its attachment form before loading slots', async () => {
  const ui = page('builder.js');
  await openBuilder(ui);
  ui.run("library.monomers = [{abbr:'NovelCap', degenerate:true, nterm_abbr:'NovelCap_', cterm_abbr:'_NovelCap'}]");
  ui.run("build.useMonomer('NovelCap')");
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
  const note = 'Converted legacy positional notation to bracket form.';
  ui.requests[0].resolve({ ...rendered, normalized_cabiln: normalized, cabiln_echo: normalized,
    normalization_note: note, warnings: [] });
  await pending;
  assert.equal(ui.element('cabiln-input').value, normalized);
  assert.equal(ui.run('drawing.cabiln'), normalized);
  assert.equal(ui.element('cabiln-status').textContent, 'oldNotation normalized · Details' + note);
  assert.equal(ui.element('cabiln-status').className, 'statusbar ok');
  assert.match(ui.element('cabiln-status').innerHTML, /<details.*Notation normalized/);
});

test('clearing the main input discards an already pending render and exports', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doRenderCabiln("A-G")');
  await ui.input('cabiln-input', '');
  ui.requests[0].resolve(rendered);
  await pending;
  assert.doesNotMatch(ui.element('render-inner').innerHTML, /OLD STRUCTURE/);
  assert.equal(ui.run('drawing.cabiln'), '');
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
  await openBuilder(ui, 'K');
  const pending = ui.run('build.selectResidue({abbr:"K", idx:0})');
  await ui.element('build-close').click();
  ui.requests[0].resolve({ svg: '<svg>OLD</svg>', rgroups: [{ slot: 4 }] });
  await pending;
  assert.equal(ui.run('build.filters.left'), null);
  assert.equal(ui.element('build-left-abbr').textContent, '—');
});

test('notation conversion cannot replace an input edited during the request', async () => {
  for (const edited of ['K-C', ' A-G ']) {
    const ui = page('builder.js');
    await ui.input('cabiln-input', 'A-G');
    const pending = ui.run('convertNotation("bracket")');
    await new Promise(setImmediate);
    assert.deepEqual(ui.requests.map(request => request.url), ['/convert_notation']);
    await ui.input('cabiln-input', edited);
    ui.requests[0].resolve({ result: 'OLD CONVERSION' });
    await pending;
    assert.equal(ui.element('cabiln-input').value, edited);
    assert.equal(ui.requests.length, 1);
    await ui.timers();
    assert.equal(ui.requests.filter(request => request.url === '/render').length, 1);
    assert.equal(JSON.parse(latestRequest(ui, '/render').options.body).cabiln, edited.trim());
  }
});

test('canonical formatting is explicit and the original spelling remains undoable', async () => {
  for (const [original, result] of [
    ['Ala-Gly', 'A-G'],
    ['D.(4,1)-G%G-A', 'D.[G(4,1).[A(2,1)]]-G'],
  ]) {
    const ui = page('builder.js');
    await ui.input('cabiln-input', original);
    ui.element('notation-policy').value = 'canonical';
    const pending = ui.element('btn-to-bracket').click();
    await new Promise(setImmediate);
    await ui.timers();
    assert.deepEqual(ui.requests.map(request => request.url), ['/convert_notation']);
    const request = latestRequest(ui, '/convert_notation');
    assert.deepEqual(JSON.parse(request.options.body), {
      cabiln: original, target: 'bracket', canonical: true,
    });
    request.resolve({ result });
    await pending;
    assert.equal(ui.element('cabiln-input').value, result);
    assert.deepEqual(ui.requests.map(request => request.url), ['/convert_notation', '/render']);
    assert.equal(JSON.parse(latestRequest(ui, '/render').options.body).cabiln, result);
    latestRequest(ui, '/render').resolve(rendered);
    await new Promise(setImmediate);
    ui.run('travelHistory("undo")');
    assert.equal(ui.element('cabiln-input').value, original);
    ui.run('travelHistory("redo")');
    assert.equal(ui.element('cabiln-input').value, result);
  }
});

function notationDrawing() {
  return { ...rendered, context: binding, cabiln_echo: 'A%G',
    residue_map: { 0: [0, 1], 1: [2, 3, 4] },
    residues: [{ idx: 0, abbr: 'A' }, { idx: 1, abbr: 'G' }],
    layout: { segments: [{ roots: [0], members: [0] }, { roots: [1], members: [1] }],
      groups: [], markers: [] },
  };
}

function notationResult() {
  return { result: 'G%A', source_echo: 'A%G', occurrence_order: [1, 0], context: binding,
    presentation: { cabiln_echo: 'G%A', warnings: ['Target warning'],
      residues: [{ idx: 0, abbr: 'G' }, { idx: 1, abbr: 'A' }],
      layout: { segments: [{ roots: [0], members: [0] }, { roots: [1], members: [1] }],
        groups: [], markers: [] }, crosslink_groups: [] },
  };
}

async function drawNotation(ui) {
  await ui.input('cabiln-input', 'A%G');
  await ui.timers();
  latestRequest(ui, '/render').resolve(notationDrawing());
  await new Promise(setImmediate);
}

test('validated notation conversion keeps the drawing, remaps owners and preserves Undo', async () => {
  const ui = page('builder.js');
  await drawNotation(ui);
  await ui.element('btn-reroll').click();
  latestRequest(ui, '/render').resolve(notationDrawing());
  await new Promise(setImmediate);
  assert.equal(ui.run('drawing.seed'), 2);
  ui.element('render-inner').style.transform = 'scale(1.5)';
  const beforeSelection = [...ui.requests];
  await openBuilder(ui, 'A%G');
  await loadSelection(ui, 'left', 'G', 1);
  const result = notationResult();
  result.canonical = { ...binding.canonical, binding: binding.library_binding };
  ui.element('notation-policy').value = 'canonical';
  const pending = ui.element('btn-to-branch').click();
  await new Promise(setImmediate);
  latestRequest(ui, '/convert_notation').resolve(result);
  await pending;
  assert.deepEqual([...beforeSelection, ...ui.requests].map(request => request.url).filter(url => !url.startsWith('/monomer_rgroups')), ['/render', '/render', '/convert_notation']);
  assert.equal(ui.element('cabiln-input').value, result.result);
  assert.equal(ui.run('drawing.cabiln'), result.result);
  assert.equal(ui.run('drawing.svg'), rendered.svg);
  assert.equal(ui.run('drawing.molBlock'), rendered.mol_block);
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(residueView.atoms)')), { 0: [2, 3, 4], 1: [0, 1] });
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(residueView.atomOwners)')), { 0: 1, 1: 1, 2: 0, 3: 0, 4: 0 });
  assert.equal(ui.run('residueView.residues[0].abbr'), 'G');
  assert.equal(ui.run('build.selection[0]'), null);
  assert.equal(ui.run('drawing.seed'), 2);
  assert.equal(ui.element('render-inner').style.transform, 'scale(1.5)');
  assert.equal(ui.element('btn-mol').disabled, false);
  assert.equal(ui.element('btn-reroll').disabled, false);
  assert.equal(ui.element('cabiln-input').className, 'ok');
  assert.match(ui.element('cabiln-status').textContent, /Target warning/);
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(editor.present.canonical)')), result.canonical);
  ui.element('notation-policy').value = 'layout';
  const second = ui.element('btn-to-bracket').click();
  await new Promise(setImmediate);
  latestRequest(ui, '/convert_notation').resolve({ ...notationResult(), result: 'A%G',
    source_echo: 'G%A', presentation: { ...notationResult().presentation,
      cabiln_echo: 'A%G', residues: notationDrawing().residues } });
  await second;
  assert.deepEqual([...beforeSelection, ...ui.requests].map(request => request.url).filter(url => !url.startsWith('/monomer_rgroups')), ['/render', '/render', '/convert_notation', '/convert_notation']);
  assert.equal(ui.run('drawing.molBlock'), rendered.mol_block);
  assert.deepEqual(JSON.parse(ui.run('JSON.stringify(residueView.atoms)')), notationDrawing().residue_map);
  assert.equal(ui.run('drawing.seed'), 2);
  ui.run('travelHistory("undo")');
  assert.equal(ui.element('cabiln-input').value, 'G%A');
  assert.equal(JSON.parse(latestRequest(ui, '/render').options.body).cabiln, 'G%A');
  assert.equal(ui.element('btn-mol').disabled, true);
});

test('notation conversion redraws when drawing correspondence is unavailable or no longer current', async () => {
  for (const scenario of ['old-server', 'source', 'target', 'mapping', 'missing-owner',
    'library', 'aliases', 'reactions', 'caps', 'convention', 'resize', 'reroll']) {
    const ui = page('builder.js');
    await drawNotation(ui);
    const pending = ui.element('btn-to-branch').click();
    await new Promise(setImmediate);
    const result = notationResult();
    if (scenario === 'old-server') delete result.presentation;
    if (scenario === 'source') result.source_echo = 'G%A';
    if (scenario === 'target') result.presentation.cabiln_echo = 'A%G';
    if (scenario === 'mapping') result.occurrence_order = [0, 0];
    if (scenario === 'missing-owner') ui.run('delete residueView.atoms[1]');
    if (['library', 'aliases', 'reactions', 'caps'].includes(scenario)) result.context = {
      ...binding, library_binding: { ...binding.library_binding,
        [scenario === 'library' ? 'monomers' : scenario]: 'changed' },
    };
    if (scenario === 'convention') result.context = { ...binding, canonical: { ...binding.canonical, format: 'changed' } };
    if (scenario === 'resize') ui.element('render-canvas').clientWidth = 900;
    if (scenario === 'reroll') {
      await ui.element('btn-reroll').click();
      latestRequest(ui, '/render').resolve({ ...notationDrawing(), svg: '<svg>NEW LAYOUT</svg>' });
      await new Promise(setImmediate);
    }
    latestRequest(ui, '/convert_notation').resolve(result);
    await pending;
    assert.equal(ui.requests.at(-1).url, '/render', scenario);
    assert.equal(JSON.parse(ui.requests.at(-1).options.body).cabiln, result.result, scenario);
    assert.equal(ui.element('btn-mol').disabled, true, scenario);
  }
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

test('layout errors retain the document without persisting as scientific warnings', async () => {
  const ui = page('builder.js');
  const original = 'K.{G(4,2).ac(1,2)}-A';
  await ui.input('cabiln-input', original);
  await ui.timers();
  const history = ui.run('editor.past.length');
  let pending = ui.element('btn-to-branch').click();
  assert.equal(latestRequest(ui, '/convert_notation'), undefined);
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
  const request = ui.requests.find(request => request.url === '/convert_notation');
  assert.equal(JSON.parse(request.options.body).canonical, false);
  request.resolve({ error: 'Cannot preserve the requested attachment' }, false);
  await pending;
  assert.equal(ui.element('cabiln-input').value, original);
  assert.equal(ui.element('cabiln-status').textContent, 'Cannot preserve the requested attachment');
  assert.equal(ui.element('cabiln-status').title, 'Cannot preserve the requested attachment');
  assert.equal(ui.element('conversion-status').hidden, true);
  assert.equal(ui.run('editor.past.length'), history);
  ui.run('saveDraft()');
  const saved = JSON.parse(ui.storage.get('cabiln.draft.v1'));
  assert.equal(saved.document.warning, '');
  assert.equal(ui.run('projectSnapshot().document.warning'), '');
  const restored = page('builder.js', { storedDraft: saved });
  const restoring = restored.element('btn-restore-draft').click();
  const validation = latestRequest(restored, '/prepare_project');
  validation.resolve({ project: JSON.parse(validation.options.body).project });
  await restoring;
  assert.equal(restored.element('cabiln-input').value, original);
  assert.equal(restored.element('conversion-status').hidden, true);
  assert.doesNotMatch(restored.element('cabiln-status').textContent, /Cannot preserve/);

  pending = ui.element('btn-to-branch').click();
  await new Promise(setImmediate);
  latestRequest(ui, '/convert_notation').resolve({ result: original, context: binding });
  await pending;
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
  assert.doesNotMatch(ui.element('cabiln-status').textContent, /Cannot preserve/);
  assert.doesNotMatch(ui.element('cabiln-status').title, /Cannot preserve/);
  assert.equal(ui.element('conversion-status').hidden, true);
  ui.run('saveDraft()');
  assert.equal(JSON.parse(ui.storage.get('cabiln.draft.v1')).document.warning, '');
});

test('format conversion uses a source normalized by its pending drawing', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'Ala-Gly');
  await ui.timers();
  const history = ui.run('editor.past.length');
  const converting = ui.element('btn-to-bracket').click();
  assert.equal(latestRequest(ui, '/convert_notation'), undefined);
  latestRequest(ui, '/render').resolve({ ...rendered, normalized_cabiln: 'A-G' });
  await new Promise(setImmediate);
  const request = latestRequest(ui, '/convert_notation');
  assert.deepEqual(JSON.parse(request.options.body), { cabiln: 'A-G', target: 'bracket', canonical: false });
  request.resolve({ result: 'A-G' });
  await converting;
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  assert.equal(ui.run('editor.past.length'), history);
});

test('failed immediate formatting draws pending input without hiding its error', async () => {
  for (const failure of ['server', 'network', 'replacement']) {
    const ui = page('builder.js');
    await ui.input('cabiln-input', 'G');
    await ui.timers();
    latestRequest(ui, '/render').resolve(rendered);
    await new Promise(setImmediate);
    await ui.element('btn-verify').click();
    ui.run("lastSmiles = 'NCC(=O)O'");
    const original = 'D.(4,1)-G%G-A';
    const normalized = 'D.[G(4,1).A(2,1)]-G';
    await ui.input('cabiln-input', original);
    const history = ui.run('editor.past.length');
    let converting = ui.element('btn-to-bracket').click();
    await new Promise(setImmediate);
    assert.equal(ui.requests.at(-1).url, '/convert_notation');
    if (failure === 'replacement') {
      const obsolete = converting;
      const oldRequest = latestRequest(ui, '/convert_notation');
      converting = ui.element('btn-to-branch').click();
      await new Promise(setImmediate);
      oldRequest.resolve({ result: 'OLD CONVERSION' });
      await obsolete;
      assert.equal(ui.requests.filter(request => request.url === '/render').length, 1);
      assert.equal(ui.element('cabiln-input').value, original);
    }
    const request = latestRequest(ui, '/convert_notation');
    const error = failure === 'network' ? 'Notation conversion failed. Try again.' : 'Cannot preserve this layout';
    if (failure === 'network') request.reject(new Error('Connection lost'));
    else request.resolve({ error }, false);
    await new Promise(setImmediate);
    const drawing = latestRequest(ui, '/render');
    assert.equal(JSON.parse(drawing.options.body).cabiln, original);
    assert.match(ui.element('render-inner').innerHTML, /OLD STRUCTURE/);
    assert.equal(ui.element('residue-chips').getAttribute('aria-disabled'), 'true');
    assert.equal(ui.element('btn-mol').disabled, true);
    drawing.resolve({ ...rendered, normalized_cabiln: normalized });
    await converting;
    assert.equal(ui.element('cabiln-input').value, normalized);
    assert.equal(ui.element('cabiln-input').className, 'ok');
    assert.equal(ui.element('cabiln-status').textContent, error);
    assert.equal(ui.element('cabiln-status').title, error);
    assert.equal(ui.element('btn-mol').disabled, false);
    assert.equal(ui.run('editor.past.length'), history);
    assert.equal(ui.run('editor.present.warning'), '');
    assert.equal(JSON.parse(latestRequest(ui, '/verify').options.body).cabiln, normalized);
    latestRequest(ui, '/verify').resolve({ match: false });
  }
});

test('immediate formatting waits for the reference and verifies the final drawing', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'Ala-Gly');
  await ui.element('btn-verify').click();
  await ui.input('smiles-input', 'NCC(=O)O');
  const converting = ui.element('btn-to-bracket').click();
  await new Promise(setImmediate);
  assert.deepEqual(ui.requests.map(request => request.url), ['/render_reference']);
  ui.requests[0].resolve(glycineDrawing);
  await new Promise(setImmediate);
  assert.deepEqual(ui.requests.map(request => request.url), ['/render_reference', '/convert_notation']);
  latestRequest(ui, '/convert_notation').resolve({ result: 'A-G' });
  await converting;
  assert.equal(ui.requests.filter(request => request.url === '/render').length, 1);
  latestRequest(ui, '/render').resolve(rendered);
  await new Promise(setImmediate);
  assert.deepEqual(JSON.parse(latestRequest(ui, '/verify').options.body), { smiles: 'NCC(=O)O', cabiln: 'A-G' });
  latestRequest(ui, '/verify').resolve({ match: false });
});

test('editing only the reference restores a cancelled formatter\'s deferred main drawing', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'A-G');
  const converting = ui.element('btn-to-bracket').click();
  await new Promise(setImmediate);
  const obsolete = latestRequest(ui, '/convert_notation');
  await ui.input('smiles-input', 'NCC(=O)O');
  assert.equal(obsolete.options.signal.aborted, true);
  await ui.timers();
  const drawing = latestRequest(ui, '/render');
  assert.ok(drawing, 'The unchanged main input still needs its drawing');
  assert.equal(JSON.parse(drawing.options.body).cabiln, 'A-G');
  obsolete.resolve({ result: 'OLD CONVERSION' });
  await converting;
  drawing.resolve(rendered);
  latestRequest(ui, '/render_reference').resolve(glycineDrawing);
  await new Promise(setImmediate);
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  assert.equal(ui.element('cabiln-input').className, 'ok');
});

test('a failed format request preserves scientific warnings and recognition evidence', async () => {
  const ui = page('builder.js');
  setImported(ui);
  const converting = ui.element('btn-to-bracket').click();
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
  const before = ui.run('JSON.stringify(editor.present)');
  const history = ui.run('editor.past.length');
  latestRequest(ui, '/convert_notation').reject(new Error('Connection lost'));
  await converting;
  assert.equal(ui.element('cabiln-status').textContent, 'Notation conversion failed. Try again.');
  assert.equal(ui.element('conversion-status').textContent, imported.warnings[0]);
  assert.equal(ui.run('JSON.stringify(editor.present)'), before);
  assert.equal(ui.run('editor.past.length'), history);
  ui.run('saveDraft()');
  assert.equal(JSON.parse(ui.storage.get('cabiln.draft.v1')).document.warning, imported.warnings[0]);
});

test('render transport failures preserve source without reporting invalid syntax', async () => {
  for (const mode of ['cabiln', 'smiles', 'reference']) {
    for (const failure of ['network', 503, 400]) {
      const ui = page('builder.js');
      const reference = mode === 'reference';
      const input = reference ? 'smiles-input' : 'cabiln-input';
      const status = reference ? 'smiles-status' : 'cabiln-status';
      const url = mode === 'cabiln' ? '/render' : '/render_reference';
      if (mode === 'smiles') ui.element('notation-select').value = 'smiles';
      await ui.input(input, mode === 'cabiln' ? 'A-G' : 'NCC(=O)O');
      await ui.timers();
      latestRequest(ui, url).resolve(mode === 'cabiln' ? rendered : glycineDrawing);
      await new Promise(setImmediate);
      const drawing = ui.element('render-inner').innerHTML;
      const source = mode === 'cabiln' ? 'G-A' : 'CC(=O)O';
      await ui.input(input, source);
      await ui.timers();
      const request = latestRequest(ui, url);
      if (failure === 'network') request.reject(new Error('Connection closed'));
      else request.resolve({ error: failure === 400 ? 'Unsupported input' : 'Chemistry capacity is busy' }, false, failure);
      await new Promise(setImmediate);
      const label = `${mode}/${failure}`;
      assert.equal(ui.element(input).value, source, label);
      assert.equal(ui.element(input).className, failure === 400 ? 'err' : '', label);
      assert.match(ui.element(status).textContent, failure === 'network' ? /input is preserved/ : failure === 400 ? /Unsupported input/ : /capacity is busy/, label);
      assert.equal(ui.element(status).title, ui.element(status).textContent.replace(/ Retry$/, ''), label);
      if (!reference) {
        assert.equal(ui.element('render-inner').innerHTML, drawing, label);
        assert.equal(ui.element('render-pane').getAttribute('aria-busy'), 'false', label);
        assert.match(ui.element('render-progress-label').textContent,
          failure === 400 ? /correct the input/ : /drawing unavailable; try again/, label);
        assert.equal(ui.element('btn-mol').disabled, true, label);
      }
    }
  }
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
  ui.run('var hovered = null; residueView.highlightGroup = ids => { hovered = ids; }');
  const brackets = chips.filter(chip => chip.textContent === '[');
  for (const [index, members] of [[0, [2, 3, 4, 5]], [1, [3, 4]], [2, [4]], [3, [5]]]) {
    await brackets[index].dispatchEvent({ type: 'mouseenter' });
    assert.equal(ui.run('JSON.stringify(hovered)'), JSON.stringify(members));
  }
});

test('editing monomer SMILES invalidates detected chemistry and registration', async () => {
  const ui = page('register.js');
  ui.element('smiles-in').value = 'NCC(=O)O';
  const pending = ui.element('preview-form').dispatchEvent({ type: 'submit' });
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
  const pending = ui.element('preview-form').dispatchEvent({ type: 'submit' });
  await ui.input('smiles-in', 'CC(=O)O');
  ui.requests[0].resolve(preview);
  await pending;
  assert.equal(ui.run('detectedData'), null);
  assert.equal(ui.element('btn-register').disabled, true);
});

test('a busy monomer preview retries, but editing cancels its retry', async () => {
  for (const edited of [false, true]) {
    const ui = page('register.js');
    ui.element('smiles-in').value = 'NCC(=O)O';
    const pending = ui.element('preview-form').dispatchEvent({ type: 'submit' });
    ui.requests[0].resolve({ error: 'Chemistry capacity is busy' }, false, 503, { 'Retry-After': '1' });
    await new Promise(setImmediate);
    assert.equal(ui.element('btn-preview').disabled, true);
    if (edited) await ui.input('smiles-in', 'CC(=O)O');
    await ui.timers();
    if (edited) {
      await pending;
      assert.equal(ui.requests.length, 1);
      assert.equal(ui.run('detectedData'), null);
      assert.equal(ui.element('preview-section').hidden, true);
    } else {
      assert.equal(ui.requests.length, 2);
      assert.equal(ui.requests[1].options.body, ui.requests[0].options.body);
      ui.requests[1].resolve(preview);
      await pending;
      assert.equal(ui.element('smiles-in').className, 'ok');
      assert.equal(ui.element('preview-canvas').innerHTML, preview.svg);
    }
  }
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
  assert.equal(ui.run('drawing.cabiln'), 'G-A');
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
    await ui.element('btn-verify').click();
    ui.run('drawing.accept({}, "A-G", "cabiln"); lastSmiles = "NCC(=O)O"');
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
  for (const side of ['left', 'right']) {
    const ui = page('builder.js');
    await openBuilder(ui, 'G-K-C');
    if (side === 'right') await loadSelection(ui, 'left', 'G', 0);
    ui.requests.length = 0;
    const first = ui.run('build.selectResidue({abbr:"K", idx:1})');
    const second = ui.run('build.selectResidue({abbr:"C", idx:2})');
    ui.requests[1].resolve({ svg: '<svg>C</svg>', rgroups: [{ slot: 4 }] });
    await second;
    ui.requests[0].resolve({ svg: '<svg>K</svg>', rgroups: [{ slot: 5 }] });
    await first;
    assert.equal(ui.element(`build-${side}-abbr`).textContent, 'C');
    assert.match(ui.element(`build-${side}-rgroups`).children[0].textContent, /^R4 /);
    assert.deepEqual(JSON.parse(ui.run('JSON.stringify(build.selection)')), side === 'left' ? [2, null] : [0, 2]);
  }
});

test('a cleared R-group selection cannot be enabled by old bond validation', async () => {
  const ui = page('builder.js');
  await prepareConnection(ui, {source: 'K', left: 'K', right: 'D', leftSlot: 4, rightSlot: 4});
  await chooseSite(ui, 'left', 4);
  ui.requests[0].resolve({ valid: true, reaction: 'amide' });
  await flush();
  assert.equal(ui.element('build-connect').disabled, true);
  assert.equal(ui.element('build-status').textContent, 'Choose a site on the current residue');
});

test('superseded connection previews cannot publish either a queued proposal or drawing', async () => {
  for (const stage of ['proposal', 'drawing']) {
    for (const action of ['site', 'residue', 'close', 'edit', 'library', 'hide']) {
      const ui = page('builder.js');
      await prepareConnection(ui);
      ui.requests[0].resolve({ valid: true, reaction: 'amide' });
      await flush();
      const pending = ui.element('build-preview-button').click();
      let response = ui.requests[1];
      if (stage === 'drawing') {
        response.resolve({ result: 'G-A' });
        await new Promise(setImmediate);
        response = ui.requests[2];
      }
      if (action === 'site') await chooseSite(ui, 'left', 2);
      if (action === 'residue') ui.run('build.useMonomer("C")');
      if (action === 'close') await ui.element('build-close').click();
      if (action === 'edit') await ui.input('cabiln-input', 'G-K');
      if (action === 'hide') await ui.element('build-preview-close').click();
      if (action === 'library') {
        const refresh = ui.run('library.load()');
        ui.requests.at(-1).resolve([], true, 200, { 'X-Library-Version': 'new' });
        await refresh;
      }
      const count = ui.requests.length;
      response.resolve(stage === 'proposal' ? { result: 'G-A' } : { ...rendered, svg: '<svg>STALE PREVIEW</svg>' });
      await pending;
      assert.equal(ui.requests.length, count, `${stage}/${action}: no follow-up render`);
      assert.equal(ui.element('build-preview').hidden, true, `${stage}/${action}`);
      assert.equal(ui.element('build-preview-inner').innerHTML, '');
      assert.equal(ui.element('build-preview').getAttribute('aria-busy'), 'false');
      assert.equal(ui.element('build-connect').disabled, action !== 'hide');
      assert.equal(ui.element('cabiln-input').value, action === 'edit' ? 'G-K' : 'G');
    }
  }
});

test('preview waits for active drawings without cancelling them or accepting their library context', async () => {
  for (const dismiss of [false, true]) {
    const ui = page('builder.js');
    await prepareConnection(ui);
    ui.run("replaceDocument({context: {project_version: 1, library_binding: 'original', canonical: {version: 1}}})");
    ui.requests[0].resolve({ valid: true, reaction: 'amide' });
    await flush();
    await ui.input('smiles-input', 'NCC(=O)O');
    await ui.timers();
    const reference = ui.requests[1];
    assert.equal(reference.url, '/render_reference');
    const pending = ui.element('build-preview-button').click();
    assert.equal(ui.requests.length, 2);
    assert.match(ui.element('build-preview-status').textContent, /Waiting/);
    if (dismiss) await ui.element('build-preview-close').click();
    assert.equal(reference.options.signal.aborted, false);
    reference.resolve(glycineDrawing);
    await new Promise(setImmediate);
    if (!dismiss) {
      assert.equal(ui.requests[2].url, '/insert_bond');
      ui.requests[2].resolve({ result: 'G-A' });
      await new Promise(setImmediate);
      ui.requests[3].resolve({ ...rendered, context: {
        project_version: 1, library_binding: 'changed', canonical: {version: 1},
      }, warnings: ['Preview warning'] });
    }
    await pending;
    if (dismiss) assert.equal(ui.requests.length, 2);
    else {
      assert.match(ui.element('build-preview-status').textContent, /Library changed.*Preview warning/);
      assert.equal(ui.run('editor.present.context.library_binding'), 'original');
    }
    assert.equal(ui.element('build-connect').disabled, false);
  }
});

const swapContext = { project_version: 1, library_binding: 'library', canonical: {version: 1} };
async function swapPage() {
  const ui = page('builder.js');
  ui.element('build-action').value = 'swap';
  await openBuilder(ui, 'A-G');
  await loadSelection(ui, 'left', 'A', 0, [{slot: 2, used: true}]);
  latestRequest(ui, '/replacement_options').resolve({ source_echo: 'A-G', residue_idx: 0,
    context: swapContext, requirements: [{slot: 2}],
    candidates: [{abbr: 'S', mapping: {2: 2}, choices: {2: [2]}, sites: [{slot: 2}]}] });
  await flush();
  await loadSelection(ui, 'right', 'S', null, [{slot: 2}]);
  ui.requests.length = 0;
  return ui;
}

test('Swap requires a matching product preview and rechecks that exact proposal before applying', async () => {
  for (const failure of [null, 'preview-context', 'apply-context', 'apply-product', 'apply-error']) {
    const ui = await swapPage();
    const context = swapContext;
    const changed = {...context, library_binding: 'changed'};
    assert.equal(ui.element('build-connect').disabled, true);
    const previewing = ui.element('build-preview-button').click();
    assert.equal(ui.requests[0].url, '/replace_monomer');
    ui.requests[0].resolve({ result: 'S-G', context });
    await new Promise(setImmediate);
    ui.requests[1].resolve({...rendered, context: failure === 'preview-context' ? changed : context});
    await previewing;
    assert.equal(ui.element('cabiln-input').value, 'A-G');
    if (failure === 'preview-context') {
      assert.equal(ui.element('build-connect').disabled, true);
      assert.match(ui.element('build-preview-status').textContent, /library changed/i);
      continue;
    }
    assert.equal(ui.element('build-connect').disabled, false);
    const applying = ui.element('build-connect').click();
    assert.deepEqual(JSON.parse(ui.requests[2].options.body), JSON.parse(ui.requests[0].options.body));
    ui.requests[2].resolve(failure === 'apply-error' ? {error: 'Library changed'} : {
      result: failure === 'apply-product' ? 'K-G' : 'S-G',
      context: failure === 'apply-context' ? changed : context,
    });
    await applying;
    assert.equal(ui.element('cabiln-input').value, failure ? 'A-G' : 'S-G');
    if (failure) {
      assert.equal(ui.element('build-connect').disabled, true);
      assert.equal(ui.element('build-preview-button').disabled, false);
    }
  }
});

test('cancelled Swap responses cannot restore a mapping, preview or apply an old source', async () => {
  for (const stage of ['options', 'proposal', 'drawing', 'apply']) {
    for (const action of ['edit', 'mode', 'close', 'library', 'mapping']) {
      const ui = await swapPage();
      const context = swapContext;
      let pending;
      if (stage === 'options') {
        pending = ui.run('build.selectResidue({abbr:"A", idx:0})');
        ui.requests[0].resolve({svg: '<svg/>', rgroups: [{slot: 2, used: true}]});
        await pending;
      }
      else {
        pending = ui.element('build-preview-button').click();
        if (stage !== 'proposal') {
          ui.requests[0].resolve({result:'S-G', context});
          await new Promise(setImmediate);
        }
        if (stage === 'apply') {
          ui.requests[1].resolve({...rendered, context});
          await pending;
          pending = ui.element('build-connect').click();
        }
      }
      const response = ui.requests.at(-1);
      if (action === 'edit') await ui.input('cabiln-input', 'K-G');
      if (action === 'mode') {
        ui.element('build-action').value = 'connect';
        await ui.element('build-action').dispatchEvent({ type: 'change' });
      }
      if (action === 'close') await ui.element('build-close').click();
      if (action === 'library') {
        const refresh = ui.run('library.load()');
        ui.requests.at(-1).resolve([], true, 200, { 'X-Library-Version': 'new' });
        await refresh;
      }
      if (action === 'mapping') {
        // Choosing another occurrence cancels both candidate discovery and the proposal.
        ui.run('build.selectResidue({abbr:"G", idx:1})');
      }
      const count = ui.requests.length;
      response.resolve(stage === 'options' ? {source_echo:'A-G', residue_idx:0, requirements:[], candidates:[], context}
        : stage === 'drawing' ? {...rendered, context} : {result:'S-G',context});
      await pending;
      await flush();
      assert.equal(ui.requests.length, count, `${stage}/${action}`);
      assert.equal(ui.element('cabiln-input').value, action === 'edit' ? 'K-G' : 'A-G');
      assert.equal(ui.element('build-connect').disabled, true);
      assert.equal(ui.element('build-preview').hidden, true);
      assert.equal(ui.element('build-preview-source').textContent, '');
    }
  }
});

test('cancelled SMILES conversion restores the button without replacing main input', async () => {
  const ui = page('builder.js');
  ui.element('btn-s2c').textContent = '→ %';
  ui.element('smiles-input').value = 'NCC(=O)O';
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doS2c("percent")');
  await new Promise(setImmediate);
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
  await new Promise(setImmediate);
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

test('an immediate conversion waits for a restored draft drawing and preserves its history', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'K-G');
  await ui.timers();
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
  ui.element('notation-select').value = 'smiles';
  await ui.element('notation-select').dispatchEvent({ type: 'change' });
  await ui.input('cabiln-input', 'NCC(=O)O');
  await ui.timers();
  latestRequest(ui, '/render_reference').resolve({ ...glycineDrawing, context: binding });
  await new Promise(setImmediate);
  ui.element('notation-select').value = 'cabiln';
  await ui.element('notation-select').dispatchEvent({ type: 'change' });
  const obsolete = latestRequest(ui, '/render');
  ui.element('notation-select').value = 'smiles';
  await ui.element('notation-select').dispatchEvent({ type: 'change' });
  const history = ui.run('editor.past.length');
  const converting = ui.element('btn-to-cabiln-pct').click();
  assert.equal(ui.element('conversion-progress').hidden, false);
  for (let attempt = 0; attempt < 2; attempt++) {
    latestRequest(ui, '/render_reference').resolve({ error: 'Chemistry capacity is busy' }, false, 503, { 'Retry-After': '1' });
    await new Promise(setImmediate);
    assert.equal(latestRequest(ui, '/to_cabiln'), undefined);
    assert.equal(ui.element('cabiln-input').value, 'NCC(=O)O');
    await ui.timers();
  }
  latestRequest(ui, '/render_reference').resolve({ ...glycineDrawing, context: binding });
  await new Promise(setImmediate);
  const request = latestRequest(ui, '/to_cabiln');
  assert.deepEqual(JSON.parse(request.options.body), {
    input: 'NCC(=O)O', input_format: 'smiles', notation: 'percent',
  });
  request.resolve(recognizedGlycine);
  await converting;
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  obsolete.resolve({ ...rendered, info: 'OBSOLETE DRAWING' });
  await new Promise(setImmediate);
  assert.equal(ui.element('cabiln-input').value, 'G');
  assert.equal(ui.element('notation-select').value, 'cabiln');
  assert.equal(ui.element('conversion-progress').hidden, true);
  assert.equal(ui.element('btn-to-cabiln-pct').disabled, false);
  assert.equal(ui.run('editor.present.quality.recognition_status'), 'complete');
  assert.equal(ui.run('editor.past.length'), history + 1);
  assert.doesNotMatch(ui.element('cabiln-status').textContent, /busy|OBSOLETE/);
  await ui.element('btn-undo').click();
  assert.equal(ui.element('cabiln-input').value, 'NCC(=O)O');
  assert.equal(ui.element('notation-select').value, 'smiles');
  await ui.element('btn-redo').click();
  assert.equal(ui.element('cabiln-input').value, 'G');
  assert.equal(ui.run('editor.present.quality.recognition_status'), 'complete');
});

test('editing cancels a conversion waiting for an independent reference drawing', async () => {
  const ui = page('builder.js');
  ui.element('notation-select').value = 'smiles';
  await ui.input('cabiln-input', 'NCC(=O)O');
  await ui.input('smiles-input', 'CC(=O)O');
  const converting = ui.element('btn-to-cabiln-pct').click();
  latestRequest(ui, '/render_reference').resolve(glycineDrawing);
  await new Promise(setImmediate);
  const reference = latestRequest(ui, '/render_reference');
  assert.equal(JSON.parse(reference.options.body).input, 'CC(=O)O');
  assert.equal(latestRequest(ui, '/to_cabiln'), undefined);
  await ui.input('cabiln-input', 'CC');
  await converting;
  assert.equal(ui.element('conversion-progress').hidden, true);
  assert.equal(ui.element('btn-to-cabiln-pct').disabled, false);
  assert.equal(reference.options.signal.aborted, false);
  reference.resolve({ svg: '<svg>REFERENCE</svg>', smiles: 'CC(=O)O', format: 'SMILES' });
  await new Promise(setImmediate);
  await ui.timers();
  assert.equal(latestRequest(ui, '/to_cabiln'), undefined);
  assert.equal(ui.element('cabiln-input').value, 'CC');
});

test('a terminal conversion error stays visible after its drawing has settled', async () => {
  const ui = page('builder.js');
  ui.element('notation-select').value = 'smiles';
  await ui.input('cabiln-input', 'NCC(=O)O');
  const converting = ui.element('btn-to-cabiln-bracket').click();
  latestRequest(ui, '/render_reference').resolve({ ...glycineDrawing, context: binding });
  await new Promise(setImmediate);
  for (let attempt = 0; attempt < 5; attempt++) {
    latestRequest(ui, '/to_cabiln').resolve({ error: 'Chemistry capacity is busy; retry shortly' }, false, 503, { 'Retry-After': '1' });
    await new Promise(setImmediate);
    if (attempt < 4) await ui.timers();
  }
  await converting;
  ui.run(`acceptDocumentContext(${JSON.stringify(binding)}); saveDraft()`);
  assert.match(ui.element('cabiln-status').textContent, /capacity is busy/);
  assert.equal(ui.element('cabiln-status').className, 'statusbar');
  assert.equal(ui.element('cabiln-status').title, ui.element('cabiln-status').textContent);
  assert.equal(ui.element('cabiln-input').className, 'ok');
  assert.equal(ui.element('cabiln-input').value, 'NCC(=O)O');
  assert.equal(ui.element('notation-select').value, 'smiles');
  assert.equal(ui.element('conversion-progress').hidden, true);
  assert.equal(ui.element('btn-to-cabiln-bracket').disabled, false);
  assert.equal(JSON.parse(ui.storage.get('cabiln.draft.v1')).document.warning, '');
});

test('reference conversion waits for original MOL rendering and its automatic Verify', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'G');
  await ui.timers();
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
  await ui.element('btn-verify').click();
  const original = 'ORIGINAL MOL\r\nwith exact newlines\n';
  const uploading = ui.element('mol-upload').dispatchEvent({ type: 'change', target: { files: [{
    name: 'original.mol', text: async () => original,
  }] } });
  await new Promise(setImmediate);
  const converting = ui.element('btn-s2c').click();
  assert.equal(latestRequest(ui, '/smiles_to_cabiln'), undefined);
  latestRequest(ui, '/render_mol').resolve({ ...glycineDrawing, context: binding });
  await uploading;
  await new Promise(setImmediate);
  assert.equal(latestRequest(ui, '/smiles_to_cabiln'), undefined);
  const verification = latestRequest(ui, '/verify');
  assert.deepEqual(JSON.parse(verification.options.body), { smiles: 'NCC(=O)O', cabiln: 'G' });
  verification.resolve({ match: true });
  await new Promise(setImmediate);
  const request = latestRequest(ui, '/smiles_to_cabiln');
  assert.deepEqual(JSON.parse(request.options.body), { smiles: 'NCC(=O)O', notation: 'percent' });
  request.resolve(recognizedGlycine);
  await converting;
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
  latestRequest(ui, '/verify').resolve({ match: true });
  await new Promise(setImmediate);
  assert.equal(ui.element('cabiln-input').value, 'G');
  assert.equal(ui.run('referenceOriginal.content'), original);
  assert.equal(ui.run('referenceOriginal.name'), 'original.mol');
  assert.equal(ui.run('referenceContext.library_binding.monomers'), 'library-one');
  ui.run('saveDraft()');
  assert.equal(JSON.parse(ui.storage.get('cabiln.draft.v1')).reference_original.content, original);
  assert.equal(ui.run('editor.present.quality.recognition_status'), 'complete');
});

test('conversion transport errors offer retry and preserve the current sources', async () => {
  for (const reference of [false, true]) {
    const ui = page('builder.js');
    ui.element('notation-select').value = 'smiles';
    ui.element('cabiln-input').value = 'NCC(=O)O';
    ui.element('smiles-input').value = 'CC(=O)O';
    const button = reference ? 'btn-s2c' : 'btn-to-cabiln-pct';
    const converting = ui.element(button).click();
    await new Promise(setImmediate);
    latestRequest(ui, reference ? '/smiles_to_cabiln' : '/to_cabiln').reject(new Error('Connection lost'));
    await converting;
    assert.equal(ui.element(reference ? 'smiles-status' : 'cabiln-status').textContent,
      reference ? 'Could not convert the reference. Try again.' : 'Could not convert the input. Try again.');
    assert.equal(ui.element('cabiln-input').value, 'NCC(=O)O');
    assert.equal(ui.element('smiles-input').value, 'CC(=O)O');
    assert.equal(ui.element(button).disabled, false);
    assert.equal(ui.element('conversion-progress').hidden, true);
    assert.equal(ui.run('editor.present.warning'), '');
  }
});

test('an unreadable MOL reports that peptide input is preserved', async () => {
  const ui = page('builder.js');
  await ui.input('cabiln-input', 'A-G');
  const before = ui.run('JSON.stringify(editor.present)');
  await ui.element('mol-upload').dispatchEvent({ type: 'change', target: { files: [{
    name: 'unreadable.mol', text: async () => { throw new Error('File unavailable'); },
  }] } });
  assert.equal(ui.run('JSON.stringify(editor.present)'), before);
  assert.equal(ui.element('cabiln-input').value, 'A-G');
  assert.equal(ui.element('smiles-status').textContent,
    'Could not read the MOL/SDF file. Your peptide input is preserved.');
});

test('a pending connection cannot overwrite manual sequence edits', async () => {
  const ui = page('builder.js');
  await prepareConnection(ui, {source: 'A-K', left: 'K', right: 'D', leftSlot: 4, rightSlot: 4, index: 1});
  ui.requests[0].resolve({valid: true});
  await flush();
  ui.requests.length = 0;
  const pending = ui.element('build-connect').dispatchEvent({ type: 'click' });
  await ui.input('cabiln-input', 'A-C');
  ui.requests[0].resolve({ result: 'OLD' });
  await pending;
  assert.equal(ui.element('cabiln-input').value, 'A-C');
});

test('cancelling backbone insertion discards the pending insertion result', async () => {
  const ui = page('builder.js');
  await openBuilder(ui, 'A-K');
  ui.run('residueView.render({layout: {segments: [{roots: [0, 1]}]}})');
  await loadSelection(ui, 'left', 'A', 0);
  await loadSelection(ui, 'right', 'K', 1);
  await ui.element('build-insert-btn').click();
  ui.requests.length = 0;
  ui.run('build.useMonomer("G")');
  await ui.element('build-insert-btn').click();
  ui.requests[0].resolve({ result: 'A-G-K' });
  await flush();
  assert.equal(ui.element('cabiln-input').value, 'A-K');
  assert.equal(ui.run('build.filters.insertBetween'), false);
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
  let pending = ui.element('preview-form').dispatchEvent({ type: 'submit' });
  ui.requests[0].resolve(preview);
  await pending;
  pending = ui.element('preview-form').dispatchEvent({ type: 'submit' });
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
  let pending = ui.element('preview-form').dispatchEvent({ type: 'submit' });
  assert.equal(ui.element('preview-section').hidden, false);
  ui.requests[0].resolve(preview);
  await pending;
  assert.equal(ui.element('btn-register').disabled, false);
  ui.element('chuckles-out').value = 'UNRELATED EDIT';
  pending = ui.element('registration-form').dispatchEvent({ type: 'submit' });
  assert.equal(JSON.parse(ui.requests[1].options.body).chuckles, preview.chuckles);
  ui.requests[1].resolve({ total: 1000 });
  await pending;
  assert.equal(ui.element('btn-register').disabled, true);
  assert.equal(ui.element('btn-register').textContent, '✓ Registered');
  await ui.input('smiles-in', 'CC(=O)O');
  await ui.input('abbr-in', 'TestAc');
  pending = ui.element('preview-form').dispatchEvent({ type: 'submit' });
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
const recognizedGlycine = {
  cabiln: 'G', from: 'SMILES', context: binding, recognition_status: 'complete',
  search_complete: true, inferred_stereo: false, synthetic_components: [], warnings: [],
  assignments: [{ symbol: 'G', recognized: true, residue_index: 0,
    source_atoms: [0, 1, 2, 3, 4], attachments: [] }],
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
  latestRequest(ui, '/render_reference').resolve({ ...rendered, format: 'SMILES' });
  await new Promise(setImmediate);
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
  latestRequest(ui, '/render').resolve({ ...rendered, context: binding });
  await new Promise(setImmediate);
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
  ui.run(`residueView.render({residues: ${JSON.stringify(residues)}, layout: {segments: [{roots: [0, 1, 2], members: [0, 1, 2]}], groups: [], markers: []}}, editor.present.quality)`);
  ui.run('let selected = []; build.selectResidue = residue => selected.push(residue.idx)');
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
  const pending = ui.run('library.load()');
  latestRequest(ui, '/monomers').resolve([monomer]);
  await pending;
  assert.match(ui.element('lib-list').innerHTML, /lib-quality.*Review &lt;atom&gt; mapping/);
  assert.match(ui.element('lib-list').innerHTML, /aria-label="Choose G: Glycine. Library quality: Review/);
  previewRow(ui);
  ui.run(`library.showPreview(${JSON.stringify({ ...preview, quality: monomer.quality })}, previewRow)`);
  assert.match(ui.element('lib-preview').innerHTML, /Library quality: Review &lt;atom&gt; mapping/);
  assert.doesNotMatch(ui.element('lib-preview').innerHTML, /<atom>/);
});

test('informational library notes remain escaped and separate from amber warnings', async () => {
  const ui = page('builder.js');
  const monomer = { abbr: 'G', name: 'Glycine', type: 'aa', subtype: 'natural', chem_types: '',
    quality: { issues: [{ severity: 'info', code: 'legacy_numbering', message: 'Legacy <atom> numbering & labels' }] } };
  const pending = ui.run('library.load()');
  latestRequest(ui, '/monomers').resolve([monomer]);
  await pending;
  assert.doesNotMatch(ui.element('lib-list').innerHTML, /lib-quality|Library quality:/);
  assert.match(ui.element('lib-list').innerHTML, /aria-label="Choose G: Glycine"/);
  previewRow(ui);
  ui.run(`library.showPreview(${JSON.stringify({ ...preview, quality: monomer.quality })}, previewRow)`);
  assert.match(ui.element('lib-preview').innerHTML,
    /class="prev-meta">Library notes: Legacy &lt;atom&gt; numbering &amp; labels/);
  assert.doesNotMatch(ui.element('lib-preview').innerHTML, /prev-warn|<atom>/);

  monomer.quality.issues.push({ severity: 'warning', code: 'uncertain_stereo', message: 'Review <stereo> assignment' });
  ui.run(`library.monomers = ${JSON.stringify([monomer])}; library.render('')`);
  assert.match(ui.element('lib-list').innerHTML, /lib-quality.*Review &lt;stereo&gt; assignment/);
  assert.match(ui.element('lib-list').innerHTML, /aria-label="Choose G: Glycine. Library quality: Review &lt;stereo&gt; assignment"/);
  assert.doesNotMatch(ui.element('lib-list').innerHTML, /Legacy/);
  ui.run(`library.showPreview(${JSON.stringify({ ...preview, quality: monomer.quality })}, previewRow)`);
  assert.match(ui.element('lib-preview').innerHTML, /class="prev-meta">Library notes: Legacy &lt;atom&gt; numbering &amp; labels/);
  assert.match(ui.element('lib-preview').innerHTML, /class="prev-meta prev-warn">Library quality: Review &lt;stereo&gt; assignment/);
  assert.doesNotMatch(ui.element('lib-preview').innerHTML, /<atom>|<stereo>/);
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

test('worker renewal retries after Retry-After and uses the same current document', async () => {
  const ui = page('builder.js');
  ui.element('cabiln-input').value = 'A-G';
  const pending = ui.run('doRenderCabiln("A-G")');
  for (let attempt = 0; attempt < 4; attempt++) {
    latestRequest(ui, '/render').resolve({ error: 'Chemistry workers are busy' }, false, 503, { 'Retry-After': '1' });
    await new Promise(setImmediate);
    assert.equal(ui.requests.length, attempt + 1);
    await ui.timers();
    assert.equal(ui.requests.length, attempt + 2);
    assert.equal(ui.requests.at(-1).options.body, ui.requests[0].options.body);
  }
  ui.requests.at(-1).resolve(rendered);
  await pending;
  assert.equal(ui.element('cabiln-input').className, 'ok');
});

test('drawing during Build startup survives library admission retries', async () => {
  for (const [mode, source, endpoint, data] of [
    ['cabiln', 'C-A-C', '/render', rendered],
    ['smiles', 'NCC(=O)O', '/render_reference', glycineDrawing],
  ]) {
    const ui = page('builder.js');
    await ui.element('btn-build').click();
    latestRequest(ui, '/monomers').resolve([]);
    await flush();
    ui.element('notation-select').value = mode;
    await ui.input('cabiln-input', source);
    await ui.timers();
    const rejected = new Set();
    for (let attempt = 0; attempt < 3; attempt++) {
      // Reproduce recovery followed by the reaction request occupying the
      // worker during the drawing's final retry, as observed in browser CI.
      for (const request of ui.requests.filter(request => request.url === endpoint && !rejected.has(request))) {
        request.resolve({ error: 'Chemistry capacity is busy; retry shortly' }, false, 503, { 'Retry-After': '1' });
        rejected.add(request);
      }
      latestRequest(ui, '/reactions').resolve(attempt < 2 ? { error: 'Busy' } : [],
        attempt === 2, attempt < 2 ? 503 : 200, { 'Retry-After': '1' });
      await flush();
      await ui.timers();
    }
    await flush();
    for (const request of ui.requests.filter(request => request.url === endpoint && !rejected.has(request))) {
      request.resolve(data);
    }
    await flush();
    assert.equal(ui.element('cabiln-input').className, 'ok', mode);
    assert.equal(ui.element('cabiln-input').value, source);
  }
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

test('clearing a drawing waiting for library metadata cancels its submission', async () => {
  for (const [mode, source] of [['cabiln', 'C-A-C'], ['smiles', 'NCC(=O)O']]) {
    const ui = page('builder.js');
    await ui.element('btn-build').click();
    ui.element('notation-select').value = mode;
    await ui.input('cabiln-input', source);
    await ui.timers();
    await ui.input('cabiln-input', '');
    latestRequest(ui, '/monomers').resolve([]);
    await flush();
    latestRequest(ui, '/reactions').resolve([]);
    await flush();
    assert.deepEqual(ui.requests.map(request => request.url), ['/monomers', '/reactions']);
    assert.equal(ui.element('cabiln-input').value, '');
  }
});

test('persistent overload is visible after four retries and other errors are never retried', async () => {
  for (const [status, header, attempts] of [[503, '1', 5], [503, '30', 1], [503, '', 1], [500, '1', 1]]) {
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
  let pending = ui.run('library.loadReactions()');
  latestRequest(ui, '/reactions').resolve({ error: 'Server remains busy' }, false, 503);
  await pending;
  assert.equal(ui.run('library.reactionPairs'), null);
  assert.match(ui.element('lib-filter-status').textContent, /Server remains busy/);
  pending = ui.run('library.loadReactions()');
  latestRequest(ui, '/reactions').resolve([['backbone_n', 'backbone_c']]);
  await pending;
  assert.equal(ui.run('library.reactionPairs.length'), 1);
  assert.doesNotMatch(ui.element('lib-filter-status').textContent, /unavailable|Server remains busy/);
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
