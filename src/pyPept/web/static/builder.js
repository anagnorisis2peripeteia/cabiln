// ─── state ────────────────────────────────────────────────────────────────────
let darkMode   = true;
let verifyMode = false;
let hlEnabled  = true;
let libLoaded  = false;
let allMonomers = [];
let cabilnTimer = null;
let smilesTimer = null;
let lastCabiln  = '';
let lastSmiles  = '';
let lastSvg     = '';
let lastMolBlock = '';
let residueMap  = {};
let atomToRes   = {};
let residueList = [];
let previewCache = {};
let previewTimer = null;
let buildMode    = false;
let buildLeft    = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
let buildRight   = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
let buildLeftRIdx = null;
let buildRightRIdx = null;
let insertBetweenActive = false;
let rerollSeed   = 0;
let reactionPairs = null;  // lazy-loaded list of [ct_a, ct_b] pairs
let rxnFilterActive = false;

// A response may already be queued when abort() runs. Only the current request
// in each group may change the UI, even if a cancelled fetch still resolves.
const activeRequests = new Map();
function cancelRequests(...keys) {
  for (const key of keys) {
    const pending = activeRequests.get(key);
    if (!pending) continue;
    activeRequests.delete(key);
    pending.controller.abort();
    pending.onEnd();
  }
}
function startRequest(key, onEnd = () => {}) {
  cancelRequests(key);
  const pending = { controller: new AbortController(), onEnd };
  activeRequests.set(key, pending);
  return {
    signal: pending.controller.signal,
    current: () => activeRequests.get(key) === pending,
    finish() {
      if (activeRequests.get(key) !== pending) return;
      activeRequests.delete(key);
      onEnd();
    },
  };
}

async function readResponse(response) {
  const data = await response.json();
  if (!response.ok && !data.error) {
    const detail = data.detail;
    data.error = Array.isArray(detail)
      ? detail.map(item => item.msg).join('; ')
      : detail || 'The request could not be completed.';
  }
  return data;
}

const RES_COLORS = [
  '#2a5080','#2a8050','#802a50','#806a2a','#502a80',
  '#2a6080','#80502a','#2a8070','#6a2a80','#80802a',
  '#3a6080','#3a8060','#603a50','#706a3a','#403a70',
];

// ─── elements ─────────────────────────────────────────────────────────────────
const cabilnInput   = document.getElementById('cabiln-input');
const cabilnStatus  = document.getElementById('cabiln-status');
const conversionStatus = document.getElementById('conversion-status');
const renderInner   = document.getElementById('render-inner');
const renderCanvas  = document.getElementById('render-canvas');
const molUpload     = document.getElementById('mol-upload');
const smilesInput   = document.getElementById('smiles-input');
const smilesStatus  = document.getElementById('smiles-status');
const smilesInner   = document.getElementById('smiles-inner');
const compareBar    = document.getElementById('compare-bar');
const btnDark       = document.getElementById('btn-dark');
const btnVerify     = document.getElementById('btn-verify');
const btnLib        = document.getElementById('btn-lib');
const btnHl         = document.getElementById('btn-hl');
const btnPng        = document.getElementById('btn-png');
const btnMol        = document.getElementById('btn-mol');
const btnToBracket  = document.getElementById('btn-to-bracket');
const btnToBranch   = document.getElementById('btn-to-branch');
const libPanel      = document.getElementById('lib-panel');
const libSearch     = document.getElementById('lib-search');
const libClose      = document.getElementById('lib-close');
const libList       = document.getElementById('lib-list');
const libCount      = document.getElementById('lib-count');
const libPreview    = document.getElementById('lib-preview');
const btnExamples    = document.getElementById('btn-examples');
const examplesPanel  = document.getElementById('examples-panel');
const examplesClose  = document.getElementById('examples-close');
const examplesList   = document.getElementById('examples-list');
const resChips      = document.getElementById('residue-chips');
const buildPanel    = document.getElementById('build-panel');
const btnBuild      = document.getElementById('btn-build');
const buildClose    = document.getElementById('build-close');
const buildConnect  = document.getElementById('build-connect');
const buildStatus   = document.getElementById('build-status');
const buildHint     = document.getElementById('build-hint');
const buildLeftAbbr = document.getElementById('build-left-abbr');
const buildLeftSvg  = document.getElementById('build-left-svg');
const buildLeftRg   = document.getElementById('build-left-rgroups');
const buildRightAbbr= document.getElementById('build-right-abbr');
const buildRightSvg = document.getElementById('build-right-svg');
const buildRightRg  = document.getElementById('build-right-rgroups');
const buildInsertRow = document.getElementById('build-insert-row');
const buildInsertInfo = document.getElementById('build-insert-info');
const buildInsertBtn = document.getElementById('build-insert-btn');
const btnReroll     = document.getElementById('btn-reroll');
const btnS2c        = document.getElementById('btn-s2c');
const btnS2cBracket = document.getElementById('btn-s2c-bracket');
const btnRxnFilter  = document.getElementById('btn-rxn-filter');
const notationSelect = document.getElementById('notation-select');
const btnToCabilnPct     = document.getElementById('btn-to-cabiln-pct');
const btnToCabilnBracket = document.getElementById('btn-to-cabiln-bracket');

async function loadCapabilities() {
  const link = document.getElementById('register-link');
  link.hidden = true;
  try {
    const response = await fetch('/capabilities');
    const capabilities = await response.json();
    link.hidden = !response.ok || capabilities.registration !== true;
  } catch (error) {
    link.hidden = true;
  }
}
loadCapabilities();

// ─── dark mode ────────────────────────────────────────────────────────────────
btnDark.addEventListener('click', () => {
  darkMode = !darkMode;
  btnDark.classList.toggle('active', darkMode);
  for (const el of document.querySelectorAll('.canvas-wrap')) {
    el.classList.toggle('dark', darkMode);
  }
  libPreview.classList.toggle('dark', darkMode);
  document.querySelectorAll('.build-box').forEach(el => el.classList.toggle('dark', darkMode));
});
// apply dark mode on load
btnDark.classList.add('active');
document.querySelectorAll('.canvas-wrap').forEach(el => el.classList.add('dark'));
libPreview.classList.add('dark');
document.querySelectorAll('.build-box').forEach(el => el.classList.add('dark'));

// ─── highlight toggle ─────────────────────────────────────────────────────────
btnHl.addEventListener('click', () => {
  hlEnabled = !hlEnabled;
  btnHl.classList.toggle('active', hlEnabled);
  if (!hlEnabled) clearHighlight();
});

// ─── verify mode ──────────────────────────────────────────────────────────────
btnVerify.addEventListener('click', () => {
  verifyMode = !verifyMode;
  btnVerify.classList.toggle('active', verifyMode);
  document.getElementById('verify-pane').style.display = verifyMode ? '' : 'none';
  compareBar.style.display = verifyMode ? '' : 'none';
  clearComparison();
  if (verifyMode) triggerVerify();
});

// ─── library sidebar (persistent) ────────────────────────────────────────────
function openLib() {
  libPanel.classList.add('open');
  btnLib.classList.add('active');
  loadMonomers();
  libSearch.focus();
  loadReactions();
}
function closeLib() {
  libPanel.classList.remove('open');
  btnLib.classList.remove('active');
  hidePreview();
}

btnLib.addEventListener('click', () =>
  libPanel.classList.contains('open') ? closeLib() : openLib());
libClose.addEventListener('click', closeLib);

libSearch.addEventListener('input', () => renderLibList(libSearch.value.trim().toLowerCase()));

btnRxnFilter.addEventListener('click', () => {
  rxnFilterActive = !rxnFilterActive;
  btnRxnFilter.classList.toggle('active', rxnFilterActive);
  renderLibList(libSearch.value.trim().toLowerCase());
});

// ─── example peptide sidebar ──────────────────────────────────────────────────
let examplesLoaded = false;

function openExamples() {
  examplesPanel.classList.add('open');
  btnExamples.classList.add('active');
  if (!examplesLoaded) loadExamples();
}
function closeExamples() {
  examplesPanel.classList.remove('open');
  btnExamples.classList.remove('active');
}

btnExamples.addEventListener('click', () =>
  examplesPanel.classList.contains('open') ? closeExamples() : openExamples());
examplesClose.addEventListener('click', closeExamples);

async function loadExamples() {
  try {
    const res = await fetch('/examples');
    const data = await readResponse(res);
    renderExamples(data);
    examplesLoaded = true;
  } catch (e) {
    examplesList.innerHTML = '<div class="placeholder err">Failed to load examples</div>';
  }
}

function renderExamples(categories) {
  const rows = [];
  for (const cat of categories) {
    rows.push(`<div class="example-cat">${escHtml(cat.category)}</div>`);
    for (const item of cat.items) {
      const preview = item.cabiln.length > 55
        ? item.cabiln.slice(0, 52) + '…'
        : item.cabiln;
      rows.push(`<div class="example-row" data-cabiln="${escAttr(item.cabiln)}">
        <div class="example-name">${escHtml(item.name)}</div>
        <div class="example-desc">${escHtml(item.description)}</div>
        <div class="example-seq" title="${escAttr(item.cabiln)}">${escHtml(preview)}</div>
      </div>`);
    }
  }
  examplesList.innerHTML = rows.join('');
  examplesList.querySelectorAll('.example-row').forEach(row => {
    row.addEventListener('click', () => {
      cabilnInput.value = row.dataset.cabiln;
      cabilnInput.dispatchEvent(new Event('input'));
      cabilnInput.focus();
    });
  });
}

async function loadMonomers() {
  const request = startRequest('library');
  if (!libLoaded) libList.innerHTML = '<div class="placeholder">Loading…</div>';
  try {
    const res = await fetch('/monomers', { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current()) return;
    if (!Array.isArray(data)) throw new Error(data.error || 'Invalid monomer library');
    allMonomers = data;
    libLoaded = true;
    renderLibList(libSearch.value.trim().toLowerCase());
  } catch (e) {
    if (!request.current()) return;
    libList.innerHTML = '<div class="placeholder err">Failed to load monomers</div>';
  } finally {
    request.finish();
  }
}

window.addEventListener('focus', () => {
  if (libPanel.classList.contains('open')) loadMonomers();
});

async function loadReactions() {
  if (reactionPairs !== null) return;
  try {
    const res = await fetch('/reactions');
    reactionPairs = await res.json();
  } catch (e) { reactionPairs = []; }
}

function parseCts(cts) {
  if (!cts) return [];
  return cts.split(',').map(p => p.includes(':') ? p.split(':')[1].trim() : p.trim()).filter(Boolean);
}

function renderLibList(q) {
  let filtered = q
    ? allMonomers.filter(m =>
        m.abbr.toLowerCase().includes(q) ||
        m.name.toLowerCase().includes(q) ||
        m.type.toLowerCase().includes(q) ||
        (m.chem_types || '').toLowerCase().includes(q)
      )
    : allMonomers;

  if (rxnFilterActive && buildLeft && Array.isArray(reactionPairs)) {
    const pairSet = new Set(reactionPairs.map(([a, b]) => a + '|' + b));
    let lcts;
    if (buildLeft.selectedSlot !== null) {
      const selRg = buildLeft.rgroups.find(r => r.slot === buildLeft.selectedSlot);
      lcts = selRg && !selRg.used ? [selRg.chem_type] : [];
    } else {
      lcts = buildLeft.rgroups.filter(r => !r.used).map(r => r.chem_type);
    }
    filtered = filtered.filter(m => {
      const mcts = parseCts(m.chem_types);
      return mcts.some(mct => lcts.some(lct =>
        pairSet.has(lct + '|' + mct) || pairSet.has(mct + '|' + lct)
      ));
    });
  }

  if (insertBetweenActive) {
    filtered = filtered.filter(m => m.backbone_insertable);
  }

  libCount.textContent = `${filtered.length} / ${allMonomers.length} monomers`;

  if (!filtered.length) {
    libList.innerHTML = '<div class="placeholder">No matches</div>';
    return;
  }

  const rows = filtered.map(m => {
    const badge = m.degenerate
      ? `<span class="lib-badge cap">N/C cap</span>`
      : m.subtype === 'modified' || m.subtype === 'natural'
      ? `<span class="lib-badge aa">${m.type}</span>`
      : m.type === 'cap' && m.subtype === 'protecting'
      ? `<span class="lib-badge protect">cap</span>`
      : `<span class="lib-badge cap">${m.type}</span>`;

    let lg = m.leaving ? `  LG: ${escHtml(m.leaving)}` : '';
    if (m.degenerate) {
      const parts = [];
      if (m.nterm_abbr) parts.push(`N: ${escHtml(m.nterm_abbr)} (${escHtml(m.nterm_leaving)})`);
      if (m.cterm_abbr) parts.push(`C: ${escHtml(m.cterm_abbr)} (${escHtml(m.cterm_leaving)})`);
      lg = '  ' + parts.join(' | ');
    }
    return `<div class="lib-row" data-abbr="${escAttr(m.abbr)}">
      <div class="lib-abbr">${escHtml(m.abbr)}</div>
      <div class="lib-info">
        <div class="lib-name" title="${escAttr(m.name)}">${escHtml(m.name)}</div>
        <div class="lib-meta">${escHtml(m.chem_types || '')}${lg}</div>
      </div>
      ${badge}
    </div>`;
  });
  libList.innerHTML = rows.join('');

  libList.querySelectorAll('.lib-row').forEach(row => {
    row.addEventListener('click', () => {
      if (buildMode && insertBetweenActive) {
        doInsertBetween(row.dataset.abbr);
      } else if (buildMode) {
        if (!cabilnInput.value.trim()) {
          // No sequence yet — insert as first monomer
          cabilnInput.value = row.dataset.abbr;
          cabilnInput.dispatchEvent(new Event('input'));
          buildHint.textContent = 'First monomer added — click its chip, then right-click another monomer';
        } else if (!buildLeft) {
          // Sequence exists but no chip selected — prompt
          buildHint.textContent = 'Click a chip on the sequence first, then click a library monomer';
        } else {
          loadBuildRight(row.dataset.abbr);
        }
      } else {
        insertAbbr(row.dataset.abbr);
      }
    });
    row.addEventListener('contextmenu', e => {
      if (buildMode) {
        e.preventDefault();
        loadBuildRight(row.dataset.abbr);
      }
    });
    row.addEventListener('mouseenter', e => startPreview(row.dataset.abbr, row));
    row.addEventListener('mouseleave', () => hidePreview());
  });
}

function insertAbbr(abbr) {
  const ta = cabilnInput;
  const start = ta.selectionStart;
  const end   = ta.selectionEnd;
  const val   = ta.value;
  const before = val.slice(0, start);
  const after  = val.slice(end);
  const needDash = before.length > 0 && !before.endsWith('-') && !before.endsWith('\n');
  const insert = (needDash ? '-' : '') + abbr;
  ta.value = before + insert + after;
  const newPos = start + insert.length;
  ta.setSelectionRange(newPos, newPos);
  ta.focus();
  ta.dispatchEvent(new Event('input'));
}

// ─── monomer preview tooltip ──────────────────────────────────────────────────
function startPreview(abbr, row) {
  clearTimeout(previewTimer);
  previewTimer = setTimeout(async () => {
    if (previewCache[abbr]) {
      showPreview(previewCache[abbr], row);
      return;
    }
    try {
      const res = await fetch(`/monomer_svg?abbr=${encodeURIComponent(abbr)}`);
      const data = await readResponse(res);
      if (data.svg) {
        previewCache[abbr] = data;
        showPreview(data, row);
      }
    } catch (e) { /* silent */ }
  }, 200);
}

function showPreview(data, row) {
  let html;
  if (data.degenerate && data.variants) {
    // Degenerate: variant panels with optional reagent info
    let panels = data.variants.map(v => {
      let pane = `<div class="prev-pane"><span class="prev-label">${v.label}</span>${v.svg}`;
      if (v.reagent) {
        pane += `<span class="prev-rxn">${v.reagent.reaction} (LG: ${v.reagent.reagent_lg})</span>`;
      }
      pane += `</div>`;
      return pane;
    }).join('');
    // Show reagent form from first variant that has one
    const withReagent = data.variants.find(v => v.svg_reagent);
    if (withReagent) {
      panels += `<div class="prev-pane"><span class="prev-label">Reagent</span>${withReagent.svg_reagent}</div>`;
    }
    panels += `<div class="prev-pane"><span class="prev-label">R-groups</span>${data.svg}</div>`;
    html = `<div class="prev-row">${panels}</div>`;
  } else if (data.degenerate) {
    // Legacy format fallback
    html = `<div class="prev-row">
      <div class="prev-pane"><span class="prev-label">N-term</span>${data.svg_nterm || data.svg}</div>
      <div class="prev-pane"><span class="prev-label">C-term</span>${data.svg_cterm || data.svg}</div>
    </div>`;
  } else {
    const restored = data.svg_restored || data.svg;
    let reagentPane = '';
    let metaLine = '';
    if (data.svg_reagent) {
      reagentPane = `<div class="prev-pane"><span class="prev-label">Reagent</span>${data.svg_reagent}</div>`;
      const r = data.reagent;
      metaLine = `<div class="prev-meta">${r.reaction}` +
        (r.reagent_note ? ` — ${r.reagent_note}` : '') +
        (r.issue ? ` <span class="prev-warn">⚠ ${r.issue}</span>` : '') +
        `</div>`;
    }
    html = `<div class="prev-row">
      <div class="prev-pane"><span class="prev-label">Monomer</span>${restored}</div>
      ${reagentPane}
      <div class="prev-pane"><span class="prev-label">R-groups</span>${data.svg}</div>
    </div>${metaLine}`;
  }
  libPreview.innerHTML = html;
  const hasReagent = !!(data.svg_reagent || (data.variants && data.variants.some(v => v.svg_reagent)));
  libPreview.classList.toggle('has-reagent', hasReagent);
  const rect = row.getBoundingClientRect();
  const previewH = hasReagent ? 240 : 202;
  let top = Math.max(8, rect.top - 40);
  if (top + previewH > window.innerHeight - 8) {
    top = window.innerHeight - 8 - previewH;
  }
  top = Math.max(8, top);
  libPreview.style.left = (rect.right + 8) + 'px';
  libPreview.style.top  = top + 'px';
  libPreview.style.display = 'block';
}

function hidePreview() {
  clearTimeout(previewTimer);
  libPreview.style.display = 'none';
}

// ─── residue chips + bidirectional highlighting ───────────────────────────────
let chainData = [];
let currentBranchSet = new Set();
let xlinkByRes = {};

function highlightGroup(idxList) {
  if (!hlEnabled) return;
  clearHighlight();
  activeRIdx = -999;
  const svg = document.querySelector('#render-inner svg');
  if (!svg) return;
  svg.classList.add('has-highlight');
  for (const rIdx of idxList) {
    const atoms = residueMap[rIdx] || [];
    for (const aidx of atoms) {
      svg.querySelectorAll(`.atom-${aidx}`).forEach(el =>
        el.classList.add('res-hl'));
    }
  }
  resChips.classList.add('dimmed');
  resChips.querySelectorAll('.res-chip').forEach(c => {
    const ri = parseInt(c.dataset.residue);
    c.classList.toggle('hover', idxList.includes(ri));
  });
  resChips.querySelectorAll('.branch-chip').forEach(c => {
    const cm = JSON.parse(c.dataset.members || '[]');
    // Highlight a branch chip only when ALL its members are in the hover set,
    // not merely "any overlap". Otherwise hovering on one !2 lights up the
    // sibling [!1] and [!3] brackets too because they share the scaffold
    // monomer (e.g. TBMB) — overzealous and confusing for tri-arm scaffolds.
    c.classList.toggle('hover', cm.length > 0 && cm.every(m => idxList.includes(m)));
  });
}

function buildResidueUI(resMap, residues, layout, crosslinkGroups) {
  residueMap = resMap || {};
  residueList = residues || [];
  const segments = layout?.segments || [];
  chainData = segments.map(segment => ({ residues: segment.roots }));
  atomToRes = {};
  for (const [rIdx, atoms] of Object.entries(residueMap)) {
    for (const aidx of atoms) atomToRes[aidx] = parseInt(rIdx);
  }

  resChips.innerHTML = '';
  resChips.style.position = '';
  resChips.style.paddingLeft = '';
  if (!residueList.length) return;

  const resById = {};
  residueList.forEach((r, i) => { resById[r.idx] = { ...r, colorIdx: i }; });

  const xlinkByMember = {};
  (crosslinkGroups || []).forEach(g => {
    g.members.forEach(mIdx => {
      if (!xlinkByMember[mIdx]) xlinkByMember[mIdx] = [];
      xlinkByMember[mIdx].push(g);
    });
  });
  xlinkByRes = xlinkByMember;

  function makeChip(rIdx, simpleHover) {
    const r = resById[rIdx];
    if (!r) return null;
    const chip = document.createElement('span');
    chip.className = 'res-chip';
    chip.textContent = r.abbr;
    chip.dataset.residue = r.idx;
    chip.style.background = RES_COLORS[r.colorIdx % RES_COLORS.length];
    const xlinks = xlinkByMember[r.idx];
    if (!simpleHover && xlinks && xlinks.length) {
      const allMembers = [...new Set(xlinks.flatMap(g => g.members))];
      chip.addEventListener('mouseenter', () => highlightGroup(allMembers));
    } else {
      chip.addEventListener('mouseenter', () => highlightResidue(r.idx));
    }
    chip.addEventListener('mouseleave', clearHighlight);
    chip.addEventListener('click', () => {
      if (buildMode) {
        if (buildLeft && buildLeftRIdx !== r.idx) {
          loadBuildRight(r.abbr, r.idx);
          chip.style.outline = '2px solid #e0a05a';
          resChips.querySelectorAll('.res-chip').forEach(c => {
            if (c !== chip && parseInt(c.dataset.residue) !== buildLeftRIdx)
              c.style.outline = '';
          });
        } else {
          loadBuildLeft(r.abbr, r.idx);
          chip.style.outline = '2px solid #5a9ae0';
          resChips.querySelectorAll('.res-chip').forEach(c => {
            if (c !== chip) c.style.outline = '';
          });
        }
      }
    });
    return chip;
  }

  function makeSeparator(text, memberIdxs) {
    const el = document.createElement('span');
    el.className = 'res-chip branch-chip';
    el.style.background = '#3a3a50';
    el.style.fontWeight = '700';
    el.textContent = text;
    el.dataset.members = JSON.stringify(memberIdxs || []);
    if (memberIdxs && memberIdxs.length) {
      el.addEventListener('mouseenter', () => highlightGroup(memberIdxs));
      el.addEventListener('mouseleave', clearHighlight);
    }
    return el;
  }

  function makeXlinkChip(tag, members) {
    const el = document.createElement('span');
    el.className = 'res-chip branch-chip xlink-chip';
    el.style.background = '#503a4a';
    el.style.fontWeight = '700';
    el.style.fontSize = '0.8em';
    el.textContent = tag;
    el.dataset.members = JSON.stringify(members);
    el.addEventListener('mouseenter', () => highlightGroup(members));
    el.addEventListener('mouseleave', clearHighlight);
    return el;
  }

  const groups = layout?.groups || [];
  const markers = layout?.markers || [];
  const groupsByHost = new Map();
  for (const group of groups) {
    if (!groupsByHost.has(group.host)) groupsByHost.set(group.host, []);
    groupsByHost.get(group.host).push(group);
  }
  currentBranchSet = new Set(groups.flatMap(group => group.members));

  function expandWithXlinks(ids) {
    const members = new Set(ids);
    for (const id of ids) {
      for (const link of xlinkByMember[id] || []) {
        link.members.forEach(member => members.add(member));
      }
    }
    return [...members];
  }

  function groupMembers(group) {
    return expandWithXlinks([
      ...group.members,
      ...markers.filter(marker => marker.group === group.id).flatMap(marker => marker.members),
    ]);
  }

  function appendGroup(group) {
    const members = groupMembers(group);
    if (group.opening) resChips.appendChild(makeSeparator(group.opening, members));
    // A marker-only bracket belongs to its host but contains no monomer.
    if (!group.roots.length) {
      markers.filter(marker => marker.group === group.id).forEach(marker => {
        resChips.appendChild(makeXlinkChip(marker.tag, marker.members));
      });
    }
    group.roots.forEach(id => appendOccurrence(id, group.id));
    if (group.closing) resChips.appendChild(makeSeparator(group.closing, members));
  }

  function appendOccurrence(id, context = null) {
    const ownMarkers = markers.filter(marker => marker.residue === id && marker.group === context);
    ownMarkers.filter(marker => marker.before).forEach(marker => {
      resChips.appendChild(makeXlinkChip(marker.tag, marker.members));
    });
    const links = xlinkByMember[id] || [];
    const chip = makeChip(id, context !== null || links.length > 1);
    if (chip) resChips.appendChild(chip);
    const children = (groupsByHost.get(id) || []).filter(group => group.parent === context);
    if (children.length || links.length > 1) {
      const members = expandWithXlinks([id, ...children.flatMap(groupMembers)]);
      const whole = makeSeparator('$', members);
      whole.style.background = '#3a5050';
      resChips.appendChild(whole);
    }
    children.forEach(appendGroup);
    ownMarkers.filter(marker => !marker.before).forEach(marker => {
      resChips.appendChild(makeXlinkChip(marker.tag, marker.members));
    });
  }

  for (const segment of segments) {
    if (segments.length > 1) {
      resChips.appendChild(makeSeparator('%', expandWithXlinks(segment.members)));
    }
    segment.roots.forEach(id => appendOccurrence(id));
  }

  wireUpSvgHover();
}

let activeRIdx = null;

function highlightResidue(rIdx) {
  if (!hlEnabled) return;
  if (rIdx === activeRIdx) return;
  clearHighlight();
  activeRIdx = rIdx;
  const atoms = residueMap[rIdx] || [];
  const svg = document.querySelector('#render-inner svg');
  if (!svg) return;

  svg.classList.add('has-highlight');
  for (const aidx of atoms) {
    svg.querySelectorAll(`.atom-${aidx}`).forEach(el =>
      el.classList.add('res-hl'));
  }

  resChips.classList.add('dimmed');
  resChips.querySelectorAll('.res-chip').forEach(c =>
    c.classList.toggle('hover', parseInt(c.dataset.residue) === rIdx));
}

function clearHighlight() {
  activeRIdx = null;
  const svg = document.querySelector('#render-inner svg');
  if (svg) {
    svg.classList.remove('has-highlight');
    svg.querySelectorAll('.res-hl').forEach(el => el.classList.remove('res-hl'));
  }
  resChips.classList.remove('dimmed');
  resChips.querySelectorAll('.res-chip.hover').forEach(c => c.classList.remove('hover'));
}

function wireUpSvgHover() {
  const svg = document.querySelector('#render-inner svg');
  if (!svg) return;
  svg.addEventListener('mousemove', e => {
    if (!hlEnabled) return;
    let el = e.target;
    while (el && el !== svg) {
      const cls = (el.getAttribute('class') || '');
      const m = cls.match(/atom-(\d+)/);
      if (m) {
        const rIdx = atomToRes[parseInt(m[1])];
        if (rIdx !== undefined) {
          const xlinks = xlinkByRes[rIdx];
          if (xlinks && xlinks.length) {
            const allMembers = [...new Set(xlinks.flatMap(g => g.members))];
            highlightGroup(allMembers);
          } else {
            highlightResidue(rIdx);
          }
          return;
        }
      }
      el = el.parentElement;
    }
    clearHighlight();
  });
  svg.addEventListener('mouseleave', clearHighlight);
  svg.addEventListener('click', e => {
    if (!buildMode) return;
    let el = e.target;
    while (el && el !== svg) {
      const cls = (el.getAttribute('class') || '');
      const m = cls.match(/atom-(\d+)/);
      if (m) {
        const rIdx = atomToRes[parseInt(m[1])];
        if (rIdx !== undefined) {
          const r = residueList.find(r => r.idx === rIdx);
          if (r) {
            if (buildLeft && buildLeftRIdx !== rIdx) {
              loadBuildRight(r.abbr, r.idx);
              resChips.querySelectorAll('.res-chip').forEach(c => {
                const cIdx = parseInt(c.dataset.residue);
                if (cIdx === rIdx) c.style.outline = '2px solid #e0a05a';
                else if (cIdx !== buildLeftRIdx) c.style.outline = '';
              });
            } else {
              loadBuildLeft(r.abbr, r.idx);
              resChips.querySelectorAll('.res-chip').forEach(c => {
                c.style.outline = parseInt(c.dataset.residue) === rIdx
                  ? '2px solid #5a9ae0' : '';
              });
            }
          }
          return;
        }
      }
      el = el.parentElement;
    }
  });
}

// ─── zoom / pan (shared, wired per canvas) ────────────────────────────────────
function makeZoomable(canvas, inner) {
  let scale = 1, tx = 0, ty = 0;
  let dragging = false, startX, startY, startTx, startTy;

  function applyTransform() {
    inner.style.transform = `translate(${tx}px,${ty}px) scale(${scale})`;
  }

  canvas.addEventListener('wheel', e => {
    e.preventDefault();
    const factor = e.deltaY < 0 ? 1.12 : 1 / 1.12;
    scale = Math.min(Math.max(scale * factor, 0.2), 20);
    applyTransform();
  }, { passive: false });

  canvas.addEventListener('mousedown', e => {
    if (e.button !== 0) return;
    dragging = true;
    startX = e.clientX; startY = e.clientY;
    startTx = tx; startTy = ty;
    canvas.style.cursor = 'grabbing';
  });
  window.addEventListener('mousemove', e => {
    if (!dragging) return;
    tx = startTx + e.clientX - startX;
    ty = startTy + e.clientY - startY;
    applyTransform();
  });
  window.addEventListener('mouseup', () => {
    dragging = false;
    canvas.style.cursor = 'grab';
  });
  canvas.addEventListener('dblclick', () => {
    scale = 1; tx = 0; ty = 0;
    applyTransform();
  });
}

makeZoomable(renderCanvas, renderInner);
makeZoomable(document.getElementById('smiles-canvas'), smilesInner);

// ─── PNG download (client-side SVG → canvas → PNG) ────────────────────────────
btnPng.addEventListener('click', () => {
  if (!lastSvg) return;
  const blob = new Blob([lastSvg], { type: 'image/svg+xml;charset=utf-8' });
  const url  = URL.createObjectURL(blob);
  const img  = new Image();
  img.onload = () => {
    const w = img.naturalWidth  || 1200;
    const h = img.naturalHeight || 900;
    const cv = document.createElement('canvas');
    cv.width = w; cv.height = h;
    const ctx = cv.getContext('2d');
    ctx.fillStyle = '#fff';
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0);
    const a = document.createElement('a');
    a.download = 'structure.png';
    a.href = cv.toDataURL('image/png');
    a.click();
    URL.revokeObjectURL(url);
  };
  img.onerror = () => URL.revokeObjectURL(url);
  img.src = url;
});

// ─── MOL download (server-side mol block) ────────────────────────────────────
btnMol.addEventListener('click', async () => {
  if (!lastMolBlock) return;
  const blob = new Blob([lastMolBlock], { type: 'chemical/x-mdl-molfile' });
  const a = document.createElement('a');
  a.download = 'structure.mol';
  a.href = URL.createObjectURL(blob);
  a.click();
});

// ─── build mode ──────────────────────────────────────────────────────────────
function openBuild() {
  buildMode = true;
  buildPanel.classList.add('open');
  btnBuild.classList.add('active');
  if (!libPanel.classList.contains('open')) openLib();
  clearBuild();
}
function closeBuild() {
  buildMode = false;
  buildPanel.classList.remove('open');
  btnBuild.classList.remove('active');
  clearBuild();
}
function clearBuild() {
  cancelRequests('build-left', 'build-right', 'bond-check', 'sequence-edit');
  buildLeft = null; buildRight = null; buildLeftRIdx = null; buildRightRIdx = null;
  insertBetweenActive = false;
  buildInsertRow.style.display = 'none';
  buildInsertBtn.textContent = '⊕ Insert Between';
  buildLeftAbbr.textContent = '—';
  buildLeftSvg.innerHTML = '<div class="box-placeholder">Click a chip above</div>';
  buildLeftRg.innerHTML = '';
  buildRightAbbr.textContent = '—';
  buildRightSvg.innerHTML = '<div class="box-placeholder">Right-click from library</div>';
  buildRightRg.innerHTML = '';
  buildConnect.disabled = true;
  buildStatus.textContent = '';
  buildStatus.className = 'build-status';
  buildHint.textContent = 'Select a chip on the sequence, then right-click a monomer in the library';
  resChips.querySelectorAll('.res-chip').forEach(c => c.style.outline = '');
  btnRxnFilter.disabled = true;
  if (rxnFilterActive) {
    rxnFilterActive = false;
    btnRxnFilter.classList.remove('active');
    if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  }
}

btnBuild.addEventListener('click', () => buildMode ? closeBuild() : openBuild());
buildClose.addEventListener('click', closeBuild);

function selectedBackbone() {
  return chainData.find(chain =>
    chain.residues.includes(buildLeftRIdx) && chain.residues.includes(buildRightRIdx)
  )?.residues;
}

function checkAdjacentBackbone() {
  if (buildLeftRIdx == null || buildRightRIdx == null) return false;
  if (buildLeftRIdx === buildRightRIdx) return false;
  const main = selectedBackbone();
  if (!main) return false;
  if (currentBranchSet.has(buildLeftRIdx) || currentBranchSet.has(buildRightRIdx)) return false;
  const posL = main.indexOf(buildLeftRIdx);
  const posR = main.indexOf(buildRightRIdx);
  return Math.abs(posL - posR) === 1;
}

function updateInsertBetweenUI() {
  if (buildMode && checkAdjacentBackbone()) {
    const la = (residueList.find(r => r.idx === buildLeftRIdx) || {}).abbr || '?';
    const ra = (residueList.find(r => r.idx === buildRightRIdx) || {}).abbr || '?';
    buildInsertInfo.textContent = `${la} and ${ra} are adjacent on the backbone`;
    buildInsertRow.style.display = 'flex';
  } else {
    buildInsertRow.style.display = 'none';
    if (insertBetweenActive) {
      insertBetweenActive = false;
      buildInsertBtn.textContent = '⊕ Insert Between';
      if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
    }
  }
}

buildInsertBtn.addEventListener('click', () => {
  cancelRequests('sequence-edit');
  if (insertBetweenActive) {
    insertBetweenActive = false;
    buildInsertBtn.textContent = '⊕ Insert Between';
    buildHint.textContent = 'Select a chip on the sequence, then right-click a monomer in the library';
    if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
    return;
  }
  if (!checkAdjacentBackbone()) return;
  insertBetweenActive = true;
  buildInsertBtn.textContent = '✕ Cancel';
  buildHint.textContent = 'Click a backbone monomer in the library to insert between the selected residues';
  if (!libPanel.classList.contains('open')) openLib();
  if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
});

async function doInsertBetween(abbr) {
  const main = selectedBackbone();
  if (!main) return;
  const posL = main.indexOf(buildLeftRIdx);
  const posR = main.indexOf(buildRightRIdx);
  const after_idx = posL < posR ? buildLeftRIdx : buildRightRIdx;
  const val = cabilnInput.value.trim();
  const request = startRequest('sequence-edit');
  buildHint.textContent = 'Inserting…';
  try {
    const res = await fetch('/insert_backbone', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ cabiln: val, after_idx, new_abbr: abbr }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== val) return;
    if (data.error) {
      buildHint.textContent = 'Insert failed: ' + data.error;
      return;
    }
    cabilnInput.value = data.result;
    cabilnInput.dispatchEvent(new Event('input'));
    // Reset build state for next operation
    buildLeft = null; buildRight = null; buildLeftRIdx = null; buildRightRIdx = null;
    insertBetweenActive = false;
    buildInsertRow.style.display = 'none';
    buildInsertBtn.textContent = '⊕ Insert Between';
    buildLeftAbbr.textContent = '—';
    buildLeftSvg.innerHTML = '<div class="box-placeholder">Click a chip above</div>';
    buildLeftRg.innerHTML = '';
    buildRightAbbr.textContent = '—';
    buildRightSvg.innerHTML = '<div class="box-placeholder">Right-click from library</div>';
    buildRightRg.innerHTML = '';
    buildConnect.disabled = true;
    buildHint.textContent = `${abbr} inserted — select chips to continue building`;
    resChips.querySelectorAll('.res-chip').forEach(c => c.style.outline = '');
    if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  } catch (e) {
    if (!request.current()) return;
    buildHint.textContent = 'Insert failed';
  } finally {
    request.finish();
  }
}
document.getElementById('build-left-change').addEventListener('click', clearBuild);
document.getElementById('build-right-change').addEventListener('click', clearBuild);

async function loadBuildLeft(abbr, rIdx) {
  cancelRequests('bond-check', 'sequence-edit');
  const request = startRequest('build-left');
  const sequence = cabilnInput.value.trim();
  buildLeftRIdx = rIdx;
  buildLeftAbbr.textContent = abbr;
  buildLeftSvg.innerHTML = '<div class="spinner"></div>';
  buildLeftRg.innerHTML = '';
  buildLeft = null;
  buildConnect.disabled = true;
  buildStatus.textContent = '';

  try {
    const res = await fetch(`/monomer_rgroups?abbr=${encodeURIComponent(abbr)}&residue_idx=${rIdx}&cabiln=${encodeURIComponent(sequence)}`, { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== sequence) return;
    if (data.error) {
      buildLeftSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
      return;
    }
    buildLeftSvg.innerHTML = data.svg || '';
    buildLeft = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
    renderRgroupButtons(buildLeftRg, buildLeft, 'left');
    buildHint.textContent = buildRight ? 'Select R-groups to connect' : 'Now right-click a monomer in the library';
    btnRxnFilter.disabled = false;
    if (rxnFilterActive && libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
    updateInsertBetweenUI();
  } catch (e) {
    if (!request.current()) return;
    buildLeftSvg.innerHTML = '<div class="box-placeholder">Error loading monomer</div>';
  } finally {
    request.finish();
  }
}

async function loadBuildRight(abbr, rIdx) {
  cancelRequests('bond-check', 'sequence-edit');
  const request = startRequest('build-right');
  const sequence = cabilnInput.value.trim();
  buildRightRIdx = rIdx !== undefined ? rIdx : null;
  buildRightAbbr.textContent = abbr;
  buildRightSvg.innerHTML = '<div class="spinner"></div>';
  buildRightRg.innerHTML = '';
  buildRight = null;
  buildConnect.disabled = true;
  buildStatus.textContent = '';

  const family = rIdx === undefined && allMonomers.find(m => m.abbr === abbr && m.degenerate);
  if (family) {
    buildRightSvg.innerHTML = '<div class="box-placeholder">Choose the attachment form</div>';
    for (const [label, symbol] of [['N-terminal', family.nterm_abbr], ['C-terminal', family.cterm_abbr]]) {
      if (!symbol) continue;
      const button = document.createElement('button');
      button.className = 'rgroup-btn';
      button.textContent = `${label}: ${symbol}`;
      button.addEventListener('click', () => loadBuildRight(symbol));
      buildRightRg.appendChild(button);
    }
    buildHint.textContent = 'Choose a monomer form, then its attachment site';
    request.finish();
    return;
  }

  try {
    let url = `/monomer_rgroups?abbr=${encodeURIComponent(abbr)}`;
    if (rIdx !== undefined) {
      url += `&residue_idx=${rIdx}&cabiln=${encodeURIComponent(sequence)}`;
    }
    const res = await fetch(url, { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== sequence) return;
    if (data.error) {
      buildRightSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
      return;
    }
    buildRightSvg.innerHTML = data.svg || '';
    buildRight = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
    renderRgroupButtons(buildRightRg, buildRight, 'right');
    buildHint.textContent = 'Select R-groups to connect';
    updateInsertBetweenUI();
  } catch (e) {
    if (!request.current()) return;
    buildRightSvg.innerHTML = '<div class="box-placeholder">Error loading monomer</div>';
  } finally {
    request.finish();
  }
}

function renderRgroupButtons(container, state, side) {
  container.innerHTML = '';
  state.rgroups.forEach(rg => {
    const btn = document.createElement('button');
    btn.className = 'rgroup-btn';
    if (rg.used) btn.classList.add('used');
    if (state.selectedSlot === rg.slot) btn.classList.add('selected');
    btn.textContent = `R${rg.slot} ${rg.chem_type || ''}`;
    btn.title = `R${rg.slot}: ${rg.chem_type || 'unknown'}${rg.leaving ? ' (LG: ' + rg.leaving + ')' : ''}`;
    if (!rg.used) {
      btn.addEventListener('click', () => selectRgroup(side, rg.slot));
    }
    container.appendChild(btn);
  });
}

function selectRgroup(side, slot) {
  cancelRequests('sequence-edit');
  if (side === 'left' && buildLeft) {
    buildLeft.selectedSlot = buildLeft.selectedSlot === slot ? null : slot;
    renderRgroupButtons(buildLeftRg, buildLeft, 'left');
    if (rxnFilterActive && libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  } else if (side === 'right' && buildRight) {
    buildRight.selectedSlot = buildRight.selectedSlot === slot ? null : slot;
    renderRgroupButtons(buildRightRg, buildRight, 'right');
  }
  checkBuildValidity();
}

async function checkBuildValidity() {
  cancelRequests('bond-check');
  buildConnect.disabled = true;
  if (!buildLeft?.selectedSlot || !buildRight?.selectedSlot) {
    buildConnect.disabled = true;
    buildStatus.textContent = '';
    buildStatus.className = 'build-status';
    return;
  }

  const leftRg = buildLeft.rgroups.find(r => r.slot === buildLeft.selectedSlot);
  const rightRg = buildRight.rgroups.find(r => r.slot === buildRight.selectedSlot);
  if (!leftRg || !rightRg) return;

  const request = startRequest('bond-check');

  buildStatus.textContent = 'Checking bond...';
  buildStatus.className = 'build-status';

  try {
    const res = await fetch('/validate_bond', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        chem_type_a: leftRg.chem_type,
        chem_type_b: rightRg.chem_type,
        abbr_a: buildLeft.abbr,
        slot_a: buildLeft.selectedSlot,
        abbr_b: buildRight.abbr,
        slot_b: buildRight.selectedSlot
      }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current()) return;
    if (data.valid) {
      buildStatus.textContent = `Valid: ${data.reaction || 'bond'} (R${buildLeft.selectedSlot}↔R${buildRight.selectedSlot})`;
      buildStatus.className = 'build-status valid';
      buildConnect.disabled = false;
    } else {
      buildStatus.textContent = data.error || data.reason || 'No compatible reaction found';
      buildStatus.className = 'build-status invalid';
      buildConnect.disabled = true;
    }
  } catch (e) {
    if (!request.current()) return;
    buildStatus.textContent = 'Validation error';
    buildStatus.className = 'build-status invalid';
    buildConnect.disabled = true;
  } finally {
    request.finish();
  }
}

buildConnect.addEventListener('click', async () => {
  if (!buildLeft || !buildRight || !buildLeft.selectedSlot || !buildRight.selectedSlot) return;

  const val = cabilnInput.value.trim();
  const rHost = buildLeft.selectedSlot;
  const rNew = buildRight.selectedSlot;
  const newAbbr = buildRight.abbr;
  const request = startRequest('sequence-edit');

  buildConnect.disabled = true;
  buildStatus.textContent = 'Inserting...';
  buildStatus.className = 'build-status';

  try {
    const res = await fetch('/insert_bond', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        cabiln: val,
        host_residue_idx: buildLeftRIdx ?? 0,
        new_abbr: newAbbr,
        r_host: rHost,
        r_new: rNew,
        target_residue_idx: (buildRightRIdx !== null && buildRightRIdx !== buildLeftRIdx)
                             ? buildRightRIdx : -1
      }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== val) return;
    if (data.error) {
      buildStatus.textContent = data.error;
      buildStatus.className = 'build-status invalid';
      return;
    }
    cabilnInput.value = data.result;
    cabilnInput.dispatchEvent(new Event('input'));

    // Clear right side for next addition
    buildRight = null;
    buildRightAbbr.textContent = '—';
    buildRightSvg.innerHTML = '<div class="box-placeholder">Right-click from library</div>';
    buildRightRg.innerHTML = '';
    buildConnect.disabled = true;
    buildStatus.textContent = '';
    buildStatus.className = 'build-status';
    buildHint.textContent = 'Connection added — select a chip and right-click another monomer';
  } catch (e) {
    if (!request.current()) return;
    buildStatus.textContent = 'Insert failed';
    buildStatus.className = 'build-status invalid';
  } finally {
    request.finish();
  }
});

// ─── notation conversion ──────────────────────────────────────────────────────
async function convertNotation(target) {
  const val = cabilnInput.value.trim();
  if (!val) return;
  const request = startRequest('sequence-edit');
  try {
    const res = await fetch('/convert_notation', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ cabiln: val, target }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== val) return;
    if (data.error) { console.error('convert error:', data.error); return; }
    if (data.result) {
      cabilnInput.value = data.result;
      cabilnInput.dispatchEvent(new Event('input'));
    }
  } catch (e) {
    if (request.current()) console.error('convertNotation error:', e);
  } finally {
    request.finish();
  }
}
btnToBracket.addEventListener('click', () => convertNotation('bracket'));
btnToBranch.addEventListener('click', () => convertNotation('branch'));

btnReroll.addEventListener('click', () => {
  if (!lastCabiln) return;
  rerollSeed++;
  btnReroll.textContent = rerollSeed % 2 === 1 ? '⟳ Indigo' : '⟳ CoordGen';
  showSpinner(renderInner);
  doRenderCabiln(lastCabiln);
});

// ─── SMILES → CABILN conversion ───────────────────────────────────────────────
function startConversion(button) {
  cancelRequests('sequence-edit');
  const label = button.textContent;
  const request = startRequest('sequence-edit', () => {
    button.textContent = label;
    button.disabled = false;
  });
  button.textContent = '…';
  button.disabled = true;
  return request;
}

function useConvertedCabiln(sequence, warning = '') {
  notationSelect.value = 'cabiln';
  btnToCabilnPct.style.display = 'none';
  btnToCabilnBracket.style.display = 'none';
  cabilnInput.placeholder = NOTATION_PLACEHOLDER.cabiln;
  cabilnInput.value = sequence;
  cabilnInput.dispatchEvent(new Event('input'));
  conversionStatus.textContent = warning;
  conversionStatus.title = warning;
  conversionStatus.hidden = !warning;
}

async function doS2c(notation) {
  const smiles = smilesInput.value.trim();
  if (!smiles) return;
  const originalMain = cabilnInput.value.trim();
  const btn = notation === 'bracket' ? btnS2cBracket : btnS2c;
  const request = startConversion(btn);
  try {
    const res = await fetch('/smiles_to_cabiln', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles, notation }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || smilesInput.value.trim() !== smiles ||
        cabilnInput.value.trim() !== originalMain) return;
    if (data.error) {
      smilesStatus.textContent = 'S2C: ' + data.error;
      smilesStatus.className = 'statusbar';
    } else {
      useConvertedCabiln(data.cabiln, data.warning);
      if (data.warning) {
        smilesStatus.textContent = `⚠ ${data.warning}`;
        smilesStatus.className = 'statusbar warn';
      } else {
        const count = data.assignments?.length ?? data.details.length;
        smilesStatus.textContent = `Converted (${notation}): ${count} monomer(s)`;
        smilesStatus.className = 'statusbar ok';
      }
    }
  } catch (e) {
    if (!request.current()) return;
    smilesStatus.textContent = 'S2C error — is server running?';
    smilesStatus.className = 'statusbar';
  } finally {
    request.finish();
  }
}
btnS2c.addEventListener('click', () => doS2c('percent'));
btnS2cBracket.addEventListener('click', () => doS2c('bracket'));

// ─── main input → CABILN convert buttons ─────────────────────────────────────
async function doToCabiln(notation) {
  const txt = cabilnInput.value.trim();
  if (!txt) return;
  const btn = notation === 'bracket' ? btnToCabilnBracket : btnToCabilnPct;
  const request = startConversion(btn);
  try {
    const res = await fetch('/to_cabiln', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt, notation }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== txt) return;
    if (data.error) {
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
    } else {
      useConvertedCabiln(data.cabiln, data.warning);
      cabilnStatus.textContent = `Converted from ${data.from}: ${data.cabiln}`;
      cabilnStatus.className = 'statusbar ok';
    }
  } catch (e) {
    if (!request.current()) return;
    cabilnStatus.textContent = '→ CABILN error — is server running?';
    cabilnStatus.className = 'statusbar';
  } finally {
    request.finish();
  }
}
btnToCabilnPct.addEventListener('click',     () => doToCabiln('percent'));
btnToCabilnBracket.addEventListener('click', () => doToCabiln('bracket'));

function setExportReady(svg, molBlock) {
  lastSvg      = svg || '';
  lastMolBlock = molBlock || '';
  btnPng.disabled = !lastSvg;
  btnMol.disabled = !lastMolBlock;
}

function clearExports() {
  lastSvg = lastMolBlock = '';
  btnPng.disabled = true;
  btnMol.disabled = true;
}

// ─── render helpers ───────────────────────────────────────────────────────────
function canvasSize(el) {
  return { w: Math.max(el.clientWidth || 600, 400),
           h: Math.max(el.clientHeight || 500, 300) };
}

function setInner(inner, html) {
  inner.innerHTML = html;
  inner.style.transform = '';
}

function showSpinner(inner) {
  setInner(inner, '<div class="spinner"></div>');
}

// ─── notation selector ────────────────────────────────────────────────────────
const NOTATION_PLACEHOLDER = {
  cabiln: 'e.g.  fmoc-A-G-L-am\nfmoc-C.trt(4,2)-A-K.boc(4,2)-am\nfmoc-K.!1(4,4)-G-G-E.!1-am',
  smiles: 'Paste SMILES here… e.g. O=C1CNC(=O)[C@@H](C)N1',
  biln:   'Paste BILN here… e.g. fmoc-A-G-L-am  (use Token(bid,rg) for crosslinks)',
  helm:   'Paste HELM here… e.g. PEPTIDE1{A.G.L}$$$$',
};

notationSelect.addEventListener('change', () => {
  const mode = notationSelect.value;
  cabilnInput.placeholder = NOTATION_PLACEHOLDER[mode] || '';
  const showConvert = mode !== 'cabiln';
  btnToCabilnPct.style.display     = showConvert ? '' : 'none';
  btnToCabilnBracket.style.display = showConvert ? '' : 'none';
  cabilnInput.value = '';
  cabilnInput.className = '';
  resetCabiln();
});

// ─── CABILN render ────────────────────────────────────────────────────────────
cabilnInput.addEventListener('input', () => {
  const seq = cabilnInput.value.trim();
  const mode = notationSelect.value;
  resetCabiln();
  if (!seq) return;
  if (mode !== 'cabiln') {
    // Non-CABILN: render via render_reference for a live preview
    showSpinner(renderInner);
    cabilnTimer = setTimeout(() => doRenderForeign(seq), 400);
    return;
  }
  showSpinner(renderInner);
  cabilnTimer = setTimeout(() => doRenderCabiln(seq), 300);
});

async function doRenderForeign(txt) {
  const request = startRequest('main-render');
  const mode = notationSelect.value;
  lastCabiln = '';
  clearExports();
  clearComparison();
  btnReroll.disabled = true;
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res = await fetch('/render_reference', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt, width: w, height: h }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== txt ||
        notationSelect.value !== mode) return;
    if (data.error) {
      setInner(renderInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
      cabilnInput.className = 'err';
    } else {
      setInner(renderInner, data.svg);
      cabilnStatus.textContent = `${data.format}: ${data.info || ''}`;
      cabilnStatus.className = 'statusbar ok';
      cabilnInput.className = 'ok';
    }
  } catch (e) {
    if (!request.current()) return;
    cabilnStatus.textContent = 'Server error';
    cabilnStatus.className = 'statusbar';
  } finally {
    request.finish();
  }
}

function resetCabiln() {
  clearTimeout(cabilnTimer);
  cancelRequests('main-render', 'sequence-edit');
  clearComparison();
  clearBuild();
  lastCabiln = '';
  rerollSeed = 0;
  btnReroll.disabled = true;
  btnReroll.textContent = '⟳ Layout';
  cabilnInput.className = '';
  cabilnStatus.textContent = '';
  cabilnStatus.className = 'statusbar';
  conversionStatus.textContent = '';
  conversionStatus.hidden = true;
  setInner(renderInner, '<div class="placeholder">Start typing a sequence…</div>');
  clearExports();
  resChips.innerHTML = '';
  residueMap = {}; atomToRes = {}; residueList = [];
  chainData = [];
  currentBranchSet = new Set();
}

async function doRenderCabiln(seq) {
  const request = startRequest('main-render');
  const _seqChanged = seq !== lastCabiln;
  lastCabiln = '';
  clearExports();
  clearComparison();
  btnReroll.disabled = true;
  if (buildMode && _seqChanged) clearBuild();
  resChips.innerHTML = '';
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res  = await fetch('/render', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: seq, width: w, height: h, seed: rerollSeed }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== seq ||
        notationSelect.value !== 'cabiln') return;
    if (data.error) {
      const html = `<div class="placeholder err">${escHtml(data.error)}</div>`;
      setInner(renderInner, html);
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
      cabilnInput.className = 'err';
      clearExports();
      btnReroll.disabled = true;
      resChips.innerHTML = '';
    } else {
      const displayedSequence = data.normalized_cabiln || seq;
      if (data.normalized_cabiln) cabilnInput.value = displayedSequence;
      lastCabiln = displayedSequence;
      setInner(renderInner, data.svg);
      cabilnStatus.textContent = [data.info, ...(data.warnings || [])].filter(Boolean).join(' · ');
      cabilnStatus.title = cabilnStatus.textContent;
      cabilnStatus.className = data.warnings?.length ? 'statusbar warn' : 'statusbar ok';
      cabilnInput.className = 'ok';
      btnReroll.disabled = false;
      setExportReady(data.svg, data.mol_block);
      buildResidueUI(data.residue_map, data.residues, data.layout, data.crosslink_groups);
      if (verifyMode && lastSmiles) triggerVerify();
    }
  } catch (e) {
    if (!request.current()) return;
    cabilnStatus.textContent = 'Server error — is the renderer running?';
    cabilnStatus.className = 'statusbar';
  } finally {
    request.finish();
  }
}

// ─── reference render (verify mode) — auto-detects SMILES / BILN / HELM ─────
function clearReference() {
  clearTimeout(smilesTimer);
  cancelRequests('reference-render', 'sequence-edit');
  lastSmiles = '';
  smilesInput.className = '';
  smilesStatus.textContent = '';
  smilesStatus.className = 'statusbar';
  clearComparison();
}

smilesInput.addEventListener('input', () => {
  clearReference();
  const txt = smilesInput.value.trim();
  if (!txt) {
    setInner(smilesInner, '<div class="placeholder">Paste SMILES, BILN, or HELM…</div>');
    compareBar.innerHTML = '';
    return;
  }
  showSpinner(smilesInner);
  smilesTimer = setTimeout(() => doRenderRef(txt), 300);
});

molUpload.addEventListener('change', async (e) => {
  const file = e.target.files[0];
  if (!file) return;
  // The selected File is retained locally; allow choosing it again after cancel.
  molUpload.value = '';
  clearReference();
  const request = startRequest('reference-render');
  showSpinner(smilesInner);
  smilesStatus.textContent = `Loaded: ${file.name}`;
  smilesStatus.className = 'statusbar ok';
  try {
    const text = await file.text();
    if (!request.current()) return;
    const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
    const res = await fetch('/render_mol', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ mol_block: text, width: w, height: h }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current()) return;
    if (data.error) {
      setInner(smilesInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      smilesStatus.textContent = data.error;
      smilesStatus.className = 'statusbar';
    } else {
      setInner(smilesInner, data.svg);
      lastSmiles = data.smiles || '';
      smilesInput.value = lastSmiles;
      smilesInput.className = 'ok';
      if (lastCabiln) triggerVerify();
    }
  } catch (err) {
    if (!request.current()) return;
    smilesStatus.textContent = 'Failed to render .mol file';
    smilesStatus.className = 'statusbar';
  } finally {
    request.finish();
  }
});

async function doRenderRef(txt) {
  const request = startRequest('reference-render');
  lastSmiles = '';
  clearComparison();
  const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
  try {
    const res = await fetch('/render_reference', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt, width: w, height: h }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || smilesInput.value.trim() !== txt) return;
    if (data.error) {
      setInner(smilesInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      smilesStatus.textContent = data.error;
      smilesStatus.className = 'statusbar';
      smilesInput.className = 'err';
    } else {
      setInner(smilesInner, data.svg);
      lastSmiles = data.smiles || '';
      smilesStatus.textContent = `${data.format}: ${data.info || ''}`;
      smilesStatus.className = 'statusbar ok';
      smilesInput.className = 'ok';
      if (lastCabiln) triggerVerify();
    }
  } catch (e) {
    if (!request.current()) return;
    smilesStatus.textContent = 'Server error';
    smilesStatus.className = 'statusbar';
  } finally {
    request.finish();
  }
}

// ─── verify comparison ────────────────────────────────────────────────────────
function clearComparison() {
  cancelRequests('verify');
  compareBar.innerHTML = '';
}

async function triggerVerify() {
  clearComparison();
  if (!verifyMode || !lastSmiles || !lastCabiln) return;
  const request = startRequest('verify');
  const smiles = lastSmiles;
  const cabiln = lastCabiln;
  try {
    const res  = await fetch('/verify', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles, cabiln }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || lastSmiles !== smiles || lastCabiln !== cabiln) return;
    if (data.error) {
      compareBar.innerHTML = `<span class="nomatch">Error: ${escHtml(data.error)}</span>`;
      return;
    }
    const badge = data.match
      ? '<span class="match">✓ EXACT MATCH</span>'
      : data.stereo_compatible
        ? '<span class="warn-badge">Compatible; unspecified stereo</span>'
        : '<span class="nomatch">✗ MISMATCH</span>';
    const warnBadge = data.warning
      ? `<span class="warn-badge" title="${escAttr(data.warning)}">⚠ ${escHtml(data.warning)}</span>`
      : '';
    compareBar.innerHTML =
      badge + warnBadge +
      `<span class="canon" title="SMILES canonical: ${escAttr(data.smiles_canonical)}">` +
      `SMILES: ${escHtml(data.smiles_canonical.slice(0, 80))}${data.smiles_canonical.length > 80 ? '…' : ''}</span>` +
      `<span class="canon" title="CABILN canonical: ${escAttr(data.cabiln_canonical)}">` +
      `CABILN: ${escHtml(data.cabiln_canonical.slice(0, 80))}${data.cabiln_canonical.length > 80 ? '…' : ''}</span>`;
  } catch (e) {
    if (!request.current()) return;
    compareBar.innerHTML = '<span class="nomatch">Verify error</span>';
  } finally {
    request.finish();
  }
}

function escHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
function escAttr(s) { return escHtml(s); }

// Auto-reload on server restart
(function(){
  let sid = null;
  setInterval(async () => {
    try {
      const r = await fetch('/server_id');
      const id = await r.text();
      if (sid === null) { sid = id; return; }
      if (id !== sid) location.reload();
    } catch(e) {}
  }, 5000);
})();
