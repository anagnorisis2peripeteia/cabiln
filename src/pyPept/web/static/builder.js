// ─── state ────────────────────────────────────────────────────────────────────
let darkMode   = true;
let verifyMode = false;
let hlEnabled  = true;
let libLoaded  = false;
let libraryVersion = null;
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
let mainStale = false;
let hasMainDrawing = false;
let displayedSource = '';
let displayedNotation = '';
let projectContext = null;
let projectRevision = 0;
let draftCleared = false;
let referenceOriginal = null;
let referenceContext = null;

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
const notationPolicy = document.getElementById('notation-policy');
const libPanel      = document.getElementById('lib-panel');
const libSearch     = document.getElementById('lib-search');
const libClose      = document.getElementById('lib-close');
const libList       = document.getElementById('lib-list');
const libCount      = document.getElementById('lib-count');
const libStatus     = document.getElementById('lib-status');
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
const btnUndo = document.getElementById('btn-undo');
const btnRedo = document.getElementById('btn-redo');
const draftNotice = document.getElementById('draft-notice');
const draftStatus = document.getElementById('draft-status');
const btnRestoreDraft = document.getElementById('btn-restore-draft');
const btnDismissDraft = document.getElementById('btn-dismiss-draft');
const renderPane = document.getElementById('render-pane');
const renderProgress = document.getElementById('render-progress');
const renderProgressLabel = document.getElementById('render-progress-label');
const conversionProgress = document.getElementById('conversion-progress');
const conversionProgressLabel = document.getElementById('conversion-progress-label');
const btnProjectSave = document.getElementById('btn-project-save');
const projectUpload = document.getElementById('project-upload');
const projectStatus = document.getElementById('project-status');
const importQuality = document.getElementById('import-quality');
const qualitySummary = document.getElementById('quality-summary');
const qualityDetails = document.getElementById('quality-details');
const canonicalStatus = document.getElementById('canonical-status');
const btnHelp = document.getElementById('btn-help');
const helpPanel = document.getElementById('help-panel');

// Document data belongs to the editor; these timers belong to browser storage.
const DRAFT_KEY = 'cabiln.draft.v1';
const editor = new CabilnDocument();
let saveDraftTimer = null;
let savedDraft = null;

function displayDocument() {
  const state = editor.present;
  if (cabilnInput.value !== state.text) cabilnInput.value = state.text;
  notationSelect.value = state.notation;
  updateNotationControls();
  setConversionWarning(state.warning);
  displayDocumentEvidence();
}

function replaceDocument(patch, drafts = editor.drafts) {
  editor.replace(patch, drafts);
  displayDocument();
}

function scheduleDraftSave() {
  clearTimeout(saveDraftTimer);
  saveDraftTimer = setTimeout(saveDraft, 200);
}

function updateHistoryControls() {
  btnUndo.disabled = !editor.past.length;
  btnRedo.disabled = !editor.future.length;
}

function saveDraft() {
  clearTimeout(saveDraftTimer);
  if (savedDraft || draftCleared) return; // Recovery and explicit clearing stay intact.
  try {
    window.localStorage.setItem(DRAFT_KEY, JSON.stringify({ version: 1,
      document: editor.present, drafts: editor.drafts, reference: smilesInput.value,
      reference_original: referenceOriginal, reference_context: referenceContext,
      context: projectContext || editor.present.context }));
    draftStatus.textContent = 'Draft saved in this browser';
  } catch (error) {
    draftStatus.textContent = 'Draft storage is unavailable; this session still supports Undo';
  }
  draftNotice.hidden = false;
}

function recordDocument(next, typing = false) {
  if (editor.commit(next, typing)) {
    projectChanged();
    if (next.text) beginDraftEdit();
    scheduleDraftSave();
  }
  displayDocument();
  updateHistoryControls();
}

function updateNotationControls() {
  const mode = notationSelect.value;
  cabilnInput.placeholder = NOTATION_PLACEHOLDER[mode] || '';
  btnToCabilnPct.style.display = mode === 'cabiln' ? 'none' : '';
  btnToCabilnBracket.style.display = mode === 'cabiln' ? 'none' : '';
  btnToBracket.disabled = btnToBranch.disabled = mode !== 'cabiln';
  notationPolicy.disabled = mode !== 'cabiln';
}

function setConversionWarning(warning = '') {
  conversionStatus.textContent = warning;
  conversionStatus.title = warning;
  conversionStatus.hidden = !warning;
}

function commitDocument(text, notation = editor.present.notation, warning = '', evidence = {}) {
  recordDocument({ text, notation, warning, quality: evidence.quality || null,
    canonical: evidence.canonical || null, context: evidence.context || null });
  renderDocument(true);
}

function travelHistory(direction) {
  if (!editor.travel(direction)) return;
  projectChanged();
  displayDocument();
  updateHistoryControls();
  saveDraft();
  renderDocument(true);
}

btnUndo.addEventListener('click', () => travelHistory('undo'));
btnRedo.addEventListener('click', () => travelHistory('redo'));
window.addEventListener('keydown', event => {
  if (!(event.ctrlKey || event.metaKey) || event.altKey) return;
  const target = event.target;
  if (target !== cabilnInput && (target?.matches?.('input, textarea') || target?.isContentEditable)) return;
  const key = event.key.toLowerCase();
  if (key !== 'z' && !(event.ctrlKey && key === 'y')) return;
  event.preventDefault();
  travelHistory(key === 'y' || event.shiftKey ? 'redo' : 'undo');
});
window.addEventListener('pagehide', saveDraft);

function finishDraftRecovery() {
  savedDraft = null;
  btnRestoreDraft.hidden = btnDismissDraft.hidden = true;
}
function beginDraftEdit() {
  draftCleared = false;
  if (!savedDraft) return;
  finishDraftRecovery();
  draftStatus.textContent = 'New draft started in this browser';
}
btnRestoreDraft.addEventListener('click', async () => {
  if (!savedDraft) return;
  const saved = savedDraft;
  const project = projectSnapshot();
  project.document = saved.document;
  project.drafts = saved.drafts;
  project.reference = saved.reference;
  project.context = saved.context;
  if (project.context) {
    // Autosaved drafts can include edits since their last server signature.
    // Preparation checks their stored binding before stamping the current text.
    await restoreBoundDraft(project);
  } else {
    // Earlier releases did not bind drafts. Never imply their definitions were checked.
    applyProject(project);
    setProjectStatus('Older draft restored without a library binding. Check any custom monomer definitions before using it.', true);
  }
});
btnDismissDraft.addEventListener('click', () => {
  cancelRequests('project-open');
  finishDraftRecovery();
  saveDraft();
});

function offerSavedDraft() {
  btnRestoreDraft.hidden = btnDismissDraft.hidden = true;
  try {
    const saved = JSON.parse(window.localStorage.getItem(DRAFT_KEY));
    if (!saved || saved.version !== 1 || !saved.document ||
        typeof saved.document.text !== 'string' ||
        !Object.hasOwn(NOTATION_PLACEHOLDER, saved.document.notation)) return;
    const drafts = CabilnProject.drafts(saved.drafts);
    const reference = CabilnProject.reference({ text: saved.reference,
      original: saved.reference_original, context: saved.reference_context }, saved.context);
    if (!saved.document.text && !Object.values(drafts).some(draft => draft.text) &&
        !reference.text && !reference.original) return;
    savedDraft = { document: CabilnProject.document(saved.document), drafts, reference,
      context: saved.context || null };
    draftStatus.textContent = 'A saved draft is available in this browser';
    draftNotice.hidden = btnRestoreDraft.hidden = btnDismissDraft.hidden = false;
  } catch (error) { /* Storage can be disabled, full, or contain an older format. */ }
}

function contextHeader(context) {
  return context ? CabilnProject.clone({ project_version: context.project_version,
    library_binding: context.library_binding, canonical: context.canonical }) : null;
}

function displayDocumentEvidence() {
  const quality = editor.present.quality;
  importQuality.hidden = !quality;
  qualitySummary.textContent = CabilnProject.qualitySummary(quality);
  importQuality.dataset.status = quality && (quality.recognition_status !== 'complete' ||
    quality.search_complete === false || quality.inferred_stereo === true) ? 'review' : 'complete';
  if (quality) {
    const notes = [...(quality.warnings || [])];
    if (quality.search_complete === false) notes.push('The search reached a resource limit; another valid monomer interpretation may exist.');
    if (quality.inferred_stereo === true) notes.push('Some stereochemistry was inferred from the library rather than specified by the input.');
    if (quality.assignments?.some(item => item.recognized === false) || quality.recognition_status === 'unresolved') {
      notes.push('Synthetic or opaque regions preserve input chemistry without claiming a library match. Available attachment sites remain editable.');
    }
    qualityDetails.innerHTML = '<p>Import evidence for this source. Editing or reordering the sequence clears this report; Undo restores it.</p>' +
      (notes.length ? '<ul>' + [...new Set(notes)].map(note => `<li>${escHtml(note)}</li>`).join('') + '</ul>' :
        '<p>All reported monomer assignments matched the library. Search completeness describes this recognition run, not a unique chemical interpretation.</p>');
  } else qualityDetails.innerHTML = '';
  const canonical = editor.present.canonical;
  canonicalStatus.hidden = !canonical;
  canonicalStatus.textContent = canonical
    ? `Canonical export: ${canonical.format || 'recorded convention'} · RDKit ${canonical.rdkit || 'recorded version'} · library binding retained in projects`
    : '';
}

function acceptDocumentContext(context) {
  if (!context?.library_binding) return;
  retainReferenceContext();
  const current = editor.present;
  const changed = !CabilnProject.sameContext(current.context, context);
  const patch = { context: changed ? contextHeader(context) : current.context };
  if (changed) projectChanged(false);
  if (changed && (current.quality || current.canonical)) {
    Object.assign(patch, { warning: '', quality: null, canonical: null });
    setProjectStatus('The library binding or convention changed. Import evidence was cleared; convert the original input again to refresh it.', true);
  }
  if (!CabilnProject.sameContext(projectContext, context)) {
    projectContext = contextHeader(context);
  }
  replaceDocument(patch);
  saveDraft();
}

function retainReferenceContext() {
  if (!referenceContext && (smilesInput.value || referenceOriginal) && projectContext) {
    referenceContext = CabilnProject.clone(projectContext);
  }
}

function acceptReferenceContext(context) {
  if (!context?.library_binding) return;
  if (!CabilnProject.sameContext(referenceContext, context)) {
    projectChanged(false);
    referenceContext = contextHeader(context);
  }
  if (!CabilnProject.sameContext(projectContext, context)) projectContext = contextHeader(context);
  saveDraft();
}

function projectChanged(edited = true) {
  projectRevision++;
  if (edited) draftCleared = false;
  if (activeRequests.has('project-open') || activeRequests.has('project-save')) {
    setProjectStatus('The document changed during the project check. Save or open again when ready.', true);
  }
  cancelRequests('project-open', 'project-save');
}

function setProjectStatus(message, error = false) {
  projectStatus.textContent = message;
  projectStatus.className = error ? 'statusbar warn' : 'statusbar ok';
  projectStatus.hidden = !message;
}

function projectSnapshot() {
  return { format: CabilnProject.FORMAT, version: CabilnProject.VERSION,
    document: CabilnProject.clone(editor.present), drafts: CabilnProject.clone(editor.drafts),
    reference: { text: smilesInput.value, original: CabilnProject.clone(referenceOriginal),
      context: CabilnProject.clone(referenceContext) },
    context: CabilnProject.clone(projectContext || editor.present.context),
    saved_at: new Date().toISOString() };
}

function downloadProject(project) {
  const blob = new Blob([JSON.stringify(project)], { type: 'application/json' });
  if (blob.size > CabilnProject.MAX_FILE_BYTES) {
    throw new Error('The project exceeds 2 MiB. Shorten unused drafts or the reference before saving.');
  }
  const link = document.createElement('a');
  const url = URL.createObjectURL(blob);
  link.href = url;
  link.download = 'peptide.cabiln.json';
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

async function restoreBoundDraft(project) {
  const revision = projectRevision;
  const request = startRequest('project-open');
  setProjectStatus('Checking the saved draft library binding and definitions…');
  try {
    const response = await fetchCalculation('/prepare_project', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ project }), signal: request.signal,
    });
    const data = await readResponse(response);
    if (!request.current() || revision !== projectRevision) return;
    if (data.error) throw new Error(data.error);
    const prepared = CabilnProject.read(data.project);
    request.finish();
    applyProject(prepared);
    setProjectStatus('Saved draft restored with its source, reference and import details.');
  } catch (error) {
    if (request.current()) setProjectStatus((error.message || 'Could not restore the saved draft.') +
      ' Your current work is unchanged.', true);
  } finally { request.finish(); }
}

async function saveProject() {
  const revision = projectRevision;
  const project = projectSnapshot();
  const request = startRequest('project-save', () => { btnProjectSave.disabled = false; });
  btnProjectSave.disabled = true;
  setProjectStatus('Checking project definitions before saving…');
  try {
    const response = await fetchCalculation('/prepare_project', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ project }), signal: request.signal,
    });
    const data = await readResponse(response);
    if (!request.current() || revision !== projectRevision) return;
    if (data.error) throw new Error(data.error);
    const prepared = CabilnProject.read(data.project);
    downloadProject(prepared);
    request.finish();
    projectContext = prepared.context;
    referenceContext = prepared.reference.context || prepared.context;
    replaceDocument({ context: prepared.document.context || prepared.context }, prepared.drafts);
    saveDraft();
    setProjectStatus('Project saved. It includes source, notation drafts, original reference and import details.');
  } catch (error) {
    if (request.current()) setProjectStatus(error.message || 'Could not save the project. Your work is unchanged.', true);
  } finally { request.finish(); }
}
btnProjectSave.addEventListener('click', saveProject);

async function validateAndOpenProject(project, successMessage = 'Project opened') {
  const revision = projectRevision;
  const request = startRequest('project-open');
  setProjectStatus('Checking the project library binding and definitions…');
  try {
    const response = await fetchCalculation('/validate_project', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ project }), signal: request.signal,
    });
    const data = await readResponse(response);
    if (!request.current() || revision !== projectRevision) return;
    if (data.error || data.valid !== true) throw new Error(data.error || 'The project could not be validated.');
    request.finish();
    // Validation proves every saved document against these exact definitions.
    const context = data.context || project.context;
    applyProject({ ...project, context, reference: { ...project.reference, context },
      document: { ...project.document, context },
      drafts: Object.fromEntries(Object.entries(project.drafts).map(([mode, state]) =>
        [mode, { ...state, context }])) });
    setProjectStatus(successMessage + '. Source and reference are preserved; use Undo to return to the previous sequence.');
  } catch (error) {
    if (request.current()) setProjectStatus((error.message || 'Could not open the project.') + ' Your current work is unchanged.', true);
  } finally { request.finish(); }
}

projectUpload.addEventListener('change', async event => {
  const file = event.target.files?.[0];
  projectUpload.value = '';
  if (!file) return;
  const revision = projectRevision;
  const request = startRequest('project-open');
  setProjectStatus('Reading project…');
  try {
    if (file.size > CabilnProject.MAX_FILE_BYTES) throw new Error('Project files must be no larger than 2 MiB.');
    const text = await file.text();
    if (!request.current() || revision !== projectRevision) return;
    const project = CabilnProject.read(JSON.parse(text));
    request.finish();
    await validateAndOpenProject(project);
  } catch (error) {
    if (request.current()) setProjectStatus((error instanceof SyntaxError ? 'This file is not valid project JSON.' : error.message) +
      ' Your current work is unchanged.', true);
  } finally { request.finish(); }
});

function applyProject(project) {
  finishDraftRecovery();
  editor.restoreDrafts(CabilnProject.drafts(project.drafts));
  projectContext = CabilnProject.clone(project.context);
  clearReference();
  const reference = CabilnProject.reference(project.reference, project.context);
  smilesInput.value = reference.original?.kind === 'text' ? reference.original.content : reference.text;
  referenceOriginal = reference.original;
  referenceContext = reference.context;
  setInner(smilesInner, '<div class="placeholder">Reference restored. Open Verify to compare.</div>');
  const state = CabilnProject.document(project.document);
  commitDocument(state.text, state.notation, state.warning, state);
  saveDraft();
  if (verifyMode) restoreReferenceDrawing();
}

function restoreReferenceDrawing() {
  if (referenceOriginal?.kind === 'mol') {
    renderMolReference(referenceOriginal.content, referenceOriginal.name);
  } else {
    if (referenceOriginal?.kind === 'text' && smilesInput.value !== referenceOriginal.content) {
      projectChanged(false);
      smilesInput.value = referenceOriginal.content;
      referenceContext = contextHeader(referenceContext);
      saveDraft();
    }
    if (smilesInput.value.trim()) doRenderRef(smilesInput.value.trim());
  }
}

function showHelp(open) {
  helpPanel.hidden = !open;
  btnHelp.setAttribute('aria-expanded', String(open));
  if (!open) btnHelp.focus();
}
btnHelp.addEventListener('click', () => showHelp(helpPanel.hidden));
document.getElementById('help-close').addEventListener('click', () => showHelp(false));
window.addEventListener('keydown', event => {
  if (event.key === 'Escape' && !helpPanel.hidden) showHelp(false);
});
document.getElementById('btn-clear-draft').addEventListener('click', () => {
  cancelRequests('project-open');
  clearTimeout(saveDraftTimer);
  finishDraftRecovery();
  draftCleared = true;
  try {
    window.localStorage.removeItem(DRAFT_KEY);
    draftStatus.textContent = 'Saved browser draft cleared. Current work and Undo are unchanged; the next edit starts a new draft.';
  } catch (error) { draftStatus.textContent = 'Browser storage is unavailable. No draft could be removed.'; }
  draftNotice.hidden = false;
});

async function loadCapabilities() {
  const link = document.getElementById('register-link');
  link.hidden = true;
  try {
    const response = await fetchCalculation('/capabilities');
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
  btnDark.setAttribute('aria-pressed', String(darkMode));
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
  btnHl.setAttribute('aria-pressed', String(hlEnabled));
  if (!hlEnabled) clearHighlight();
});

// ─── verify mode ──────────────────────────────────────────────────────────────
btnVerify.addEventListener('click', () => {
  verifyMode = !verifyMode;
  btnVerify.classList.toggle('active', verifyMode);
  btnVerify.setAttribute('aria-expanded', String(verifyMode));
  document.getElementById('verify-pane').style.display = verifyMode ? '' : 'none';
  compareBar.style.display = verifyMode ? '' : 'none';
  clearComparison();
  if (verifyMode && (smilesInput.value.trim() || referenceOriginal) && !lastSmiles) restoreReferenceDrawing();
  else if (verifyMode) triggerVerify();
});

// ─── library sidebar (persistent) ────────────────────────────────────────────
function openLib() {
  libPanel.classList.add('open');
  btnLib.classList.add('active');
  btnLib.setAttribute('aria-expanded', 'true');
  loadMonomers();
  libSearch.focus();
  loadReactions();
}
function closeLib() {
  libPanel.classList.remove('open');
  btnLib.classList.remove('active');
  btnLib.setAttribute('aria-expanded', 'false');
  hidePreview();
}

btnLib.addEventListener('click', () =>
  libPanel.classList.contains('open') ? closeLib() : openLib());
libClose.addEventListener('click', closeLib);

libSearch.addEventListener('input', () => renderLibList(libSearch.value.trim().toLowerCase()));

btnRxnFilter.addEventListener('click', async () => {
  if (!Array.isArray(reactionPairs)) {
    await loadReactions();
    if (!Array.isArray(reactionPairs) || !buildLeft) return;
  }
  rxnFilterActive = !rxnFilterActive;
  btnRxnFilter.classList.toggle('active', rxnFilterActive);
  btnRxnFilter.setAttribute('aria-pressed', String(rxnFilterActive));
  renderLibList(libSearch.value.trim().toLowerCase());
});

// ─── example peptide sidebar ──────────────────────────────────────────────────
let examplesLoaded = false;

function openExamples() {
  examplesPanel.classList.add('open');
  btnExamples.classList.add('active');
  btnExamples.setAttribute('aria-expanded', 'true');
  if (!examplesLoaded) loadExamples();
}
function closeExamples() {
  examplesPanel.classList.remove('open');
  btnExamples.classList.remove('active');
  btnExamples.setAttribute('aria-expanded', 'false');
}

btnExamples.addEventListener('click', () =>
  examplesPanel.classList.contains('open') ? closeExamples() : openExamples());
examplesClose.addEventListener('click', closeExamples);

async function loadExamples() {
  try {
    const res = await fetchCalculation('/examples');
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
      commitDocument(row.dataset.cabiln, 'cabiln');
      cabilnInput.focus();
    });
  });
}

async function loadMonomers() {
  const request = startRequest('library');
  if (!libLoaded) libList.innerHTML = '<div class="placeholder">Loading…</div>';
  try {
    const res = await fetchCalculation('/monomers', { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current()) return;
    if (!Array.isArray(data)) throw new Error(data.error || 'Invalid monomer library');
    const version = res.headers?.get('X-Library-Version') || null;
    if (libLoaded && version && version === libraryVersion) return;
    libraryVersion = version;
    allMonomers = data.map(monomer => ({ ...monomer,
      searchText: [monomer.abbr, monomer.name, monomer.type, monomer.chem_types].join(' ').toLowerCase(),
      attachmentTypes: parseCts(monomer.chem_types),
    }));
    previewCache = {};
    hidePreview();
    libLoaded = true;
    renderLibList(libSearch.value.trim().toLowerCase());
  } catch (e) {
    if (!request.current()) return;
    if (!libLoaded) libList.innerHTML = '<div class="placeholder err">Failed to load monomers. Close and reopen the library to retry.</div>';
  } finally {
    request.finish();
  }
}

window.addEventListener('focus', () => {
  if (libPanel.classList.contains('open')) loadMonomers();
});

async function loadReactions() {
  if (reactionPairs !== null) return;
  const request = startRequest('reactions');
  libStatus.hidden = true;
  try {
    const res = await fetchCalculation('/reactions', { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current()) return;
    if (!Array.isArray(data)) throw new Error(data.error || 'Reaction data is unavailable');
    reactionPairs = data;
    if (rxnFilterActive && libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  } catch (error) {
    if (!request.current()) return;
    reactionPairs = null;
    libStatus.textContent = 'Reaction filter unavailable. Reopen Library or click Filter to retry. ' + error.message;
    libStatus.hidden = false;
  } finally { request.finish(); }
}

function parseCts(cts) {
  if (!cts) return [];
  return cts.split(',').map(p => p.includes(':') ? p.split(':')[1].trim() : p.trim()).filter(Boolean);
}

function renderLibList(q) {
  let filtered = q
    ? allMonomers.filter(m => m.searchText.includes(q))
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
      const mcts = m.attachmentTypes;
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
      ? `<span class="lib-badge aa">${escHtml(m.type)}</span>`
      : m.type === 'cap' && m.subtype === 'protecting'
      ? `<span class="lib-badge protect">cap</span>`
      : `<span class="lib-badge cap">${escHtml(m.type)}</span>`;

    const issues = Array.isArray(m.quality?.issues) ? m.quality.issues : [];
    const qualityText = issues.map(issue => issue.message).filter(Boolean).join(' · ');
    let lg = m.leaving ? `  LG: ${escHtml(m.leaving)}` : '';
    if (m.degenerate) {
      const parts = [];
      if (m.nterm_abbr) parts.push(`N: ${escHtml(m.nterm_abbr)} (${escHtml(m.nterm_leaving)})`);
      if (m.cterm_abbr) parts.push(`C: ${escHtml(m.cterm_abbr)} (${escHtml(m.cterm_leaving)})`);
      lg = '  ' + parts.join(' | ');
    }
    const label = m.abbr + ': ' + m.name + (qualityText ? '. Library quality: ' + qualityText : '');
    return `<div class="lib-row" data-abbr="${escAttr(m.abbr)}" tabindex="0" aria-label="${escAttr(label)}">
      <div class="lib-abbr">${escHtml(m.abbr)}</div>
      <div class="lib-info">
        <div class="lib-name" title="${escAttr(m.name)}">${escHtml(m.name)}</div>
        <div class="lib-meta">${escHtml(m.chem_types || '')}${lg}</div>
        ${qualityText ? `<div class="lib-quality" title="${escAttr(qualityText)}">${escHtml(qualityText)}</div>` : ''}
      </div>
      ${badge}
      <button type="button" class="lib-use" aria-label="Use ${escAttr(m.abbr)} in builder" title="Choose this monomer in the builder">Use</button>
    </div>`;
  });
  libList.innerHTML = rows.join('');

}

function useLibraryMonomer(abbr, explicit = false) {
  if (explicit && !buildMode) openBuild();
  if (buildMode && notationSelect.value !== 'cabiln') {
    buildHint.textContent = 'Convert this input to CABILN before building';
  } else if (buildMode && insertBetweenActive) {
    doInsertBetween(abbr);
  } else if (buildMode && !cabilnInput.value.trim()) {
    commitDocument(abbr, 'cabiln');
    buildHint.textContent = 'First monomer added — select its residue, then choose another monomer';
  } else if (buildMode && mainStale) {
    buildHint.textContent = 'Wait for a valid drawing before choosing attachment sites';
  } else if (buildMode) {
    loadBuildRight(abbr);
    if (!buildLeft) buildHint.textContent = 'Now select a residue in the sequence';
  } else {
    insertAbbr(abbr);
  }
}

// Keep one set of listeners while search replaces the rows.
libList.addEventListener('click', event => {
  const row = event.target.closest('.lib-row');
  if (row) useLibraryMonomer(row.dataset.abbr, !!event.target.closest('.lib-use'));
});
libList.addEventListener('contextmenu', event => {
  const row = event.target.closest('.lib-row');
  if (row && buildMode) {
    event.preventDefault();
    useLibraryMonomer(row.dataset.abbr, true);
  }
});
libList.addEventListener('keydown', event => {
  if (event.target.matches('.lib-row') && ['Enter', ' '].includes(event.key)) {
    event.preventDefault();
    useLibraryMonomer(event.target.dataset.abbr);
  }
});
for (const type of ['mouseover', 'focusin']) {
  libList.addEventListener(type, event => {
    const row = event.target.closest('.lib-row');
    if (row && !row.contains(event.relatedTarget)) startPreview(row.dataset.abbr, row);
  });
}
for (const type of ['mouseout', 'focusout']) {
  libList.addEventListener(type, event => {
    const row = event.target.closest('.lib-row');
    if (row && !row.contains(event.relatedTarget)) hidePreview();
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
  commitDocument(before + insert + after);
  const newPos = start + insert.length;
  ta.setSelectionRange(newPos, newPos);
  ta.focus();
}

// ─── monomer preview tooltip ──────────────────────────────────────────────────
function startPreview(abbr, row) {
  clearTimeout(previewTimer);
  const request = startRequest('monomer-preview');
  previewTimer = setTimeout(async () => {
    try {
      const data = previewCache[abbr] || await readResponse(await fetchCalculation(
        `/monomer_svg?abbr=${encodeURIComponent(abbr)}`, { signal: request.signal }
      ));
      if (!request.current()) return;
      if (data.svg) {
        previewCache[abbr] = data;
        showPreview(data, row);
      }
    } catch (e) { /* silent */ }
    finally { request.finish(); }
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
  const quality = data.quality || allMonomers.find(item => item.abbr === row.dataset?.abbr)?.quality;
  const issues = Array.isArray(quality?.issues) ? quality.issues : [];
  if (issues.length) html += `<div class="prev-meta prev-warn">Library quality: ${issues.map(issue => escHtml(issue.message || issue.code)).join(' · ')}</div>`;
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
  cancelRequests('monomer-preview');
  libPreview.style.display = 'none';
}

// ─── residue chips + bidirectional highlighting ───────────────────────────────
let chainData = [];
let currentBranchSet = new Set();
let xlinkByRes = {};

function highlightGroup(idxList) {
  if (!hlEnabled || mainStale) return;
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
    const chip = document.createElement('button');
    chip.type = 'button';
    chip.className = 'res-chip';
    chip.textContent = r.abbr;
    chip.dataset.residue = r.idx;
    chip.title = `Select residue ${r.idx + 1}: ${r.abbr}`;
    const assignment = editor.present.quality?.assignments?.find(item => item.residue_index === r.idx);
    const kind = r.kind || r.quality_kind;
    const preserved = ['opaque', 'synthetic'].includes(kind) ? kind :
      !kind && (r.is_synthetic || assignment?.recognized === false) ? 'synthetic' : '';
    if (preserved) {
      chip.dataset.quality = preserved;
      chip.title += ` · ${preserved === 'opaque' ? 'Opaque preserved fragment' : 'Synthetic preserved region'}; available sites remain editable`;
      chip.setAttribute('aria-label', chip.title);
    }
    chip.style.background = RES_COLORS[r.colorIdx % RES_COLORS.length];
    const xlinks = xlinkByMember[r.idx];
    if (!simpleHover && xlinks && xlinks.length) {
      const allMembers = [...new Set(xlinks.flatMap(g => g.members))];
      chip.addEventListener('mouseenter', () => highlightGroup(allMembers));
    } else {
      chip.addEventListener('mouseenter', () => highlightResidue(r.idx));
    }
    chip.addEventListener('mouseleave', clearHighlight);
    chip.addEventListener('focus', () => highlightResidue(r.idx));
    chip.addEventListener('blur', clearHighlight);
    chip.addEventListener('click', () => selectBuildResidue(r));
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
      el.tabIndex = 0;
      el.setAttribute('aria-label', `Highlight ${text} group`);
      el.addEventListener('mouseenter', () => highlightGroup(memberIdxs));
      el.addEventListener('mouseleave', clearHighlight);
      el.addEventListener('focus', () => highlightGroup(memberIdxs));
      el.addEventListener('blur', clearHighlight);
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
    el.tabIndex = 0;
    el.setAttribute('aria-label', `Highlight connection ${tag}`);
    el.addEventListener('mouseenter', () => highlightGroup(members));
    el.addEventListener('mouseleave', clearHighlight);
    el.addEventListener('focus', () => highlightGroup(members));
    el.addEventListener('blur', clearHighlight);
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
  if (!hlEnabled || mainStale) return;
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

function svgResidueIndex(target, svg) {
  for (let el = target; el && el !== svg; el = el.parentElement) {
    const atom = (el.getAttribute('class') || '').match(/atom-(\d+)/);
    if (atom) {
      const rIdx = atomToRes[parseInt(atom[1])];
      if (rIdx !== undefined) return rIdx;
    }
  }
}

function wireUpSvgHover() {
  const svg = document.querySelector('#render-inner svg');
  if (!svg) return;
  svg.addEventListener('mousemove', e => {
    if (!hlEnabled || mainStale) return;
    const rIdx = svgResidueIndex(e.target, svg);
    if (rIdx === undefined) return clearHighlight();
    const xlinks = xlinkByRes[rIdx];
    if (xlinks && xlinks.length) {
      highlightGroup([...new Set(xlinks.flatMap(g => g.members))]);
    } else {
      highlightResidue(rIdx);
    }
  });
  svg.addEventListener('mouseleave', clearHighlight);
  svg.addEventListener('click', e => {
    if (!buildMode || mainStale) return;
    const rIdx = svgResidueIndex(e.target, svg);
    const residue = residueList.find(r => r.idx === rIdx);
    if (residue) selectBuildResidue(residue);
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
  return { reset() { scale = 1; tx = 0; ty = 0; applyTransform(); } };
}

const mainViewport = makeZoomable(renderCanvas, renderInner);
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
function selectBuildResidue(residue) {
  if (!buildMode || mainStale) return;
  const right = buildLeft && buildLeftRIdx !== residue.idx;
  const pending = right
    ? loadBuildRight(residue.abbr, residue.idx)
    : loadBuildLeft(residue.abbr, residue.idx);
  resChips.querySelectorAll('.res-chip').forEach(chip => {
    const idx = parseInt(chip.dataset.residue);
    if (idx === residue.idx) {
      chip.style.outline = right ? '2px solid #e0a05a' : '2px solid #5a9ae0';
    } else if (!right || idx !== buildLeftRIdx) {
      chip.style.outline = '';
    }
  });
  return pending;
}

function openBuild() {
  buildMode = true;
  buildPanel.classList.add('open');
  btnBuild.classList.add('active');
  btnBuild.setAttribute('aria-expanded', 'true');
  if (!libPanel.classList.contains('open')) openLib();
  clearBuild();
}
function closeBuild() {
  buildMode = false;
  buildPanel.classList.remove('open');
  btnBuild.classList.remove('active');
  btnBuild.setAttribute('aria-expanded', 'false');
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
  buildRightSvg.innerHTML = '<div class="box-placeholder">Choose Use in the library</div>';
  buildRightRg.innerHTML = '';
  buildConnect.disabled = true;
  buildStatus.textContent = 'Choose a residue and a monomer to connect';
  buildStatus.className = 'build-status';
  buildHint.textContent = cabilnInput.value.trim()
    ? 'Select a residue, then choose Use beside a library monomer'
    : 'Choose Use beside a monomer to start a peptide';
  resChips.querySelectorAll('.res-chip').forEach(c => c.style.outline = '');
  btnRxnFilter.disabled = true;
  if (rxnFilterActive) {
    rxnFilterActive = false;
    btnRxnFilter.classList.remove('active');
    btnRxnFilter.setAttribute('aria-pressed', 'false');
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
    buildHint.textContent = 'Select a residue, then choose Use beside a library monomer';
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
    const res = await fetchCalculation('/insert_backbone', {
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
    commitDocument(data.result, 'cabiln');
    buildHint.textContent = `${abbr} inserted — select chips to continue building`;
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
    const res = await fetchCalculation(`/monomer_rgroups?abbr=${encodeURIComponent(abbr)}&residue_idx=${rIdx}&cabiln=${encodeURIComponent(sequence)}`, { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== sequence) return;
    if (data.error) {
      buildLeftSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
      return;
    }
    buildLeftSvg.innerHTML = data.svg || '';
    buildLeft = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
    renderRgroupButtons(buildLeftRg, buildLeft, 'left');
    buildHint.textContent = buildRight ? 'Choose an attachment site on each side' : 'Choose Use beside a library monomer';
    btnRxnFilter.disabled = false;
    if (rxnFilterActive && libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
    updateInsertBetweenUI();
    checkBuildValidity();
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
    const res = await fetchCalculation(url, { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== sequence) return;
    if (data.error) {
      buildRightSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
      return;
    }
    buildRightSvg.innerHTML = data.svg || '';
    buildRight = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
    renderRgroupButtons(buildRightRg, buildRight, 'right');
    buildHint.textContent = buildLeft ? 'Choose an attachment site on each side' : 'Select a residue in the sequence';
    updateInsertBetweenUI();
    checkBuildValidity();
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
    btn.disabled = !!rg.used;
    btn.setAttribute('aria-pressed', String(state.selectedSlot === rg.slot));
    if (state.selectedSlot === rg.slot) btn.classList.add('selected');
    btn.textContent = `R${rg.slot} ${(rg.chem_type || '').replaceAll('_', ' ')}${rg.used ? ' · used' : ''}`;
    btn.title = `R${rg.slot}: ${rg.chem_type || 'unknown'}${rg.leaving ? ' (LG: ' + rg.leaving + ')' : ''}${rg.used ? ' — already connected' : ''}`;
    if (!rg.used) {
      btn.addEventListener('click', () => selectRgroup(side, rg.slot));
    }
    container.appendChild(btn);
  });
  const drawing = side === 'left' ? buildLeftSvg : buildRightSvg;
  drawing.querySelectorAll('.site-selected').forEach(path => path.classList.remove('site-selected'));
  const selected = state.rgroups.find(rg => rg.slot === state.selectedSlot);
  if (Number.isInteger(selected?.atom_idx)) {
    drawing.querySelectorAll(`.atom-${selected.atom_idx}`).forEach(path => path.classList.add('site-selected'));
  }
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
    buildStatus.textContent = !buildLeft ? 'Select a residue in the sequence'
      : !buildRight ? 'Choose a monomer from the library'
      : !buildLeft.selectedSlot ? 'Choose a site on the current residue'
      : 'Choose a site on the other monomer';
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
    const res = await fetchCalculation('/validate_bond', {
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
    const res = await fetchCalculation('/insert_bond', {
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
    commitDocument(data.result, 'cabiln');
    buildHint.textContent = 'Connection added — select a residue to continue building';
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
  const canonical = notationPolicy.value === 'canonical';
  const request = startRequest('sequence-edit');
  try {
    const res = await fetchCalculation('/convert_notation', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ cabiln: val, target, canonical }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== val) return;
    if (data.error) { replaceDocument({ warning: data.error }); return; }
    if (data.result) {
      commitDocument(data.result, 'cabiln', editor.present.warning, {
        canonical: data.canonical || null, context: data.context || editor.present.context,
      });
    }
  } catch (e) {
    if (request.current()) replaceDocument({ warning: 'Notation conversion failed. Try again.' });
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
  setMainProgress(true);
  doRenderCabiln(lastCabiln);
});

// ─── SMILES → CABILN conversion ───────────────────────────────────────────────
function startConversion(button) {
  cancelRequests('sequence-edit');
  const label = button.textContent;
  const request = startRequest('sequence-edit', () => {
    button.textContent = label;
    button.disabled = false;
    conversionProgress.hidden = true;
  });
  button.textContent = '…';
  button.disabled = true;
  conversionProgress.hidden = false;
  conversionProgressLabel.textContent = 'Recognizing monomers and checking the structure. Large peptides can take several seconds.';
  return request;
}

function useConvertedCabiln(sequence, warning = '', result = {}) {
  commitDocument(sequence, 'cabiln', warning, CabilnProject.fromConversion(result));
}

async function doS2c(notation) {
  const smiles = smilesInput.value.trim();
  if (!smiles) return;
  const originalMain = cabilnInput.value.trim();
  const btn = notation === 'bracket' ? btnS2cBracket : btnS2c;
  const request = startConversion(btn);
  try {
    const res = await fetchCalculation('/smiles_to_cabiln', {
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
      useConvertedCabiln(data.cabiln, data.warning, data);
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
  const inputFormat = notationSelect.value;
  if (!txt) return;
  const btn = notation === 'bracket' ? btnToCabilnBracket : btnToCabilnPct;
  const request = startConversion(btn);
  try {
    const res = await fetchCalculation('/to_cabiln', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt, input_format: inputFormat, notation }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== txt ||
        notationSelect.value !== inputFormat) return;
    if (data.error) {
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
    } else {
      useConvertedCabiln(data.cabiln, data.warning, data);
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
  const draft = editor.drafts[notationSelect.value];
  commitDocument(draft?.text || '', notationSelect.value, draft?.warning || '', draft || {});
});

// ─── CABILN render ────────────────────────────────────────────────────────────
cabilnInput.addEventListener('input', () => {
  recordDocument({ ...editor.present, text: cabilnInput.value,
    notation: notationSelect.value, warning: '', quality: null, canonical: null }, true);
  renderDocument();
});

function setMainProgress(pending) {
  renderPane.setAttribute('aria-busy', String(pending));
  renderProgress.hidden = !pending && !(mainStale && hasMainDrawing);
  renderProgressLabel.textContent = pending
    ? (hasMainDrawing ? 'Updating — previous drawing shown' : 'Drawing structure…')
    : 'Previous drawing — correct the input to update';
}

function invalidateDocument() {
  clearTimeout(cabilnTimer);
  cancelRequests('main-render', 'sequence-edit');
  clearComparison();
  clearBuild();
  clearHighlight();
  clearExports();
  lastCabiln = '';
  mainStale = true;
  renderCanvas.classList.add('stale');
  resChips.classList.add('stale');
  resChips.setAttribute('aria-disabled', 'true');
  resChips.querySelectorAll('button').forEach(button => { button.disabled = true; });
  rerollSeed = 0;
  btnReroll.disabled = true;
  btnReroll.textContent = '⟳ Layout';
  cabilnInput.className = '';
  cabilnStatus.textContent = '';
  cabilnStatus.className = 'statusbar';
}

function renderDocument(immediate = false) {
  const seq = cabilnInput.value.trim();
  const mode = notationSelect.value;
  if (!seq) {
    resetCabiln();
    return;
  }
  invalidateDocument();
  if (!hasMainDrawing) showSpinner(renderInner);
  setMainProgress(true);
  const render = () => mode === 'cabiln' ? doRenderCabiln(seq) : doRenderForeign(seq);
  if (immediate) render();
  else cabilnTimer = setTimeout(render, 180);
}

function acceptMainDrawing(svg, source, notation) {
  setInner(renderInner, svg);
  hasMainDrawing = true;
  mainStale = false;
  displayedSource = source;
  displayedNotation = notation;
  renderCanvas.classList.remove('stale');
  resChips.classList.remove('stale');
  resChips.setAttribute('aria-disabled', 'false');
}

function mainRenderError(message) {
  mainStale = true;
  if (!hasMainDrawing) setInner(renderInner, `<div class="placeholder err">${escHtml(message)}</div>`);
  cabilnStatus.textContent = message;
  cabilnStatus.className = 'statusbar';
  cabilnInput.className = 'err';
}

async function doRenderForeign(txt) {
  const request = startRequest('main-render');
  const mode = notationSelect.value;
  lastCabiln = '';
  clearExports();
  clearComparison();
  btnReroll.disabled = true;
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res = await fetchCalculation('/render_reference', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt, input_format: mode, width: w, height: h }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== txt ||
        notationSelect.value !== mode) return;
    if (data.error) {
      mainRenderError(data.error);
    } else {
      acceptDocumentContext(data.context);
      acceptMainDrawing(data.svg, txt, mode);
      resChips.innerHTML = '';
      residueMap = {}; atomToRes = {}; residueList = [];
      cabilnStatus.textContent = `${data.format}: ${data.info || ''}`;
      cabilnStatus.className = 'statusbar ok';
      cabilnInput.className = 'ok';
    }
  } catch (e) {
    if (!request.current()) return;
    mainRenderError('Could not reach the renderer. Your input is preserved.');
  } finally {
    if (request.current()) setMainProgress(false);
    request.finish();
  }
}

function resetCabiln() {
  invalidateDocument();
  hasMainDrawing = false;
  mainStale = false;
  displayedSource = displayedNotation = '';
  renderCanvas.classList.remove('stale');
  resChips.classList.remove('stale');
  resChips.setAttribute('aria-disabled', 'false');
  mainViewport.reset();
  setInner(renderInner, '<div class="placeholder">Start typing a sequence…</div>');
  resChips.innerHTML = '';
  residueMap = {}; atomToRes = {}; residueList = [];
  chainData = [];
  currentBranchSet = new Set();
  setMainProgress(false);
}

async function doRenderCabiln(seq) {
  const request = startRequest('main-render');
  const sameDocument = seq === displayedSource && displayedNotation === 'cabiln';
  lastCabiln = '';
  clearExports();
  clearComparison();
  btnReroll.disabled = true;
  if (buildMode && !sameDocument) clearBuild();
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res  = await fetchCalculation('/render', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: seq, width: w, height: h, seed: rerollSeed }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== seq ||
        notationSelect.value !== 'cabiln') return;
    if (data.error) {
      mainRenderError(data.error);
      clearExports();
      btnReroll.disabled = true;
    } else {
      const displayedSequence = data.normalized_cabiln || seq;
      if (data.normalized_cabiln) {
        projectChanged();
        replaceDocument({ text: displayedSequence, quality: null, canonical: null });
        saveDraft();
      }
      acceptDocumentContext(data.context);
      lastCabiln = displayedSequence;
      acceptMainDrawing(data.svg, displayedSequence, 'cabiln');
      cabilnStatus.textContent = [data.info, ...(data.warnings || [])].filter(Boolean).join(' · ');
      cabilnStatus.title = cabilnStatus.textContent;
      cabilnStatus.className = data.warnings?.length ? 'statusbar warn' : 'statusbar ok';
      cabilnInput.className = 'ok';
      btnReroll.disabled = false;
      setExportReady(data.svg, data.mol_block);
      buildResidueUI(data.residue_map, data.residues, data.layout, data.crosslink_groups);
      if (sameDocument) {
        resChips.querySelectorAll('.res-chip').forEach(chip => {
          const idx = parseInt(chip.dataset.residue);
          if (idx === buildLeftRIdx) chip.style.outline = '2px solid #5a9ae0';
          else if (idx === buildRightRIdx) chip.style.outline = '2px solid #e0a05a';
        });
      }
      if (verifyMode && lastSmiles) triggerVerify();
    }
  } catch (e) {
    if (!request.current()) return;
    mainRenderError('Could not reach the renderer. Your input is preserved.');
  } finally {
    if (request.current()) setMainProgress(false);
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
  projectChanged();
  referenceContext = smilesInput.value ? contextHeader(referenceContext) : null;
  referenceOriginal = smilesInput.value ? { kind: 'text', content: smilesInput.value, name: '' } : null;
  clearReference();
  if (smilesInput.value) beginDraftEdit();
  scheduleDraftSave();
  const txt = smilesInput.value.trim();
  if (!txt) {
    setInner(smilesInner, '<div class="placeholder">Paste SMILES, BILN, or HELM…</div>');
    compareBar.innerHTML = '';
    return;
  }
  smilesStatus.textContent = 'Updating reference…';
  smilesTimer = setTimeout(() => doRenderRef(txt), 180);
});

molUpload.addEventListener('change', async (e) => {
  const file = e.target.files[0];
  if (!file) return;
  beginDraftEdit();
  projectChanged();
  referenceContext = contextHeader(referenceContext);
  referenceOriginal = null;
  smilesInput.value = '';
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
    await renderMolReference(text, file.name, request);
  } catch (err) {
    if (!request.current()) return;
    smilesStatus.textContent = 'Could not read the MOL/SDF file. The previous document is unchanged.';
    smilesStatus.className = 'statusbar';
  } finally { request.finish(); }
});

async function renderMolReference(text, name, request = startRequest('reference-render')) {
  if (referenceOriginal?.kind !== 'mol' || referenceOriginal.content !== text ||
      referenceOriginal.name !== name) {
    projectChanged(false);
    referenceContext = contextHeader(referenceContext);
  }
  referenceOriginal = { kind: 'mol', content: text, name };
  saveDraft();
  showSpinner(smilesInner);
  smilesStatus.textContent = `Reference: ${name || 'uploaded structure'}`;
  try {
    const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
    const res = await fetchCalculation('/render_mol', {
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
      if (smilesInput.value !== lastSmiles) projectChanged(false);
      smilesInput.value = lastSmiles;
      acceptReferenceContext(data.context);
      saveDraft();
      smilesInput.className = 'ok';
      if (lastCabiln) triggerVerify();
    }
  } catch (err) {
    if (!request.current()) return;
    smilesStatus.textContent = 'Failed to render the original MOL/SDF reference';
    smilesStatus.className = 'statusbar';
  } finally { request.finish(); }
}

async function doRenderRef(txt) {
  const request = startRequest('reference-render');
  lastSmiles = '';
  clearComparison();
  const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
  try {
    const res = await fetchCalculation('/render_reference', {
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
      acceptReferenceContext(data.context);
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
    const res  = await fetchCalculation('/verify', {
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

offerSavedDraft();
updateHistoryControls();
