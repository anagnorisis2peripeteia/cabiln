// state
let darkMode   = true;
let verifyMode = false;
let hlEnabled  = true;
let cabilnTimer = null;
let smilesTimer = null;
let lastCabiln  = '';
let lastSmiles  = '';
let lastSvg     = '';
let lastMolBlock = '';
let buildMode    = false;
let buildLeft    = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
let buildRight   = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
let buildLeftRIdx = null;
let buildRightRIdx = null;
let buildReaction = '';
let insertBetweenActive = false;
let swapState = null;  // Current source, eligible definitions, reviewed mapping and preview.
let rerollSeed   = 0;
let mainStale = false;
let hasMainDrawing = false;
let displayedSource = '';
let displayedNotation = '';
let cabilnDrawing = null;
let projectContext = null;
let projectRevision = 0;
let draftCleared = false;
let referenceOriginal = null;
let referenceContext = null;

function setStatus(element, message, kind = '') {
  element.textContent = message;
  element.title = message;
  element.className = kind ? `statusbar ${kind}` : 'statusbar';
}

function invalidInputResponse(response) {
  return response.status === 400 || response.status === 422;
}

// elements
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
const btnHl         = document.getElementById('btn-hl');
const btnPng        = document.getElementById('btn-png');
const btnMol        = document.getElementById('btn-mol');
const btnToBracket  = document.getElementById('btn-to-bracket');
const btnToBranch   = document.getElementById('btn-to-branch');
const notationPolicy = document.getElementById('notation-policy');
const btnExamples    = document.getElementById('btn-examples');
const examplesPanel  = document.getElementById('examples-panel');
const examplesClose  = document.getElementById('examples-close');
const examplesList   = document.getElementById('examples-list');
const resChips      = document.getElementById('residue-chips');
const buildPanel    = document.getElementById('build-panel');
const btnBuild      = document.getElementById('btn-build');
const buildClose    = document.getElementById('build-close');
const buildConnect  = document.getElementById('build-connect');
const buildAction = document.getElementById('build-action');
const swapMapping = document.getElementById('build-swap-mapping');
const swapSites = document.getElementById('build-swap-sites');
const buildPreviewButton = document.getElementById('build-preview-button');
const buildPreviewPanel = document.getElementById('build-preview');
const buildPreviewStatus = document.getElementById('build-preview-status');
const buildPreviewReaction = document.getElementById('build-preview-reaction');
const buildPreviewInner = document.getElementById('build-preview-inner');
const buildPreviewSource = document.getElementById('build-preview-source');
const buildPreviewNotation = document.getElementById('build-preview-notation');
const buildStatus   = document.getElementById('build-status');
const buildHint     = document.getElementById('build-hint');
const buildLeftAbbr = document.getElementById('build-left-abbr');
const buildLeftSvg  = document.getElementById('build-left-svg');
const buildLeftRg   = document.getElementById('build-left-rgroups');
const buildLeftSite = document.getElementById('build-left-site');
const buildRightAbbr= document.getElementById('build-right-abbr');
const buildRightSvg = document.getElementById('build-right-svg');
const buildRightRg  = document.getElementById('build-right-rgroups');
const buildRightSite = document.getElementById('build-right-site');
const buildInsertRow = document.getElementById('build-insert-row');
const buildInsertInfo = document.getElementById('build-insert-info');
const buildInsertBtn = document.getElementById('build-insert-btn');
const btnReroll     = document.getElementById('btn-reroll');
const btnS2c        = document.getElementById('btn-s2c');
const btnS2cBracket = document.getElementById('btn-s2c-bracket');
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

const residueView = createResidueView({
  inner: renderInner, chips: resChips,
  canHighlight: () => hlEnabled && !mainStale,
  onSelect: residue => selectBuildResidue(residue),
});

const library = new MonomerLibrary({
  getFilters: () => ({
    building: buildMode, left: buildLeft, insertBetween: insertBetweenActive,
    replacements: buildMode && isSwapMode() ? (swapState?.candidates || []) : null,
    awaitingSelection: !swapState,
  }),
  onUse: (abbr, explicit) => useLibraryMonomer(abbr, explicit),
  onChanged: () => { if (buildMode) clearBuild(); },
});

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
    // Drafts without a library binding require verification.
    applyProject(project);
    setProjectStatus('Older draft restored without a library binding. Check any custom monomer definitions before using it.', true);
  }
});
btnDismissDraft.addEventListener('click', () => {
  requests.cancel('project-open');
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
  if (requests.has('project-open') || requests.has('project-save')) {
    setProjectStatus('The document changed during the project check. Save or open again when ready.', true);
  }
  requests.cancel('project-open', 'project-save');
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
  const request = requests.start('project-open');
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
  const request = requests.start('project-save', () => { btnProjectSave.disabled = false; });
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
  const request = requests.start('project-open');
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
  const request = requests.start('project-open');
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
  requests.cancel('project-open');
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

// dark mode
btnDark.addEventListener('click', () => {
  darkMode = !darkMode;
  btnDark.classList.toggle('active', darkMode);
  btnDark.setAttribute('aria-pressed', String(darkMode));
  for (const el of document.querySelectorAll('.canvas-wrap')) {
    el.classList.toggle('dark', darkMode);
  }
  library.setDark(darkMode);
  document.querySelectorAll('.build-box').forEach(el => el.classList.toggle('dark', darkMode));
});
btnDark.classList.add('active');
document.querySelectorAll('.canvas-wrap').forEach(el => el.classList.add('dark'));
library.setDark(true);
document.querySelectorAll('.build-box').forEach(el => el.classList.add('dark'));

// highlight toggle
btnHl.addEventListener('click', () => {
  hlEnabled = !hlEnabled;
  btnHl.classList.toggle('active', hlEnabled);
  btnHl.setAttribute('aria-pressed', String(hlEnabled));
  if (!hlEnabled) residueView.clearHighlight();
});

// verify mode
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

// example peptide sidebar
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
  } else if (buildMode && isSwapMode()) {
    chooseSwapMonomer(abbr);
  } else if (buildMode) {
    loadBuildRight(abbr);
    if (!buildLeft) buildHint.textContent = 'Now select a residue in the sequence';
  } else {
    insertAbbr(abbr);
  }
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

// zoom / pan (shared, wired per canvas)
const mainViewport = makeZoomable(renderCanvas, renderInner);
makeZoomable(document.getElementById('smiles-canvas'), smilesInner);
const buildPreviewViewport = makeZoomable(document.getElementById('build-preview-canvas'), buildPreviewInner);

// PNG download (client-side SVG → canvas → PNG)
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

// MOL download (server-side mol block)
btnMol.addEventListener('click', async () => {
  if (!lastMolBlock) return;
  const blob = new Blob([lastMolBlock], { type: 'chemical/x-mdl-molfile' });
  const a = document.createElement('a');
  a.download = 'structure.mol';
  a.href = URL.createObjectURL(blob);
  a.click();
});

// build mode
function selectBuildResidue(residue) {
  if (!buildMode || mainStale) return;
  const right = !isSwapMode() && buildLeft && buildLeftRIdx !== residue.idx;
  const pending = right
    ? loadBuildRight(residue.abbr, residue.idx)
    : loadBuildLeft(residue.abbr, residue.idx);
  residueView.select(right ? buildLeftRIdx : residue.idx, right ? residue.idx : null);
  return pending;
}

function openBuild() {
  buildMode = true;
  buildPanel.classList.add('open');
  btnBuild.classList.add('active');
  btnBuild.setAttribute('aria-expanded', 'true');
  if (!library.isOpen) library.open();
  clearBuild();
}
function closeBuild() {
  buildMode = false;
  buildPanel.classList.remove('open');
  btnBuild.classList.remove('active');
  btnBuild.setAttribute('aria-expanded', 'false');
  clearBuild();
  if (library.loaded && isSwapMode()) library.render();
}
function clearBuild() {
  clearBuildPreview();
  requests.cancel('build-left', 'build-right', 'bond-check', 'sequence-edit', 'swap-options');
  swapState = null;
  swapMapping.hidden = true;
  swapSites.innerHTML = '';
  buildLeft = null; buildRight = null; buildLeftRIdx = null; buildRightRIdx = null;
  buildReaction = '';
  insertBetweenActive = false;
  buildInsertRow.style.display = 'none';
  buildInsertBtn.textContent = '⊕ Insert Between';
  buildLeftAbbr.textContent = '—';
  buildLeftSvg.innerHTML = '<div class="box-placeholder">Click a chip above</div>';
  buildLeftRg.innerHTML = '';
  buildLeftSite.hidden = true;
  buildRightAbbr.textContent = '—';
  document.getElementById('build-right-label').textContent = isSwapMode() ? 'Replacement monomer' : 'New monomer';
  document.getElementById('build-title').textContent = isSwapMode() ? 'Swap Monomer' : 'Build by Connection';
  buildPanel.classList.toggle('swapping', isSwapMode());
  buildConnect.textContent = isSwapMode() ? 'Apply swap' : 'Connect';
  buildRightSvg.innerHTML = '<div class="box-placeholder">Choose Use in the library</div>';
  buildRightRg.innerHTML = '';
  buildRightSite.hidden = true;
  setBuildReady(false);
  buildStatus.textContent = 'Choose a residue and a monomer to connect';
  buildStatus.className = 'build-status';
  buildHint.textContent = isSwapMode() ? 'Select the residue to replace; every existing connection will be kept'
    : cabilnInput.value.trim()
    ? 'Select a residue, then choose Use beside a library monomer'
    : 'Choose Use beside a monomer to start a peptide';
  residueView.select();
  const filtered = library.resetFilter();
  if (library.loaded && (filtered || (buildMode && isSwapMode()))) library.render();
}

function isSwapMode() { return buildAction.value === 'swap'; }

buildAction.addEventListener('change', () => {
  const selected = residueView.residue(buildLeftRIdx);
  clearBuild();
  if (library.loaded) library.render();
  if (selected && !mainStale) selectBuildResidue(selected);
});

async function loadSwapOptions() {
  const request = requests.start('swap-options');
  const source = cabilnInput.value.trim();
  const index = buildLeftRIdx;
  buildStatus.textContent = 'Finding replacements for all connected sites…';
  try {
    while (requests.has('library') || requests.has('reactions')) {
      const pending = requests.pending('library', 'reactions');
      await Promise.race([pending, request.done]);
      if (!request.current()) return;
    }
    const data = await readResponse(await fetchCalculation('/replacement_options', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: source, residue_idx: index }), signal: request.signal,
    }));
    if (!request.current() || source !== cabilnInput.value.trim() || index !== buildLeftRIdx) return;
    if (data.error) throw new Error(data.error);
    if (data.source_echo !== source || data.residue_idx !== index || !data.context) {
      throw new Error('The replacement list is out of date. Select the residue again.');
    }
    const quality = new Map(library.monomers.map(m => [m.abbr, m.quality]));
    const candidates = data.candidates.map(m => ({ ...m, quality: quality.get(m.abbr),
      searchText: [m.abbr, m.name, m.type, m.chem_types].join(' ').toLowerCase() }));
    swapState = { ...data, candidates, candidate: null, mapping: {}, preview: null };
    buildStatus.textContent = 'Choose a replacement from the filtered library';
    buildHint.textContent = 'Only monomers with compatible sites for every connection are shown';
    library.render();
  } catch (error) {
    if (request.current()) {
      buildStatus.textContent = error.message || 'Could not find replacements. Select the residue to retry.';
      buildStatus.className = 'build-status invalid';
    }
  } finally { request.finish(); }
}

function chooseSwapMonomer(abbr) {
  const candidate = swapState?.candidates.find(item => item.abbr === abbr);
  if (!candidate) return;
  swapState.candidate = candidate;
  swapState.mapping = { ...candidate.mapping };
  swapState.preview = null;
  renderSwapMapping();
  loadBuildRight(abbr);
}

function renderSwapMapping() {
  swapSites.innerHTML = '';
  swapMapping.hidden = !swapState?.candidate;
  if (swapMapping.hidden) return;
  const candidate = swapState.candidate;
  for (const need of swapState.requirements) {
    const label = document.createElement('label');
    label.className = 'swap-site';
    const text = document.createElement('span');
    text.textContent = `Current R${need.slot} →`;
    const select = document.createElement('select');
    select.className = 'hbtn';
    select.dataset.slot = need.slot;
    select.setAttribute('aria-label', `Replacement site for R${need.slot}`);
    for (const slot of candidate.choices[need.slot]) {
      const option = document.createElement('option');
      option.value = slot;
      const site = candidate.sites.find(item => item.slot === slot);
      option.textContent = `R${slot} ${(site?.chem_type || '').replaceAll('_', ' ')}`;
      select.appendChild(option);
    }
    select.value = swapState.mapping[need.slot];
    select.addEventListener('change', () => {
      requests.cancel('sequence-edit');
      swapState.mapping[need.slot] = Number(select.value);
      checkBuildValidity();
      if (buildRight) renderRgroupButtons(buildRightRg, buildRight, 'right');
    });
    const partner = document.createElement('small');
    partner.textContent = need.internal ? `joins the replacement site mapped from R${need.partner_slot}`
      : `keeps ${need.partner_abbr} (residue ${need.partner_idx + 1}) R${need.partner_slot}`;
    label.appendChild(text);
    label.appendChild(select);
    label.appendChild(partner);
    swapSites.appendChild(label);
  }
  if (!swapState.requirements.length) swapSites.textContent = 'This monomer has no existing connections.';
}

function swapRequest() {
  if (!swapState?.candidate || !buildLeft || !buildRight || mainStale ||
      swapState.source_echo !== cabilnInput.value.trim() || swapState.residue_idx !== buildLeftRIdx ||
      swapState.candidate.abbr !== buildRight.abbr) return null;
  const slots = Object.values(swapState.mapping);
  if (slots.length !== swapState.requirements.length || new Set(slots).size !== slots.length) return null;
  return { cabiln: swapState.source_echo, residue_idx: swapState.residue_idx,
    new_abbr: swapState.candidate.abbr, slot_map: { ...swapState.mapping }, context: swapState.context };
}

btnBuild.addEventListener('click', () => buildMode ? closeBuild() : openBuild());
buildClose.addEventListener('click', closeBuild);

function checkAdjacentBackbone() {
  return residueView.insertionAnchor(buildLeftRIdx, buildRightRIdx) !== null;
}

function updateInsertBetweenUI() {
  if (buildMode && !isSwapMode() && checkAdjacentBackbone()) {
    const la = (residueView.residue(buildLeftRIdx) || {}).abbr || '?';
    const ra = (residueView.residue(buildRightRIdx) || {}).abbr || '?';
    buildInsertInfo.textContent = `${la} and ${ra} are adjacent on the backbone`;
    buildInsertRow.style.display = 'flex';
  } else {
    buildInsertRow.style.display = 'none';
    if (insertBetweenActive) {
      insertBetweenActive = false;
      buildInsertBtn.textContent = '⊕ Insert Between';
      if (library.loaded) library.render();
    }
  }
}

buildInsertBtn.addEventListener('click', () => {
  clearBuildPreview();
  requests.cancel('sequence-edit');
  if (insertBetweenActive) {
    insertBetweenActive = false;
    buildInsertBtn.textContent = '⊕ Insert Between';
    buildHint.textContent = 'Select a residue, then choose Use beside a library monomer';
    if (library.loaded) library.render();
    return;
  }
  if (!checkAdjacentBackbone()) return;
  insertBetweenActive = true;
  buildInsertBtn.textContent = '✕ Cancel';
  buildHint.textContent = 'Click a backbone monomer in the library to insert between the selected residues';
  if (!library.isOpen) library.open();
  if (library.loaded) library.render();
});

async function doInsertBetween(abbr) {
  clearBuildPreview();
  const after_idx = residueView.insertionAnchor(buildLeftRIdx, buildRightRIdx);
  if (after_idx === null) return;
  const val = cabilnInput.value.trim();
  const request = requests.start('sequence-edit');
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
    if (library.loaded) library.render();
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
  clearBuildPreview();
  requests.cancel('bond-check', 'sequence-edit', 'swap-options');
  if (isSwapMode()) {
    requests.cancel('build-right');
    swapState = null;
    swapMapping.hidden = true;
    buildRight = null;
    buildRightRIdx = null;
    buildRightAbbr.textContent = '—';
    buildRightSvg.innerHTML = '<div class="box-placeholder">Choose a replacement in the library</div>';
    buildRightRg.innerHTML = '';
    buildRightSite.hidden = true;
    library.render();
  }
  const request = requests.start('build-left');
  const sequence = cabilnInput.value.trim();
  buildLeftRIdx = rIdx;
  buildLeftAbbr.textContent = abbr;
  buildLeftSvg.innerHTML = '<div class="spinner"></div>';
  buildLeftRg.innerHTML = '';
  buildLeftSite.hidden = true;
  buildLeft = null;
  setBuildReady(false);
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
    library.setFilterEnabled(!isSwapMode());
    if (library.filterActive && library.loaded) library.render();
    updateInsertBetweenUI();
    if (isSwapMode()) loadSwapOptions();
    else checkBuildValidity();
  } catch (e) {
    if (!request.current()) return;
    buildLeftSvg.innerHTML = '<div class="box-placeholder">Error loading monomer</div>';
  } finally {
    request.finish();
  }
}

async function loadBuildRight(abbr, rIdx) {
  clearBuildPreview();
  requests.cancel('bond-check', 'sequence-edit');
  const request = requests.start('build-right');
  const sequence = cabilnInput.value.trim();
  buildRightRIdx = rIdx !== undefined ? rIdx : null;
  buildRightAbbr.textContent = abbr;
  document.getElementById('build-right-label').textContent = isSwapMode() ? 'Replacement monomer'
    : rIdx !== undefined ? 'Current residue' : 'New monomer';
  buildRightSvg.innerHTML = '<div class="spinner"></div>';
  buildRightRg.innerHTML = '';
  buildRightSite.hidden = true;
  buildRight = null;
  setBuildReady(false);
  buildStatus.textContent = '';

  const family = !isSwapMode() && rIdx === undefined && library.monomers.find(m => m.abbr === abbr && m.degenerate);
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
    buildHint.textContent = isSwapMode() ? 'Review each site mapping, then preview the replacement'
      : buildLeft ? 'Choose an attachment site on each side' : 'Select a residue in the sequence';
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
    btn.disabled = isSwapMode() || !!rg.used;
    btn.setAttribute('aria-pressed', String(state.selectedSlot === rg.slot));
    if (state.selectedSlot === rg.slot) btn.classList.add('selected');
    btn.textContent = `R${rg.slot} ${(rg.chem_type || '').replaceAll('_', ' ')}${rg.used ? (isSwapMode() ? ' · connected' : ' · used') : ''}`;
    btn.title = `R${rg.slot}: ${rg.chem_type || 'unknown'} · Free-site group: ${rg.leaving || '[H] (implicit)'}${rg.used ? ' — already connected' : ''}`;
    if (!isSwapMode() && !rg.used) {
      btn.addEventListener('click', () => selectRgroup(side, rg.slot));
    }
    container.appendChild(btn);
  });
  const drawing = side === 'left' ? buildLeftSvg : buildRightSvg;
  drawing.querySelectorAll('.site-selected').forEach(path => path.classList.remove('site-selected'));
  const selected = state.rgroups.find(rg => rg.slot === state.selectedSlot);
  const detail = side === 'left' ? buildLeftSite : buildRightSite;
  detail.hidden = !selected;
  if (selected) {
    detail.textContent = `R${selected.slot} · Free-site group: ${selected.leaving || '[H] (implicit)'}`;
    detail.title = 'Group present at this site when it is unconnected. The product preview shows the selected reaction.';
  }
  if (Number.isInteger(selected?.atom_idx)) {
    drawing.querySelectorAll(`.atom-${selected.atom_idx}`).forEach(path => path.classList.add('site-selected'));
  }
  if (isSwapMode()) {
    const occupied = side === 'left' ? state.rgroups.filter(rg => rg.used)
      : state.rgroups.filter(rg => Object.values(swapState?.mapping || {}).includes(rg.slot));
    for (const site of occupied) {
      drawing.querySelectorAll(`.atom-${site.atom_idx}`).forEach(path => path.classList.add('site-selected'));
    }
  }
}

function selectRgroup(side, slot) {
  requests.cancel('sequence-edit');
  if (side === 'left' && buildLeft) {
    buildLeft.selectedSlot = buildLeft.selectedSlot === slot ? null : slot;
    renderRgroupButtons(buildLeftRg, buildLeft, 'left');
    if (library.filterActive && library.loaded) library.render();
  } else if (side === 'right' && buildRight) {
    buildRight.selectedSlot = buildRight.selectedSlot === slot ? null : slot;
    renderRgroupButtons(buildRightRg, buildRight, 'right');
  }
  checkBuildValidity();
}

async function checkBuildValidity() {
  clearBuildPreview();
  requests.cancel('bond-check');
  buildReaction = '';
  setBuildReady(false);
  if (isSwapMode()) {
    const valid = !!swapRequest();
    buildStatus.textContent = valid ? 'Preview to validate the complete replacement product'
      : swapState?.candidate && buildRight ? 'Choose a different replacement site for each connection'
      : 'Choose a replacement from the filtered library';
    buildStatus.className = 'build-status';
    setBuildReady(valid);
    return;
  }
  if (!buildLeft?.selectedSlot || !buildRight?.selectedSlot) {
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

  const request = requests.start('bond-check');

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
      buildReaction = data.reaction ? `Reaction: ${data.reaction.replaceAll('_', ' ')}` : 'Bond formation';
      buildStatus.textContent = `Valid: ${data.reaction || 'bond'} (R${buildLeft.selectedSlot}↔R${buildRight.selectedSlot})`;
      buildStatus.className = 'build-status valid';
      setBuildReady(true);
    } else {
      buildStatus.textContent = data.error || data.reason || 'No compatible reaction found';
      buildStatus.className = 'build-status invalid';
      setBuildReady(false);
    }
  } catch (e) {
    if (!request.current()) return;
    buildStatus.textContent = 'Validation error';
    buildStatus.className = 'build-status invalid';
    setBuildReady(false);
  } finally {
    request.finish();
  }
}

function setBuildReady(ready) {
  buildConnect.disabled = !ready || (isSwapMode() && !swapState?.preview);
  buildPreviewButton.disabled = !ready;
}

function buildConnection() {
  if (!buildLeft?.selectedSlot || !buildRight?.selectedSlot) return null;
  return {
    cabiln: cabilnInput.value.trim(),
    host_residue_idx: buildLeftRIdx ?? 0,
    new_abbr: buildRight.abbr,
    r_host: buildLeft.selectedSlot,
    r_new: buildRight.selectedSlot,
    target_residue_idx: (buildRightRIdx !== null && buildRightRIdx !== buildLeftRIdx)
      ? buildRightRIdx : -1,
  };
}

function clearBuildPreview() {
  if (swapState) swapState.preview = null;
  if (isSwapMode()) buildConnect.disabled = true;
  const pending = requests.has('build-preview');
  requests.cancel('build-preview');
  if (pending) setBuildReady(true);
  buildPreviewPanel.hidden = true;
  buildPreviewInner.innerHTML = '';
  buildPreviewStatus.textContent = '';
  buildPreviewReaction.textContent = '';
  buildPreviewSource.textContent = '';
  buildPreviewNotation.hidden = true;
  buildPreviewButton.setAttribute('aria-expanded', 'false');
}

document.getElementById('build-preview-close').addEventListener('click', () => {
  clearBuildPreview();
  buildPreviewButton.focus();
});
buildPreviewButton.addEventListener('click', async () => {
  const swapping = isSwapMode();
  const connection = swapping ? swapRequest() : buildConnection();
  if (!connection || buildPreviewButton.disabled) return;
  clearBuildPreview();
  const request = requests.start('build-preview', () => buildPreviewPanel.setAttribute('aria-busy', 'false'));
  const left = swapping ? `${buildLeft.abbr} (residue ${connection.residue_idx + 1})`
    : `${buildLeft.abbr} (residue ${connection.host_residue_idx + 1}) R${connection.r_host}`;
  const right = swapping ? buildRight.abbr
    : `${buildRight.abbr} (${connection.target_residue_idx < 0 ? 'new' : `residue ${connection.target_residue_idx + 1}`}) R${connection.r_new}`;
  const pair = `${left} ↔ ${right}`;
  setBuildReady(false);
  buildPreviewPanel.hidden = false;
  buildPreviewPanel.setAttribute('aria-busy', 'true');
  buildPreviewButton.setAttribute('aria-expanded', 'true');
  buildPreviewStatus.textContent = `${pair} · Preparing preview…`;
  buildPreviewStatus.className = '';
  buildPreviewReaction.textContent = swapping ? Object.entries(connection.slot_map)
    .map(([old, next]) => `R${old} → R${next}`).join(' · ') || 'No existing connections' : buildReaction;
  buildPreviewInner.innerHTML = '<div class="spinner"></div>';
  buildPreviewViewport.reset();
  buildPreviewPanel.scrollIntoView({ block: 'nearest' });
  try {
    // Drawing and verification share one worker. Do not spend preview retries
    // competing with work that is already running for this page.
    while (request.current()) {
      const pending = requests.pending('main-render', 'reference-render', 'verify');
      if (!pending) break;
      buildPreviewStatus.textContent = `${pair} · Waiting for the current drawing or verification…`;
      await Promise.race([pending, request.done]);
    }
    if (!request.current()) return;
    buildPreviewStatus.textContent = `${pair} · Preparing preview…`;
    // These endpoints calculate a candidate; only Connect commits the document.
    const proposal = await readResponse(await fetchCalculation(swapping ? '/replace_monomer' : '/insert_bond', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(connection), signal: request.signal,
    }));
    if (!request.current()) return;
    if (proposal.error) throw new Error(proposal.error);
    const drawing = await readResponse(await fetchCalculation('/render', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: proposal.result, width: 1000, height: 320 }),
      signal: request.signal,
    }));
    if (!request.current()) return;
    if (drawing.error) throw new Error(drawing.error);
    if (swapping && (!CabilnProject.sameContext(connection.context, proposal.context) ||
        !CabilnProject.sameContext(connection.context, drawing.context))) {
      throw new Error('The library changed. Select the residue again to refresh replacements');
    }
    const changedLibrary = editor.present.context?.library_binding && drawing.context?.library_binding
      && !CabilnProject.sameContext(editor.present.context, drawing.context);
    const notice = changedLibrary ? 'Library changed since the current drawing; preview uses current definitions.' : '';
    buildPreviewInner.innerHTML = drawing.svg;
    buildPreviewStatus.textContent = [pair, drawing.info, notice, drawing.normalization_note, ...(drawing.warnings || [])].filter(Boolean).join(' · ');
    buildPreviewSource.textContent = proposal.result;
    buildPreviewNotation.hidden = false;
    if (swapping) {
      swapState.preview = { result: proposal.result, request: connection };
      buildStatus.textContent = 'Product validated. Apply swap keeps every mapped connection.';
      buildStatus.className = 'build-status valid';
      if (Object.entries(connection.slot_map).some(([old, next]) => Number(old) !== next)) {
        buildPreviewStatus.textContent += ' · Site renumbering can reformat the notation; review Proposed CABILN.';
      }
    }
  } catch (error) {
    if (!request.current()) return;
    buildPreviewInner.innerHTML = '';
    buildPreviewStatus.textContent = `Preview unavailable: ${error.message}. Your sequence is unchanged.`;
    buildPreviewStatus.className = 'error';
  } finally {
    if (request.current()) setBuildReady(true);
    request.finish();
  }
});

buildConnect.addEventListener('click', async () => {
  const swapping = isSwapMode();
  const reviewed = swapping ? swapState?.preview : null;
  const connection = swapping ? reviewed?.request : buildConnection();
  if (!connection) return;
  clearBuildPreview();
  const request = requests.start('sequence-edit');

  setBuildReady(false);
  buildStatus.textContent = swapping ? 'Checking and applying swap…' : 'Inserting...';
  buildStatus.className = 'build-status';

  try {
    const res = await fetchCalculation(swapping ? '/replace_monomer' : '/insert_bond', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(connection),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== connection.cabiln) return;
    if (data.error) {
      buildStatus.textContent = data.error;
      buildStatus.className = 'build-status invalid';
      return;
    }
    if (swapping && (data.result !== reviewed.result ||
        !CabilnProject.sameContext(connection.context, data.context))) {
      buildStatus.textContent = 'The replacement changed. Preview it again before applying.';
      buildStatus.className = 'build-status invalid';
      return;
    }
    commitDocument(data.result, 'cabiln');
    buildHint.textContent = swapping ? 'Monomer replaced — Undo restores the original peptide'
      : 'Connection added — select a residue to continue building';
  } catch (e) {
    if (!request.current()) return;
    buildStatus.textContent = swapping ? 'Swap failed. Your sequence is unchanged; preview to retry.' : 'Insert failed';
    buildStatus.className = 'build-status invalid';
  } finally {
    if (swapping && request.current()) setBuildReady(!!swapRequest());
    request.finish();
  }
});

// notation conversion
async function settleDrawingsForConversion(request) {
  // Aborting a running drawing retires the server's chemistry worker. Let the
  // current drawing (and its automatic Verify) finish before competing with it.
  while (request.current()) {
    const pending = requests.pending('main-render', 'reference-render', 'verify');
    if (pending) {
      await Promise.race([pending, request.done]);
      continue;
    }
    if (cabilnTimer !== null) {
      clearTimeout(cabilnTimer);
      cabilnTimer = null;
      const text = cabilnInput.value.trim();
      if (text) {
        if (notationSelect.value === 'cabiln') doRenderCabiln(text);
        else doRenderForeign(text);
      }
      continue;
    }
    if (smilesTimer !== null) {
      clearTimeout(smilesTimer);
      smilesTimer = null;
      const text = smilesInput.value.trim();
      if (text) doRenderRef(text);
      continue;
    }
    return;
  }
}

function drawingForNotation(snapshot, data) {
  const view = canvasSize(renderCanvas);
  const presentation = data.presentation;
  const order = data.occurrence_order;
  if (!snapshot || snapshot !== cabilnDrawing || mainStale ||
      requests.has('main-render') || displayedNotation !== 'cabiln' ||
      displayedSource !== snapshot.source || data.source_echo !== snapshot.source ||
      !lastSvg || !lastMolBlock || view.w !== snapshot.canvas.w || view.h !== snapshot.canvas.h ||
      !CabilnProject.sameContext(snapshot.data.context, data.context) ||
      !CabilnProject.sameContext(editor.present.context, data.context) ||
      presentation?.cabiln_echo !== data.result || !presentation?.layout ||
      !Array.isArray(presentation.residues) || !Array.isArray(order)) return null;
  const count = snapshot.data.residues.length;
  if (order.length !== count || presentation.residues.length !== count ||
      new Set(order).size !== count || !order.every(index =>
        Number.isInteger(index) && index >= 0 && index < count &&
        Array.isArray(snapshot.data.residue_map[index]))) return null;
  // Keep the SVG and MOL atom order together. Only occurrence IDs change;
  // molecular equality alone cannot identify the atom owners of a depiction.
  return { ...presentation, context: data.context,
    svg: snapshot.data.svg, mol_block: snapshot.data.mol_block, info: snapshot.data.info,
    residue_map: Object.fromEntries(order.map((previous, index) =>
      [index, snapshot.data.residue_map[previous]])),
  };
}

async function convertNotation(target) {
  const original = cabilnInput.value.trim();
  if (!original) return;
  const canonical = notationPolicy.value === 'canonical';
  let deferredDrawing = false;
  const request = requests.start('sequence-edit', () => {
    // Let the cancelling action schedule its own drawing or formatter first.
    Promise.resolve().then(() => {
      if (deferredDrawing && mainStale && cabilnTimer === null &&
          !requests.has('main-render') && !requests.has('sequence-edit') &&
          notationSelect.value === 'cabiln' && cabilnInput.value.trim() === original) renderDocument();
    });
  });
  // The formatter validates and normalizes unsent input. Draw its result once,
  // or draw the original input if formatting fails. Stale input also covers a
  // second formatter click taking over the first click's deferred drawing.
  deferredDrawing = cabilnTimer !== null || (mainStale && !requests.has('main-render'));
  clearTimeout(cabilnTimer);
  cabilnTimer = null;
  let conversionError = '';
  try {
    await settleDrawingsForConversion(request);
    if (!request.current()) return;
    const val = cabilnInput.value.trim();
    const snapshot = cabilnDrawing;
    const res = await fetchCalculation('/convert_notation', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ cabiln: val, target, canonical }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== val) return;
    if (data.error) { conversionError = data.error; return; }
    if (data.result) {
      deferredDrawing = false;
      const drawing = drawingForNotation(snapshot, data);
      const document = { text: data.result, notation: 'cabiln', warning: editor.present.warning,
        quality: null, canonical: data.canonical || null, context: data.context || editor.present.context };
      recordDocument(document);
      if (drawing) {
        clearComparison();
        clearBuild();
        residueView.clearHighlight();
        acceptCabilnDrawing(drawing, data.result, snapshot);
        setMainProgress(false);
      } else renderDocument(true);
    }
  } catch (e) {
    conversionError = 'Notation conversion failed. Try again.';
  } finally {
    if (request.current() && deferredDrawing) {
      deferredDrawing = false;
      await doRenderCabiln(cabilnInput.value.trim());
    }
    if (request.current() && conversionError) setStatus(cabilnStatus, conversionError);
    request.finish();
  }
}
btnToBracket.addEventListener('click', () => convertNotation('bracket'));
btnToBranch.addEventListener('click', () => convertNotation('branch'));

btnReroll.addEventListener('click', () => {
  if (!lastCabiln) return;
  // The default already uses Indigo; start with the alternate CoordGen layout.
  rerollSeed = rerollSeed === 0 ? 2 : rerollSeed + 1;
  btnReroll.textContent = rerollSeed % 2 === 1 ? '⟳ Indigo' : '⟳ CoordGen';
  setMainProgress(true);
  doRenderCabiln(lastCabiln);
});

// SMILES → CABILN conversion
function startConversion(button) {
  requests.cancel('sequence-edit');
  const label = button.textContent;
  const request = requests.start('sequence-edit', () => {
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
  if (!smilesInput.value.trim() && !requests.has('reference-render')) return;
  const btn = notation === 'bracket' ? btnS2cBracket : btnS2c;
  const request = startConversion(btn);
  setStatus(smilesStatus, '');
  try {
    await settleDrawingsForConversion(request);
    if (!request.current()) return;
    const smiles = smilesInput.value.trim();
    if (!smiles) return;
    const originalMain = cabilnInput.value.trim();
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
      setStatus(smilesStatus, 'S2C: ' + data.error);
    } else {
      useConvertedCabiln(data.cabiln, data.warning, data);
      if (data.warning) {
        setStatus(smilesStatus, `⚠ ${data.warning}`, 'warn');
      } else {
        const count = data.assignments?.length ?? data.details.length;
        setStatus(smilesStatus, `Converted (${notation}): ${count} monomer(s)`, 'ok');
      }
    }
  } catch (e) {
    if (!request.current()) return;
    setStatus(smilesStatus, 'Could not convert the reference. Try again.');
  } finally {
    request.finish();
  }
}
btnS2c.addEventListener('click', () => doS2c('percent'));
btnS2cBracket.addEventListener('click', () => doS2c('bracket'));

// main input → CABILN convert buttons
async function doToCabiln(notation) {
  if (!cabilnInput.value.trim()) return;
  const btn = notation === 'bracket' ? btnToCabilnBracket : btnToCabilnPct;
  const request = startConversion(btn);
  setStatus(cabilnStatus, '');
  try {
    await settleDrawingsForConversion(request);
    if (!request.current()) return;
    const txt = cabilnInput.value.trim();
    const inputFormat = notationSelect.value;
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
      setStatus(cabilnStatus, data.error);
    } else {
      useConvertedCabiln(data.cabiln, data.warning, data);
      setStatus(cabilnStatus, `Converted from ${data.from}: ${data.cabiln}`, 'ok');
    }
  } catch (e) {
    if (!request.current()) return;
    setStatus(cabilnStatus, 'Could not convert the input. Try again.');
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

// render helpers
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

// notation selector
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

// CABILN render
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
    : cabilnInput.className === 'err'
      ? 'Previous drawing — correct the input to update'
      : 'Previous drawing — drawing unavailable; try again';
}

function invalidateDocument() {
  clearTimeout(cabilnTimer);
  cabilnTimer = null;
  requests.cancel('main-render', 'sequence-edit');
  clearComparison();
  clearBuild();
  residueView.clearHighlight();
  clearExports();
  cabilnDrawing = null;
  lastCabiln = '';
  mainStale = true;
  renderCanvas.classList.add('stale');
  residueView.setStale(true);
  rerollSeed = 0;
  btnReroll.disabled = true;
  btnReroll.textContent = '⟳ Layout';
  cabilnInput.className = '';
  setStatus(cabilnStatus, '');
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
  const render = () => {
    cabilnTimer = null;
    return mode === 'cabiln' ? doRenderCabiln(seq) : doRenderForeign(seq);
  };
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
  residueView.setStale(false);
}

function mainRenderError(message, invalidInput = false) {
  mainStale = true;
  if (!hasMainDrawing) setInner(renderInner, `<div class="placeholder err">${escHtml(message)}</div>`);
  setStatus(cabilnStatus, message);
  cabilnInput.className = invalidInput ? 'err' : '';
}

async function doRenderForeign(txt) {
  const request = requests.start('main-render');
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
      mainRenderError(data.error, invalidInputResponse(res));
    } else {
      acceptDocumentContext(data.context);
      acceptMainDrawing(data.svg, txt, mode);
      residueView.clear();
      setStatus(cabilnStatus, `${data.format}: ${data.info || ''}`, 'ok');
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
  residueView.setStale(false);
  mainViewport.reset();
  setInner(renderInner, '<div class="placeholder">Start typing a sequence…</div>');
  residueView.clear();
  setMainProgress(false);
}

function acceptCabilnDrawing(data, source, view, sameDocument = false) {
  const displayedSequence = data.normalized_cabiln || source;
  if (data.normalized_cabiln) {
    projectChanged();
    replaceDocument({ text: displayedSequence, quality: null, canonical: null });
    saveDraft();
  }
  acceptDocumentContext(data.context);
  lastCabiln = displayedSequence;
  acceptMainDrawing(data.svg, displayedSequence, 'cabiln');
  rerollSeed = view.seed;
  setStatus(cabilnStatus, [data.info, data.normalization_note, ...(data.warnings || [])].filter(Boolean).join(' · '),
    data.warnings?.length ? 'warn' : 'ok');
  cabilnInput.className = 'ok';
  btnReroll.disabled = false;
  setExportReady(data.svg, data.mol_block);
  residueView.render(data, editor.present.quality);
  // Tabs can change the canvas height on the first render. Detect later resizes
  // against the accepted view, while retaining the dimensions of the SVG itself.
  cabilnDrawing = { source: displayedSequence, data, w: view.w, h: view.h,
    seed: view.seed, canvas: canvasSize(renderCanvas) };
  if (sameDocument) residueView.select(buildLeftRIdx, buildRightRIdx);
  if (verifyMode && lastSmiles) triggerVerify();
}

async function doRenderCabiln(seq) {
  const request = requests.start('main-render');
  const sameDocument = seq === displayedSource && displayedNotation === 'cabiln';
  lastCabiln = '';
  clearExports();
  cabilnDrawing = null;
  clearComparison();
  btnReroll.disabled = true;
  if (buildMode && !sameDocument) clearBuild();
  const view = { ...canvasSize(renderCanvas), seed: rerollSeed };
  try {
    const res  = await fetchCalculation('/render', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: seq, width: view.w, height: view.h, seed: view.seed }),
      signal: request.signal,
    });
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== seq ||
        notationSelect.value !== 'cabiln') return;
    if (data.error) {
      mainRenderError(data.error, invalidInputResponse(res));
      clearExports();
      btnReroll.disabled = true;
    } else {
      acceptCabilnDrawing(data, seq, view, sameDocument);
    }
  } catch (e) {
    if (!request.current()) return;
    mainRenderError('Could not reach the renderer. Your input is preserved.');
  } finally {
    if (request.current()) setMainProgress(false);
    request.finish();
  }
}

// reference render (verify mode) — auto-detects SMILES / BILN / HELM
function clearReference() {
  clearTimeout(smilesTimer);
  smilesTimer = null;
  requests.cancel('reference-render', 'sequence-edit');
  lastSmiles = '';
  smilesInput.className = '';
  setStatus(smilesStatus, '');
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
  setStatus(smilesStatus, 'Updating reference…');
  smilesTimer = setTimeout(() => {
    smilesTimer = null;
    doRenderRef(txt);
  }, 180);
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
  const request = requests.start('reference-render');
  showSpinner(smilesInner);
  setStatus(smilesStatus, `Loaded: ${file.name}`, 'ok');
  try {
    const text = await file.text();
    if (!request.current()) return;
    await renderMolReference(text, file.name, request);
  } catch (err) {
    if (!request.current()) return;
    setStatus(smilesStatus, 'Could not read the MOL/SDF file. Your peptide input is preserved.');
  } finally { request.finish(); }
});

async function renderMolReference(text, name, request = requests.start('reference-render')) {
  if (referenceOriginal?.kind !== 'mol' || referenceOriginal.content !== text ||
      referenceOriginal.name !== name) {
    projectChanged(false);
    referenceContext = contextHeader(referenceContext);
  }
  referenceOriginal = { kind: 'mol', content: text, name };
  saveDraft();
  showSpinner(smilesInner);
  setStatus(smilesStatus, `Reference: ${name || 'uploaded structure'}`);
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
      setStatus(smilesStatus, data.error);
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
    setStatus(smilesStatus, 'Failed to render the original MOL/SDF reference');
  } finally { request.finish(); }
}

async function doRenderRef(txt) {
  const request = requests.start('reference-render');
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
      setStatus(smilesStatus, data.error);
      smilesInput.className = invalidInputResponse(res) ? 'err' : '';
    } else {
      setInner(smilesInner, data.svg);
      lastSmiles = data.smiles || '';
      acceptReferenceContext(data.context);
      setStatus(smilesStatus, `${data.format}: ${data.info || ''}`, 'ok');
      smilesInput.className = 'ok';
      if (lastCabiln) triggerVerify();
    }
  } catch (e) {
    if (!request.current()) return;
    setStatus(smilesStatus, 'Could not reach the renderer. Your reference input is preserved.');
    smilesInput.className = '';
  } finally {
    request.finish();
  }
}

// verify comparison
function clearComparison() {
  requests.cancel('verify');
  compareBar.innerHTML = '';
}

async function triggerVerify() {
  clearComparison();
  if (!verifyMode || !lastSmiles || !lastCabiln) return;
  const request = requests.start('verify');
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

offerSavedDraft();
updateHistoryControls();
