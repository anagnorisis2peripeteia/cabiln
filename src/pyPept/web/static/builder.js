// Page preferences and reference lifetime.
const drawing = new DrawingState();
let darkMode   = true;
let hlEnabled  = true;
let cabilnTimer = null;
let smilesTimer = null;
let lastSmiles  = '';
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
  canHighlight: () => hlEnabled && !drawing.stale,
  onSelect: residue => build.selectResidue(residue),
});

const library = new MonomerLibrary({
  getFilters: () => build.filters,
  onUse: (abbr, explicit) => { if (!build.useMonomer(abbr, explicit)) insertAbbr(abbr); },
  onChanged: () => { if (build.filters.building) build.clear(); },
});

// Document data belongs to the editor; these timers belong to browser storage.
const DRAFT_KEY = 'cabiln.draft.v1';
const editor = new CabilnDocument();
const build = createBuildPanel({ library, residueView,
  getDocument: () => editor.present, commitDocument, isStale: () => drawing.stale,
});
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
  btnToBracket.hidden = btnToBranch.hidden = mode !== 'cabiln';
  btnToBracket.disabled = btnToBranch.disabled = mode !== 'cabiln';
  notationPolicy.hidden = mode !== 'cabiln';
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
    const response = await postCalculation('/prepare_project', { project }, request.signal);
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
    const response = await postCalculation('/prepare_project', { project }, request.signal);
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
    const response = await postCalculation('/validate_project', { project }, request.signal);
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
  if (verifyPanel.isOpen) restoreReferenceDrawing();
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

createPanel({ panel: helpPanel, button: btnHelp,
  closeButton: document.getElementById('help-close'),
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

const verifyPanel = createPanel({ panel: document.getElementById('verify-pane'), button: btnVerify,
  closeButton: document.getElementById('verify-close'), focus: smilesInput,
  onChange(open) {
    compareBar.hidden = !open;
    clearComparison();
    if (open && (smilesInput.value.trim() || referenceOriginal) && !lastSmiles) restoreReferenceDrawing();
    else if (open) triggerVerify();
  },
});

let examplesLoaded = false;
createPanel({ panel: examplesPanel, button: btnExamples, closeButton: examplesClose,
  onChange(open) { if (open && !examplesLoaded) loadExamples(); },
});

async function loadExamples() {
  const request = requests.start('examples');
  showLoading(examplesList, 'Loading examples…');
  try {
    const res = await fetchCalculation('/examples', { signal: request.signal });
    const data = await readResponse(res);
    if (!request.current()) return;
    if (!Array.isArray(data)) throw new Error(data.error || 'Invalid examples');
    renderExamples(data);
    examplesLoaded = true;
  } catch (e) {
    if (request.current()) showRetry(examplesList, 'Could not load examples.', loadExamples);
  } finally { request.finish(); }
}

function renderExamples(categories) {
  const rows = [];
  for (const cat of categories) {
    rows.push(`<div class="example-cat">${escHtml(cat.category)}</div>`);
    for (const item of cat.items) {
      const preview = item.cabiln.length > 55
        ? item.cabiln.slice(0, 52) + '…'
        : item.cabiln;
      rows.push(`<button type="button" class="example-row" data-cabiln="${escAttr(item.cabiln)}">
        <span class="example-name">${escHtml(item.name)}</span>
        <span class="example-desc">${escHtml(item.description)}</span>
        <span class="example-seq" title="${escAttr(item.cabiln)}">${escHtml(preview)}</span>
      </button>`);
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

// PNG download (client-side SVG → canvas → PNG)
btnPng.addEventListener('click', () => {
  if (!drawing.svg) return;
  const blob = new Blob([drawing.svg], { type: 'image/svg+xml;charset=utf-8' });
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
  if (!drawing.molBlock) return;
  const blob = new Blob([drawing.molBlock], { type: 'chemical/x-mdl-molfile' });
  const a = document.createElement('a');
  a.download = 'structure.mol';
  a.href = URL.createObjectURL(blob);
  a.click();
});

// notation conversion
async function settleDrawingsForConversion(request) {
  // Aborting a running drawing retires the server's chemistry worker. Let the
  // current drawing (and its automatic Verify) finish before competing with it.
  while (request.current()) {
    if (requests.has('main-render', 'reference-render', 'verify')) {
      if (!await request.waitFor('main-render', 'reference-render', 'verify')) return;
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
  if (!snapshot || snapshot !== drawing.snapshot || drawing.stale ||
      requests.has('main-render') || !drawing.matches(snapshot.source, 'cabiln') || data.source_echo !== snapshot.source ||
      !drawing.svg || !drawing.molBlock || view.w !== snapshot.canvas.w || view.h !== snapshot.canvas.h ||
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
  const button = target === 'bracket' ? btnToBracket : btnToBranch;
  const request = startConversion(button, () => {
    // Let the cancelling action schedule its own drawing or formatter first.
    Promise.resolve().then(() => {
      if (deferredDrawing && drawing.stale && cabilnTimer === null &&
          !requests.has('main-render') && !requests.has('sequence-edit') &&
          notationSelect.value === 'cabiln' && cabilnInput.value.trim() === original) renderDocument();
    });
  });
  conversionProgressLabel.textContent = 'Converting notation and checking the structure…';
  // The formatter validates and normalizes unsent input. Draw its result once,
  // or draw the original input if formatting fails. Stale input also covers a
  // second formatter click taking over the first click's deferred drawing.
  deferredDrawing = cabilnTimer !== null || (drawing.stale && !requests.has('main-render'));
  clearTimeout(cabilnTimer);
  cabilnTimer = null;
  let conversionError = '';
  try {
    await settleDrawingsForConversion(request);
    if (!request.current()) return;
    const val = cabilnInput.value.trim();
    const snapshot = drawing.snapshot;
    const res = await postCalculation('/convert_notation', { cabiln: val, target, canonical }, request.signal);
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== val) return;
    if (data.error) { conversionError = data.error; return; }
    if (data.result) {
      deferredDrawing = false;
      const reused = drawingForNotation(snapshot, data);
      const document = { text: data.result, notation: 'cabiln', warning: editor.present.warning,
        quality: null, canonical: data.canonical || null, context: data.context || editor.present.context };
      recordDocument(document);
      if (reused) {
        clearComparison();
        build.clear();
        residueView.clearHighlight();
        acceptCabilnDrawing(reused, data.result, snapshot);
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
  if (!drawing.cabiln) return;
  // The default already uses Indigo; start with the alternate CoordGen layout.
  drawing.nextLayout();
  btnReroll.textContent = drawing.seed % 2 === 1 ? '⟳ Indigo' : '⟳ CoordGen';
  setMainProgress(true);
  doRenderCabiln(drawing.cabiln);
});

// SMILES → CABILN conversion
function startConversion(button, onEnd = () => {}) {
  requests.cancel('sequence-edit');
  const request = requests.start('sequence-edit', () => {
    button.setAttribute('aria-busy', 'false');
    button.disabled = false;
    conversionProgress.hidden = true;
    onEnd();
    updateNotationControls();
  });
  button.setAttribute('aria-busy', 'true');
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
    const res = await postCalculation('/smiles_to_cabiln', { smiles, notation }, request.signal);
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
    const res = await postCalculation('/to_cabiln', {
      input: txt, input_format: inputFormat, notation,
    }, request.signal);
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

function updateDrawingControls() {
  btnPng.disabled = !drawing.svg;
  btnMol.disabled = !drawing.molBlock;
  btnReroll.disabled = !drawing.cabiln;
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
  showLoading(inner, 'Drawing structure…');
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
  renderProgress.hidden = !pending && !(drawing.stale && drawing.hasDrawing);
  renderProgressLabel.textContent = pending
    ? (drawing.hasDrawing ? 'Updating — previous drawing shown' : 'Drawing structure…')
    : cabilnInput.className === 'err'
      ? 'Previous drawing — correct the input to update'
      : 'Previous drawing — drawing unavailable; try again';
}

function invalidateDocument() {
  clearTimeout(cabilnTimer);
  cabilnTimer = null;
  requests.cancel('main-render', 'sequence-edit');
  clearComparison();
  build.clear();
  residueView.clearHighlight();
  drawing.invalidate();
  updateDrawingControls();
  renderCanvas.classList.add('stale');
  residueView.setStale(true);
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
  if (!drawing.hasDrawing) showSpinner(renderInner);
  setMainProgress(true);
  const render = () => {
    cabilnTimer = null;
    return mode === 'cabiln' ? doRenderCabiln(seq) : doRenderForeign(seq);
  };
  if (immediate) render();
  else cabilnTimer = setTimeout(render, 180);
}

function displayMainDrawing(svg) {
  setInner(renderInner, svg);
  renderCanvas.classList.remove('stale');
  residueView.setStale(false);
}

function mainRenderError(message, invalidInput = false) {
  drawing.invalidate();
  updateDrawingControls();
  if (!drawing.hasDrawing) setInner(renderInner, `<div class="placeholder err">${escHtml(message)}</div>`);
  setStatus(cabilnStatus, message, invalidInput ? 'error' : 'warn');
  if (!invalidInput) showRetry(cabilnStatus, message, () => renderDocument(true));
  cabilnInput.className = invalidInput ? 'err' : '';
}

async function doRenderForeign(txt) {
  const request = requests.start('main-render');
  const mode = notationSelect.value;
  drawing.begin();
  updateDrawingControls();
  clearComparison();
  const { w, h } = canvasSize(renderCanvas);
  try {
    if (requests.has('library', 'reactions') &&
        !await request.waitFor('library', 'reactions')) return;
    const res = await postCalculation('/render_reference', {
      input: txt, input_format: mode, width: w, height: h,
    }, request.signal);
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== txt ||
        notationSelect.value !== mode) return;
    if (data.error) {
      mainRenderError(data.error, invalidInputResponse(res));
    } else {
      acceptDocumentContext(data.context);
      displayMainDrawing(data.svg);
      drawing.accept(data, txt, mode);
      updateDrawingControls();
      residueView.clear();
      setStatus(cabilnStatus, '', 'ok');
      showStructureInfo(cabilnStatus, data, { lead: data.format });
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
  drawing.clear();
  updateDrawingControls();
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
  displayMainDrawing(data.svg);
  setStatus(cabilnStatus, '', data.warnings?.length ? 'warn' : 'ok');
  showStructureInfo(cabilnStatus, data);
  cabilnInput.className = 'ok';
  residueView.render(data, editor.present.quality);
  // Tabs can change the canvas height on the first render. Detect later resizes
  // against the accepted view, while retaining the dimensions of the SVG itself.
  drawing.accept(data, displayedSequence, 'cabiln', { w: view.w, h: view.h,
    seed: view.seed, canvas: canvasSize(renderCanvas) });
  updateDrawingControls();
  if (sameDocument) residueView.select(...build.selection);
  if (verifyPanel.isOpen && lastSmiles) triggerVerify();
}

async function doRenderCabiln(seq) {
  const request = requests.start('main-render');
  const sameDocument = drawing.matches(seq, 'cabiln');
  drawing.begin();
  updateDrawingControls();
  clearComparison();
  if (build.filters.building && !sameDocument) build.clear();
  const view = { ...canvasSize(renderCanvas), seed: drawing.seed };
  try {
    if (requests.has('library', 'reactions') &&
        !await request.waitFor('library', 'reactions')) return;
    const res = await postCalculation('/render', {
      cabiln: seq, width: view.w, height: view.h, seed: view.seed,
    }, request.signal);
    const data = await readResponse(res);
    if (!request.current() || cabilnInput.value.trim() !== seq ||
        notationSelect.value !== 'cabiln') return;
    if (data.error) {
      mainRenderError(data.error, invalidInputResponse(res));
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
    const res = await postCalculation('/render_mol', { mol_block: text, width: w, height: h }, request.signal);
    const data = await readResponse(res);
    if (!request.current()) return;
    if (data.error) {
      setInner(smilesInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      setStatus(smilesStatus, data.error, 'error');
      if (!invalidInputResponse(res)) showRetry(smilesStatus, data.error, restoreReferenceDrawing);
    } else {
      setInner(smilesInner, data.svg);
      lastSmiles = data.smiles || '';
      if (smilesInput.value !== lastSmiles) projectChanged(false);
      smilesInput.value = lastSmiles;
      acceptReferenceContext(data.context);
      saveDraft();
      smilesInput.className = 'ok';
      if (drawing.cabiln) triggerVerify();
    }
  } catch (err) {
    if (!request.current()) return;
    showRetry(smilesStatus, 'Could not render the original MOL/SDF reference.', restoreReferenceDrawing);
  } finally { request.finish(); }
}

async function doRenderRef(txt) {
  const request = requests.start('reference-render');
  lastSmiles = '';
  clearComparison();
  const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
  try {
    const res = await postCalculation('/render_reference', { input: txt, width: w, height: h }, request.signal);
    const data = await readResponse(res);
    if (!request.current() || smilesInput.value.trim() !== txt) return;
    if (data.error) {
      setInner(smilesInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      setStatus(smilesStatus, data.error, 'error');
      if (!invalidInputResponse(res)) showRetry(smilesStatus, data.error, restoreReferenceDrawing);
      smilesInput.className = invalidInputResponse(res) ? 'err' : '';
    } else {
      setInner(smilesInner, data.svg);
      lastSmiles = data.smiles || '';
      acceptReferenceContext(data.context);
      setStatus(smilesStatus, '', 'ok');
      showStructureInfo(smilesStatus, data, { lead: data.format });
      smilesInput.className = 'ok';
      if (drawing.cabiln) triggerVerify();
    }
  } catch (e) {
    if (!request.current()) return;
    showRetry(smilesStatus, 'Could not reach the renderer. Your reference input is preserved.', restoreReferenceDrawing);
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
  if (!verifyPanel.isOpen || !lastSmiles || !drawing.cabiln) return;
  const request = requests.start('verify');
  const smiles = lastSmiles;
  const cabiln = drawing.cabiln;
  try {
    const res = await postCalculation('/verify', { smiles, cabiln }, request.signal);
    const data = await readResponse(res);
    if (!request.current() || lastSmiles !== smiles || drawing.cabiln !== cabiln) return;
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
