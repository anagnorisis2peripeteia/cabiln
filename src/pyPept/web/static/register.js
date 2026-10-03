let detectedData = null;
let detectedSmiles = '';
let formRevision = 0;
let registering = false;
let registeredPayload = null;
let attachmentChoices = [];
let choiceSmiles = '';
let persistentAllowed = false;
let scopeReady = false;
let registeredDestination = null;
let registeredAbbr = '';

const smilesIn   = document.getElementById('smiles-in');
const btnPreview = document.getElementById('btn-preview');
const prevSec    = document.getElementById('preview-section');
const prevCanvas = document.getElementById('preview-canvas');
const detDisp    = document.getElementById('detected-display');
const chucklesOut= document.getElementById('chuckles-out');
const abbrIn     = document.getElementById('abbr-in');
const nameIn     = document.getElementById('name-in');
const typeIn     = document.getElementById('type-in');
const subtypeIn  = document.getElementById('subtype-in');
const btnRegister= document.getElementById('btn-register');
const statusMsg  = document.getElementById('status-msg');
const previewError = document.getElementById('preview-error');
const registrationGuide = document.getElementById('registration-guide');
const choiceSection = document.getElementById('attachment-choices');
const choiceList = document.getElementById('attachment-choice-list');
const doneButton = document.getElementById('registration-done');
const destinationIn = document.getElementById('destination-in');
destinationIn.value = 'session';

function destination() {
  return persistentAllowed && destinationIn.value === 'installed' ? 'installed' : 'session';
}

async function loadRegistrationScope() {
  try {
    const response = await fetchCalculation('/capabilities', { timeoutMs: 5000 });
    const capabilities = await readResponse(response);
    persistentAllowed = capabilities.registration === true;
    document.getElementById('registration-destination').hidden = !persistentAllowed;
    if (persistentAllowed) destinationIn.value = 'installed';
  } catch (_) { /* Temporary registration remains available without administrative access. */ }
  scopeReady = true;
  updateRegisterButton();
}
loadRegistrationScope();

destinationIn.addEventListener('change', () => { formRevision++; updateRegisterButton(); });

function returnToLibrary(event) {
  if (window.parent && window.parent !== window) {
    event.preventDefault();
    window.parent.postMessage({ type: 'cabiln-registration-close', abbr: registeredAbbr }, window.location.origin);
  } else if (event.currentTarget === doneButton) window.location.href = '/';
}
document.getElementById('registration-return').addEventListener('click', returnToLibrary);
doneButton.addEventListener('click', returnToLibrary);
window.addEventListener('keydown', event => {
  if (event.key === 'Escape' && window.parent && window.parent !== window) {
    event.preventDefault();
    window.parent.postMessage({ type: 'cabiln-registration-close', abbr: registeredAbbr }, window.location.origin);
  }
});

function registrationPayload() {
  if (!detectedData || detectedSmiles !== smilesIn.value.trim()) return null;
  return {
    chuckles: detectedData.chuckles,
    chem_types: detectedData.chem_types,
    leaving: detectedData.leaving,
    activation_policy: detectedData.activation_policy,
    abbr: abbrIn.value.trim(),
    name: nameIn.value.trim(),
    type: typeIn.value,
    subtype: subtypeIn.value,
  };
}

function updateRegisterButton() {
  const payload = registrationPayload();
  const installed = destination() === 'installed';
  const registered = payload && JSON.stringify(payload) === registeredPayload &&
    registeredDestination === destination();
  doneButton.hidden = !registered;
  btnRegister.disabled = !scopeReady || registering || !payload || !payload.abbr ||
    !payload.name || registered;
  btnRegister.textContent = registering ? 'Adding…' :
    registered ? '✓ Added' : installed ? 'Install monomer' : 'Add to this tab';
  registrationGuide.textContent = attachmentChoices.length && !payload ? 'Choose which numbered attachments this monomer should use.' :
    !payload ? 'Preview the current SMILES before registering.' :
    !payload.abbr || !payload.name ? 'Enter an abbreviation and full name to register.' :
    registered ? 'Available in Library when you return to the renderer.' :
    installed ? 'Permanently adds this monomer to the installed library.' :
    'Available only in this tab. Save a project to keep its definition.';
  document.getElementById('registration-scope').textContent = installed
    ? 'Installed library: additions remain available after restarting the app.'
    : 'Custom monomers stay in this tab. Save a project to keep them after closing it.';
}

function invalidatePreview() {
  requests.cancel('registration-preview');
  detectedData = null;
  detectedSmiles = '';
  registeredPayload = null;
  attachmentChoices = [];
  choiceSmiles = '';
  choiceSection.hidden = true;
  choiceList.innerHTML = '';
  chucklesOut.value = '';
  detDisp.innerHTML = '';
  prevCanvas.innerHTML = '';
  prevCanvas.hidden = false;
  prevSec.hidden = true;
  prevCanvas.setAttribute('aria-busy', 'false');
  setInputState(smilesIn, '');
  previewError.textContent = '';
  btnPreview.disabled = false;
  btnPreview.textContent = 'Preview & detect R-groups';
  updateRegisterButton();
}

smilesIn.addEventListener('input', () => {
  formRevision++;
  statusMsg.style.display = 'none';
  invalidatePreview();
});
for (const field of [abbrIn, nameIn, typeIn, subtypeIn]) {
  const event = field === typeIn || field === subtypeIn ? 'change' : 'input';
  field.addEventListener(event, () => {
    formRevision++;
    statusMsg.style.display = 'none';
    updateRegisterButton();
  });
}

document.getElementById('preview-form').addEventListener('submit', async event => {
  event.preventDefault();
  const smi = smilesIn.value.trim();
  invalidatePreview();
  if (!smi) {
    previewError.textContent = 'Enter a SMILES string first, for example NCC(=O)O for glycine.';
    setInputState(smilesIn, 'err');
    smilesIn.focus();
    return;
  }
  const request = requests.start('registration-preview');
  const current = () => request.current() && smilesIn.value.trim() === smi;
  btnPreview.disabled = true;
  btnPreview.textContent = 'Detecting…';
  showLoading(prevCanvas, 'Detecting attachment sites…');
  prevCanvas.setAttribute('aria-busy', 'true');
  detDisp.innerHTML = '';
  prevSec.hidden = false;
  try {
    const res = await postCalculation('/preview_monomer', { smiles: smi }, request.signal);
    const data = await readResponse(res);
    if (!current()) return;
    if (data.error) {
      previewError.textContent = data.error + ' Check the SMILES or try Preview again.';
      prevSec.hidden = true;
      if (res.status === 400 || res.status === 422) {
        setInputState(smilesIn, 'err');
      }
      btnRegister.disabled = true;
    } else if (data.choices) {
      attachmentChoices = data.choices;
      choiceSmiles = smi;
      choiceSection.hidden = false;
      prevCanvas.hidden = true;
      document.getElementById('attachment-choice-help').textContent = data.message;
      choiceList.innerHTML = data.choices.map((choice, index) =>
        `<label class="attachment-choice"><input type="radio" name="attachment-choice" value="${index}">` +
        `<span>Option ${index + 1}</span><span class="choice-drawing" aria-hidden="true">${choice.svg}</span></label>`
      ).join('');
    } else {
      usePreview(data, smi);
    }
  } catch (e) {
    if (!current()) return;
    detectedData = null;
    detectedSmiles = '';
    chucklesOut.value = '';
    detDisp.innerHTML = '';
    prevSec.hidden = true;
    previewError.textContent = 'Could not reach the preview service. Your input is preserved. Check your connection and try Preview again.';
  } finally {
    if (current()) {
      btnPreview.disabled = false;
      btnPreview.textContent = 'Preview & detect R-groups';
      prevCanvas.setAttribute('aria-busy', 'false');
      updateRegisterButton();
    }
    request.finish();
  }
});

function usePreview(data, smi) {
  detectedData = data;
  detectedSmiles = smi;
  prevCanvas.hidden = false;
  prevCanvas.innerHTML = data.svg;
  chucklesOut.value = data.chuckles;
  setInputState(smilesIn, 'ok');
  renderDetected(data.chem_types, data.leaving);
  updateRegisterButton();
}

choiceList.addEventListener('change', event => {
  if (choiceSmiles !== smilesIn.value.trim() || event.target.name !== 'attachment-choice') return;
  const choice = attachmentChoices[Number(event.target.value)];
  if (!choice) return;
  formRevision++;
  statusMsg.style.display = 'none';
  usePreview(choice, choiceSmiles);
});

function renderDetected(chem_types, leaving) {
  const lines = Object.entries(chem_types).sort(([a],[b]) => +a - +b).map(([slot, ct]) => {
    const lg = leaving[slot] ? `<span class="lg"> · leaving group: ${escHtml(leaving[slot])}</span>` : '';
    return `<div><span class="slot">R${slot}</span> ${escHtml(ct.replaceAll('_', ' '))}${lg}</div>`;
  });
  detDisp.innerHTML = lines.join('');
}

document.getElementById('registration-form').addEventListener('submit', async event => {
  event.preventDefault();
  if (registering) return;
  const payload = registrationPayload();
  if (!payload || !payload.abbr || !payload.name) {
    showStatus('err', 'Preview the current SMILES and fill in abbreviation and name.');
    return;
  }
  const revision = formRevision;
  const body = JSON.stringify(payload);
  const target = destination();
  if (body === registeredPayload && registeredDestination === target) return;
  registering = true;
  updateRegisterButton();
  try {
    let data;
    if (target === 'installed') {
      const res = await fetch('/register_monomer', {
        method: 'POST', headers: { 'Content-Type': 'application/json' }, body,
      });
      data = await readResponse(res);
    } else data = await CabilnLibrary.add(payload);
    if (data.error) {
      if (revision === formRevision) showStatus('err', data.error);
    } else {
      registeredPayload = body;
      registeredDestination = target;
      registeredAbbr = payload.abbr;
      if (target === 'installed' && window.parent && window.parent !== window) {
        window.parent.postMessage({ type: 'cabiln-monomer-added', abbr: payload.abbr }, window.location.origin);
      }
      if (revision === formRevision) {
        const message = target === 'installed'
          ? `${payload.abbr} installed. Return to Library to use it.`
          : `${payload.abbr} added to this tab. Return to Library to use it. Save a project to keep it.`;
        showStatus('ok', message + (target === 'session' && !CabilnLibrary.retained
          ? ' Tab storage is unavailable; save before reloading.' : ''));
      }
    }
  } catch (e) {
    if (revision === formRevision) showStatus('err', e.message || 'Could not confirm registration. Your entries are preserved. Check the Library before trying again.');
  } finally {
    registering = false;
    updateRegisterButton();
    if (revision === formRevision && !doneButton.hidden) doneButton.focus();
  }
});

function showStatus(cls, msg) {
  statusMsg.className = cls;
  statusMsg.textContent = msg;
  statusMsg.style.display = '';
}
