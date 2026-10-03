let detectedData = null;
let detectedSmiles = '';
let formRevision = 0;
let registering = false;
let registeredPayload = null;
let attachmentChoices = [];
let choiceSmiles = '';

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
  const registered = payload && JSON.stringify(payload) === registeredPayload;
  btnRegister.disabled = registering || !payload || !payload.abbr ||
    !payload.name || registered;
  btnRegister.textContent = registering ? 'Registering…' :
    registered ? '✓ Registered' : 'Register monomer';
  registrationGuide.textContent = attachmentChoices.length && !payload ? 'Choose which numbered attachments this monomer should use.' :
    !payload ? 'Preview the current SMILES before registering.' :
    !payload.abbr || !payload.name ? 'Enter an abbreviation and full name to register.' :
    registered ? 'Available in Library when you return to the renderer.' :
    'Adds this monomer and its detected sites to this server’s library.';
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
  if (body === registeredPayload) return;
  registering = true;
  updateRegisterButton();
  try {
    const res  = await fetch('/register_monomer', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body,
    });
    const data = await readResponse(res);
    if (data.error) {
      if (revision === formRevision) showStatus('err', data.error);
    } else {
      registeredPayload = body;
      if (revision === formRevision) {
        showStatus('ok', `${payload.abbr} registered successfully. Monomer count: ${data.total}`);
      }
    }
  } catch (e) {
    if (revision === formRevision) showStatus('err', 'Could not confirm registration. Your entries are preserved. Check the Library for this abbreviation before trying again.');
  } finally {
    registering = false;
    updateRegisterButton();
  }
});

function showStatus(cls, msg) {
  statusMsg.className = cls;
  statusMsg.textContent = msg;
  statusMsg.style.display = '';
}
