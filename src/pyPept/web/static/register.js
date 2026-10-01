let detectedData = null;
let detectedSmiles = '';
let formRevision = 0;
let registering = false;
let registeredPayload = null;

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

function registrationPayload() {
  if (!detectedData || detectedSmiles !== smilesIn.value.trim()) return null;
  return {
    chuckles: detectedData.chuckles,
    chem_types: detectedData.chem_types,
    leaving: detectedData.leaving,
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
}

function invalidatePreview() {
  requests.cancel('registration-preview');
  detectedData = null;
  detectedSmiles = '';
  registeredPayload = null;
  chucklesOut.value = '';
  detDisp.innerHTML = '';
  prevCanvas.innerHTML = '';
  prevSec.style.display = 'none';
  smilesIn.className = '';
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

btnPreview.addEventListener('click', async () => {
  const smi = smilesIn.value.trim();
  invalidatePreview();
  if (!smi) return;
  const request = requests.start('registration-preview');
  const current = () => request.current() && smilesIn.value.trim() === smi;
  btnPreview.disabled = true;
  btnPreview.textContent = 'Detecting…';
  prevCanvas.innerHTML = '<div class="preview-message">Analysing…</div>';
  detDisp.innerHTML = '';
  prevSec.style.display = 'block';
  try {
    const res = await postCalculation('/preview_monomer', { smiles: smi }, request.signal);
    const data = await readResponse(res);
    if (!current()) return;
    if (data.error) {
      prevCanvas.innerHTML = `<div class="preview-message err">${escHtml(data.error)}</div>`;
      smilesIn.className = 'err';
      btnRegister.disabled = true;
    } else {
      detectedData = data;
      detectedSmiles = smi;
      prevCanvas.innerHTML = data.svg;
      chucklesOut.value = data.chuckles;
      smilesIn.className = 'ok';
      renderDetected(data.chem_types, data.leaving);
      updateRegisterButton();
    }
  } catch (e) {
    if (!current()) return;
    detectedData = null;
    detectedSmiles = '';
    chucklesOut.value = '';
    detDisp.innerHTML = '';
    prevCanvas.innerHTML = '<div class="preview-message err">Server error</div>';
  } finally {
    if (current()) {
      btnPreview.disabled = false;
      btnPreview.textContent = 'Preview & detect R-groups';
      updateRegisterButton();
    }
    request.finish();
  }
});

function renderDetected(chem_types, leaving) {
  const lines = Object.entries(chem_types).sort(([a],[b]) => +a - +b).map(([slot, ct]) => {
    const lg = leaving[slot] ? `<span class="lg">  LG: ${escHtml(leaving[slot])}</span>` : '';
    return `<div><span class="slot">R${slot}</span>  ${escHtml(ct)}${lg}</div>`;
  });
  detDisp.innerHTML = lines.join('');
}

btnRegister.addEventListener('click', async () => {
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
    if (revision === formRevision) showStatus('err', 'Server error');
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
