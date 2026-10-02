// Shared DOM presentation helpers; no document or chemistry state.
function makeZoomable(canvas, inner) {
  let scale = 1, tx = 0, ty = 0;
  const pointers = new Map();
  let moved = false;

  function applyTransform() {
    inner.style.transform = `translate(${tx}px,${ty}px) scale(${scale})`;
  }

  function zoom(factor) {
    scale = Math.min(Math.max(scale * factor, 0.2), 20);
    applyTransform();
  }

  function reset() {
    scale = 1; tx = 0; ty = 0;
    applyTransform();
  }

  function focus(elements) {
    if (!elements.length) return;
    reset();
    const boxes = [...elements].map(element => element.getBoundingClientRect());
    const left = Math.min(...boxes.map(box => box.left)), right = Math.max(...boxes.map(box => box.right));
    const top = Math.min(...boxes.map(box => box.top)), bottom = Math.max(...boxes.map(box => box.bottom));
    const viewport = canvas.getBoundingClientRect();
    scale = Math.max(0.2, Math.min(20, (viewport.width - 48) / Math.max(1, right - left),
      (viewport.height - 48) / Math.max(1, bottom - top)));
    tx = (viewport.left + viewport.width / 2 - (left + right) / 2) * scale;
    ty = (viewport.top + viewport.height / 2 - (top + bottom) / 2) * scale;
    applyTransform();
    canvas.scrollIntoView({ block: 'center' });
  }

  canvas.addEventListener('wheel', e => {
    e.preventDefault();
    zoom(e.deltaY < 0 ? 1.12 : 1 / 1.12);
  }, { passive: false });

  canvas.addEventListener('pointerdown', e => {
    if (e.button !== 0) return;
    if (!pointers.size) moved = false;
    pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  });
  window.addEventListener('pointermove', e => {
    const previous = pointers.get(e.pointerId);
    if (!previous) return;
    const next = { x: e.clientX, y: e.clientY };
    const dx = next.x - previous.x, dy = next.y - previous.y;
    if (!moved && Math.hypot(dx, dy) < 3) return;
    moved = true;
    canvas.setPointerCapture(e.pointerId);
    canvas.style.cursor = 'grabbing';
    const other = [...pointers.entries()].find(([id]) => id !== e.pointerId)?.[1];
    if (other) {
      const distance = Math.hypot(previous.x - other.x, previous.y - other.y);
      const nextDistance = Math.hypot(next.x - other.x, next.y - other.y);
      const rect = canvas.getBoundingClientRect();
      const x = (previous.x + other.x) / 2 - rect.left - rect.width / 2;
      const y = (previous.y + other.y) / 2 - rect.top - rect.height / 2;
      const nextScale = Math.min(Math.max(scale * nextDistance / (distance || 1), 0.2), 20);
      tx = x - (x - tx) * nextScale / scale + dx / 2;
      ty = y - (y - ty) * nextScale / scale + dy / 2;
      scale = nextScale;
    } else {
      tx += dx; ty += dy;
    }
    pointers.set(e.pointerId, next);
    applyTransform();
  });
  for (const type of ['pointerup', 'pointercancel', 'lostpointercapture']) {
    window.addEventListener(type, e => {
      pointers.delete(e.pointerId);
      if (!pointers.size) canvas.style.cursor = 'grab';
    });
  }
  window.addEventListener('blur', () => {
    for (const id of pointers.keys()) {
      if (canvas.hasPointerCapture(id)) canvas.releasePointerCapture(id);
    }
    pointers.clear();
    canvas.style.cursor = 'grab';
  });
  canvas.addEventListener('click', e => {
    // A completed pan must not also select the atom under the pointer.
    if (moved) { e.preventDefault(); e.stopImmediatePropagation(); moved = false; }
  }, true);
  canvas.addEventListener('dblclick', reset);
  canvas.tabIndex = 0;
  canvas.style.touchAction = 'none';
  canvas.addEventListener('keydown', e => {
    if (e.target !== canvas || e.ctrlKey || e.metaKey || e.altKey) return;
    const actions = {
      '+': () => zoom(1.2), '=': () => zoom(1.2), '-': () => zoom(1 / 1.2),
      '0': reset, Home: reset,
      ArrowLeft: () => { tx -= 30; }, ArrowRight: () => { tx += 30; },
      ArrowUp: () => { ty -= 30; }, ArrowDown: () => { ty += 30; },
    };
    if (!actions[e.key]) return;
    e.preventDefault();
    actions[e.key]();
    applyTransform();
  });

  document.querySelectorAll(`[data-viewport="${canvas.id}"]`).forEach(controls => {
    controls.addEventListener('click', e => {
      const action = e.target.closest('[data-zoom]')?.dataset.zoom;
      if (action === 'in') zoom(1.2);
      if (action === 'out') zoom(1 / 1.2);
      if (action === 'reset') reset();
    });
  });
  return { reset, focus };
}

// Nonmodal panels retain independent lifetimes; Escape dismisses the panel
// containing focus, or the most recently opened panel if focus is elsewhere.
const createPanel = (() => {
  const openPanels = [];
  window.addEventListener('keydown', event => {
    if (event.key !== 'Escape' || event.defaultPrevented) return;
    const panel = [...openPanels].reverse().find(item => item.containsFocus()) || openPanels.at(-1);
    if (panel) { event.preventDefault(); panel.close({ focus: true }); }
  });
  return ({ panel, button, closeButton, focus = closeButton, onChange = () => {}, toggle = true }) => {
    panel.hidden = true;
    const control = {
      get isOpen() { return !panel.hidden; },
      containsFocus: () => panel.contains(document.activeElement) || button === document.activeElement,
      open({ focus: moveFocus = true } = {}) {
        if (control.isOpen) return;
        panel.hidden = false;
        panel.classList.add('open');
        button.classList.add('active');
        button.setAttribute('aria-expanded', 'true');
        openPanels.push(control);
        onChange(true);
        if (moveFocus) focus?.focus();
      },
      close({ focus: moveFocus = control.containsFocus() } = {}) {
        if (!control.isOpen) return;
        panel.hidden = true;
        panel.classList.remove('open');
        button.classList.remove('active');
        button.setAttribute('aria-expanded', 'false');
        openPanels.splice(openPanels.indexOf(control), 1);
        onChange(false);
        if (moveFocus) button.focus();
      },
    };
    if (toggle) button.addEventListener('click', () => control.isOpen ? control.close() : control.open());
    closeButton?.addEventListener('click', () => control.close({ focus: true }));
    return control;
  };
})();

function showLoading(element, message) {
  element.innerHTML = `<div class="loading-state" role="status"><span class="spinner" aria-hidden="true"></span><span>${escHtml(message)}</span></div>`;
}

function setInputState(element, state = '') {
  element.className = state;
  element.setAttribute('aria-invalid', String(state === 'err'));
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.download = filename;
  link.href = url;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function showRetry(element, message, retry) {
  element.textContent = message + ' ';
  element.title = message;
  const button = document.createElement('button');
  button.type = 'button';
  button.className = 'hbtn retry-button';
  button.textContent = 'Retry';
  button.addEventListener('click', retry);
  element.appendChild(button);
}

function showStructureInfo(element, data, { lead = '', notices = [] } = {}) {
  const warnings = [...notices, ...(data.warnings || [])].filter(Boolean);
  const info = [lead, data.info].filter(Boolean).join(' · ');
  element.innerHTML = `<span class="structure-info">${escHtml(info)}</span>` +
    warnings.map(warning => `<p class="structure-warning">${escHtml(warning)}</p>`).join('') +
    (data.normalization_note ? `<details class="structure-details"><summary>Notation normalized · Details</summary><p>${escHtml(data.normalization_note)}</p></details>` : '');
}

function escHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
function escAttr(s) { return escHtml(s); }
