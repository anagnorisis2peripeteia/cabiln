// Shared DOM presentation helpers; no document or chemistry state.
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

function escHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
function escAttr(s) { return escHtml(s); }
