#!/usr/bin/env python3
"""
Live CABILN -> 2D structure renderer.

    pip install fastapi uvicorn
    python tools/live_renderer.py

Modes
-----
Render  — compact input bar at top, full-width structure canvas below.
Verify  — side-by-side SMILES (left) vs CABILN (right); canonical SMILES
          comparison tells you whether the structures are identical.

Controls
--------
  Dark mode toggle  — moon icon in header; inverts the canvas colours.
  Zoom / pan        — scroll wheel to zoom; click-drag to pan; double-click to reset.
"""

__credits__ = ["Cameron Beesley"]
__license__ = "MIT"

import threading
import time
import webbrowser

try:
    import uvicorn
    from fastapi import FastAPI
    from fastapi.responses import HTMLResponse, JSONResponse
    from pydantic import BaseModel
except ImportError:
    raise SystemExit("pip install fastapi uvicorn  (then re-run)")

_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>CABILN Live Renderer</title>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  html, body {
    height: 100%;
    font-family: system-ui, -apple-system, sans-serif;
    background: #12192b;
    color: #d0dce8;
    display: flex;
    flex-direction: column;
    overflow: hidden;
  }

  /* ── header ── */
  header {
    flex-shrink: 0;
    padding: 7px 14px;
    background: #0d1422;
    border-bottom: 1px solid #1e3050;
    font-size: 13px;
    letter-spacing: .07em;
    color: #7aaeff;
    display: flex;
    align-items: center;
    gap: 10px;
  }
  header .title { flex: 1; }
  header span   { color: #3a5580; font-weight: 400; }
  .hbtn {
    background: #1a2640;
    border: 1px solid #2a4070;
    border-radius: 5px;
    color: #8ab4e8;
    font-size: 12px;
    padding: 3px 10px;
    cursor: pointer;
    transition: background .15s, color .15s;
    white-space: nowrap;
  }
  .hbtn:hover       { background: #263550; color: #c8daf0; }
  .hbtn.active      { background: #2a4a80; border-color: #4a6faa; color: #d8e8ff; }

  /* ── input bar ── */
  #input-bar {
    flex-shrink: 0;
    background: #0d1422;
    border-bottom: 2px solid #1e3050;
    padding: 7px 14px 5px;
    display: flex;
    flex-direction: column;
    gap: 3px;
  }
  #input-bar-row { display: flex; align-items: center; gap: 9px; }
  .seq-label {
    font-size: 10px;
    text-transform: uppercase;
    letter-spacing: .12em;
    color: #3a5580;
    white-space: nowrap;
    flex-shrink: 0;
  }
  textarea {
    width: 100%;
    background: #0a1018;
    color: #c8daf0;
    border: 1px solid #1e3050;
    border-radius: 5px;
    padding: 6px 11px;
    font-family: "Cascadia Code", "Fira Mono", "Courier New", monospace;
    font-size: 13px;
    line-height: 1.6;
    resize: vertical;
    min-height: 32px;
    height: 62px;
    max-height: 220px;
    outline: none;
    transition: border-color .15s;
  }
  textarea:focus { border-color: #3a6fd8; }
  textarea.err   { border-color: #d9534f; }
  textarea.ok    { border-color: #28a745; }
  .statusbar {
    font-size: 11px;
    font-family: "Cascadia Code", "Fira Mono", monospace;
    min-height: 15px;
    color: #d9534f;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    padding-left: 2px;
  }
  .statusbar.ok { color: #3dbe6c; }

  /* ── main area ── */
  #main { flex: 1; display: flex; min-height: 0; }

  /* ── render mode ── */
  #render-pane {
    flex: 1;
    display: flex;
    min-height: 0;
    position: relative;
  }

  /* ── verify mode ── */
  #verify-pane {
    flex: 1;
    display: flex;
    min-height: 0;
  }
  .vpanel {
    flex: 1;
    display: flex;
    flex-direction: column;
    padding: 8px 12px;
    gap: 5px;
    min-width: 0;
    border-right: 1px solid #1e3050;
  }
  .vpanel:last-child { border-right: none; }
  .vpanel textarea   { height: 52px; }

  /* ── canvas (shared) ── */
  .canvas-wrap {
    flex: 1;
    background: #fff;
    display: flex;
    align-items: center;
    justify-content: center;
    overflow: hidden;
    min-height: 0;
    cursor: grab;
    position: relative;
    border-radius: 2px;
  }
  .canvas-wrap:active { cursor: grabbing; }
  .canvas-wrap.dark   { filter: invert(1); }

  .canvas-inner {
    transform-origin: center center;
    display: flex;
    align-items: center;
    justify-content: center;
    width: 100%;
    height: 100%;
  }
  .canvas-inner svg { max-width: 100%; max-height: 100%; display: block; }

  .placeholder {
    color: #b0bec5;
    font-size: 13px;
    text-align: center;
    padding: 20px;
    user-select: none;
    pointer-events: none;
  }
  .placeholder.err { color: #d9534f; }

  .spinner {
    width: 28px; height: 28px;
    border: 3px solid #e0e0e0;
    border-top-color: #3a6fd8;
    border-radius: 50%;
    animation: spin .7s linear infinite;
    pointer-events: none;
  }
  @keyframes spin { to { transform: rotate(360deg); } }

  /* ── verify comparison bar ── */
  #compare-bar {
    flex-shrink: 0;
    padding: 5px 16px;
    background: #0d1422;
    border-top: 1px solid #1e3050;
    font-size: 12px;
    font-family: "Cascadia Code", "Fira Mono", monospace;
    min-height: 24px;
    display: flex;
    align-items: center;
    gap: 12px;
  }
  #compare-bar .match   { color: #3dbe6c; font-weight: 600; }
  #compare-bar .nomatch { color: #d9534f; font-weight: 600; }
  #compare-bar .canon   { color: #7aaeff; font-size: 10.5px; flex: 1;
                          overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
</style>
</head>
<body>

<header>
  <span class="title">CABILN Live Renderer <span>— pyPept</span></span>
  <button id="btn-dark"   class="hbtn" title="Toggle dark canvas">🌙 Dark</button>
  <button id="btn-verify" class="hbtn" title="SMILES vs CABILN comparison">⚖ Verify</button>
</header>

<!-- input bar (always visible) -->
<div id="input-bar">
  <div id="input-bar-row">
    <span class="seq-label" id="cabiln-label">Sequence</span>
    <textarea id="cabiln-input"
              placeholder="e.g.  fmoc-A-G-L-am&#10;fmoc-C.trt(4,1)-A-K.boc(4,1)-am&#10;fmoc-K.!1(4,1)-G-G-E.!1-am"
              spellcheck="false" autocomplete="off"></textarea>
  </div>
  <div id="cabiln-status" class="statusbar"></div>
</div>

<!-- render mode -->
<div id="main">
  <div id="render-pane">
    <div id="render-canvas" class="canvas-wrap">
      <div class="canvas-inner" id="render-inner">
        <div class="placeholder">Start typing a sequence…</div>
      </div>
    </div>
  </div>

  <!-- verify mode (hidden by default) -->
  <div id="verify-pane" style="display:none;">
    <div class="vpanel">
      <span class="seq-label">SMILES</span>
      <textarea id="smiles-input" placeholder="Paste reference SMILES here…"
                spellcheck="false" autocomplete="off"></textarea>
      <div id="smiles-status" class="statusbar"></div>
      <div id="smiles-canvas" class="canvas-wrap">
        <div class="canvas-inner" id="smiles-inner">
          <div class="placeholder">Enter a SMILES string…</div>
        </div>
      </div>
    </div>
    <div class="vpanel">
      <span class="seq-label">CABILN (generated)</span>
      <div id="verify-canvas" class="canvas-wrap">
        <div class="canvas-inner" id="verify-inner">
          <div class="placeholder">Waiting for CABILN…</div>
        </div>
      </div>
    </div>
  </div>
</div>

<!-- compare bar (verify mode only) -->
<div id="compare-bar" style="display:none;"></div>

<script>
// ─── state ────────────────────────────────────────────────────────────────────
let darkMode   = false;
let verifyMode = false;
let cabilnTimer = null;
let smilesTimer = null;
let lastCabiln  = '';
let lastSmiles  = '';

// ─── elements ─────────────────────────────────────────────────────────────────
const cabilnInput   = document.getElementById('cabiln-input');
const cabilnStatus  = document.getElementById('cabiln-status');
const renderInner   = document.getElementById('render-inner');
const renderCanvas  = document.getElementById('render-canvas');
const verifyInner   = document.getElementById('verify-inner');
const smilesInput   = document.getElementById('smiles-input');
const smilesStatus  = document.getElementById('smiles-status');
const smilesInner   = document.getElementById('smiles-inner');
const compareBar    = document.getElementById('compare-bar');
const btnDark       = document.getElementById('btn-dark');
const btnVerify     = document.getElementById('btn-verify');

// ─── dark mode ────────────────────────────────────────────────────────────────
btnDark.addEventListener('click', () => {
  darkMode = !darkMode;
  btnDark.classList.toggle('active', darkMode);
  for (const el of document.querySelectorAll('.canvas-wrap')) {
    el.classList.toggle('dark', darkMode);
  }
});

// ─── verify mode ──────────────────────────────────────────────────────────────
btnVerify.addEventListener('click', () => {
  verifyMode = !verifyMode;
  btnVerify.classList.toggle('active', verifyMode);
  document.getElementById('render-pane').style.display = verifyMode ? 'none' : '';
  document.getElementById('verify-pane').style.display = verifyMode ? '' : 'none';
  compareBar.style.display = verifyMode ? '' : 'none';
  if (verifyMode) triggerVerify();
});

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
makeZoomable(document.getElementById('verify-canvas'), verifyInner);

// ─── render helpers ───────────────────────────────────────────────────────────
function canvasSize(el) {
  return { w: Math.max(el.clientWidth || 600, 400),
           h: Math.max(el.clientHeight || 500, 300) };
}

function setInner(inner, html) {
  inner.innerHTML = html;
  // reset transform on content change
  inner.style.transform = '';
}

function showSpinner(inner) {
  setInner(inner, '<div class="spinner"></div>');
}

// ─── CABILN render ────────────────────────────────────────────────────────────
cabilnInput.addEventListener('input', () => {
  clearTimeout(cabilnTimer);
  const seq = cabilnInput.value.trim();
  if (seq === lastCabiln) return;
  if (!seq) { resetCabiln(); return; }
  showSpinner(renderInner);
  if (verifyMode) showSpinner(verifyInner);
  cabilnTimer = setTimeout(() => doRenderCabiln(seq), 300);
});

function resetCabiln() {
  lastCabiln = '';
  cabilnInput.className = '';
  cabilnStatus.textContent = '';
  cabilnStatus.className = 'statusbar';
  setInner(renderInner, '<div class="placeholder">Start typing a sequence…</div>');
  if (verifyMode) {
    setInner(verifyInner, '<div class="placeholder">Waiting for CABILN…</div>');
    compareBar.innerHTML = '';
  }
}

async function doRenderCabiln(seq) {
  lastCabiln = seq;
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res  = await fetch('/render', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: seq, width: w, height: h })
    });
    const data = await res.json();
    if (data.error) {
      const html = `<div class="placeholder err">${escHtml(data.error)}</div>`;
      setInner(renderInner, html);
      if (verifyMode) setInner(verifyInner, html);
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
      cabilnInput.className = 'err';
    } else {
      setInner(renderInner, data.svg);
      if (verifyMode) setInner(verifyInner, data.svg);
      cabilnStatus.textContent = data.info || '';
      cabilnStatus.className = 'statusbar ok';
      cabilnInput.className = 'ok';
      if (verifyMode && lastSmiles) triggerVerify();
    }
  } catch (e) {
    cabilnStatus.textContent = 'Server error — is the renderer running?';
    cabilnStatus.className = 'statusbar';
  }
}

// ─── SMILES render (verify mode) ─────────────────────────────────────────────
smilesInput.addEventListener('input', () => {
  clearTimeout(smilesTimer);
  const smi = smilesInput.value.trim();
  if (smi === lastSmiles) return;
  if (!smi) {
    lastSmiles = '';
    setInner(smilesInner, '<div class="placeholder">Enter a SMILES string…</div>');
    compareBar.innerHTML = '';
    return;
  }
  showSpinner(smilesInner);
  smilesTimer = setTimeout(() => doRenderSmiles(smi), 300);
});

async function doRenderSmiles(smi) {
  lastSmiles = smi;
  const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
  try {
    const res  = await fetch('/render_smiles', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles: smi, width: w, height: h })
    });
    const data = await res.json();
    if (data.error) {
      setInner(smilesInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      smilesStatus.textContent = data.error;
      smilesStatus.className = 'statusbar';
      smilesInput.className = 'err';
    } else {
      setInner(smilesInner, data.svg);
      smilesStatus.textContent = data.info || '';
      smilesStatus.className = 'statusbar ok';
      smilesInput.className = 'ok';
      if (lastCabiln) triggerVerify();
    }
  } catch (e) {
    smilesStatus.textContent = 'Server error';
    smilesStatus.className = 'statusbar';
  }
}

// ─── verify comparison ────────────────────────────────────────────────────────
async function triggerVerify() {
  if (!lastSmiles || !lastCabiln) return;
  try {
    const res  = await fetch('/verify', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles: lastSmiles, cabiln: lastCabiln })
    });
    const data = await res.json();
    if (data.error) {
      compareBar.innerHTML = `<span class="nomatch">Error: ${escHtml(data.error)}</span>`;
      return;
    }
    const badge = data.match
      ? '<span class="match">✓ EXACT MATCH</span>'
      : '<span class="nomatch">✗ MISMATCH</span>';
    compareBar.innerHTML =
      badge +
      `<span class="canon" title="SMILES canonical: ${escHtml(data.smiles_canonical)}">` +
      `SMILES: ${escHtml(data.smiles_canonical.slice(0, 80))}${data.smiles_canonical.length > 80 ? '…' : ''}</span>` +
      `<span class="canon" title="CABILN canonical: ${escHtml(data.cabiln_canonical)}">` +
      `CABILN: ${escHtml(data.cabiln_canonical.slice(0, 80))}${data.cabiln_canonical.length > 80 ? '…' : ''}</span>`;
  } catch (e) {
    compareBar.innerHTML = '<span class="nomatch">Verify error</span>';
  }
}

function escHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}
</script>
</body>
</html>
"""


app = FastAPI()


# ── shared drawing helper ─────────────────────────────────────────────────────

def _draw_mol(romol, width: int, height: int) -> str:
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(romol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.drawOptions().padding = 0.12
    drawer.DrawMolecule(romol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _canon(romol) -> str:
    from rdkit.Chem import MolToSmiles
    return MolToSmiles(romol, canonical=True)


# ── request models ────────────────────────────────────────────────────────────

class _CabilnReq(BaseModel):
    cabiln: str
    width:  int = 960
    height: int = 680

class _SmilesReq(BaseModel):
    smiles: str
    width:  int = 960
    height: int = 680

class _VerifyReq(BaseModel):
    smiles: str
    cabiln: str


# ── routes ────────────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def index():
    return _HTML


@app.post("/render")
async def render(req: _CabilnReq):
    try:
        from rdkit.Chem.Descriptors import ExactMolWt
        from pyPept.molecule import Molecule
        from pyPept.sequence import Sequence

        seq   = Sequence(req.cabiln)
        mol   = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        if romol is None:
            return JSONResponse({"error": "Assembly produced no molecule"}, status_code=400)

        w, h = max(400, req.width), max(300, req.height)
        svg  = _draw_mol(romol, w, h)
        return {"svg": svg, "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}"}

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


@app.post("/render_smiles")
async def render_smiles(req: _SmilesReq):
    try:
        from rdkit import Chem
        from rdkit.Chem.Descriptors import ExactMolWt

        romol = Chem.MolFromSmiles(req.smiles)
        if romol is None:
            return JSONResponse({"error": "Invalid SMILES"}, status_code=400)

        w, h = max(400, req.width), max(300, req.height)
        svg  = _draw_mol(romol, w, h)
        return {"svg": svg, "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}"}

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


@app.post("/verify")
async def verify(req: _VerifyReq):
    try:
        from rdkit import Chem
        from pyPept.molecule import Molecule
        from pyPept.sequence import Sequence

        smiles_mol = Chem.MolFromSmiles(req.smiles)
        if smiles_mol is None:
            return JSONResponse({"error": "Invalid SMILES"}, status_code=400)

        seq   = Sequence(req.cabiln)
        mol   = Molecule(seq)
        cabiln_mol = mol.get_molecule(fmt='ROMol')
        if cabiln_mol is None:
            return JSONResponse({"error": "CABILN assembly produced no molecule"}, status_code=400)

        smi_canon    = _canon(smiles_mol)
        cabiln_canon = _canon(cabiln_mol)
        return {
            "match":           smi_canon == cabiln_canon,
            "smiles_canonical": smi_canon,
            "cabiln_canonical": cabiln_canon,
        }

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


# ── launch ────────────────────────────────────────────────────────────────────

def _open_browser():
    time.sleep(0.9)
    webbrowser.open("http://localhost:8732")


if __name__ == "__main__":
    print("CABILN Live Renderer")
    print("  Local     ->  http://localhost:8732")
    print("  Tailscale ->  http://100.119.0.78:8732")
    threading.Thread(target=_open_browser, daemon=True).start()
    uvicorn.run(app, host="0.0.0.0", port=8732, log_level="warning")
