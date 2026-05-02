#!/usr/bin/env python3
"""
Live CABILN -> 2D structure renderer.

    pip install fastapi uvicorn
    python tools/live_renderer.py

Input bar at top (Excel-style, vertically resizable).
Structure canvas fills the remaining space.
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

_HTML = """\
<!DOCTYPE html>
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
    padding: 8px 18px;
    background: #0d1422;
    border-bottom: 1px solid #1e3050;
    font-size: 13px;
    letter-spacing: .08em;
    color: #7aaeff;
    display: flex;
    align-items: center;
    gap: 8px;
  }
  header span { color: #3a5580; font-weight: 400; }

  /* ── input bar (Excel formula-bar style) ── */
  #input-bar {
    flex-shrink: 0;
    display: flex;
    flex-direction: column;
    gap: 4px;
    padding: 8px 16px 6px;
    background: #0d1422;
    border-bottom: 2px solid #1e3050;
  }
  #input-bar-top {
    display: flex;
    align-items: center;
    gap: 10px;
  }
  #seq-label {
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
    padding: 7px 12px;
    font-family: "Cascadia Code", "Fira Mono", "Courier New", monospace;
    font-size: 13.5px;
    line-height: 1.6;
    resize: vertical;
    min-height: 34px;
    height: 68px;
    max-height: 240px;
    outline: none;
    transition: border-color .15s;
  }
  textarea:focus { border-color: #3a6fd8; }
  textarea.err   { border-color: #d9534f; }
  textarea.ok    { border-color: #28a745; }
  #statusbar {
    font-size: 11.5px;
    font-family: "Cascadia Code", "Fira Mono", monospace;
    min-height: 16px;
    color: #d9534f;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    padding-left: 2px;
  }
  #statusbar.ok { color: #3dbe6c; }

  /* ── structure canvas ── */
  #svg-wrap {
    flex: 1;
    background: #ffffff;
    display: flex;
    align-items: center;
    justify-content: center;
    overflow: hidden;
    min-height: 0;
  }
  #svg-wrap svg {
    max-width: 100%;
    max-height: 100%;
    display: block;
  }
  #placeholder {
    color: #b0bec5;
    font-size: 13px;
    text-align: center;
    padding: 20px;
    user-select: none;
  }
  #spinner {
    display: none;
    width: 28px; height: 28px;
    border: 3px solid #e0e0e0;
    border-top-color: #3a6fd8;
    border-radius: 50%;
    animation: spin .7s linear infinite;
  }
  @keyframes spin { to { transform: rotate(360deg); } }
</style>
</head>
<body>

<header>
  CABILN Live Renderer
  <span>— pyPept</span>
</header>

<div id="input-bar">
  <div id="input-bar-top">
    <span id="seq-label">Sequence</span>
    <textarea id="input"
              placeholder="e.g.  fmoc-A-G-L-am&#10;fmoc-C.trt(4,1)-A-K.boc(4,1)-am&#10;fmoc-K.!1(4,1)-G-G-E.!1-am"
              spellcheck="false" autocomplete="off"></textarea>
  </div>
  <div id="statusbar"></div>
</div>

<div id="svg-wrap">
  <div id="placeholder">Start typing a sequence&hellip;</div>
  <div id="spinner"></div>
</div>

<script>
  const input   = document.getElementById('input');
  const status  = document.getElementById('statusbar');
  const svgWrap = document.getElementById('svg-wrap');
  const spinner = document.getElementById('spinner');
  let timer = null;
  let lastSeq = '';

  input.addEventListener('input', () => {
    clearTimeout(timer);
    const seq = input.value.trim();
    if (seq === lastSeq) return;
    if (!seq) { reset(); return; }
    showSpinner();
    timer = setTimeout(() => doRender(seq), 300);
  });

  function reset() {
    lastSeq = '';
    input.className = '';
    status.textContent = '';
    status.className = '';
    svgWrap.innerHTML = '<div id="placeholder">Start typing a sequence…</div>';
  }

  function showSpinner() {
    spinner.style.display = 'block';
    svgWrap.innerHTML = '';
    svgWrap.appendChild(spinner);
  }

  async function doRender(seq) {
    lastSeq = seq;
    // pass canvas size so server can render at native resolution
    const w = svgWrap.clientWidth  || 960;
    const h = svgWrap.clientHeight || 680;
    try {
      const res  = await fetch('/render', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({cabiln: seq, width: w, height: h})
      });
      const data = await res.json();
      if (data.error) {
        svgWrap.innerHTML = '<div id="placeholder" style="color:#d9534f">' +
          escHtml(data.error) + '</div>';
        status.textContent = data.error;
        status.className   = '';
        input.className    = 'err';
      } else {
        svgWrap.innerHTML  = data.svg;
        status.textContent = data.info || '';
        status.className   = 'ok';
        input.className    = 'ok';
      }
    } catch (e) {
      status.textContent = 'Server error — is the renderer running?';
      status.className   = '';
    }
  }

  function escHtml(s) {
    return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
  }
</script>
</body>
</html>
"""


app = FastAPI()


class _Req(BaseModel):
    cabiln: str
    width: int = 960
    height: int = 680


@app.get("/", response_class=HTMLResponse)
async def index():
    return _HTML


@app.post("/render")
async def render(req: _Req):
    try:
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Descriptors import ExactMolWt
        from rdkit.Chem.Draw import rdMolDraw2D

        from pyPept.molecule import Molecule
        from pyPept.sequence import Sequence

        seq   = Sequence(req.cabiln)
        mol   = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        if romol is None:
            return JSONResponse({"error": "Assembly produced no molecule"}, status_code=400)

        # CoordGen produces ChemDraw-quality layouts for macrocycles and
        # complex ring systems — avoids the bond-overlap issues that the
        # default RDKit depiction algorithm produces on larger peptides.
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(romol)

        w = max(400, req.width)
        h = max(300, req.height)
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().padding = 0.12
        drawer.DrawMolecule(romol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        n_atoms = romol.GetNumAtoms()
        mw      = ExactMolWt(romol)
        return {"svg": svg, "info": f"{n_atoms} heavy atoms · MW {mw:.2f}"}

    except Exception as exc:
        msg = str(exc).split("\n")[0]
        return JSONResponse({"error": msg}, status_code=400)


def _open_browser():
    time.sleep(0.9)
    webbrowser.open("http://localhost:8732")


if __name__ == "__main__":
    print("CABILN Live Renderer")
    print("  Local     ->  http://localhost:8732")
    print("  Tailscale ->  http://100.119.0.78:8732")
    threading.Thread(target=_open_browser, daemon=True).start()
    uvicorn.run(app, host="0.0.0.0", port=8732, log_level="warning")
