#!/usr/bin/env python3
"""
Live CABILN -> 2D structure renderer.

    pip install fastapi uvicorn
    python tools/live_renderer.py

Left pane: type a CABILN sequence.
Right pane: molecule updates ~300 ms after you stop typing.
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
  html, body { height: 100%; }
  body {
    font-family: system-ui, -apple-system, sans-serif;
    background: #12192b;
    color: #d0dce8;
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
  }
  header {
    flex-shrink: 0;
    padding: 10px 20px;
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
  .panes {
    display: flex;
    flex: 1;
    min-height: 0;
  }
  .pane {
    flex: 1;
    display: flex;
    flex-direction: column;
    padding: 14px 16px;
    gap: 8px;
    min-width: 0;
  }
  .pane + .pane { border-left: 1px solid #1e3050; }
  .pane-label {
    font-size: 10px;
    text-transform: uppercase;
    letter-spacing: .12em;
    color: #3a5580;
    flex-shrink: 0;
  }
  textarea {
    flex: 1;
    background: #0a1018;
    color: #c8daf0;
    border: 1px solid #1e3050;
    border-radius: 6px;
    padding: 12px 14px;
    font-family: "Cascadia Code", "Fira Mono", "Courier New", monospace;
    font-size: 14px;
    line-height: 1.7;
    resize: none;
    outline: none;
    transition: border-color .15s;
  }
  textarea:focus  { border-color: #3a6fd8; }
  textarea.err    { border-color: #d9534f; }
  textarea.ok     { border-color: #28a745; }
  #statusbar {
    flex-shrink: 0;
    font-size: 11.5px;
    font-family: "Cascadia Code", "Fira Mono", monospace;
    min-height: 18px;
    color: #d9534f;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  #statusbar.ok { color: #3dbe6c; }
  #svg-wrap {
    flex: 1;
    background: #ffffff;
    border-radius: 6px;
    border: 1px solid #1e3050;
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
<div class="panes">

  <div class="pane">
    <div class="pane-label">Sequence</div>
    <textarea id="input"
              placeholder="e.g.  fmoc-A-G-L-am
fmoc-C.trt(4,1)-A-K.boc(4,1)-am
fmoc-K.!1(4,1)-G-G-E.!1-am"
              spellcheck="false" autocomplete="off"></textarea>
    <div id="statusbar"></div>
  </div>

  <div class="pane">
    <div class="pane-label">Structure</div>
    <div id="svg-wrap">
      <div id="placeholder">Start typing a sequence&hellip;</div>
      <div id="spinner"></div>
    </div>
  </div>

</div>
<script>
  const input    = document.getElementById('input');
  const status   = document.getElementById('statusbar');
  const svgWrap  = document.getElementById('svg-wrap');
  const spinner  = document.getElementById('spinner');
  let timer = null;
  let lastSeq = '';

  input.addEventListener('input', () => {
    clearTimeout(timer);
    const seq = input.value.trim();
    if (seq === lastSeq) return;
    if (!seq) { reset(); return; }
    showSpinner();
    timer = setTimeout(() => doRender(seq), 280);
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
    try {
      const res  = await fetch('/render', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({cabiln: seq})
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

        seq  = Sequence(req.cabiln)
        mol  = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        if romol is None:
            return JSONResponse({"error": "Assembly produced no molecule"}, status_code=400)

        rdDepictor.Compute2DCoords(romol)
        drawer = rdMolDraw2D.MolDraw2DSVG(680, 480)
        drawer.drawOptions().addStereoAnnotation = True
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
    print("CABILN Live Renderer  →  http://localhost:8732")
    threading.Thread(target=_open_browser, daemon=True).start()
    uvicorn.run(app, host="127.0.0.1", port=8732, log_level="warning")
