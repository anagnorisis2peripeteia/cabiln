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
  Library overlay   — searchable monomer browser; click abbreviation to insert.
  Exports           — PNG (client-side) and MOL (server-side) download buttons.
  Register page     — /register to add new monomers from plain SMILES.
"""

__credits__ = ["Cameron Beesley"]
__license__ = "MIT"

import threading
import time
import uuid
import webbrowser

_SERVER_ID = uuid.uuid4().hex[:8]

# ── in-memory SDF cache ──────────────────────────────────────────────────────
import pathlib as _pathlib

_SDF_PATH = _pathlib.Path(__file__).resolve().parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
_sdf_cache: dict = {'mols': None, 'by_abbr': None, 'mtime': 0.0}
_sdf_lock = threading.Lock()


def _load_sdf():
    """Return (mols_list, by_abbr_dict) from a cached SDF load.  Reloads only when the file mtime changes."""
    from rdkit import Chem
    mtime = _SDF_PATH.stat().st_mtime
    with _sdf_lock:
        if _sdf_cache['mols'] is not None and _sdf_cache['mtime'] == mtime:
            return _sdf_cache['mols'], _sdf_cache['by_abbr']
        suppl = Chem.SDMolSupplier(str(_SDF_PATH), removeHs=False)
        mols, by_abbr = [], {}
        for mol in suppl:
            if mol is None:
                continue
            mols.append(mol)
            abbr = mol.GetPropsAsDict().get('m_abbr', '')
            if abbr:
                by_abbr[abbr] = mol
        _sdf_cache['mols'] = mols
        _sdf_cache['by_abbr'] = by_abbr
        _sdf_cache['mtime'] = mtime
        return mols, by_abbr


def _invalidate_sdf():
    with _sdf_lock:
        _sdf_cache['mols'] = None
        _sdf_cache['by_abbr'] = None
        _sdf_cache['mtime'] = 0.0

# ──────────────────────────────────────────────────────────────────────────────

try:
    import uvicorn
    from fastapi import FastAPI
    from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
    from pydantic import BaseModel
except ImportError:
    raise SystemExit("pip install fastapi uvicorn  (then re-run)")

EXAMPLES = [
    {
        "category": "GLP-1 / Incretin Agonists",
        "items": [
            {
                "name": "Retatrutide",
                "description": "GIP/GLP-1/glucagon triple agonist · 39 AA · C20 lipid conjugate",
                "cabiln": "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am",
            },
            {
                "name": "Semaglutide-like",
                "description": "GLP-1 scaffold · γGlu–AEEA–C20 lipid linker on K8",
                "cabiln": "His-Aib-Glu-Gly-Thr-Phe-Thr-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am",
            },
        ],
    },
    {
        "category": "Bicyclic Peptides (Bicycle Therapeutics)",
        "items": [
            {
                "name": "TBMB Bicycle — minimal",
                "description": "Three-Cys scaffold cyclised with TBMB · Bicycle Therapeutics platform",
                "cabiln": "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3",
            },
            {
                "name": "TBMB Bicycle — asymmetric loops",
                "description": "Variable loop lengths for epitope shape mimicry",
                "cabiln": "ac-C.!1(4,4)-G-A-K-C.!2(4,5)-E-L-F-C.!3(4,6)-am%TBMB.!1.!2.!3",
            },
        ],
    },
    {
        "category": "Cyclic Peptides",
        "items": [
            {
                "name": "Head-to-tail cyclic with lipid",
                "description": "N→C cyclised backbone with C20 lipid branch on K",
                "cabiln": "!1-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1",
            },
            {
                "name": "Head-to-tail cyclic",
                "description": "Simple N→C cyclised hexapeptide",
                "cabiln": "!1-A-K-G-E-L-F-!1",
            },
        ],
    },
    {
        "category": "Lipid Conjugates",
        "items": [
            {
                "name": "Dual-lipid backbone",
                "description": "Two independent γGlu–AEEA–C20 arms on a 3-residue backbone",
                "cabiln": "ac-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am",
            },
            {
                "name": "AEEA PEG spacer",
                "description": "Tandem AEEA mini-PEG linker — common in half-life extension",
                "cabiln": "ac-K.[AEEA(4,2).AEEA(1,2)]-G-am",
            },
        ],
    },
    {
        "category": "Branched Peptides",
        "items": [
            {
                "name": "K-branched dipeptide",
                "description": "Lysine with G–A pendant arm via R4→R1 bond",
                "cabiln": "ac-A-K.[G(4,1).A(2,1).am(2,1)]-G-am",
            },
            {
                "name": "Triple-branched backbone",
                "description": "Three positional branches — each K carries a distinct arm",
                "cabiln": "ac-K.[G(4,1).am(2,1)]-A-K.[A(4,1).am(2,1)]-A-K.[G(4,1).G(2,1).am(2,1)]-am",
            },
        ],
    },
]

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
    flex-wrap: wrap;
  }
  header .title { flex: 1; min-width: 120px; }
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
  .hbtn:hover        { background: #263550; color: #c8daf0; }
  .hbtn.active       { background: #2a4a80; border-color: #4a6faa; color: #d8e8ff; }
  .hbtn:disabled     { opacity: .35; cursor: default; pointer-events: none; }
  .hbtn.green        { border-color: #1e5030; color: #5dba7a; }
  .hbtn.green:hover  { background: #1a3828; color: #7dd098; }

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
  #main { flex: 1; display: flex; min-height: 0; position: relative; }

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
    border-left: 1px solid #1e3050;
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

  /* ── monomer library sidebar ── */
  #lib-panel {
    width: 420px;
    max-width: 50%;
    background: #0d1422;
    border-right: 1px solid #1e3050;
    display: none;
    flex-direction: column;
    flex-shrink: 0;
    overflow: hidden;
  }
  #lib-panel.open { display: flex; }

  /* ── example peptide sidebar (right) ── */
  #examples-panel {
    width: 340px;
    max-width: 50%;
    background: #0d1422;
    border-left: 1px solid #1e3050;
    display: none;
    flex-direction: column;
    flex-shrink: 0;
    overflow: hidden;
    order: 99;
  }
  #examples-panel.open { display: flex; }

  #examples-header {
    flex-shrink: 0;
    padding: 10px 14px 8px;
    background: #0a1018;
    border-bottom: 1px solid #1e3050;
    display: flex;
    align-items: center;
    gap: 8px;
  }
  #examples-header span { flex: 1; font-size: 12px; color: #7aaeff; letter-spacing: .08em; }
  #examples-close {
    background: none;
    border: none;
    color: #4a6a9a;
    font-size: 16px;
    cursor: pointer;
    padding: 2px 4px;
    line-height: 1;
  }
  #examples-close:hover { color: #d8e8ff; }

  #examples-list {
    flex: 1;
    overflow-y: auto;
    padding: 0 0 8px;
  }
  #examples-list::-webkit-scrollbar { width: 6px; }
  #examples-list::-webkit-scrollbar-track { background: #0a1018; }
  #examples-list::-webkit-scrollbar-thumb { background: #1e3050; border-radius: 3px; }

  .example-cat {
    padding: 8px 14px 4px;
    font-size: 9.5px;
    text-transform: uppercase;
    letter-spacing: .1em;
    color: #3a5070;
    border-bottom: 1px solid #0d1830;
    margin-top: 6px;
  }
  .example-cat:first-child { margin-top: 0; }
  .example-row {
    padding: 8px 14px 7px;
    border-bottom: 1px solid #0d1830;
    cursor: pointer;
    transition: background .1s;
  }
  .example-row:hover { background: #131f35; }
  .example-name {
    font-size: 12px;
    font-weight: 600;
    color: #8ab8f0;
    margin-bottom: 2px;
  }
  .example-desc {
    font-size: 10.5px;
    color: #4a6a8a;
    margin-bottom: 5px;
    line-height: 1.4;
  }
  .example-seq {
    font-family: "Cascadia Code", "Fira Mono", monospace;
    font-size: 9px;
    color: #2e4560;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }

  #lib-header {
    flex-shrink: 0;
    padding: 10px 14px 8px;
    background: #0a1018;
    border-bottom: 1px solid #1e3050;
    display: flex;
    align-items: center;
    gap: 8px;
  }
  #lib-header span { flex: 1; font-size: 12px; color: #7aaeff; letter-spacing: .08em; }
  #lib-search {
    flex: 2;
    background: #12192b;
    border: 1px solid #1e3050;
    border-radius: 4px;
    color: #c8daf0;
    padding: 4px 9px;
    font-size: 12px;
    font-family: "Cascadia Code", "Fira Mono", monospace;
    outline: none;
    transition: border-color .15s;
  }
  #lib-search:focus { border-color: #3a6fd8; }
  #lib-close {
    background: none;
    border: none;
    color: #4a6a9a;
    font-size: 16px;
    cursor: pointer;
    padding: 2px 4px;
    line-height: 1;
  }
  #lib-close:hover { color: #d8e8ff; }

  #lib-count {
    flex-shrink: 0;
    padding: 4px 14px;
    font-size: 10.5px;
    color: #3a5580;
    border-bottom: 1px solid #1e3050;
  }

  #lib-list {
    flex: 1;
    overflow-y: auto;
    padding: 4px 0;
  }
  #lib-list::-webkit-scrollbar { width: 6px; }
  #lib-list::-webkit-scrollbar-track { background: #0a1018; }
  #lib-list::-webkit-scrollbar-thumb { background: #1e3050; border-radius: 3px; }

  .lib-row {
    display: grid;
    grid-template-columns: 64px 1fr auto;
    align-items: start;
    gap: 6px;
    padding: 6px 14px;
    border-bottom: 1px solid #0d1830;
    cursor: pointer;
    transition: background .1s;
  }
  .lib-row:hover { background: #131f35; }
  .lib-row.hidden { display: none; }
  .lib-abbr {
    font-family: "Cascadia Code", "Fira Mono", monospace;
    font-size: 12px;
    font-weight: 700;
    color: #7aaeff;
    word-break: break-all;
  }
  .lib-info { min-width: 0; }
  .lib-name {
    font-size: 11.5px;
    color: #b0c8e0;
    margin-bottom: 2px;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  .lib-meta {
    font-size: 10px;
    color: #3a5580;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  .lib-badge {
    font-size: 9.5px;
    padding: 1px 5px;
    border-radius: 3px;
    background: #1a2848;
    color: #4a7ab0;
    white-space: nowrap;
    align-self: start;
    margin-top: 1px;
  }
  .lib-badge.aa      { background: #1a3028; color: #4ab070; }
  .lib-badge.cap     { background: #281a1a; color: #b07040; }
  .lib-badge.protect { background: #28281a; color: #a8a040; }

  /* ── monomer preview tooltip ── */
  #lib-preview {
    position: fixed;
    z-index: 200;
    background: #0a1018;
    border: 1px solid #2a4070;
    border-radius: 6px;
    padding: 0;
    box-shadow: 0 8px 32px rgba(0,0,0,.6);
    display: none;
    pointer-events: none;
    width: 480px;
    max-width: 680px;
    overflow: hidden;
  }
  #lib-preview.has-reagent {
    width: 640px;
  }
  #lib-preview .prev-row {
    display: flex;
  }
  #lib-preview .prev-pane {
    flex: 1;
    height: 200px;
    position: relative;
    overflow: hidden;
  }
  #lib-preview .prev-pane + .prev-pane {
    border-left: 1px solid #1e3050;
  }
  #lib-preview .prev-label {
    position: absolute;
    top: 4px; left: 6px;
    font-size: 9px;
    color: #5a7aa0;
    letter-spacing: .06em;
    text-transform: uppercase;
    z-index: 1;
  }
  #lib-preview .prev-pane svg { width: 100%; height: 100%; display: block; }
  #lib-preview.dark .prev-pane svg { filter: invert(1); }
  #lib-preview .prev-meta {
    padding: 4px 8px;
    font-size: 10px;
    color: #8ab4e8;
    border-top: 1px solid #1e3050;
    background: #0d1422;
  }
  #lib-preview .prev-warn {
    color: #e8a84a;
  }
  #lib-preview .prev-rxn {
    position: absolute;
    bottom: 4px; left: 6px;
    font-size: 9px;
    color: #6a9a6a;
    z-index: 1;
  }

  /* ── build panel ── */
  #build-panel {
    display: none;
    flex-shrink: 0;
    background: #0d1422;
    border-top: 2px solid #2a4070;
    padding: 8px 14px;
  }
  #build-panel.open { display: block; }
  #build-panel .build-header {
    display: flex;
    align-items: center;
    gap: 10px;
    margin-bottom: 6px;
  }
  #build-panel .build-header span {
    font-size: 11px;
    color: #7aaeff;
    letter-spacing: .08em;
    text-transform: uppercase;
  }
  #build-panel .build-header .build-hint {
    font-size: 10px;
    color: #3a5580;
    text-transform: none;
    letter-spacing: normal;
    flex: 1;
  }
  .build-boxes {
    display: flex;
    gap: 12px;
    align-items: stretch;
    min-height: 180px;
  }
  .build-box {
    flex: 1;
    background: #0a1018;
    border: 1px solid #1e3050;
    border-radius: 6px;
    display: flex;
    flex-direction: column;
    min-width: 0;
    position: relative;
  }
  .build-box .box-label {
    font-size: 9px;
    text-transform: uppercase;
    letter-spacing: .1em;
    color: #3a5580;
    padding: 5px 10px 0;
    display: flex;
    align-items: center;
    justify-content: space-between;
  }
  .box-clear-btn {
    font-size: 9px;
    background: none;
    border: 1px solid #2a4060;
    border-radius: 3px;
    color: #3a6090;
    cursor: pointer;
    padding: 0 4px;
    line-height: 14px;
    letter-spacing: normal;
    text-transform: none;
  }
  .box-clear-btn:hover { color: #d9534f; border-color: #d9534f; }
  .build-box .box-abbr {
    font-family: "Cascadia Code", "Fira Mono", monospace;
    font-size: 14px;
    font-weight: 700;
    color: #7aaeff;
    padding: 2px 10px;
  }
  .build-box .box-svg {
    flex: 1;
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 100px;
    overflow: hidden;
  }
  .build-box .box-svg svg { max-width: 100%; max-height: 100%; display: block; }
  .build-box.dark .box-svg svg { filter: invert(1); }
  .build-box .box-placeholder {
    color: #3a5580;
    font-size: 11px;
    text-align: center;
    padding: 12px;
  }
  .rgroup-buttons {
    display: flex;
    flex-wrap: wrap;
    gap: 4px;
    padding: 6px 10px 8px;
    border-top: 1px solid #1e3050;
  }
  .rgroup-btn {
    font-family: "Cascadia Code", "Fira Mono", monospace;
    font-size: 10px;
    padding: 2px 8px;
    border-radius: 4px;
    border: 1px solid #2a4070;
    background: #1a2640;
    color: #8ab4e8;
    cursor: pointer;
    transition: all .15s;
  }
  .rgroup-btn:hover { background: #263550; color: #c8daf0; border-color: #4a6faa; }
  .rgroup-btn.selected { background: #2a5a90; border-color: #5a9ae0; color: #fff; box-shadow: 0 0 6px rgba(90,154,224,.4); }
  .rgroup-btn.used { opacity: .25; cursor: default; pointer-events: none; text-decoration: line-through; }
  .build-arrow {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    min-width: 60px;
    gap: 6px;
  }
  .build-arrow .arrow-icon {
    font-size: 24px;
    color: #3a5580;
  }
  .build-connect-btn {
    font-size: 10px;
    padding: 4px 12px;
    border-radius: 4px;
    border: 1px solid #1e5030;
    background: #1a3828;
    color: #5dba7a;
    cursor: pointer;
    transition: all .15s;
  }
  .build-connect-btn:hover { background: #264a38; color: #7dd098; }
  .build-connect-btn:disabled { opacity: .3; cursor: default; pointer-events: none; }
  .build-status {
    font-size: 10px;
    color: #3a5580;
    text-align: center;
    min-height: 14px;
  }
  .build-status.valid { color: #3dbe6c; }
  .build-status.invalid { color: #d9534f; }

  /* ── residue chips ── */
  #residue-chips {
    display: flex;
    flex-wrap: wrap;
    gap: 3px;
    padding: 2px 0;
    min-height: 0;
  }
  .res-chip {
    font-family: "Cascadia Code", "Fira Mono", monospace;
    font-size: 11px;
    padding: 1px 7px;
    border-radius: 3px;
    cursor: pointer;
    transition: opacity .15s, border-color .15s;
    border: 1px solid transparent;
    color: #d0dce8;
  }
  .res-chip:hover, .res-chip.hover {
    border-color: rgba(255,255,255,.5);
  }
  #residue-chips.dimmed .res-chip:not(.hover) { opacity: .35; }

  /* ── SVG residue highlighting ── */
  #render-inner svg.has-highlight path:not(.res-hl),
  #render-inner svg.has-highlight text:not(.res-hl) {
    opacity: .12;
    transition: opacity .15s;
  }
  #render-inner svg .res-hl { transition: opacity .15s; }
</style>
</head>
<body>

<header>
  <span class="title">CABILN Live Renderer <span>— pyPept</span></span>
  <button id="btn-lib"      class="hbtn green"  title="Browse monomer library">📋 Library</button>
  <button id="btn-examples" class="hbtn green"  title="Browse example peptides">🧬 Examples</button>
  <button id="btn-hl"     class="hbtn active" title="Toggle residue highlighting">🔗 Highlight</button>
  <button id="btn-dark"   class="hbtn"        title="Toggle dark canvas">🌙 Dark</button>
  <button id="btn-verify" class="hbtn"        title="SMILES vs CABILN comparison">⚖ Verify</button>
  <button id="btn-build"  class="hbtn green" title="Visual peptide builder">🔧 Build</button>
  <button id="btn-to-bracket" class="hbtn" title="Convert to bracket notation">→[ ]</button>
  <button id="btn-to-branch"  class="hbtn" title="Convert to branch (%) notation">→%</button>
  <button id="btn-reroll" class="hbtn"        title="New 2D layout seed" disabled>⟳ Reroll</button>
  <button id="btn-png"    class="hbtn"        title="Download PNG" disabled>⬇ PNG</button>
  <button id="btn-mol"    class="hbtn"        title="Download MOL file" disabled>⬇ MOL</button>
  <a href="/register" target="_blank" style="text-decoration:none;">
    <button class="hbtn green" title="Register a new monomer">＋ Register</button>
  </a>
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
  <div id="residue-chips"></div>
</div>

<!-- build panel -->
<div id="build-panel">
  <div class="build-header">
    <span>Build by Connection</span>
    <span class="build-hint" id="build-hint">Select a chip on the sequence, then right-click a monomer in the library</span>
    <button id="build-close" style="background:none;border:none;color:#4a6a9a;font-size:16px;cursor:pointer;padding:2px 4px;">✕</button>
  </div>
  <div class="build-boxes">
    <div class="build-box" id="build-left">
      <div class="box-label"><span>Current residue</span><button class="box-clear-btn" id="build-left-change" title="Clear selection">Clear</button></div>
      <div class="box-abbr" id="build-left-abbr">—</div>
      <div class="box-svg" id="build-left-svg">
        <div class="box-placeholder">Click a chip above</div>
      </div>
      <div class="rgroup-buttons" id="build-left-rgroups"></div>
    </div>
    <div class="build-arrow">
      <div class="arrow-icon">⟷</div>
      <button class="build-connect-btn" id="build-connect" disabled>Connect</button>
      <div class="build-status" id="build-status"></div>
    </div>
    <div class="build-box" id="build-right">
      <div class="box-label"><span id="build-right-label">New monomer</span><button class="box-clear-btn" id="build-right-change" title="Clear selection">Clear</button></div>
      <div class="box-abbr" id="build-right-abbr">—</div>
      <div class="box-svg" id="build-right-svg">
        <div class="box-placeholder">Right-click from library</div>
      </div>
      <div class="rgroup-buttons" id="build-right-rgroups"></div>
    </div>
  </div>
</div>

<!-- main area -->
<div id="main">

  <!-- library sidebar (persistent until dismissed) -->
  <div id="lib-panel">
    <div id="lib-header">
      <span>Monomer Library</span>
      <input id="lib-search" type="text" placeholder="Search abbr / name / type…" autocomplete="off">
      <button id="btn-rxn-filter" class="hbtn" title="Filter to monomers reactive with selected residue" disabled>⚗ Filter</button>
      <button id="lib-close" title="Close">✕</button>
    </div>
    <div id="lib-count"></div>
    <div id="lib-list"><div class="placeholder">Loading…</div></div>
  </div>
  <div id="lib-preview"></div>

  <!-- render mode -->
  <div id="render-pane">
    <div id="render-canvas" class="canvas-wrap">
      <div class="canvas-inner" id="render-inner">
        <div class="placeholder">Start typing a sequence…</div>
      </div>
    </div>
  </div>

  <!-- example peptide sidebar (right, persistent) -->
  <div id="examples-panel">
    <div id="examples-header">
      <span>Example Peptides</span>
      <button id="examples-close" title="Close">✕</button>
    </div>
    <div id="examples-list"><div class="placeholder">Loading…</div></div>
  </div>

  <!-- verify reference pane (hidden by default, shown beside render-pane) -->
  <div id="verify-pane" style="display:none;">
    <div class="vpanel" style="border-right:none;">
      <div style="display:flex;gap:6px;align-items:center;">
        <span class="seq-label" style="flex:none;">Reference</span>
        <label class="hbtn" style="font-size:11px;padding:2px 8px;cursor:pointer;margin:0;">
          📂 .mol
          <input type="file" id="mol-upload" accept=".mol,.sdf" style="display:none;">
        </label>
      </div>
      <textarea id="smiles-input" placeholder="Paste SMILES, BILN, or HELM here…"
                spellcheck="false" autocomplete="off"></textarea>
      <div id="smiles-status" class="statusbar"></div>
      <div id="smiles-canvas" class="canvas-wrap">
        <div class="canvas-inner" id="smiles-inner">
          <div class="placeholder">Paste SMILES or upload .mol to compare…</div>
        </div>
      </div>
    </div>
  </div>
</div>

<!-- compare bar (verify mode only) -->
<div id="compare-bar" style="display:none;"></div>

<script>
// ─── state ────────────────────────────────────────────────────────────────────
let darkMode   = true;
let verifyMode = false;
let hlEnabled  = true;
let libLoaded  = false;
let allMonomers = [];
let cabilnTimer = null;
let smilesTimer = null;
let lastCabiln  = '';
let lastSmiles  = '';
let lastSvg     = '';
let lastMolBlock = '';
let residueMap  = {};
let atomToRes   = {};
let residueList = [];
let previewCache = {};
let previewTimer = null;
let buildMode    = false;
let buildLeft    = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
let buildRight   = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
let buildLeftRIdx = null;
let buildRightRIdx = null;
let rerollSeed   = 0;
let reactionPairs = null;  // lazy-loaded list of [ct_a, ct_b] pairs
let rxnFilterActive = false;

const RES_COLORS = [
  '#2a5080','#2a8050','#802a50','#806a2a','#502a80',
  '#2a6080','#80502a','#2a8070','#6a2a80','#80802a',
  '#3a6080','#3a8060','#603a50','#706a3a','#403a70',
];

// ─── elements ─────────────────────────────────────────────────────────────────
const cabilnInput   = document.getElementById('cabiln-input');
const cabilnStatus  = document.getElementById('cabiln-status');
const renderInner   = document.getElementById('render-inner');
const renderCanvas  = document.getElementById('render-canvas');
const molUpload     = document.getElementById('mol-upload');
const smilesInput   = document.getElementById('smiles-input');
const smilesStatus  = document.getElementById('smiles-status');
const smilesInner   = document.getElementById('smiles-inner');
const compareBar    = document.getElementById('compare-bar');
const btnDark       = document.getElementById('btn-dark');
const btnVerify     = document.getElementById('btn-verify');
const btnLib        = document.getElementById('btn-lib');
const btnHl         = document.getElementById('btn-hl');
const btnPng        = document.getElementById('btn-png');
const btnMol        = document.getElementById('btn-mol');
const btnToBracket  = document.getElementById('btn-to-bracket');
const btnToBranch   = document.getElementById('btn-to-branch');
const libPanel      = document.getElementById('lib-panel');
const libSearch     = document.getElementById('lib-search');
const libClose      = document.getElementById('lib-close');
const libList       = document.getElementById('lib-list');
const libCount      = document.getElementById('lib-count');
const libPreview    = document.getElementById('lib-preview');
const btnExamples    = document.getElementById('btn-examples');
const examplesPanel  = document.getElementById('examples-panel');
const examplesClose  = document.getElementById('examples-close');
const examplesList   = document.getElementById('examples-list');
const resChips      = document.getElementById('residue-chips');
const buildPanel    = document.getElementById('build-panel');
const btnBuild      = document.getElementById('btn-build');
const buildClose    = document.getElementById('build-close');
const buildConnect  = document.getElementById('build-connect');
const buildStatus   = document.getElementById('build-status');
const buildHint     = document.getElementById('build-hint');
const buildLeftAbbr = document.getElementById('build-left-abbr');
const buildLeftSvg  = document.getElementById('build-left-svg');
const buildLeftRg   = document.getElementById('build-left-rgroups');
const buildRightAbbr= document.getElementById('build-right-abbr');
const buildRightSvg = document.getElementById('build-right-svg');
const buildRightRg  = document.getElementById('build-right-rgroups');
const btnReroll     = document.getElementById('btn-reroll');
const btnRxnFilter  = document.getElementById('btn-rxn-filter');

// ─── dark mode ────────────────────────────────────────────────────────────────
btnDark.addEventListener('click', () => {
  darkMode = !darkMode;
  btnDark.classList.toggle('active', darkMode);
  for (const el of document.querySelectorAll('.canvas-wrap')) {
    el.classList.toggle('dark', darkMode);
  }
  libPreview.classList.toggle('dark', darkMode);
  document.querySelectorAll('.build-box').forEach(el => el.classList.toggle('dark', darkMode));
});
// apply dark mode on load
btnDark.classList.add('active');
document.querySelectorAll('.canvas-wrap').forEach(el => el.classList.add('dark'));
libPreview.classList.add('dark');
document.querySelectorAll('.build-box').forEach(el => el.classList.add('dark'));

// ─── highlight toggle ─────────────────────────────────────────────────────────
btnHl.addEventListener('click', () => {
  hlEnabled = !hlEnabled;
  btnHl.classList.toggle('active', hlEnabled);
  if (!hlEnabled) clearHighlight();
});

// ─── verify mode ──────────────────────────────────────────────────────────────
btnVerify.addEventListener('click', () => {
  verifyMode = !verifyMode;
  btnVerify.classList.toggle('active', verifyMode);
  document.getElementById('verify-pane').style.display = verifyMode ? '' : 'none';
  compareBar.style.display = verifyMode ? '' : 'none';
  if (verifyMode) triggerVerify();
});

// ─── library sidebar (persistent) ────────────────────────────────────────────
function openLib() {
  libPanel.classList.add('open');
  btnLib.classList.add('active');
  if (!libLoaded) loadMonomers();
  else libSearch.focus();
  loadReactions();
}
function closeLib() {
  libPanel.classList.remove('open');
  btnLib.classList.remove('active');
  hidePreview();
}

btnLib.addEventListener('click', () =>
  libPanel.classList.contains('open') ? closeLib() : openLib());
libClose.addEventListener('click', closeLib);

libSearch.addEventListener('input', () => renderLibList(libSearch.value.trim().toLowerCase()));

btnRxnFilter.addEventListener('click', () => {
  rxnFilterActive = !rxnFilterActive;
  btnRxnFilter.classList.toggle('active', rxnFilterActive);
  renderLibList(libSearch.value.trim().toLowerCase());
});

// ─── example peptide sidebar ──────────────────────────────────────────────────
let examplesLoaded = false;

function openExamples() {
  examplesPanel.classList.add('open');
  btnExamples.classList.add('active');
  if (!examplesLoaded) loadExamples();
}
function closeExamples() {
  examplesPanel.classList.remove('open');
  btnExamples.classList.remove('active');
}

btnExamples.addEventListener('click', () =>
  examplesPanel.classList.contains('open') ? closeExamples() : openExamples());
examplesClose.addEventListener('click', closeExamples);

async function loadExamples() {
  try {
    const res = await fetch('/examples');
    const data = await res.json();
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
      cabilnInput.value = row.dataset.cabiln;
      cabilnInput.dispatchEvent(new Event('input'));
      cabilnInput.focus();
    });
  });
}

async function loadMonomers() {
  libList.innerHTML = '<div class="placeholder">Loading…</div>';
  try {
    const res = await fetch('/monomers');
    allMonomers = await res.json();
    libLoaded = true;
    renderLibList('');
    libSearch.focus();
  } catch (e) {
    libList.innerHTML = '<div class="placeholder err">Failed to load monomers</div>';
  }
}

async function loadReactions() {
  if (reactionPairs !== null) return;
  try {
    const res = await fetch('/reactions');
    reactionPairs = await res.json();
  } catch (e) { reactionPairs = []; }
}

function parseCts(cts) {
  if (!cts) return [];
  return cts.split(',').map(p => p.includes(':') ? p.split(':')[1].trim() : p.trim()).filter(Boolean);
}

function renderLibList(q) {
  let filtered = q
    ? allMonomers.filter(m =>
        m.abbr.toLowerCase().includes(q) ||
        m.name.toLowerCase().includes(q) ||
        m.type.toLowerCase().includes(q) ||
        (m.chem_types || '').toLowerCase().includes(q)
      )
    : allMonomers;

  if (rxnFilterActive && buildLeft && Array.isArray(reactionPairs)) {
    const pairSet = new Set(reactionPairs.map(([a, b]) => a + '|' + b));
    let lcts;
    if (buildLeft.selectedSlot !== null) {
      const selRg = buildLeft.rgroups.find(r => r.slot === buildLeft.selectedSlot);
      lcts = selRg && !selRg.used ? [selRg.chem_type] : [];
    } else {
      lcts = buildLeft.rgroups.filter(r => !r.used).map(r => r.chem_type);
    }
    filtered = filtered.filter(m => {
      const mcts = parseCts(m.chem_types);
      return mcts.some(mct => lcts.some(lct =>
        pairSet.has(lct + '|' + mct) || pairSet.has(mct + '|' + lct)
      ));
    });
  }

  libCount.textContent = `${filtered.length} / ${allMonomers.length} monomers`;

  if (!filtered.length) {
    libList.innerHTML = '<div class="placeholder">No matches</div>';
    return;
  }

  const rows = filtered.map(m => {
    const badge = m.degenerate
      ? `<span class="lib-badge cap">N/C cap</span>`
      : m.subtype === 'modified' || m.subtype === 'natural'
      ? `<span class="lib-badge aa">${m.type}</span>`
      : m.type === 'cap' && m.subtype === 'protecting'
      ? `<span class="lib-badge protect">cap</span>`
      : `<span class="lib-badge cap">${m.type}</span>`;

    let lg = m.leaving ? `  LG: ${escHtml(m.leaving)}` : '';
    if (m.degenerate) {
      const parts = [];
      if (m.nterm_abbr) parts.push(`N: ${escHtml(m.nterm_abbr)} (${escHtml(m.nterm_leaving)})`);
      if (m.cterm_abbr) parts.push(`C: ${escHtml(m.cterm_abbr)} (${escHtml(m.cterm_leaving)})`);
      lg = '  ' + parts.join(' | ');
    }
    return `<div class="lib-row" data-abbr="${escAttr(m.abbr)}">
      <div class="lib-abbr">${escHtml(m.abbr)}</div>
      <div class="lib-info">
        <div class="lib-name" title="${escAttr(m.name)}">${escHtml(m.name)}</div>
        <div class="lib-meta">${escHtml(m.chem_types || '')}${lg}</div>
      </div>
      ${badge}
    </div>`;
  });
  libList.innerHTML = rows.join('');

  libList.querySelectorAll('.lib-row').forEach(row => {
    row.addEventListener('click', () => {
      if (buildMode) {
        if (!cabilnInput.value.trim()) {
          // No sequence yet — insert as first monomer
          cabilnInput.value = row.dataset.abbr;
          cabilnInput.dispatchEvent(new Event('input'));
          buildHint.textContent = 'First monomer added — click its chip, then right-click another monomer';
        } else if (!buildLeft) {
          // Sequence exists but no chip selected — prompt
          buildHint.textContent = 'Click a chip on the sequence first, then click a library monomer';
        } else {
          loadBuildRight(row.dataset.abbr);
        }
      } else {
        insertAbbr(row.dataset.abbr);
      }
    });
    row.addEventListener('contextmenu', e => {
      if (buildMode) {
        e.preventDefault();
        loadBuildRight(row.dataset.abbr);
      }
    });
    row.addEventListener('mouseenter', e => startPreview(row.dataset.abbr, row));
    row.addEventListener('mouseleave', () => hidePreview());
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
  ta.value = before + insert + after;
  const newPos = start + insert.length;
  ta.setSelectionRange(newPos, newPos);
  ta.focus();
  ta.dispatchEvent(new Event('input'));
}

// ─── monomer preview tooltip ──────────────────────────────────────────────────
function startPreview(abbr, row) {
  clearTimeout(previewTimer);
  previewTimer = setTimeout(async () => {
    if (previewCache[abbr]) {
      showPreview(previewCache[abbr], row);
      return;
    }
    try {
      const res = await fetch(`/monomer_svg?abbr=${encodeURIComponent(abbr)}`);
      const data = await res.json();
      if (data.svg) {
        previewCache[abbr] = data;
        showPreview(data, row);
      }
    } catch (e) { /* silent */ }
  }, 200);
}

function showPreview(data, row) {
  let html;
  if (data.degenerate && data.variants) {
    // Degenerate: variant panels with optional reagent info
    let panels = data.variants.map(v => {
      let pane = `<div class="prev-pane"><span class="prev-label">${v.label}</span>${v.svg}`;
      if (v.reagent) {
        pane += `<span class="prev-rxn">${v.reagent.reaction} (LG: ${v.reagent.reagent_lg})</span>`;
      }
      pane += `</div>`;
      return pane;
    }).join('');
    // Show reagent form from first variant that has one
    const withReagent = data.variants.find(v => v.svg_reagent);
    if (withReagent) {
      panels += `<div class="prev-pane"><span class="prev-label">Reagent</span>${withReagent.svg_reagent}</div>`;
    }
    panels += `<div class="prev-pane"><span class="prev-label">R-groups</span>${data.svg}</div>`;
    html = `<div class="prev-row">${panels}</div>`;
  } else if (data.degenerate) {
    // Legacy format fallback
    html = `<div class="prev-row">
      <div class="prev-pane"><span class="prev-label">N-term</span>${data.svg_nterm || data.svg}</div>
      <div class="prev-pane"><span class="prev-label">C-term</span>${data.svg_cterm || data.svg}</div>
    </div>`;
  } else {
    const restored = data.svg_restored || data.svg;
    let reagentPane = '';
    let metaLine = '';
    if (data.svg_reagent) {
      reagentPane = `<div class="prev-pane"><span class="prev-label">Reagent</span>${data.svg_reagent}</div>`;
      const r = data.reagent;
      metaLine = `<div class="prev-meta">${r.reaction}` +
        (r.reagent_note ? ` — ${r.reagent_note}` : '') +
        (r.issue ? ` <span class="prev-warn">⚠ ${r.issue}</span>` : '') +
        `</div>`;
    }
    html = `<div class="prev-row">
      <div class="prev-pane"><span class="prev-label">Monomer</span>${restored}</div>
      ${reagentPane}
      <div class="prev-pane"><span class="prev-label">R-groups</span>${data.svg}</div>
    </div>${metaLine}`;
  }
  libPreview.innerHTML = html;
  const hasReagent = !!(data.svg_reagent || (data.variants && data.variants.some(v => v.svg_reagent)));
  libPreview.classList.toggle('has-reagent', hasReagent);
  const rect = row.getBoundingClientRect();
  const previewH = hasReagent ? 240 : 202;
  let top = Math.max(8, rect.top - 40);
  if (top + previewH > window.innerHeight - 8) {
    top = window.innerHeight - 8 - previewH;
  }
  top = Math.max(8, top);
  libPreview.style.left = (rect.right + 8) + 'px';
  libPreview.style.top  = top + 'px';
  libPreview.style.display = 'block';
}

function hidePreview() {
  clearTimeout(previewTimer);
  libPreview.style.display = 'none';
}

// ─── residue chips + bidirectional highlighting ───────────────────────────────
let chainData = [];
let xlinkByRes = {};

function highlightGroup(idxList) {
  if (!hlEnabled) return;
  clearHighlight();
  activeRIdx = -999;
  const svg = document.querySelector('#render-inner svg');
  if (!svg) return;
  svg.classList.add('has-highlight');
  for (const rIdx of idxList) {
    const atoms = residueMap[rIdx] || [];
    for (const aidx of atoms) {
      svg.querySelectorAll(`.atom-${aidx}`).forEach(el =>
        el.classList.add('res-hl'));
    }
  }
  resChips.classList.add('dimmed');
  resChips.querySelectorAll('.res-chip').forEach(c => {
    const ri = parseInt(c.dataset.residue);
    c.classList.toggle('hover', idxList.includes(ri));
  });
  resChips.querySelectorAll('.branch-chip').forEach(c => {
    const cm = JSON.parse(c.dataset.members || '[]');
    c.classList.toggle('hover', cm.some(m => idxList.includes(m)));
  });
}

function buildResidueUI(resMap, residues, chains, cabiln, bracketGroups, crosslinkGroups) {
  residueMap = resMap || {};
  residueList = residues || [];
  chainData = chains || [];
  atomToRes = {};
  for (const [rIdx, atoms] of Object.entries(residueMap)) {
    for (const aidx of atoms) atomToRes[aidx] = parseInt(rIdx);
  }

  resChips.innerHTML = '';
  resChips.style.position = '';
  resChips.style.paddingLeft = '';
  if (!residueList.length) return;

  const resById = {};
  residueList.forEach((r, i) => { resById[r.idx] = { ...r, colorIdx: i }; });

  const xlinkByMember = {};
  (crosslinkGroups || []).forEach(g => {
    g.members.forEach(mIdx => {
      if (!xlinkByMember[mIdx]) xlinkByMember[mIdx] = [];
      xlinkByMember[mIdx].push(g);
    });
  });
  xlinkByRes = xlinkByMember;

  function makeChip(rIdx, simpleHover) {
    const r = resById[rIdx];
    if (!r) return null;
    const chip = document.createElement('span');
    chip.className = 'res-chip';
    chip.textContent = r.abbr;
    chip.dataset.residue = r.idx;
    chip.style.background = RES_COLORS[r.colorIdx % RES_COLORS.length];
    const xlinks = xlinkByMember[r.idx];
    if (!simpleHover && xlinks && xlinks.length) {
      const allMembers = [...new Set(xlinks.flatMap(g => g.members))];
      chip.addEventListener('mouseenter', () => highlightGroup(allMembers));
    } else {
      chip.addEventListener('mouseenter', () => highlightResidue(r.idx));
    }
    chip.addEventListener('mouseleave', clearHighlight);
    chip.addEventListener('click', () => {
      if (buildMode) {
        if (buildLeft && buildLeftRIdx !== r.idx) {
          loadBuildRight(r.abbr, r.idx);
          chip.style.outline = '2px solid #e0a05a';
          resChips.querySelectorAll('.res-chip').forEach(c => {
            if (c !== chip && parseInt(c.dataset.residue) !== buildLeftRIdx)
              c.style.outline = '';
          });
        } else {
          loadBuildLeft(r.abbr, r.idx);
          chip.style.outline = '2px solid #5a9ae0';
          resChips.querySelectorAll('.res-chip').forEach(c => {
            if (c !== chip) c.style.outline = '';
          });
        }
      }
    });
    return chip;
  }

  function makeSeparator(text, memberIdxs) {
    const el = document.createElement('span');
    el.className = 'res-chip branch-chip';
    el.style.background = '#3a3a50';
    el.style.fontWeight = '700';
    el.textContent = text;
    el.dataset.members = JSON.stringify(memberIdxs || []);
    if (memberIdxs && memberIdxs.length) {
      el.addEventListener('mouseenter', () => highlightGroup(memberIdxs));
      el.addEventListener('mouseleave', clearHighlight);
    }
    return el;
  }

  function makeXlinkChip(tag, members) {
    const el = document.createElement('span');
    el.className = 'res-chip branch-chip xlink-chip';
    el.style.background = '#503a4a';
    el.style.fontWeight = '700';
    el.style.fontSize = '0.8em';
    el.textContent = tag;
    el.dataset.members = JSON.stringify(members);
    el.addEventListener('mouseenter', () => highlightGroup(members));
    el.addEventListener('mouseleave', clearHighlight);
    return el;
  }

  const nTermXlinkTags = new Set();
  if (cabiln) {
    for (const seg of cabiln.split(/[%\n]/)) {
      const m = seg.trim().match(/^(!\d+)-/);
      if (m) nTermXlinkTags.add(m[1]);
    }
  }

  function prependXlinks(rIdx) {
    const xlinks = xlinkByMember[rIdx];
    if (!xlinks) return;
    xlinks.forEach(g => {
      if (nTermXlinkTags.has(g.tag))
        resChips.appendChild(makeXlinkChip(g.tag, g.members));
    });
  }

  function appendXlinks(rIdx) {
    const xlinks = xlinkByMember[rIdx];
    if (!xlinks) return;
    xlinks.forEach(g => {
      if (!nTermXlinkTags.has(g.tag))
        resChips.appendChild(makeXlinkChip(g.tag, g.members));
    });
  }

  const hasBranch = cabiln && cabiln.includes('%');
  const groups = bracketGroups || [];
  const groupsByHost = {};
  groups.forEach(g => {
    if (!groupsByHost[g.host]) groupsByHost[g.host] = [];
    groupsByHost[g.host].push(g.members);
  });
  const branchSet = new Set();
  groups.forEach(g => g.members.forEach(m => branchSet.add(m)));

  function appendResidueWithBrackets(rIdx) {
    if (branchSet.has(rIdx)) return;
    prependXlinks(rIdx);
    const chip = makeChip(rIdx);
    if (chip) resChips.appendChild(chip);
    const hostGroups = groupsByHost[rIdx];
    if (hostGroups) {
      // $ chip selects residue + all its bracket branches (always, even single group)
      const allMembers = [rIdx, ...hostGroups.flat()];
      const el = document.createElement('span');
      el.className = 'res-chip branch-chip xlink-chip';
      el.style.background = '#3a5050';
      el.style.fontWeight = '700';
      el.style.fontSize = '0.8em';
      el.textContent = '$';
      el.dataset.members = JSON.stringify(allMembers);
      el.addEventListener('mouseenter', () => highlightGroup(allMembers));
      el.addEventListener('mouseleave', clearHighlight);
      resChips.appendChild(el);
      const brk = cabiln && cabiln.includes('{') ? ['{', '}'] : ['[', ']'];
      hostGroups.forEach(members => {
        resChips.appendChild(makeSeparator(brk[0], members));
        members.forEach(mIdx => {
          const mc = makeChip(mIdx);
          if (mc) resChips.appendChild(mc);
        });
        resChips.appendChild(makeSeparator(brk[1], members));
      });
    }
    appendXlinks(rIdx);
  }

  if (hasBranch) {
    chainData.forEach((chain, ci) => {
      if (chainData.length > 1) {
        // % chip highlights all chain residues + their crosslink partners
        const expanded = [...chain.residues];
        chain.residues.forEach(rIdx => {
          const xlinks = xlinkByMember[rIdx];
          if (xlinks) xlinks.forEach(g => g.members.forEach(m => expanded.push(m)));
        });
        resChips.appendChild(makeSeparator('%', [...new Set(expanded)]));
      }
      chain.residues.forEach(rIdx => {
        if (branchSet.has(rIdx)) return;
        prependXlinks(rIdx);

        // Crosslink brackets shown explicitly? Use simple hover on the chip itself
        const xlinks = xlinkByMember[rIdx];
        const hasXlinkBrackets = xlinks && xlinks.length > 1;
        const chip = makeChip(rIdx, hasXlinkBrackets);
        if (chip) resChips.appendChild(chip);

        // Bracket groups on this host (from bracket notation)
        const hostGroups = groupsByHost[rIdx];
        if (hostGroups) {
          // $ selects residue + all bracket branches (always, even single group)
          const bracketMembers = [rIdx, ...hostGroups.flat()];
          const el = document.createElement('span');
          el.className = 'res-chip branch-chip xlink-chip';
          el.style.background = '#3a5050';
          el.style.fontWeight = '700';
          el.style.fontSize = '0.8em';
          el.textContent = '$';
          el.dataset.members = JSON.stringify(bracketMembers);
          el.addEventListener('mouseenter', () => highlightGroup(bracketMembers));
          el.addEventListener('mouseleave', clearHighlight);
          resChips.appendChild(el);
          const brk = cabiln && cabiln.includes('{') ? ['{', '}'] : ['[', ']'];
          hostGroups.forEach(members => {
            resChips.appendChild(makeSeparator(brk[0], members));
            members.forEach(mIdx => {
              const mc = makeChip(mIdx);
              if (mc) resChips.appendChild(mc);
            });
            resChips.appendChild(makeSeparator(brk[1], members));
          });
        }

        // Crosslink chips: if multiple on one residue, add $ then wrap each in brackets
        if (hasXlinkBrackets) {
          const xlinkMembers = [rIdx, ...xlinks.flatMap(g => g.members)];
          const el = document.createElement('span');
          el.className = 'res-chip branch-chip xlink-chip';
          el.style.background = '#3a5050';
          el.style.fontWeight = '700';
          el.style.fontSize = '0.8em';
          el.textContent = '$';
          el.dataset.members = JSON.stringify(xlinkMembers);
          el.addEventListener('mouseenter', () => highlightGroup(xlinkMembers));
          el.addEventListener('mouseleave', clearHighlight);
          resChips.appendChild(el);
          xlinks.forEach(g => {
            const isNterm = nTermXlinkTags.has(g.tag);
            if (!isNterm) {
              resChips.appendChild(makeSeparator('[', g.members));
              resChips.appendChild(makeXlinkChip(g.tag, g.members));
              resChips.appendChild(makeSeparator(']', g.members));
            }
          });
        } else {
          appendXlinks(rIdx);
        }
      });
    });
  } else if (groups.length) {
    const mainChain = chainData.length ? chainData[0].residues : [];
    mainChain.forEach(rIdx => appendResidueWithBrackets(rIdx));
  } else {
    const allIds = chainData.flatMap(c => c.residues);
    allIds.forEach(rIdx => {
      prependXlinks(rIdx);
      const chip = makeChip(rIdx);
      if (chip) resChips.appendChild(chip);
      appendXlinks(rIdx);
    });
  }

  wireUpSvgHover();
}

let activeRIdx = null;

function highlightResidue(rIdx) {
  if (!hlEnabled) return;
  if (rIdx === activeRIdx) return;
  clearHighlight();
  activeRIdx = rIdx;
  const atoms = residueMap[rIdx] || [];
  const svg = document.querySelector('#render-inner svg');
  if (!svg) return;

  svg.classList.add('has-highlight');
  for (const aidx of atoms) {
    svg.querySelectorAll(`.atom-${aidx}`).forEach(el =>
      el.classList.add('res-hl'));
  }

  resChips.classList.add('dimmed');
  resChips.querySelectorAll('.res-chip').forEach(c =>
    c.classList.toggle('hover', parseInt(c.dataset.residue) === rIdx));
}

function clearHighlight() {
  activeRIdx = null;
  const svg = document.querySelector('#render-inner svg');
  if (svg) {
    svg.classList.remove('has-highlight');
    svg.querySelectorAll('.res-hl').forEach(el => el.classList.remove('res-hl'));
  }
  resChips.classList.remove('dimmed');
  resChips.querySelectorAll('.res-chip.hover').forEach(c => c.classList.remove('hover'));
}

function wireUpSvgHover() {
  const svg = document.querySelector('#render-inner svg');
  if (!svg) return;
  svg.addEventListener('mousemove', e => {
    if (!hlEnabled) return;
    let el = e.target;
    while (el && el !== svg) {
      const cls = (el.getAttribute('class') || '');
      const m = cls.match(/atom-(\d+)/);
      if (m) {
        const rIdx = atomToRes[parseInt(m[1])];
        if (rIdx !== undefined) {
          const xlinks = xlinkByRes[rIdx];
          if (xlinks && xlinks.length) {
            const allMembers = [...new Set(xlinks.flatMap(g => g.members))];
            highlightGroup(allMembers);
          } else {
            highlightResidue(rIdx);
          }
          return;
        }
      }
      el = el.parentElement;
    }
    clearHighlight();
  });
  svg.addEventListener('mouseleave', clearHighlight);
  svg.addEventListener('click', e => {
    if (!buildMode) return;
    let el = e.target;
    while (el && el !== svg) {
      const cls = (el.getAttribute('class') || '');
      const m = cls.match(/atom-(\d+)/);
      if (m) {
        const rIdx = atomToRes[parseInt(m[1])];
        if (rIdx !== undefined) {
          const r = residueList.find(r => r.idx === rIdx);
          if (r) {
            if (buildLeft && buildLeftRIdx !== rIdx) {
              loadBuildRight(r.abbr, r.idx);
              resChips.querySelectorAll('.res-chip').forEach(c => {
                const cIdx = parseInt(c.dataset.residue);
                if (cIdx === rIdx) c.style.outline = '2px solid #e0a05a';
                else if (cIdx !== buildLeftRIdx) c.style.outline = '';
              });
            } else {
              loadBuildLeft(r.abbr, r.idx);
              resChips.querySelectorAll('.res-chip').forEach(c => {
                c.style.outline = parseInt(c.dataset.residue) === rIdx
                  ? '2px solid #5a9ae0' : '';
              });
            }
          }
          return;
        }
      }
      el = el.parentElement;
    }
  });
}

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

// ─── PNG download (client-side SVG → canvas → PNG) ────────────────────────────
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

// ─── MOL download (server-side mol block) ────────────────────────────────────
btnMol.addEventListener('click', async () => {
  if (!lastMolBlock) return;
  const blob = new Blob([lastMolBlock], { type: 'chemical/x-mdl-molfile' });
  const a = document.createElement('a');
  a.download = 'structure.mol';
  a.href = URL.createObjectURL(blob);
  a.click();
});

// ─── build mode ──────────────────────────────────────────────────────────────
function openBuild() {
  buildMode = true;
  buildPanel.classList.add('open');
  btnBuild.classList.add('active');
  if (!libPanel.classList.contains('open')) openLib();
  clearBuild();
}
function closeBuild() {
  buildMode = false;
  buildPanel.classList.remove('open');
  btnBuild.classList.remove('active');
  clearBuild();
}
function clearBuild() {
  buildLeft = null; buildRight = null; buildLeftRIdx = null; buildRightRIdx = null;
  buildLeftAbbr.textContent = '—';
  buildLeftSvg.innerHTML = '<div class="box-placeholder">Click a chip above</div>';
  buildLeftRg.innerHTML = '';
  buildRightAbbr.textContent = '—';
  buildRightSvg.innerHTML = '<div class="box-placeholder">Right-click from library</div>';
  buildRightRg.innerHTML = '';
  buildConnect.disabled = true;
  buildStatus.textContent = '';
  buildStatus.className = 'build-status';
  buildHint.textContent = 'Select a chip on the sequence, then right-click a monomer in the library';
  btnRxnFilter.disabled = true;
  if (rxnFilterActive) {
    rxnFilterActive = false;
    btnRxnFilter.classList.remove('active');
    if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  }
}

btnBuild.addEventListener('click', () => buildMode ? closeBuild() : openBuild());
buildClose.addEventListener('click', closeBuild);
document.getElementById('build-left-change').addEventListener('click', clearBuild);
document.getElementById('build-right-change').addEventListener('click', clearBuild);

async function loadBuildLeft(abbr, rIdx) {
  buildLeftRIdx = rIdx;
  buildLeftAbbr.textContent = abbr;
  buildLeftSvg.innerHTML = '<div class="spinner"></div>';
  buildLeftRg.innerHTML = '';
  buildLeft = null;
  buildConnect.disabled = true;
  buildStatus.textContent = '';

  try {
    const res = await fetch(`/monomer_rgroups?abbr=${encodeURIComponent(abbr)}&residue_idx=${rIdx}&cabiln=${encodeURIComponent(cabilnInput.value.trim())}`);
    const data = await res.json();
    if (data.error) {
      buildLeftSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
      return;
    }
    buildLeftSvg.innerHTML = data.svg || '';
    buildLeft = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
    renderRgroupButtons(buildLeftRg, buildLeft, 'left');
    buildHint.textContent = buildRight ? 'Select R-groups to connect' : 'Now right-click a monomer in the library';
    btnRxnFilter.disabled = false;
    if (rxnFilterActive && libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  } catch (e) {
    buildLeftSvg.innerHTML = '<div class="box-placeholder">Error loading monomer</div>';
  }
}

async function loadBuildRight(abbr, rIdx) {
  buildRightRIdx = rIdx !== undefined ? rIdx : null;
  buildRightAbbr.textContent = abbr;
  buildRightSvg.innerHTML = '<div class="spinner"></div>';
  buildRightRg.innerHTML = '';
  buildRight = null;
  buildConnect.disabled = true;
  buildStatus.textContent = '';

  try {
    let url = `/monomer_rgroups?abbr=${encodeURIComponent(abbr)}`;
    if (rIdx !== undefined) {
      url += `&residue_idx=${rIdx}&cabiln=${encodeURIComponent(cabilnInput.value.trim())}`;
    }
    const res = await fetch(url);
    const data = await res.json();
    if (data.error) {
      buildRightSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
      return;
    }
    buildRightSvg.innerHTML = data.svg || '';
    buildRight = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
    renderRgroupButtons(buildRightRg, buildRight, 'right');
    buildHint.textContent = 'Select R-groups to connect';
  } catch (e) {
    buildRightSvg.innerHTML = '<div class="box-placeholder">Error loading monomer</div>';
  }
}

function renderRgroupButtons(container, state, side) {
  container.innerHTML = '';
  state.rgroups.forEach(rg => {
    const btn = document.createElement('button');
    btn.className = 'rgroup-btn';
    if (rg.used) btn.classList.add('used');
    if (state.selectedSlot === rg.slot) btn.classList.add('selected');
    btn.textContent = `R${rg.slot} ${rg.chem_type || ''}`;
    btn.title = `R${rg.slot}: ${rg.chem_type || 'unknown'}${rg.leaving ? ' (LG: ' + rg.leaving + ')' : ''}`;
    if (!rg.used) {
      btn.addEventListener('click', () => selectRgroup(side, rg.slot));
    }
    container.appendChild(btn);
  });
}

function selectRgroup(side, slot) {
  if (side === 'left' && buildLeft) {
    buildLeft.selectedSlot = buildLeft.selectedSlot === slot ? null : slot;
    renderRgroupButtons(buildLeftRg, buildLeft, 'left');
    if (rxnFilterActive && libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  } else if (side === 'right' && buildRight) {
    buildRight.selectedSlot = buildRight.selectedSlot === slot ? null : slot;
    renderRgroupButtons(buildRightRg, buildRight, 'right');
  }
  checkBuildValidity();
}

async function checkBuildValidity() {
  if (!buildLeft?.selectedSlot || !buildRight?.selectedSlot) {
    buildConnect.disabled = true;
    buildStatus.textContent = '';
    buildStatus.className = 'build-status';
    return;
  }

  const leftRg = buildLeft.rgroups.find(r => r.slot === buildLeft.selectedSlot);
  const rightRg = buildRight.rgroups.find(r => r.slot === buildRight.selectedSlot);
  if (!leftRg || !rightRg) return;

  buildStatus.textContent = 'Checking bond...';
  buildStatus.className = 'build-status';

  try {
    const res = await fetch('/validate_bond', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        chem_type_a: leftRg.chem_type,
        chem_type_b: rightRg.chem_type,
        abbr_a: buildLeft.abbr,
        slot_a: buildLeft.selectedSlot,
        abbr_b: buildRight.abbr,
        slot_b: buildRight.selectedSlot
      })
    });
    const data = await res.json();
    if (data.valid) {
      buildStatus.textContent = `Valid: ${data.reaction || 'bond'} (R${buildLeft.selectedSlot}↔R${buildRight.selectedSlot})`;
      buildStatus.className = 'build-status valid';
      buildConnect.disabled = false;
    } else {
      buildStatus.textContent = data.reason || 'No compatible reaction found';
      buildStatus.className = 'build-status invalid';
      buildConnect.disabled = true;
    }
  } catch (e) {
    buildStatus.textContent = 'Validation error';
    buildStatus.className = 'build-status invalid';
    buildConnect.disabled = true;
  }
}

buildConnect.addEventListener('click', async () => {
  if (!buildLeft || !buildRight || !buildLeft.selectedSlot || !buildRight.selectedSlot) return;

  const val = cabilnInput.value.trim();
  const rHost = buildLeft.selectedSlot;
  const rNew = buildRight.selectedSlot;
  const newAbbr = buildRight.abbr;

  buildConnect.disabled = true;
  buildStatus.textContent = 'Inserting...';
  buildStatus.className = 'build-status';

  try {
    const res = await fetch('/insert_bond', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        cabiln: val,
        host_residue_idx: buildLeftRIdx ?? 0,
        new_abbr: newAbbr,
        r_host: rHost,
        r_new: rNew,
        target_residue_idx: (buildRightRIdx !== null && buildRightRIdx !== buildLeftRIdx)
                             ? buildRightRIdx : -1
      })
    });
    const data = await res.json();
    if (data.error) {
      buildStatus.textContent = data.error;
      buildStatus.className = 'build-status invalid';
      return;
    }
    cabilnInput.value = data.result;
    cabilnInput.dispatchEvent(new Event('input'));

    // Clear right side for next addition
    buildRight = null;
    buildRightAbbr.textContent = '—';
    buildRightSvg.innerHTML = '<div class="box-placeholder">Right-click from library</div>';
    buildRightRg.innerHTML = '';
    buildConnect.disabled = true;
    buildStatus.textContent = '';
    buildStatus.className = 'build-status';
    buildHint.textContent = 'Connection added — select a chip and right-click another monomer';
  } catch (e) {
    buildStatus.textContent = 'Insert failed';
    buildStatus.className = 'build-status invalid';
  }
});

// ─── notation conversion ──────────────────────────────────────────────────────
async function convertNotation(target) {
  const val = cabilnInput.value.trim();
  if (!val) return;
  try {
    const res = await fetch('/convert_notation', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ cabiln: val, target })
    });
    const data = await res.json();
    if (data.error) { console.error('convert error:', data.error); return; }
    if (data.result) {
      cabilnInput.value = data.result;
      cabilnInput.dispatchEvent(new Event('input'));
    }
  } catch (e) { console.error('convertNotation error:', e); }
}
btnToBracket.addEventListener('click', () => convertNotation('bracket'));
btnToBranch.addEventListener('click', () => convertNotation('branch'));

btnReroll.addEventListener('click', () => {
  if (!lastCabiln) return;
  rerollSeed++;
  btnReroll.textContent = rerollSeed % 2 === 1 ? '⟳ Indigo' : '⟳ CoordGen';
  showSpinner(renderInner);
  doRenderCabiln(lastCabiln);
});

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

// ─── render helpers ───────────────────────────────────────────────────────────
function canvasSize(el) {
  return { w: Math.max(el.clientWidth || 600, 400),
           h: Math.max(el.clientHeight || 500, 300) };
}

function setInner(inner, html) {
  inner.innerHTML = html;
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
  rerollSeed = 0;
  btnReroll.textContent = '⟳ Layout';
  if (!seq) { resetCabiln(); return; }
  showSpinner(renderInner);
  cabilnTimer = setTimeout(() => doRenderCabiln(seq), 300);
});

function resetCabiln() {
  lastCabiln = '';
  rerollSeed = 0;
  btnReroll.disabled = true;
  btnReroll.textContent = '⟳ Layout';
  cabilnInput.className = '';
  cabilnStatus.textContent = '';
  cabilnStatus.className = 'statusbar';
  setInner(renderInner, '<div class="placeholder">Start typing a sequence…</div>');
  clearExports();
  resChips.innerHTML = '';
  residueMap = {}; atomToRes = {}; residueList = [];
  if (verifyMode) {
    compareBar.innerHTML = '';
  }
}

async function doRenderCabiln(seq) {
  lastCabiln = seq;
  if (buildMode) clearBuild();
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res  = await fetch('/render', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: seq, width: w, height: h, seed: rerollSeed })
    });
    const data = await res.json();
    if (data.error) {
      const html = `<div class="placeholder err">${escHtml(data.error)}</div>`;
      setInner(renderInner, html);
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
      cabilnInput.className = 'err';
      clearExports();
      btnReroll.disabled = true;
      resChips.innerHTML = '';
    } else {
      setInner(renderInner, data.svg);
      cabilnStatus.textContent = data.info || '';
      cabilnStatus.className = 'statusbar ok';
      cabilnInput.className = 'ok';
      btnReroll.disabled = false;
      setExportReady(data.svg, data.mol_block);
      buildResidueUI(data.residue_map, data.residues, data.chains, data.cabiln_echo, data.bracket_groups, data.crosslink_groups);
      if (verifyMode && lastSmiles) triggerVerify();
    }
  } catch (e) {
    cabilnStatus.textContent = 'Server error — is the renderer running?';
    cabilnStatus.className = 'statusbar';
  }
}

// ─── reference render (verify mode) — auto-detects SMILES / BILN / HELM ─────
smilesInput.addEventListener('input', () => {
  clearTimeout(smilesTimer);
  const txt = smilesInput.value.trim();
  if (txt === lastSmiles) return;
  if (!txt) {
    lastSmiles = '';
    setInner(smilesInner, '<div class="placeholder">Paste SMILES, BILN, or HELM…</div>');
    compareBar.innerHTML = '';
    return;
  }
  showSpinner(smilesInner);
  smilesTimer = setTimeout(() => doRenderRef(txt), 300);
});

molUpload.addEventListener('change', async (e) => {
  const file = e.target.files[0];
  if (!file) return;
  const text = await file.text();
  showSpinner(smilesInner);
  smilesStatus.textContent = `Loaded: ${file.name}`;
  smilesStatus.className = 'statusbar ok';
  try {
    const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
    const res = await fetch('/render_mol', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ mol_block: text, width: w, height: h })
    });
    const data = await res.json();
    if (data.error) {
      setInner(smilesInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      smilesStatus.textContent = data.error;
      smilesStatus.className = 'statusbar';
    } else {
      setInner(smilesInner, data.svg);
      lastSmiles = data.smiles || '';
      smilesInput.value = lastSmiles;
      smilesInput.className = 'ok';
      if (lastCabiln) triggerVerify();
    }
  } catch (err) {
    smilesStatus.textContent = 'Failed to render .mol file';
    smilesStatus.className = 'statusbar';
  }
  molUpload.value = '';
});

async function doRenderRef(txt) {
  const { w, h } = canvasSize(document.getElementById('smiles-canvas'));
  try {
    const res = await fetch('/render_reference', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt, width: w, height: h })
    });
    const data = await res.json();
    if (data.error) {
      setInner(smilesInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      smilesStatus.textContent = data.error;
      smilesStatus.className = 'statusbar';
      smilesInput.className = 'err';
    } else {
      setInner(smilesInner, data.svg);
      lastSmiles = data.smiles || '';
      smilesStatus.textContent = `${data.format}: ${data.info || ''}`;
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
      `<span class="canon" title="SMILES canonical: ${escAttr(data.smiles_canonical)}">` +
      `SMILES: ${escHtml(data.smiles_canonical.slice(0, 80))}${data.smiles_canonical.length > 80 ? '…' : ''}</span>` +
      `<span class="canon" title="CABILN canonical: ${escAttr(data.cabiln_canonical)}">` +
      `CABILN: ${escHtml(data.cabiln_canonical.slice(0, 80))}${data.cabiln_canonical.length > 80 ? '…' : ''}</span>`;
  } catch (e) {
    compareBar.innerHTML = '<span class="nomatch">Verify error</span>';
  }
}

function escHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
function escAttr(s) { return escHtml(s); }

// Auto-reload on server restart
(function(){
  let sid = null;
  setInterval(async () => {
    try {
      const r = await fetch('/server_id');
      const id = await r.text();
      if (sid === null) { sid = id; return; }
      if (id !== sid) location.reload();
    } catch(e) {}
  }, 5000);
})();
</script>
</body>
</html>
"""

# ── Registration page HTML ────────────────────────────────────────────────────

_REGISTER_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Register Monomer — pyPept</title>
<style>
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  html, body {
    min-height: 100%;
    font-family: system-ui, -apple-system, sans-serif;
    background: #12192b;
    color: #d0dce8;
  }
  header {
    padding: 10px 20px;
    background: #0d1422;
    border-bottom: 1px solid #1e3050;
    display: flex;
    align-items: center;
    gap: 14px;
  }
  header a { color: #4a7ab0; font-size: 12px; text-decoration: none; }
  header a:hover { color: #7aaeff; }
  header .title { flex: 1; font-size: 14px; color: #7aaeff; letter-spacing: .07em; }

  .page { max-width: 860px; margin: 0 auto; padding: 24px 20px; }

  .card {
    background: #0d1422;
    border: 1px solid #1e3050;
    border-radius: 8px;
    padding: 20px;
    margin-bottom: 20px;
  }
  .card h2 { font-size: 13px; color: #7aaeff; letter-spacing: .08em; margin-bottom: 14px; }

  .field { margin-bottom: 14px; }
  .field label { display: block; font-size: 11px; color: #4a7ab0; margin-bottom: 5px; letter-spacing: .1em; text-transform: uppercase; }
  .field input, .field select, .field textarea {
    width: 100%;
    background: #0a1018;
    border: 1px solid #1e3050;
    border-radius: 5px;
    color: #c8daf0;
    padding: 7px 11px;
    font-family: "Cascadia Code", "Fira Mono", monospace;
    font-size: 13px;
    outline: none;
    transition: border-color .15s;
  }
  .field input:focus, .field select:focus, .field textarea:focus { border-color: #3a6fd8; }
  .field input.err, .field textarea.err { border-color: #d9534f; }
  .field input.ok,  .field textarea.ok  { border-color: #28a745; }
  .field textarea { height: 62px; resize: vertical; }
  .field select option { background: #0d1422; }

  .field-row { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }

  .btn {
    padding: 8px 20px;
    border-radius: 6px;
    font-size: 13px;
    cursor: pointer;
    border: 1px solid;
    transition: background .15s, color .15s;
  }
  .btn-primary { background: #1e3a70; border-color: #3a6fd8; color: #c8daf0; }
  .btn-primary:hover { background: #2a4e90; }
  .btn-success { background: #1a3828; border-color: #28a745; color: #5dba7a; }
  .btn-success:hover { background: #1f4830; }
  .btn:disabled { opacity: .35; cursor: default; }

  #preview-section { display: none; }

  #preview-canvas {
    background: #fff;
    border-radius: 4px;
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 280px;
    margin-bottom: 14px;
    overflow: hidden;
  }
  #preview-canvas svg { max-width: 100%; display: block; }

  .detected-types {
    background: #0a1018;
    border: 1px solid #1e3050;
    border-radius: 5px;
    padding: 10px 14px;
    font-family: "Cascadia Code", "Fira Mono", monospace;
    font-size: 12px;
    color: #7aaeff;
    margin-bottom: 14px;
  }
  .detected-types .slot { color: #3dbe6c; }
  .detected-types .lg   { color: #a0b8d0; font-size: 11px; }

  #status-msg {
    margin-top: 10px;
    padding: 8px 14px;
    border-radius: 5px;
    font-size: 12px;
    display: none;
  }
  #status-msg.ok  { background: #1a3828; border: 1px solid #28a745; color: #5dba7a; }
  #status-msg.err { background: #28181a; border: 1px solid #d9534f; color: #e07080; }
</style>
</head>
<body>
<header>
  <span class="title">Register Monomer <span style="color:#3a5580">— pyPept</span></span>
  <a href="/">← Back to Renderer</a>
</header>

<div class="page">
  <div class="card">
    <h2>STEP 1 — SMILES INPUT</h2>
    <div class="field">
      <label>Full monomer SMILES (all atoms, including leaving groups)</label>
      <textarea id="smiles-in" placeholder="e.g.  NCC(=O)O  (glycine)  or  O=C(O)CN1C(=O)C=CC1=O  (Mal-Gly)" spellcheck="false" autocomplete="off"></textarea>
    </div>
    <button class="btn btn-primary" id="btn-preview">Preview &amp; detect R-groups</button>
  </div>

  <div id="preview-section">
    <div class="card">
      <h2>STEP 2 — DETECTED R-GROUPS &amp; PREVIEW</h2>
      <div id="preview-canvas"><div style="color:#b0bec5;font-size:13px">…</div></div>
      <div id="detected-display" class="detected-types"></div>
      <div class="field">
        <label>CHUCKLES (auto-generated — edit only if needed)</label>
        <textarea id="chuckles-out" spellcheck="false" autocomplete="off"></textarea>
      </div>
    </div>

    <div class="card">
      <h2>STEP 3 — METADATA</h2>
      <div class="field-row">
        <div class="field">
          <label>Abbreviation (unique token)</label>
          <input id="abbr-in" type="text" placeholder="e.g.  MalGly" autocomplete="off">
        </div>
        <div class="field">
          <label>Full name</label>
          <input id="name-in" type="text" placeholder="e.g.  Maleimidoglycine" autocomplete="off">
        </div>
      </div>
      <div class="field-row">
        <div class="field">
          <label>Type</label>
          <select id="type-in">
            <option value="aa">aa</option>
            <option value="cap">cap</option>
            <option value="linker">linker</option>
          </select>
        </div>
        <div class="field">
          <label>Subtype</label>
          <select id="subtype-in">
            <option value="natural">natural</option>
            <option value="modified">modified</option>
            <option value="cap">cap</option>
            <option value="protecting">protecting</option>
            <option value="label">label</option>
            <option value="linker">linker</option>
          </select>
        </div>
      </div>
      <button class="btn btn-success" id="btn-register" disabled>Register monomer</button>
      <div id="status-msg"></div>
    </div>
  </div>
</div>

<script>
let detectedData = null;

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

btnPreview.addEventListener('click', async () => {
  const smi = smilesIn.value.trim();
  if (!smi) return;
  btnPreview.disabled = true;
  btnPreview.textContent = 'Detecting…';
  prevCanvas.innerHTML = '<div style="color:#b0bec5;font-size:13px">Analysing…</div>';
  detDisp.innerHTML = '';
  prevSec.style.display = '';
  try {
    const res  = await fetch('/preview_monomer', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ smiles: smi })
    });
    const data = await res.json();
    if (data.error) {
      prevCanvas.innerHTML = `<div style="color:#d9534f;font-size:13px;padding:20px">${escHtml(data.error)}</div>`;
      smilesIn.className = 'err';
      btnRegister.disabled = true;
    } else {
      detectedData = data;
      prevCanvas.innerHTML = data.svg;
      chucklesOut.value = data.chuckles;
      smilesIn.className = 'ok';
      renderDetected(data.chem_types, data.leaving);
      btnRegister.disabled = false;
    }
  } catch (e) {
    prevCanvas.innerHTML = '<div style="color:#d9534f;font-size:13px;padding:20px">Server error</div>';
  } finally {
    btnPreview.disabled = false;
    btnPreview.textContent = 'Preview & detect R-groups';
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
  const abbr = abbrIn.value.trim();
  const name = nameIn.value.trim();
  if (!abbr || !name || !detectedData) {
    showStatus('err', 'Please fill in abbreviation and name.');
    return;
  }
  btnRegister.disabled = true;
  btnRegister.textContent = 'Registering…';
  try {
    const res  = await fetch('/register_monomer', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        chuckles: chucklesOut.value.trim() || detectedData.chuckles,
        chem_types: detectedData.chem_types,
        leaving: detectedData.leaving,
        abbr,
        name,
        type:    typeIn.value,
        subtype: subtypeIn.value,
      })
    });
    const data = await res.json();
    if (data.error) {
      showStatus('err', data.error);
    } else {
      showStatus('ok', `${abbr} registered successfully. Monomer count: ${data.total}`);
      btnRegister.textContent = '✓ Registered';
    }
  } catch (e) {
    showStatus('err', 'Server error');
  } finally {
    if (btnRegister.textContent !== '✓ Registered') btnRegister.disabled = false;
    if (btnRegister.textContent !== '✓ Registered') btnRegister.textContent = 'Register monomer';
  }
});

function showStatus(cls, msg) {
  statusMsg.className = cls;
  statusMsg.textContent = msg;
  statusMsg.style.display = '';
}

function escHtml(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}
</script>
</body>
</html>
"""


app = FastAPI()


# ── shared drawing helper ─────────────────────────────────────────────────────

def _indigo_layout(romol):
    """Lay out romol using Indigo's algorithm; copy coords back preserving atom indices."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.rdchem import Conformer
    try:
        from indigo import Indigo as _Indigo
    except ImportError:
        return None
    # mol block needs a conformer — compute a quick one just for the connectivity export
    tmp = Chem.RWMol(romol)
    if tmp.GetNumConformers() == 0:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(tmp)
    mol_block = Chem.MolToMolBlock(tmp)
    try:
        indigo = _Indigo()
        indigo.setOption('ignore-stereochemistry-errors', True)
        im = indigo.loadMolecule(mol_block)
        im.layout()
        result_block = im.molfile()
        result = Chem.MolFromMolBlock(result_block, removeHs=False, sanitize=False)
        if result is None or result.GetNumConformers() == 0:
            return None
        src = result.GetConformer()
        n = romol.GetNumAtoms()
        if romol.GetNumConformers() == 0:
            conf = Conformer(n)
            for i in range(n):
                p = src.GetAtomPosition(i)
                conf.SetAtomPosition(i, (p.x, p.y, 0.0))
            romol.AddConformer(conf, assignId=True)
        else:
            conf = romol.GetConformer()
            for i in range(n):
                p = src.GetAtomPosition(i)
                conf.SetAtomPosition(i, (p.x, p.y, 0.0))
        return romol
    except Exception:
        return None


# Cached flag: None = not yet checked, True/False = availability
_rdmoldraw2d_available: bool | None = None

def _get_rdmoldraw2d():
    global _rdmoldraw2d_available
    if _rdmoldraw2d_available is None:
        try:
            from rdkit.Chem.Draw import rdMolDraw2D  # noqa: F401
            _rdmoldraw2d_available = True
        except (ImportError, OSError):
            _rdmoldraw2d_available = False
    return _rdmoldraw2d_available


def _indigo_render_svg(romol, width: int, height: int) -> str:
    """SVG via Indigo renderer — no system X11/Cairo required (Vercel fallback)."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    try:
        from indigo import Indigo as _Indigo
        from indigo.renderer import IndigoRenderer as _Renderer
    except ImportError:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
            f'<text x="10" y="30" fill="red" font-size="14">Indigo not installed</text></svg>'
        )
    tmp = Chem.RWMol(romol)
    if tmp.GetNumConformers() == 0:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(tmp)
    mol_block = Chem.MolToMolBlock(tmp)
    try:
        ind = _Indigo()
        ind.setOption("ignore-stereochemistry-errors", True)
        ind.setOption("render-output-format", "svg")
        ind.setOption("render-image-width", width)
        ind.setOption("render-image-height", height)
        im = ind.loadMolecule(mol_block)
        rend = _Renderer(ind)
        return rend.renderToBuffer(im).decode("utf-8", errors="replace")
    except Exception as exc:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
            f'<text x="10" y="30" fill="red" font-size="14">Render error: {exc}</text></svg>'
        )


def _overlap_score(romol) -> int:
    """Count non-bonded atom pairs closer than 0.8 Å — lower is better."""
    conf = romol.GetConformer()
    n = romol.GetNumAtoms()
    pos = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(n)]
    bonded: set = set()
    for b in romol.GetBonds():
        u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonded.add((u, v)); bonded.add((v, u))
    thresh = 0.64  # 0.8² Å²
    score = 0
    for i in range(n):
        xi, yi = pos[i]
        for j in range(i + 1, n):
            if (i, j) not in bonded:
                dx, dy = xi - pos[j][0], yi - pos[j][1]
                if dx * dx + dy * dy < thresh:
                    score += 1
    return score


def _best_layout(romol, base_seed: int, n_tries: int = 6):
    """Try n_tries random atom orderings with CoordGen; apply the layout with fewest
    overlaps back to romol using the original atom indices (preserves residue map)."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.rdchem import Conformer
    import random

    n = romol.GetNumAtoms()
    best_positions = None
    best_score = 10 ** 9

    # Compute baseline CoordGen score so we never return something worse
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(romol)
    baseline_score = _overlap_score(romol)
    baseline_conf = romol.GetConformer()
    baseline_positions = {i: (baseline_conf.GetAtomPosition(i).x,
                               baseline_conf.GetAtomPosition(i).y) for i in range(n)}
    best_score = baseline_score
    best_positions = baseline_positions

    for i in range(n_tries):
        rng = random.Random(base_seed + i)
        new_order = list(range(n))
        rng.shuffle(new_order)
        candidate = Chem.RenumberAtoms(romol, new_order)
        try:
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(candidate)
        except Exception:
            continue
        score = _overlap_score(candidate)
        if score < best_score:
            best_score = score
            conf = candidate.GetConformer()
            # new_order[new_idx] == old_idx → translate back to original indexing
            best_positions = {}
            for new_idx in range(n):
                p = conf.GetAtomPosition(new_idx)
                best_positions[new_order[new_idx]] = (p.x, p.y)

    if best_positions is None:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(romol)
        return romol

    if romol.GetNumConformers() == 0:
        conf = Conformer(n)
        for idx, (x, y) in best_positions.items():
            conf.SetAtomPosition(idx, (x, y, 0.0))
        romol.AddConformer(conf, assignId=True)
    else:
        conf = romol.GetConformer()
        for idx, (x, y) in best_positions.items():
            conf.SetAtomPosition(idx, (x, y, 0.0))
    return romol


def _draw_mol(romol, width: int, height: int, used_slots: set | None = None, seed: int = 0) -> str:
    import re as _re
    from rdkit.Chem import rdDepictor

    if not _get_rdmoldraw2d():
        # Vercel / environments without system X11: use Indigo SVG renderer
        if seed == 0:
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(romol)
        elif seed % 2 == 1:
            _indigo_layout(romol)
        else:
            romol = _best_layout(romol, base_seed=seed * 6)
        rdDepictor.NormalizeDepiction(romol)
        rdDepictor.StraightenDepiction(romol)
        return _indigo_render_svg(romol, width, height)

    from rdkit.Chem.Draw import rdMolDraw2D
    if seed == 0:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(romol)
    elif seed % 2 == 1:
        # Odd clicks: Indigo layout (falls back to CoordGen reorder if unavailable)
        if _indigo_layout(romol) is None:
            romol = _best_layout(romol, base_seed=seed * 6)
    else:
        # Even clicks: CoordGen best-of-6 atom reorderings
        romol = _best_layout(romol, base_seed=seed * 6)
    rdDepictor.NormalizeDepiction(romol)
    rdDepictor.StraightenDepiction(romol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.padding = 0.12
    natoms = romol.GetNumAtoms()
    if natoms > 150:
        opts.minFontSize = 5
        opts.maxFontSize = 8
        opts.bondLineWidth = 0.8
        opts.additionalAtomLabelPadding = 0.0
    elif natoms > 80:
        opts.minFontSize = 7
        opts.maxFontSize = 10
        opts.bondLineWidth = 1.0
    drawer.DrawMolecule(romol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    if used_slots is not None:
        dummy_slot = {}
        for atom in romol.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() > 0:
                dummy_slot[atom.GetIdx()] = atom.GetIsotope()

        COLOR_USED = '#d9534f'
        COLOR_FREE = '#3dbe6c'
        OPACITY_USED = '0.35'

        for aidx, slot in dummy_slot.items():
            color = COLOR_USED if slot in used_slots else COLOR_FREE
            tag = f'atom-{aidx}'
            svg = _re.sub(
                rf"(<path\s+class='[^']*{tag}[^']*'[^>]*fill=')[^']+(')",
                rf"\g<1>{color}\2",
                svg,
            )
            if slot in used_slots:
                svg = _re.sub(
                    rf"(<path\s+class='[^']*{tag}[^']*'[^>]*)(/>)",
                    rf"\1 opacity='{OPACITY_USED}'\2",
                    svg,
                )
            svg = _re.sub(
                rf"(<path\s+class='bond-\d+\s+[^']*{tag}[^']*'[^>]*stroke:)#[0-9a-fA-F]{{6}}",
                rf"\g<1>{color}",
                svg,
            )
            if slot in used_slots:
                svg = _re.sub(
                    rf"(<path\s+class='bond-\d+\s+[^']*{tag}[^']*'[^>]*stroke-opacity:)\d[\d.]*",
                    rf"\g<1>{OPACITY_USED}",
                    svg,
                )
    return svg


def _mol_block(romol) -> str:
    from rdkit.Chem import MolToMolBlock
    return MolToMolBlock(romol)


def _canon(romol) -> str:
    from rdkit.Chem import MolToSmiles
    return MolToSmiles(romol, canonical=True)


# ── request models ────────────────────────────────────────────────────────────

class _CabilnReq(BaseModel):
    cabiln: str
    width:  int = 960
    height: int = 680
    seed:   int = 0

class _SmilesReq(BaseModel):
    smiles: str
    width:  int = 960
    height: int = 680

class _VerifyReq(BaseModel):
    smiles: str
    cabiln: str

class _PreviewReq(BaseModel):
    smiles: str
    width:  int = 640
    height: int = 440

class _RegisterReq(BaseModel):
    chuckles:   str
    chem_types: dict
    leaving:    dict
    abbr:       str
    name:       str
    type:       str = "aa"
    subtype:    str = "modified"


# ── routes ────────────────────────────────────────────────────────────────────

@app.get("/server_id")
async def server_id():
    return PlainTextResponse(_SERVER_ID)


@app.get("/", response_class=HTMLResponse)
async def index():
    return _HTML


@app.get("/register", response_class=HTMLResponse)
async def register_page():
    return _REGISTER_HTML


@app.get("/monomers")
async def list_monomers():
    try:
        from rdkit import Chem
        import sys, pathlib
        sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
        all_mols, mol_by_abbr = _load_sdf()
        monomers = []
        degen_nterm = {}  # base -> entry with trailing _
        degen_cterm = {}  # base -> entry with leading _
        for mol in all_mols:
            if mol is None:
                continue
            p = mol.GetPropsAsDict()
            abbr = p.get('m_abbr', '') or p.get('symbol', '')
            if not abbr:
                continue
            rgroups = p.get('m_Rgroups', '')
            slots   = [r.strip() for r in rgroups.split(',')]
            lg_parts = [f"R{i+1}:{slots[i]}" for i in range(len(slots))
                        if slots[i] not in ('None', '', 'none')]
            entry = {
                'abbr':      abbr,
                'name':      p.get('m_name', ''),
                'type':      p.get('m_type', ''),
                'subtype':   p.get('m_subtype', ''),
                'chem_types': p.get('m_chem_types', ''),
                'leaving':   ', '.join(lg_parts),
            }
            # Collect degenerate cap pairs for merging
            if abbr.endswith('_') and not abbr.startswith('_'):
                degen_nterm[abbr[:-1]] = entry
                continue
            if abbr.startswith('_') and not abbr.endswith('_'):
                degen_cterm[abbr[1:]] = entry
                continue
            monomers.append(entry)
        # Merge degenerate pairs into single entries
        all_bases = set(degen_nterm) | set(degen_cterm)
        for base in sorted(all_bases):
            nt = degen_nterm.get(base)
            ct = degen_cterm.get(base)
            primary = nt or ct
            merged = {
                'abbr': base,
                'name': primary['name'],
                'type': primary['type'],
                'subtype': primary['subtype'],
                'chem_types': primary['chem_types'],
                'leaving': primary['leaving'],
                'degenerate': True,
            }
            if nt:
                merged['nterm_abbr'] = nt['abbr']
                merged['nterm_leaving'] = nt['leaving']
            if ct:
                merged['cterm_abbr'] = ct['abbr']
                merged['cterm_leaving'] = ct['leaving']
            monomers.append(merged)
        return monomers
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)


_CAP_REACTIONS = None

def _load_cap_reactions():
    global _CAP_REACTIONS
    if _CAP_REACTIONS is None:
        import pathlib, yaml
        yaml_path = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'cap_reactions.yaml'
        if yaml_path.exists():
            with open(yaml_path, encoding="utf-8") as f:
                _CAP_REACTIONS = yaml.safe_load(f) or {}
        else:
            _CAP_REACTIONS = {}
    return _CAP_REACTIONS


_LG_SMILES = {
    'Cl': '[Cl]',
    'Br': '[Br]',
    'I': '[I]',
    'OH': '[OH]',
    'H': '[H]',
}


def _restore_reagent_form(mol, abbr):
    """Replace the bonding R-group dummy with the reagent's canonical leaving group.

    For an N-terminal cap (R2 bonding slot), this shows the actual reagent:
    e.g. ac with LG=Cl renders as CH₃C(=O)Cl (acetyl chloride).
    For C-terminal caps (R1 bonding slot), LG=H shows the free nucleophile.
    """
    from rdkit import Chem
    cap_rx = _load_cap_reactions()
    info = cap_rx.get(abbr)
    if not info:
        return None, None

    reagent_lg = info.get('reagent_lg', '')
    lg_smi = _LG_SMILES.get(reagent_lg)
    if not lg_smi:
        return None, None

    props = mol.GetPropsAsDict()
    rgroups = props.get('m_Rgroups', '')
    slots = [r.strip() for r in rgroups.split(',')]

    # Determine which slot is the bonding slot for this cap type
    # R2 caps (N-terminal/electrophilic): slot 2 bonds to the amine
    # R1 caps (C-terminal/nucleophilic): slot 1 bonds to the carbonyl
    # Sidechain caps: slot 1 bonds to the sidechain
    # We detect by checking which slot has the conventional LG ([OH] or [H])
    bonding_slot = None
    if len(slots) >= 2 and slots[1] in ('[OH]', '[H]'):
        bonding_slot = 2  # R2 cap
    elif len(slots) >= 1 and slots[0] in ('[OH]', '[H]'):
        bonding_slot = 1  # R1 cap
    else:
        bonding_slot = 2  # Default to R2

    emol = Chem.RWMol(Chem.RWMol(mol))
    to_remove = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        iso = atom.GetIsotope()
        if iso < 1 or iso > len(slots):
            continue
        slot_lg = slots[iso - 1]
        if slot_lg in ('None', 'none', ''):
            continue

        if iso == bonding_slot:
            # Replace with reagent LG
            if lg_smi == '[H]':
                nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
                if nb:
                    nb.SetNumExplicitHs(nb.GetNumExplicitHs() + 1)
                    nb.SetNoImplicit(False)
                to_remove.append(atom.GetIdx())
            elif lg_smi == '[OH]':
                new_atom = Chem.Atom(8)
                new_atom.SetNumExplicitHs(1)
                emol.ReplaceAtom(atom.GetIdx(), new_atom)
            elif lg_smi == '[Cl]':
                new_atom = Chem.Atom(17)
                emol.ReplaceAtom(atom.GetIdx(), new_atom)
            elif lg_smi == '[Br]':
                new_atom = Chem.Atom(35)
                emol.ReplaceAtom(atom.GetIdx(), new_atom)
            elif lg_smi == '[I]':
                new_atom = Chem.Atom(53)
                emol.ReplaceAtom(atom.GetIdx(), new_atom)
        else:
            # Non-bonding slots: restore normally
            if slot_lg == '[H]':
                nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
                if nb:
                    nb.SetNumExplicitHs(nb.GetNumExplicitHs() + 1)
                    nb.SetNoImplicit(False)
                to_remove.append(atom.GetIdx())
            elif slot_lg == '[OH]':
                new_atom = Chem.Atom(8)
                new_atom.SetNumExplicitHs(1)
                emol.ReplaceAtom(atom.GetIdx(), new_atom)
            else:
                lg_mol = Chem.MolFromSmiles(slot_lg, sanitize=True)
                if lg_mol:
                    new_atom = Chem.Atom(lg_mol.GetAtomWithIdx(0).GetAtomicNum())
                    emol.ReplaceAtom(atom.GetIdx(), new_atom)

    for idx in sorted(set(to_remove), reverse=True):
        emol.RemoveAtom(idx)

    try:
        Chem.SanitizeMol(emol)
    except Exception:
        pass

    reaction_type = info.get('reaction', '')
    reagent_note = info.get('reagent_note', '')
    issue = info.get('issue', '')
    meta = {
        'reaction': reaction_type,
        'reagent_lg': reagent_lg,
        'reagent_note': reagent_note,
    }
    if issue:
        meta['issue'] = issue
    return emol.GetMol(), meta


def _restore_leaving_groups(mol):
    """Replace dummy atoms with their leaving groups to produce the standalone monomer."""
    from rdkit import Chem
    props = mol.GetPropsAsDict()
    rgroups = props.get('m_Rgroups', '')
    slots = [r.strip() for r in rgroups.split(',')]

    emol = Chem.RWMol(Chem.RWMol(mol))
    to_remove = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        iso = atom.GetIsotope()
        if iso < 1 or iso > len(slots):
            continue
        lg_smi = slots[iso - 1]
        if lg_smi in ('None', 'none', ''):
            continue
        lg_mol = Chem.MolFromSmiles(lg_smi, sanitize=True)
        if lg_mol is None:
            continue
        if lg_smi == '[H]':
            nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
            if nb:
                nb.SetNumExplicitHs(nb.GetNumExplicitHs() + 1)
                nb.SetNoImplicit(False)
            to_remove.append(atom.GetIdx())
        elif lg_smi == '[OH]':
            new_atom = Chem.Atom(8)
            new_atom.SetNumExplicitHs(1)
            emol.ReplaceAtom(atom.GetIdx(), new_atom)
        else:
            new_atom = Chem.Atom(lg_mol.GetAtomWithIdx(0).GetAtomicNum())
            emol.ReplaceAtom(atom.GetIdx(), new_atom)

    for idx in sorted(set(to_remove), reverse=True):
        emol.RemoveAtom(idx)

    try:
        Chem.SanitizeMol(emol)
    except Exception:
        pass
    return emol.GetMol()


def _restore_one_slot(mol, keep_slot):
    """Restore all dummies except keep_slot; remove keep_slot's dummy (simulates bonding)."""
    from rdkit import Chem
    props = mol.GetPropsAsDict()
    rgroups = props.get('m_Rgroups', '')
    slots = [r.strip() for r in rgroups.split(',')]
    emol = Chem.RWMol(Chem.RWMol(mol))
    to_remove = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        iso = atom.GetIsotope()
        if iso < 1 or iso > len(slots):
            continue
        lg_smi = slots[iso - 1]
        if lg_smi in ('None', 'none', ''):
            continue
        if iso == keep_slot:
            nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
            if nb:
                nb.SetNoImplicit(False)
            to_remove.append(atom.GetIdx())
        elif lg_smi == '[H]':
            nb = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
            if nb:
                nb.SetNumExplicitHs(nb.GetNumExplicitHs() + 1)
                nb.SetNoImplicit(False)
            to_remove.append(atom.GetIdx())
        elif lg_smi == '[OH]':
            new_atom = Chem.Atom(8)
            new_atom.SetNumExplicitHs(1)
            emol.ReplaceAtom(atom.GetIdx(), new_atom)
        else:
            new_atom = Chem.Atom(Chem.MolFromSmiles(lg_smi).GetAtomWithIdx(0).GetAtomicNum())
            emol.ReplaceAtom(atom.GetIdx(), new_atom)
    for idx in sorted(set(to_remove), reverse=True):
        emol.RemoveAtom(idx)
    try:
        Chem.SanitizeMol(emol)
    except Exception:
        pass
    return emol.GetMol()


@app.get("/monomer_svg")
async def monomer_svg(abbr: str, width: int = 220, height: int = 180):
    try:
        from rdkit import Chem
        _all_mols, mol_by_abbr = _load_sdf()

        # Direct match
        if abbr in mol_by_abbr:
            mol = mol_by_abbr[abbr]
            svg = _draw_mol(mol, width, height)
            restored = _restore_leaving_groups(mol)
            svg_restored = _draw_mol(restored, width, height)
            result = {"svg": svg, "svg_restored": svg_restored}

            # Add reagent form if cap has reaction metadata
            reagent_mol, reagent_meta = _restore_reagent_form(mol, abbr)
            if reagent_mol is not None:
                result["svg_reagent"] = _draw_mol(reagent_mol, width, height)
                result["reagent"] = reagent_meta
            return result

        # Check if this is a degenerate base name (e.g. "Bn" -> Bn_/_Bn)
        # Detect pairs: look for abbr_ and _abbr variants
        nterm_key = abbr + '_'
        cterm_key = '_' + abbr
        variants_found = []
        if nterm_key in mol_by_abbr:
            variants_found.append(nterm_key)
        if cterm_key in mol_by_abbr:
            variants_found.append(cterm_key)

        if len(variants_found) >= 2:
            # This is a degenerate base name — render all variants
            panels = []
            for vkey in variants_found:
                vmol = mol_by_abbr[vkey]
                restored_v = _restore_leaving_groups(vmol)
                label = f"N-term ({vkey})" if vkey.endswith('_') else f"C-term ({vkey})"
                reagent_mol, reagent_meta = _restore_reagent_form(vmol, vkey)
                panel = {
                    "label": label,
                    "svg": _draw_mol(restored_v, width, height),
                }
                if reagent_mol is not None:
                    panel["svg_reagent"] = _draw_mol(reagent_mol, width, height)
                    panel["reagent"] = reagent_meta
                panels.append(panel)
            # R-group panel from first variant
            rgroup_svg = _draw_mol(mol_by_abbr[variants_found[0]], width, height)
            return {
                "degenerate": True,
                "variants": panels,
                "svg": rgroup_svg,
            }

        return JSONResponse({"error": f"Monomer '{abbr}' not found"}, status_code=404)
    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=500)


class _ConvertReq(BaseModel):
    cabiln: str
    target: str  # "bracket" or "branch"

@app.post("/convert_notation")
async def convert_notation(req: _ConvertReq):
    from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
    try:
        if req.target == 'bracket':
            return {"result": cabiln_to_bracket(req.cabiln)}
        elif req.target == 'branch':
            return {"result": cabiln_to_branch(req.cabiln)}
        else:
            return JSONResponse({"error": "target must be 'bracket' or 'branch'"},
                                status_code=400)
    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


@app.get("/monomer_rgroups")
async def monomer_rgroups(abbr: str, residue_idx: int = -1, cabiln: str = ''):
    """Return R-group info for a monomer, with used slots marked if in a sequence."""
    try:
        from rdkit import Chem
        _all_mols, mol_by_abbr = _load_sdf()

        target_mol = mol_by_abbr.get(abbr)
        if target_mol is None:
            return JSONResponse({"error": f"Monomer '{abbr}' not found"}, status_code=404)

        props = target_mol.GetPropsAsDict()
        chem_types_str = props.get('m_chem_types', '')
        rgroups_str = props.get('m_Rgroups', '')

        ct_map = {}
        for part in chem_types_str.split(','):
            part = part.strip()
            if ':' in part:
                slot_s, ct = part.split(':', 1)
                try:
                    ct_map[int(slot_s)] = ct.strip()
                except ValueError:
                    pass

        lg_list = [s.strip() for s in rgroups_str.split(',')]

        # Infer chem_types at runtime for slots missing from the SDF property
        from pyPept.interfaces.reaction_library import infer_chem_type
        from pyPept.sequence import _attachment_idx

        rgroups = []
        for i, lg in enumerate(lg_list):
            slot = i + 1
            if lg in ('None', 'none', ''):
                continue
            ct = ct_map.get(slot, '')
            if not ct:
                aidx = _attachment_idx(target_mol, slot)
                if aidx is not None:
                    ct = infer_chem_type(target_mol, aidx, slot=slot, leaving=lg)
            rgroups.append({
                'slot': slot,
                'chem_type': ct,
                'leaving': lg,
                'used': False
            })

        # If we have a sequence context, mark which R-groups are already used
        used_slots = set()
        if cabiln.strip() and residue_idx >= 0:
            try:
                from pyPept.sequence import Sequence, _slot_for_attachment
                seq = Sequence(cabiln)
                for bond in seq.s_bonds:
                    m1, atom1, m2, atom2 = bond[0], bond[1], bond[2], bond[3]
                    if len(bond) > 4:
                        s1, s2 = bond[4], bond[5]
                    else:
                        mol1 = seq.get_monomer(m1)['m_romol']
                        mol2 = seq.get_monomer(m2)['m_romol']
                        s1 = _slot_for_attachment(mol1, atom1)
                        s2 = _slot_for_attachment(mol2, atom2)
                    if m1 == residue_idx and s1 is not None:
                        used_slots.add(s1)
                    if m2 == residue_idx and s2 is not None:
                        used_slots.add(s2)
                for rg in rgroups:
                    if rg['slot'] in used_slots:
                        rg['used'] = True
            except Exception:
                pass

        svg = _draw_mol(target_mol, 180, 140, used_slots if used_slots else None)
        return {"svg": svg, "rgroups": rgroups, "abbr": abbr}

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=500)


class _InsertBondReq(BaseModel):
    cabiln: str
    host_residue_idx: int
    new_abbr: str
    r_host: int
    r_new: int
    target_residue_idx: int = -1  # -1 = add new monomer; >=0 = crosslink to existing

@app.post("/insert_bond")
async def insert_bond(req: _InsertBondReq):
    """Insert a bracket branch at the given host residue in CABILN notation."""
    import re

    def _rest_to_arms(s):
        """Rewrite rest_in as explicit [.arm] sub-brackets.

        Existing [.arm] groups pass through unchanged; any flat content between
        (or after) them is wrapped in its own [...].  This makes every arm from
        the hub unambiguous regardless of order.
        """
        parts = []
        i = 0
        flat_start = 0
        while i < len(s):
            if s[i] == '[':
                if i > flat_start:
                    parts.append('[' + s[flat_start:i] + ']')
                depth = 0
                j = i
                while j < len(s):
                    if s[j] == '[':
                        depth += 1
                    elif s[j] == ']':
                        depth -= 1
                        if depth == 0:
                            break
                    j += 1
                parts.append(s[i:j + 1])
                i = j + 1
                flat_start = i
            else:
                i += 1
        if flat_start < len(s):
            parts.append('[' + s[flat_start:] + ']')
        return ''.join(parts)
    try:
        cabiln = req.cabiln.strip()
        if not cabiln:
            return {"result": req.new_abbr}

        # ── Intramolecular crosslink: annotate two EXISTING residues ─────────
        if req.target_residue_idx >= 0 and req.target_residue_idx != req.host_residue_idx:
            from pyPept.sequence import Sequence as _Seq
            try:
                _seq2 = _Seq(cabiln)
                _chain_ids2 = _seq2.s_chains.get('s_monomerIDs', [])
                _main2 = set(_chain_ids2[0]) if _chain_ids2 else set()
            except Exception:
                _seq2 = None
                _main2 = set()

            # Pick the next unused !n tag
            _used_ns = {int(m.group(1)) for m in re.finditer(r'\.?!(\d+)', cabiln)}
            _n = 1
            while _n in _used_ns:
                _n += 1
            _tag = f'!{_n}'

            def _token_split(s):
                """Split CABILN by '-' respecting bracket depth."""
                toks, depth, cur = [], 0, ''
                for ch in s:
                    if ch == '[': depth += 1
                    elif ch == ']': depth -= 1
                    if ch == '-' and depth == 0:
                        toks.append(cur); cur = ''
                    else:
                        cur += ch
                if cur: toks.append(cur)
                return toks

            def _annotate_main(cabiln, res_idx, ann):
                """Append ann to the main-chain token for res_idx."""
                toks = _token_split(cabiln)
                midx = 0
                for i, tok in enumerate(toks):
                    if midx == res_idx:
                        toks[i] = tok + ann
                        return '-'.join(toks)
                    midx += 1
                return cabiln

            def _annotate_pendant(cabiln, res_idx, ann, seq, main_set):
                """Insert ann after the pendant residue's Abbr(r,r) in its bracket."""
                try:
                    abbr = seq.s_monomers[res_idx].get('m_abbr', '')
                except (IndexError, KeyError):
                    return cabiln
                pendant_indices = [i for i in range(len(seq.s_monomers))
                                   if i not in main_set]
                same_abbr = [i for i in pendant_indices
                             if seq.s_monomers[i].get('m_abbr', '') == abbr]
                try:
                    occ = same_abbr.index(res_idx)
                except ValueError:
                    occ = 0
                pat = re.compile(re.escape(abbr) + r'\(\d+,\d+\)')
                cnt = 0
                def _spans(s):
                    i = 0
                    while i < len(s):
                        if s[i] == '[':
                            depth = 0; j = i
                            while j < len(s):
                                if s[j] == '[': depth += 1
                                elif s[j] == ']':
                                    depth -= 1
                                    if depth == 0: break
                                j += 1
                            yield (i, j + 1); i += 1
                        else:
                            i += 1
                for bs, be in _spans(cabiln):
                    content = cabiln[bs + 1:be - 1]
                    for pm in pat.finditer(content):
                        pfix = content[:pm.start()]
                        if pfix.count('[') != pfix.count(']'):
                            continue
                        if cnt < occ:
                            cnt += 1; continue
                        new_content = content[:pm.end()] + ann + content[pm.end():]
                        return cabiln[:bs] + '[' + new_content + ']' + cabiln[be:]
                return cabiln

            ann_host   = f'.{_tag}({req.r_host},{req.r_new})'
            ann_target = f'.{_tag}({req.r_new},{req.r_host})'

            result = cabiln
            for _ridx, _ann in [(req.host_residue_idx, ann_host),
                                 (req.target_residue_idx, ann_target)]:
                if _main2 and _ridx in _main2:
                    result = _annotate_main(result, _ridx, _ann)
                elif _seq2 is not None:
                    result = _annotate_pendant(result, _ridx, _ann, _seq2, _main2)
                else:
                    result = _annotate_main(result, _ridx, _ann)

            return {"result": result}
        # ─────────────────────────────────────────────────────────────────────

        is_backbone = (req.r_host == 2 and req.r_new == 1)

        new_entry = f'.{req.new_abbr}({req.r_host},{req.r_new})'
        bracket = f'.[{req.new_abbr}({req.r_host},{req.r_new})]'

        # Detect bracket-internal host BEFORE the backbone early-return so that
        # pendant residues (e.g. a Cys arm from TBMB) can be extended correctly.
        from pyPept.sequence import Sequence
        try:
            seq = Sequence(cabiln)
            chain_ids = seq.s_chains.get('s_monomerIDs', [])
            main_set = set(chain_ids[0]) if chain_ids else set()
        except Exception:
            main_set = set()

        host_abbr = None
        target_occurrence = 0
        if req.host_residue_idx not in main_set and main_set:
            try:
                host_abbr = seq.s_monomers[req.host_residue_idx].get('m_abbr', '')
                # Which ordinal is this among all pendant monomers with the same abbr?
                # s_monomers order matches left-to-right CABILN order, so the ordinal
                # here maps directly to the ordinal we'll find during bracket scanning.
                pendant_indices = [i for i in range(len(seq.s_monomers))
                                   if i not in main_set]
                same_abbr = [i for i in pendant_indices
                             if seq.s_monomers[i].get('m_abbr', '') == host_abbr]
                target_occurrence = same_abbr.index(req.host_residue_idx)
            except (IndexError, KeyError, ValueError):
                pass

        # Backbone on a main-chain residue: simple append.
        # Backbone on a bracket-internal residue falls through to bracket logic.
        if is_backbone and not host_abbr:
            return {"result": cabiln + '-' + req.new_abbr}

        if host_abbr:
            # Host is a bracket-internal monomer — insert inside its bracket.
            # Scan ALL brackets at any nesting depth (char-by-char, not regex
            # finditer) so old-style [[hub.chain]] notation is handled by
            # finding the innermost bracket that directly contains hub.
            def _all_bracket_spans(s):
                i = 0
                while i < len(s):
                    if s[i] == '[':
                        depth = 0
                        j = i
                        while j < len(s):
                            if s[j] == '[':
                                depth += 1
                            elif s[j] == ']':
                                depth -= 1
                                if depth == 0:
                                    break
                            j += 1
                        yield (i, j + 1)
                        i += 1  # +1 so inner brackets are yielded too
                    else:
                        i += 1

            inserted = False
            occurrence_count = 0
            abbr_pat = re.compile(re.escape(host_abbr) + r'\(\d+,\d+\)')
            for br_start, br_end in _all_bracket_spans(cabiln):
                content = cabiln[br_start + 1:br_end - 1]
                # Iterate ALL top-level matches within this bracket so that two
                # same-abbr residues in the same sequential arm can both be targeted.
                for pm in abbr_pat.finditer(content):
                    prefix_in = content[:pm.start()]
                    open_count = prefix_in.count('[') - prefix_in.count(']')
                    if open_count != 0:
                        continue  # inside a nested bracket — keep looking
                    # Skip earlier occurrences until we reach the target.
                    if occurrence_count < target_occurrence:
                        occurrence_count += 1
                        continue
                    hub_in = content[pm.start():pm.end()]
                    rest_in = content[pm.end():]
                    if rest_in and not rest_in.startswith(']'):
                        # Rewrite existing content after hub as explicit arms,
                        # then add new arm: [hub[.existing...][.new]]
                        arms = _rest_to_arms(rest_in) + f'[{new_entry}]'
                        new_bracket = f'[{prefix_in}{hub_in}{arms}]'
                    else:
                        new_bracket = f'[{prefix_in}{hub_in}{new_entry}]'
                    # If the replaced bracket was the sole content of an outer
                    # [[hub...]] style bracket, unnest it now.
                    old_bracket = cabiln[br_start:br_end]
                    if (br_start > 0 and cabiln[br_start - 1] == '['
                            and br_end < len(cabiln) and cabiln[br_end] == ']'
                            and cabiln[br_start - 1:br_end + 1] == '[' + old_bracket + ']'):
                        result = cabiln[:br_start - 1] + new_bracket + cabiln[br_end + 1:]
                    else:
                        result = cabiln[:br_start] + new_bracket + cabiln[br_end:]
                    inserted = True
                    break
                if inserted:
                    break
            if inserted:
                return {"result": result}

            # Fallback: host uses inline-cap notation .HOST(r,r) without explicit
            # brackets — convert to bracket form and insert the new entry.
            inline_m = re.search(r'\.' + re.escape(host_abbr) + r'\(\d+,\d+\)', cabiln)
            if inline_m and (inline_m.start() == 0
                             or cabiln[inline_m.start() - 1] != '['):
                old_cap = inline_m.group(0)[1:]  # strip leading '.'
                result = (cabiln[:inline_m.start()]
                          + '.[' + old_cap + new_entry + ']'
                          + cabiln[inline_m.end():])
                return {"result": result}

        # Main-chain host: tokenize by '-' respecting bracket depth
        tokens = []
        depth = 0
        current = ''
        for ch in cabiln:
            if ch == '[':
                depth += 1
            elif ch == ']':
                depth -= 1
            if ch == '-' and depth == 0:
                tokens.append(current)
                current = ''
            else:
                current += ch
        if current:
            tokens.append(current)

        main_idx = 0
        insert_at = -1
        for i, tok in enumerate(tokens):
            base = re.split(r'\.\[|\.\(|\.!', tok)[0]
            if re.match(r'^!\d+$', base):
                main_idx += 1
                if main_idx - 1 == req.host_residue_idx:
                    insert_at = i
                continue
            if main_idx == req.host_residue_idx:
                insert_at = i
                break
            main_idx += 1

        if insert_at >= 0:
            tokens[insert_at] = tokens[insert_at] + bracket
        else:
            tokens[-1] = tokens[-1] + bracket

        return {"result": '-'.join(tokens)}

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=500)


class _ValidateBondReq(BaseModel):
    chem_type_a: str
    chem_type_b: str
    abbr_a: str = ''
    slot_a: int = 0
    abbr_b: str = ''
    slot_b: int = 0

@app.get("/examples")
async def get_examples():
    return JSONResponse(EXAMPLES)


@app.post("/validate_bond")
async def validate_bond(req: _ValidateBondReq):
    """Check if a bond between two R-group chemistry types is valid."""
    try:
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        entry = REACTION_INDEX.get((req.chem_type_a, req.chem_type_b))
        if entry is None:
            entry = REACTION_INDEX.get((req.chem_type_b, req.chem_type_a))
        if entry:
            return {"valid": True, "reaction": entry.get('id', 'unknown'),
                    "description": entry.get('description', '')}
        else:
            return {"valid": False,
                    "reason": f"No reaction for {req.chem_type_a} + {req.chem_type_b}"}
    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=500)


@app.get("/reactions")
async def list_reactions():
    """Return all valid (chem_type_a, chem_type_b) pairs from the reaction index."""
    try:
        from pyPept.interfaces.reaction_library import REACTION_INDEX
        pairs = [list(pair) for pair in REACTION_INDEX.keys()]
        return JSONResponse(pairs)
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)


def _build_bracket_groups(seq, chain_ids, cabiln='', crosslink_groups=None):
    """Identify bracket branch groups and their host monomers.

    Each separate attachment (``.[A.B.C]``, ``.[D]``, ``.cap(r,r)``) on a
    host residue becomes its own group.  Groups are split by connected
    components among the branch residues — branch residues that bond to
    each other stay together; those that only bond to the host get their
    own group.

    Returns list of ``{host, members}``.
    """
    if not chain_ids or len(chain_ids) < 2:
        return []
    main_set = set(chain_ids[0])
    branch_set = set()
    for ci in range(1, len(chain_ids)):
        branch_set.update(chain_ids[ci])
    if not branch_set:
        return []

    xlink_pairs = set()
    for g in (crosslink_groups or []):
        if len(g['members']) == 2:
            a, b = g['members']
            xlink_pairs.add((a, b))
            xlink_pairs.add((b, a))

    host_of = {}
    branch_adj: dict[int, set] = {idx: set() for idx in branch_set}
    for bond in seq.s_bonds:
        m1, m2 = bond[0], bond[2]
        if (m1, m2) in xlink_pairs:
            continue
        if m1 in main_set and m2 in branch_set:
            host_of[m2] = m1
        elif m2 in main_set and m1 in branch_set:
            host_of[m1] = m2
        elif m1 in branch_set and m2 in branch_set:
            branch_adj[m1].add(m2)
            branch_adj[m2].add(m1)

    # Propagate host through inter-branch bonds
    changed = True
    while changed:
        changed = False
        for idx in branch_set:
            if idx in host_of:
                continue
            for peer in branch_adj.get(idx, set()):
                if peer in host_of:
                    host_of[idx] = host_of[peer]
                    changed = True
                    break

    # Connected components among branch residues (union-find)
    parent = {idx: idx for idx in branch_set}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for idx, peers in branch_adj.items():
        for p in peers:
            union(idx, p)

    from collections import defaultdict
    components = defaultdict(list)
    for idx in branch_set:
        if idx in host_of:
            components[(host_of[idx], find(idx))].append(idx)

    return [{"host": key[0], "members": sorted(members)}
            for key, members in components.items()]


def _build_crosslink_groups(seq):
    """Extract crosslink !n pairs from the parsed BILN."""
    import re
    biln = seq.s_biln
    chains = biln.split('.')
    m_idx = 0
    tag_to_monomers = {}
    for chain in chains:
        residues = chain.split('-')
        for res in residues:
            for m in re.finditer(r'\((!\w+),\d+\)', res):
                tag = m.group(1)
                tag_to_monomers.setdefault(tag, []).append(m_idx)
            m_idx += 1
    return [{"tag": tag, "members": mems}
            for tag, mems in tag_to_monomers.items()
            if len(mems) == 2]


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

        w, h  = max(400, req.width), max(300, req.height)
        svg   = _draw_mol(romol, w, h, seed=req.seed)
        block = _mol_block(romol)

        res_map = mol.get_residue_atom_map()
        residues = [{"idx": i, "abbr": m.get("m_abbr", f"?{i}")}
                    for i, m in enumerate(seq.s_monomers)]

        chain_ids = seq.s_chains.get('s_monomerIDs', [])
        chains = [{"idx": ci, "residues": ids} for ci, ids in enumerate(chain_ids)]

        crosslink_groups = _build_crosslink_groups(seq)
        bracket_groups = _build_bracket_groups(seq, chain_ids, req.cabiln, crosslink_groups)

        return {"svg": svg, "mol_block": block,
                "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
                "residue_map": {str(k): v for k, v in res_map.items()},
                "residues": residues,
                "chains": chains,
                "bracket_groups": bracket_groups,
                "crosslink_groups": crosslink_groups,
                "cabiln_echo": req.cabiln}

    except (Exception, SystemExit) as exc:
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


class _ReferenceReq(BaseModel):
    input: str
    width: int = 960
    height: int = 680


@app.post("/render_reference")
async def render_reference(req: _ReferenceReq):
    """Auto-detect input format (SMILES, old BILN, HELM, CABILN) and render."""
    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt

    txt = req.input.strip()
    w, h = max(400, req.width), max(300, req.height)
    romol = None
    fmt = None

    romol = Chem.MolFromSmiles(txt)
    if romol is not None:
        fmt = 'SMILES'

    if romol is None and 'PEPTIDE' in txt.upper() and '$' in txt:
        try:
            import sys, pathlib
            sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
            from pyPept.converter import Converter
            from pyPept.sequence import Sequence
            from pyPept.molecule import Molecule
            conv = Converter(helm=txt)
            biln = conv.get_biln()
            seq = Sequence(biln, fmt='biln')
            mol = Molecule(seq)
            romol = mol.get_molecule(fmt='ROMol')
            fmt = 'HELM'
        except Exception:
            pass

    if romol is None:
        try:
            from pyPept.sequence import Sequence
            from pyPept.molecule import Molecule
            seq = Sequence(txt, fmt='biln')
            mol = Molecule(seq)
            romol = mol.get_molecule(fmt='ROMol')
            fmt = 'BILN'
        except Exception:
            pass

    if romol is None:
        try:
            from pyPept.sequence import Sequence
            from pyPept.molecule import Molecule
            seq = Sequence(txt)
            mol = Molecule(seq)
            romol = mol.get_molecule(fmt='ROMol')
            fmt = 'CABILN'
        except Exception:
            pass

    if romol is None:
        return JSONResponse(
            {"error": "Could not parse as SMILES, BILN, HELM, or CABILN"},
            status_code=400)

    svg = _draw_mol(romol, w, h)
    smiles = Chem.MolToSmiles(romol)
    return {"svg": svg, "smiles": smiles, "format": fmt,
            "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}"}


class _MolBlockReq(BaseModel):
    mol_block: str
    width: int = 960
    height: int = 680


@app.post("/render_mol")
async def render_mol(req: _MolBlockReq):
    try:
        from rdkit import Chem
        from rdkit.Chem.Descriptors import ExactMolWt

        romol = Chem.MolFromMolBlock(req.mol_block)
        if romol is None:
            return JSONResponse({"error": "Invalid .mol data"}, status_code=400)

        w, h = max(400, req.width), max(300, req.height)
        svg = _draw_mol(romol, w, h)
        smiles = Chem.MolToSmiles(romol)
        return {"svg": svg, "smiles": smiles,
                "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}"}

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

        seq        = Sequence(req.cabiln)
        mol        = Molecule(seq)
        cabiln_mol = mol.get_molecule(fmt='ROMol')
        if cabiln_mol is None:
            return JSONResponse({"error": "CABILN assembly produced no molecule"}, status_code=400)

        smi_canon    = _canon(smiles_mol)
        cabiln_canon = _canon(cabiln_mol)
        return {
            "match":            smi_canon == cabiln_canon,
            "smiles_canonical": smi_canon,
            "cabiln_canonical": cabiln_canon,
        }

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


@app.post("/preview_monomer")
async def preview_monomer(req: _PreviewReq):
    try:
        import sys, pathlib
        sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
        from pyPept.interfaces.monomer_pipeline import pre_activate, ActivationError

        result = pre_activate(req.smiles)

        from rdkit import Chem
        mol = Chem.MolFromSmiles(result.chuckles)
        if mol is None:
            return JSONResponse({"error": "Generated CHUCKLES is invalid"}, status_code=400)

        svg = _draw_mol(mol, req.width, req.height)
        return {
            "chuckles":   result.chuckles,
            "chem_types": {str(k): v for k, v in result.chem_types.items()},
            "leaving":    {str(k): v for k, v in result.leaving.items()},
            "svg":        svg,
        }

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


@app.post("/register_monomer")
async def register_monomer(req: _RegisterReq):
    try:
        import sys, pathlib
        sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
        from rdkit import Chem
        from rdkit.Chem import SDWriter, rdDepictor

        sdf_path = _SDF_PATH

        # Check for duplicate abbreviation
        _all_mols, mol_by_abbr = _load_sdf()
        if req.abbr in mol_by_abbr:
            return JSONResponse({"error": f"Abbreviation '{req.abbr}' already exists in library"}, status_code=400)

        mol = Chem.MolFromSmiles(req.chuckles)
        if mol is None:
            return JSONResponse({"error": "Invalid CHUCKLES SMILES"}, status_code=400)

        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(mol)

        # Format m_Rgroups (6 positions)
        rgroups_list = ['None'] * 6
        for slot_s, lg in req.leaving.items():
            slot = int(slot_s) - 1
            if 0 <= slot < 6:
                rgroups_list[slot] = lg
        rgroups_str = ','.join(rgroups_list)

        # Format m_chem_types
        chem_types_str = ','.join(
            f"{slot}:{ct}" for slot, ct in
            sorted(req.chem_types.items(), key=lambda x: int(x[0]))
        )

        mol.SetProp('m_abbr',      req.abbr)
        mol.SetProp('symbol',      req.abbr)
        mol.SetProp('m_name',      req.name)
        mol.SetProp('m_type',      req.type)
        mol.SetProp('m_subtype',   req.subtype)
        mol.SetProp('m_Rgroups',   rgroups_str)
        mol.SetProp('m_chem_types', chem_types_str)

        with open(sdf_path, 'a') as fh:
            writer = SDWriter(fh)
            writer.SetKekulize(False)
            writer.write(mol)
            writer.close()

        _invalidate_sdf()
        all_mols, _ = _load_sdf()
        return {"ok": True, "total": len(all_mols)}

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


# ── launch ────────────────────────────────────────────────────────────────────

def _open_browser():
    time.sleep(0.9)
    webbrowser.open("http://localhost:8732")


if __name__ == "__main__":
    import os as _os
    _port = int(_os.environ.get("PORT", 8732))
    print("CABILN Live Renderer")
    import time as _t; _t0 = _t.perf_counter()
    _warm_mols, _ = _load_sdf()
    print(f"  SDF cache warmed: {len(_warm_mols)} monomers in {_t.perf_counter()-_t0:.1f}s")
    print(f"  Local     ->  http://localhost:{_port}")
    print("  Tailscale ->  http://100.119.0.78:8732")
    if _port == 8732:
        threading.Thread(target=_open_browser, daemon=True).start()
    uvicorn.run(app, host="0.0.0.0", port=_port, log_level="warning")
