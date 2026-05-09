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

# ── render cache (LRU, capped at 200 entries) ─────────────────────────────────
from collections import OrderedDict as _OD
_render_cache: _OD = _OD()
_RENDER_CACHE_MAX = 200

def _to_bracket(cabiln: str) -> str:
    """Normalise positional-% branch notation to bracket form for Sequence.

    Pure crosslink branches (%TBMB.!1.!2.!3, %C20FA-AEEA-E_g.!1, etc.) are
    already parseable by Sequence and are left untouched.  Only branches where
    every segment lacks an !n marker need conversion.
    """
    if '%' not in cabiln:
        return cabiln
    branch_parts = cabiln.split('%')[1:]
    if all('!' in p for p in branch_parts):
        return cabiln          # all crosslink — Sequence handles it
    from pyPept.sequence import cabiln_to_bracket
    return cabiln_to_bracket(cabiln)


def _rc_get(key):
    if key in _render_cache:
        _render_cache.move_to_end(key)
        return _render_cache[key]
    return None

def _rc_put(key, value):
    _render_cache[key] = value
    _render_cache.move_to_end(key)
    if len(_render_cache) > _RENDER_CACHE_MAX:
        _render_cache.popitem(last=False)

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
                "name": "Semaglutide",
                "description": "GLP-1 agonist · 31 AA · C18 diacid via AEEA–AEEA–γGlu on K26 · Ozempic/Wegovy",
                "cabiln": "H-Aib-E-G-T-F-T-S-D-V-S-S-Y-L-E-G-Q-A-A-K.[AEEA(4,2).AEEA(1,2).E_g(1,2).C18FA(1,2)]-E-F-I-A-W-L-V-R-G-R-G",
            },
            {
                "name": "Liraglutide",
                "description": "GLP-1 agonist · 31 AA · C16 palmitoyl via γGlu on K26 · Victoza/Saxenda",
                "cabiln": "H-A-E-G-T-F-T-S-D-V-S-S-Y-L-E-G-Q-A-A-K.[E(4,4).Pal(1,2)]-E-F-I-A-W-L-V-R-G-R-G",
            },
            {
                "name": "Tirzepatide",
                "description": "GIP/GLP-1 dual agonist · 39 AA · C20 diacid via AEEA–AEEA–γGlu on K20 · Mounjaro/Zepbound",
                "cabiln": "Y-Aib-E-G-T-F-T-S-D-Y-S-I-Aib-L-D-K-I-A-Q-K.[AEEA(4,2).AEEA(1,2).E_g(1,2).C20FA(1,2)]-A-F-V-Q-W-L-I-A-G-G-P-S-S-G-A-P-P-P-S-am",
            },
            {
                "name": "Retatrutide",
                "description": "GIP/GLP-1/GCG triple agonist · 39 AA · C20 diacid via AEEA–γGlu on K17 · LY3437943",
                "cabiln": "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K-K.[AEEA(4,2).E_g(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am",
            },
            {
                "name": "Exenatide",
                "description": "GLP-1 agonist · 39 AA · exendin-4 · amidated C-terminus · no lipid · Byetta",
                "cabiln": "H-G-E-G-T-F-T-S-D-L-S-K-Q-M-E-E-E-A-V-R-L-F-I-E-W-L-K-N-G-G-P-S-S-G-A-P-P-P-S-am",
            },
            {
                "name": "Lixisenatide",
                "description": "GLP-1 agonist · 44 AA · exendin-4 + 6×Lys C-terminal extension · amidated · Adlyxin",
                "cabiln": "H-G-E-G-T-F-T-S-D-L-S-K-Q-M-E-E-E-A-V-R-L-F-I-E-W-L-K-N-G-G-P-S-S-G-A-P-P-S-K-K-K-K-K-K-am",
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
                "cabiln": "!1-A-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1",
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
                "cabiln": "ac-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[E_g(4,4).AEEA(1,2).C20FA(1,2)]-am",
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
                "name": "RGD pendant arm",
                "description": "Aspartate isopeptide — RGD integrin-targeting arm via β-carboxyl",
                "cabiln": "ac-G-D.[R(4,1).G(2,1).D(2,1).am(2,1)]-S-am",
            },
            {
                "name": "Dual-E isopeptide scaffold",
                "description": "Two glutamate γ-carboxyl isopeptide arms — asymmetric dual-branch scaffold",
                "cabiln": "ac-E.[G(4,1).A(2,1).am(2,1)]-G-E.[A(4,1).G(2,1).am(2,1)]-am",
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
  #notation-select {
    flex-shrink: 0;
    background: #0d1e38;
    border: 1px solid #2a4070;
    color: #7aaeff;
    font-size: 10px;
    padding: 2px 5px;
    border-radius: 4px;
    cursor: pointer;
    outline: none;
  }
  #notation-select:focus { border-color: #5a9ae0; }
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
  .statusbar.ok   { color: #3dbe6c; }
  .statusbar.warn { color: #f0a030; }

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
  #compare-bar .match      { color: #3dbe6c; font-weight: 600; }
  #compare-bar .nomatch    { color: #d9534f; font-weight: 600; }
  #compare-bar .warn-badge { color: #f0a030; font-weight: 600; font-size: 11px; cursor: help; }
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
    <select id="notation-select" title="Input notation">
      <option value="cabiln" selected>CABILN</option>
      <option value="smiles">SMILES</option>
      <option value="biln">BILN</option>
      <option value="helm">HELM</option>
    </select>
    <textarea id="cabiln-input"
              placeholder="e.g.  fmoc-A-G-L-am&#10;fmoc-C.trt(4,1)-A-K.boc(4,1)-am&#10;fmoc-K.!1(4,1)-G-G-E.!1-am"
              spellcheck="false" autocomplete="off"></textarea>
    <button id="btn-to-cabiln" class="hbtn green" style="display:none;font-size:10px;padding:2px 10px;white-space:nowrap;" title="Convert to CABILN">→ CABILN</button>
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
  <div id="build-insert-row" style="display:none;align-items:center;gap:10px;padding:6px 2px 0;">
    <span id="build-insert-info" style="font-size:10px;color:#5a8ab0;flex:1;"></span>
    <button id="build-insert-btn" class="hbtn green" style="font-size:10px;padding:2px 10px;">⊕ Insert Between</button>
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
        <button id="btn-s2c" class="hbtn green" style="font-size:10px;padding:2px 10px;" title="Convert SMILES → % CABILN (branch form)">→ %</button>
        <button id="btn-s2c-bracket" class="hbtn green" style="font-size:10px;padding:2px 10px;" title="Convert SMILES → [] CABILN (bracket form)">→ []</button>
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
let cabilnAbort = null;
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
let insertBetweenActive = false;
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
const buildInsertRow = document.getElementById('build-insert-row');
const buildInsertInfo = document.getElementById('build-insert-info');
const buildInsertBtn = document.getElementById('build-insert-btn');
const btnReroll     = document.getElementById('btn-reroll');
const btnS2c        = document.getElementById('btn-s2c');
const btnS2cBracket = document.getElementById('btn-s2c-bracket');
const btnRxnFilter  = document.getElementById('btn-rxn-filter');
const notationSelect = document.getElementById('notation-select');
const btnToCabiln   = document.getElementById('btn-to-cabiln');

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

  if (insertBetweenActive) {
    filtered = filtered.filter(m => {
      const raw = m.chem_types || '';
      return raw.includes('backbone_n') && raw.includes('backbone_c');
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
      if (buildMode && insertBetweenActive) {
        doInsertBetween(row.dataset.abbr);
      } else if (buildMode) {
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
let currentBranchSet = new Set();
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
      if (nTermXlinkTags.has(g.tag) && g.members[0] === rIdx)
        resChips.appendChild(makeXlinkChip(g.tag, g.members));
    });
  }

  function appendXlinks(rIdx) {
    const xlinks = xlinkByMember[rIdx];
    if (!xlinks) return;
    xlinks.forEach(g => {
      if (!nTermXlinkTags.has(g.tag) || g.members[0] !== rIdx)
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
  currentBranchSet = new Set(branchSet);

  function expandWithXlinks(idxList) {
    const out = new Set(idxList);
    for (const m of idxList) {
      const xls = xlinkByMember[m];
      if (xls) xls.forEach(g => g.members.forEach(x => out.add(x)));
    }
    return [...out];
  }

  function appendResidueWithBrackets(rIdx) {
    if (branchSet.has(rIdx)) return;
    prependXlinks(rIdx);
    const chip = makeChip(rIdx);
    if (chip) resChips.appendChild(chip);
    const hostGroups = groupsByHost[rIdx];
    if (hostGroups) {
      // $ chip: host + bracket members + their xlink partners
      const allMembers = expandWithXlinks([rIdx, ...hostGroups.flat()]);
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
        // bracket [ ] highlight: members + their xlink partners
        const brkMembers = expandWithXlinks(members);
        resChips.appendChild(makeSeparator(brk[0], brkMembers));
        members.forEach(mIdx => {
          // bracket member chip: simple hover (self only)
          const mc = makeChip(mIdx, true);
          if (mc) resChips.appendChild(mc);
          appendXlinks(mIdx);
        });
        resChips.appendChild(makeSeparator(brk[1], brkMembers));
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
          // $ selects host + bracket members + their xlink partners
          const bracketMembers = expandWithXlinks([rIdx, ...hostGroups.flat()]);
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
            const brkMembers = expandWithXlinks(members);
            resChips.appendChild(makeSeparator(brk[0], brkMembers));
            members.forEach(mIdx => {
              const mc = makeChip(mIdx, true);
              if (mc) resChips.appendChild(mc);
              appendXlinks(mIdx);
            });
            resChips.appendChild(makeSeparator(brk[1], brkMembers));
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
  insertBetweenActive = false;
  buildInsertRow.style.display = 'none';
  buildInsertBtn.textContent = '⊕ Insert Between';
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
  resChips.querySelectorAll('.res-chip').forEach(c => c.style.outline = '');
  btnRxnFilter.disabled = true;
  if (rxnFilterActive) {
    rxnFilterActive = false;
    btnRxnFilter.classList.remove('active');
    if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  }
}

btnBuild.addEventListener('click', () => buildMode ? closeBuild() : openBuild());
buildClose.addEventListener('click', closeBuild);

function checkAdjacentMainChain() {
  if (buildLeftRIdx == null || buildRightRIdx == null) return false;
  if (buildLeftRIdx === buildRightRIdx) return false;
  const main = chainData.length ? chainData[0].residues : [];
  if (!main.includes(buildLeftRIdx) || !main.includes(buildRightRIdx)) return false;
  if (currentBranchSet.has(buildLeftRIdx) || currentBranchSet.has(buildRightRIdx)) return false;
  const posL = main.indexOf(buildLeftRIdx);
  const posR = main.indexOf(buildRightRIdx);
  return Math.abs(posL - posR) === 1;
}

function updateInsertBetweenUI() {
  if (buildMode && checkAdjacentMainChain()) {
    const la = (residueList.find(r => r.idx === buildLeftRIdx) || {}).abbr || '?';
    const ra = (residueList.find(r => r.idx === buildRightRIdx) || {}).abbr || '?';
    buildInsertInfo.textContent = `${la} and ${ra} are adjacent on the backbone`;
    buildInsertRow.style.display = 'flex';
  } else {
    buildInsertRow.style.display = 'none';
    if (insertBetweenActive) {
      insertBetweenActive = false;
      buildInsertBtn.textContent = '⊕ Insert Between';
      if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
    }
  }
}

buildInsertBtn.addEventListener('click', () => {
  if (insertBetweenActive) {
    insertBetweenActive = false;
    buildInsertBtn.textContent = '⊕ Insert Between';
    buildHint.textContent = 'Select a chip on the sequence, then right-click a monomer in the library';
    if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
    return;
  }
  if (!checkAdjacentMainChain()) return;
  insertBetweenActive = true;
  buildInsertBtn.textContent = '✕ Cancel';
  buildHint.textContent = 'Click a backbone monomer in the library to insert between the selected residues';
  if (!libPanel.classList.contains('open')) openLib();
  if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
});

async function doInsertBetween(abbr) {
  const main = chainData.length ? chainData[0].residues : [];
  const posL = main.indexOf(buildLeftRIdx);
  const posR = main.indexOf(buildRightRIdx);
  const after_idx = posL < posR ? buildLeftRIdx : buildRightRIdx;
  const val = cabilnInput.value.trim();
  buildHint.textContent = 'Inserting…';
  try {
    const res = await fetch('/insert_backbone', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({ cabiln: val, after_idx, new_abbr: abbr })
    });
    const data = await res.json();
    if (data.error) {
      buildHint.textContent = 'Insert failed: ' + data.error;
      return;
    }
    cabilnInput.value = data.result;
    cabilnInput.dispatchEvent(new Event('input'));
    // Reset build state for next operation
    buildLeft = null; buildRight = null; buildLeftRIdx = null; buildRightRIdx = null;
    insertBetweenActive = false;
    buildInsertRow.style.display = 'none';
    buildInsertBtn.textContent = '⊕ Insert Between';
    buildLeftAbbr.textContent = '—';
    buildLeftSvg.innerHTML = '<div class="box-placeholder">Click a chip above</div>';
    buildLeftRg.innerHTML = '';
    buildRightAbbr.textContent = '—';
    buildRightSvg.innerHTML = '<div class="box-placeholder">Right-click from library</div>';
    buildRightRg.innerHTML = '';
    buildConnect.disabled = true;
    buildHint.textContent = `${abbr} inserted — select chips to continue building`;
    resChips.querySelectorAll('.res-chip').forEach(c => c.style.outline = '');
    if (libLoaded) renderLibList(libSearch.value.trim().toLowerCase());
  } catch (e) {
    buildHint.textContent = 'Insert failed';
  }
}
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
    updateInsertBetweenUI();
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
    updateInsertBetweenUI();
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

// ─── SMILES → CABILN conversion ───────────────────────────────────────────────
async function doS2c(notation) {
  const smiles = smilesInput.value.trim();
  if (!smiles) return;
  const btn = notation === 'bracket' ? btnS2cBracket : btnS2c;
  const origText = btn.textContent;
  btn.textContent = '…';
  btn.disabled = true;
  try {
    const res = await fetch('/smiles_to_cabiln', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles, notation })
    });
    const data = await res.json();
    if (data.error) {
      smilesStatus.textContent = 'S2C: ' + data.error;
      smilesStatus.className = 'statusbar';
    } else {
      cabilnInput.value = data.cabiln;
      cabilnInput.dispatchEvent(new Event('input'));
      if (data.warning) {
        smilesStatus.textContent = `⚠ ${data.warning}`;
        smilesStatus.className = 'statusbar warn';
      } else {
        smilesStatus.textContent = `Converted (${notation}): ${data.details.length} residue(s)`;
        smilesStatus.className = 'statusbar ok';
      }
    }
  } catch (e) {
    smilesStatus.textContent = 'S2C error — is server running?';
    smilesStatus.className = 'statusbar';
  } finally {
    btn.textContent = origText;
    btn.disabled = false;
  }
}
btnS2c.addEventListener('click', () => doS2c('percent'));
btnS2cBracket.addEventListener('click', () => doS2c('bracket'));

// ─── main input → CABILN convert button ──────────────────────────────────────
btnToCabiln.addEventListener('click', async () => {
  const txt = cabilnInput.value.trim();
  if (!txt) return;
  const origText = btnToCabiln.textContent;
  btnToCabiln.textContent = '…';
  btnToCabiln.disabled = true;
  try {
    const res = await fetch('/to_cabiln', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt })
    });
    const data = await res.json();
    if (data.error) {
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
    } else {
      notationSelect.value = 'cabiln';
      btnToCabiln.style.display = 'none';
      cabilnInput.placeholder = NOTATION_PLACEHOLDER.cabiln;
      cabilnInput.value = data.cabiln;
      cabilnInput.dispatchEvent(new Event('input'));
      cabilnStatus.textContent = `Converted from ${data.from}: ${data.cabiln}`;
      cabilnStatus.className = 'statusbar ok';
    }
  } catch (e) {
    cabilnStatus.textContent = '→ CABILN error — is server running?';
    cabilnStatus.className = 'statusbar';
  } finally {
    btnToCabiln.textContent = origText;
    btnToCabiln.disabled = false;
  }
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

// ─── notation selector ────────────────────────────────────────────────────────
const NOTATION_PLACEHOLDER = {
  cabiln: 'e.g.  fmoc-A-G-L-am\nfmoc-C.trt(4,1)-A-K.boc(4,1)-am\nfmoc-K.!1(4,1)-G-G-E.!1-am',
  smiles: 'Paste SMILES here… e.g. O=C1CNC(=O)[C@@H](C)N1',
  biln:   'Paste BILN here… e.g. fmoc-A-G-L-am  (use Token(bid,rg) for crosslinks)',
  helm:   'Paste HELM here… e.g. PEPTIDE1{A.G.L}$$$$',
};

notationSelect.addEventListener('change', () => {
  const mode = notationSelect.value;
  cabilnInput.placeholder = NOTATION_PLACEHOLDER[mode] || '';
  btnToCabiln.style.display = mode === 'cabiln' ? 'none' : '';
  cabilnInput.value = '';
  cabilnInput.className = '';
  resetCabiln();
});

// ─── CABILN render ────────────────────────────────────────────────────────────
cabilnInput.addEventListener('input', () => {
  clearTimeout(cabilnTimer);
  const seq = cabilnInput.value.trim();
  const mode = notationSelect.value;
  if (mode !== 'cabiln') {
    // Non-CABILN: render via render_reference for a live preview
    if (!seq) { resetCabiln(); return; }
    showSpinner(renderInner);
    cabilnTimer = setTimeout(() => doRenderForeign(seq), 400);
    return;
  }
  if (seq === lastCabiln) return;
  rerollSeed = 0;
  btnReroll.textContent = '⟳ Layout';
  if (!seq) { resetCabiln(); return; }
  showSpinner(renderInner);
  cabilnTimer = setTimeout(() => doRenderCabiln(seq), 300);
});

async function doRenderForeign(txt) {
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res = await fetch('/render_reference', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input: txt, width: w, height: h })
    });
    const data = await res.json();
    if (data.error) {
      setInner(renderInner, `<div class="placeholder err">${escHtml(data.error)}</div>`);
      cabilnStatus.textContent = data.error;
      cabilnStatus.className = 'statusbar';
      cabilnInput.className = 'err';
    } else {
      setInner(renderInner, data.svg);
      cabilnStatus.textContent = `${data.format}: ${data.info || ''}`;
      cabilnStatus.className = 'statusbar ok';
      cabilnInput.className = 'ok';
    }
  } catch (e) {
    cabilnStatus.textContent = 'Server error';
    cabilnStatus.className = 'statusbar';
  }
}

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
  // Abort any in-flight render so stale responses never overwrite newer ones.
  if (cabilnAbort) cabilnAbort.abort();
  cabilnAbort = new AbortController();
  const signal = cabilnAbort.signal;

  const _seqChanged = seq !== lastCabiln;
  lastCabiln = seq;
  if (buildMode && _seqChanged) clearBuild();
  resChips.innerHTML = '';
  const { w, h } = canvasSize(renderCanvas);
  try {
    const res  = await fetch('/render', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cabiln: seq, width: w, height: h, seed: rerollSeed }),
      signal
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
    if (e.name === 'AbortError') return;  // superseded by newer input — discard silently
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
      ? '<span class="match">✓ MATCH</span>'
      : '<span class="nomatch">✗ MISMATCH</span>';
    const warnBadge = data.warning
      ? `<span class="warn-badge" title="${escAttr(data.warning)}">⚠ achiral input</span>`
      : '';
    compareBar.innerHTML =
      badge + warnBadge +
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
    from rdkit.Chem.Draw import rdMolDraw2D
    if seed == 0:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(romol)
    elif seed % 2 == 1:
        if _indigo_layout(romol) is None:
            romol = _best_layout(romol, base_seed=seed * 6)
    else:
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


def _count_defined_stereo(mol) -> int:
    """Count atoms with explicitly defined chirality (not CHI_UNSPECIFIED)."""
    from rdkit import Chem
    return sum(
        1 for a in mol.GetAtoms()
        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    )


def _canon_flat(romol) -> str:
    """Canonical SMILES with all stereo stripped."""
    from rdkit import Chem
    rw = Chem.RWMol(romol)
    Chem.RemoveStereochemistry(rw)
    return Chem.MolToSmiles(rw, canonical=True)


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


# ── SMILES → CABILN algorithm (ported from cyclicpeptide, MIT) ───────────────
#
# Detects the peptide backbone (NCC=O repeating unit), isolates each residue by
# cutting backbone peptide bonds and capping the C-terminus with -OH, then
# substructure-matches the isolated residues against the pyPept SDF monomer
# library to produce a CABILN sequence string.

import re as _re

def _s2c_normalize(mol):
    """Normalize resonance forms (e.g. guanidinium) via InChI round-trip."""
    try:
        from rdkit.Chem.inchi import MolToInchi, MolFromInchi
        from rdkit import Chem as _Chem
        inchi = MolToInchi(mol)
        if inchi:
            m = MolFromInchi(inchi)
            if m:
                _Chem.SanitizeMol(m)
                return m
    except Exception:
        pass
    return mol


def _s2c_cap_smiles(smi, rgroups_str):
    """Replace [n*] dummy atoms with caps: [OH] slot → O atom; everything else → [H]."""
    rgroups = [r.strip() for r in rgroups_str.split(',')]
    def _rep(m):
        n = int(m.group(1))
        rg = rgroups[n - 1] if n - 1 < len(rgroups) else 'None'
        return 'O' if rg == '[OH]' else '[H]'
    return _re.sub(r'\[(\d+)\*\]', _rep, smi)


def _s2c_detect_backbone(mol):
    from rdkit import Chem as _C
    n = len(mol.GetSubstructMatches(_C.MolFromSmiles('NCC=O')))
    if n == 0:
        return '', None
    # Branched peptides (e.g. K with E_g lipid linker) have extra NCC=O units
    # in side-chain amino acids; scan from n down to 1 for longest real backbone.
    for length in range(n, 0, -1):
        matches = mol.GetSubstructMatches(_C.MolFromSmiles('C(=O)CN' * length))
        if matches:
            bb = matches[0]
            bb_idx = list(bb)
            bb_idx.reverse()
            return bb, bb_idx
    return '', None


def _s2c_order_backbone(m, bb_idx):
    from rdkit import Chem as _C
    id_list = bb_idx[:]
    for idx in [a.GetIdx() for a in m.GetAtoms()]:
        if idx not in id_list:
            id_list.append(idx)
    m_renum = _C.RenumberAtoms(m, newOrder=id_list)
    _, new_bb_idx = _s2c_detect_backbone(m_renum)
    return m_renum, new_bb_idx


def _s2c_side_chain_neighbors(m, atoms, origin):
    out = set()
    for ai in atoms:
        for nb in m.GetAtomWithIdx(ai).GetNeighbors():
            ni = nb.GetIdx()
            if ni not in origin:
                out.add(ni)
    return [i for i in out if i not in atoms]


def _s2c_expand_aa(m, aa_set):
    origin = [idx for aa in aa_set for idx in aa]
    expanded = []
    for aa in [set(a) for a in aa_set]:
        nbs = _s2c_side_chain_neighbors(m, aa, origin)
        aa.update(nbs)
        while nbs:
            nbs = _s2c_side_chain_neighbors(m, aa, origin)
            aa.update(nbs)
        expanded.append(list(aa))
    return expanded


def _s2c_split_residues(m):
    from rdkit import Chem as _C
    raw = m.GetSubstructMatches(_C.MolFromSmiles('NCC=O'))
    return _s2c_expand_aa(m, raw)


def _s2c_isolate_residues(m, aa_units):
    """Cut each residue from the peptide ring; cap C=O termini with -OH."""
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol
    aas, mappings = [], []
    for atom_idxs in aa_units:
        rw = RWMol()
        amap = {}
        for idx in atom_idxs:
            amap[idx] = rw.AddAtom(m.GetAtomWithIdx(idx))
        mappings.append(amap)
        # First pass: add all bonds (no stereo yet)
        for bond in m.GetBonds():
            bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in atom_idxs and ei in atom_idxs:
                rw.AddBond(amap[bi], amap[ei], bond.GetBondType())
            elif bi in atom_idxs or ei in atom_idxs:
                inside = bi if bi in atom_idxs else ei
                outside = ei if inside == bi else bi
                ia = m.GetAtomWithIdx(inside)
                oa = m.GetAtomWithIdx(outside)
                if ia.GetSymbol() == 'C' and oa.GetSymbol() == 'N':
                    has_dbl_o = any(
                        nb.GetSymbol() == 'O' and
                        m.GetBondBetweenAtoms(inside, nb.GetIdx()).GetBondType() == _C.BondType.DOUBLE
                        for nb in ia.GetNeighbors()
                    )
                    if has_dbl_o:
                        o_idx = rw.AddAtom(_C.Atom('O'))
                        rw.AddBond(amap[inside], o_idx, _C.BondType.SINGLE)
        # Second pass: restore E/Z stereo on double bonds (stereoAtoms must be bonded first)
        from rdkit.Chem.rdchem import BondStereo as _BS
        for bond in m.GetBonds():
            bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in atom_idxs and ei in atom_idxs and bond.GetBondTypeAsDouble() == 2.0:
                db_stereo = bond.GetStereo()
                if db_stereo not in (_BS.STEREONONE, _BS.STEREOANY):
                    sa = list(bond.GetStereoAtoms())
                    if len(sa) == 2 and sa[0] in atom_idxs and sa[1] in atom_idxs:
                        try:
                            nb = rw.GetBondBetweenAtoms(amap[bi], amap[ei])
                            nb.SetStereo(db_stereo)
                            nb.SetStereoAtoms(amap[sa[0]], amap[sa[1]])
                        except Exception:
                            pass
        mol_out = rw.GetMol()
        _C.SanitizeMol(mol_out)
        # Fix chirality: AddAtom copies ChiralTags, but they are relative to
        # the original molecule's neighbor ordering.  In the fragment the
        # neighbor order may differ, flipping the effective CW/CCW meaning.
        # Correct each chiral centre so its CIP code matches the original.
        from rdkit.Chem import AllChem as _AC2, rdchem as _RC2
        _AC2.AssignStereochemistry(mol_out, cleanIt=True, force=True)
        inv_map = {v: k for k, v in amap.items()}
        needs_fix = False
        for frag_idx, orig_idx in inv_map.items():
            orig_cip = m.GetAtomWithIdx(orig_idx).GetPropsAsDict().get('_CIPCode')
            if not orig_cip:
                continue
            frag_cip = mol_out.GetAtomWithIdx(frag_idx).GetPropsAsDict().get('_CIPCode')
            if frag_cip and frag_cip != orig_cip:
                needs_fix = True
                break
        if needs_fix:
            rw2 = _C.RWMol(mol_out)
            for frag_idx, orig_idx in inv_map.items():
                orig_cip = m.GetAtomWithIdx(orig_idx).GetPropsAsDict().get('_CIPCode')
                if not orig_cip:
                    continue
                frag_cip = mol_out.GetAtomWithIdx(frag_idx).GetPropsAsDict().get('_CIPCode')
                if frag_cip and frag_cip != orig_cip:
                    ct = rw2.GetAtomWithIdx(frag_idx).GetChiralTag()
                    if ct == _RC2.ChiralType.CHI_TETRAHEDRAL_CW:
                        rw2.GetAtomWithIdx(frag_idx).SetChiralTag(
                            _RC2.ChiralType.CHI_TETRAHEDRAL_CCW)
                    elif ct == _RC2.ChiralType.CHI_TETRAHEDRAL_CCW:
                        rw2.GetAtomWithIdx(frag_idx).SetChiralTag(
                            _RC2.ChiralType.CHI_TETRAHEDRAL_CW)
            mol_out = rw2.GetMol()
        aas.append(mol_out)
    return aas, mappings


def _s2c_connected_pairs(m, mappings):
    from rdkit import Chem as _C
    # Backbone alpha-N atoms: the N in every NCC=O backbone unit.
    # Isopeptide bonds (e.g. K epsilon-N → E_g gamma-CO) use a non-alpha N and
    # must not be counted as backbone peptide bonds — they are 'side chain'.
    backbone_n_set = {match[0] for match in m.GetSubstructMatches(_C.MolFromSmiles('NCC=O'))}
    aa_idxs = [set(mp.keys()) for mp in mappings]
    pairs = []
    for ai in range(len(aa_idxs) - 1):
        for aj in range(ai + 1, len(aa_idxs)):
            fi, fj = aa_idxs[ai], aa_idxs[aj]
            is_link = is_peptide = False
            link_count = 0
            for bond in m.GetBonds():
                bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if (bi in fi and ei in fj) or (bi in fj and ei in fi):
                    is_link = True; link_count += 1
                    for c_idx, n_idx in [(bi, ei), (ei, bi)]:
                        ca = m.GetAtomWithIdx(c_idx); na = m.GetAtomWithIdx(n_idx)
                        if ca.GetSymbol() == 'C' and na.GetSymbol() == 'N':
                            # Only a backbone peptide bond if the N is a backbone alpha-N
                            if n_idx in backbone_n_set:
                                if any(nb.GetSymbol() == 'O' and
                                       m.GetBondBetweenAtoms(c_idx, nb.GetIdx()).GetBondType() == _C.BondType.DOUBLE
                                       for nb in ca.GetNeighbors()):
                                    is_peptide = True
            if is_link:
                if is_peptide:
                    pairs.append((ai, aj, 'peptide bond'))
                if link_count > 1 or not is_peptide:
                    pairs.append((ai, aj, 'side chain'))
    return pairs


def _s2c_search_chain(qi, search, N, pairs):
    r, chain = 0, [qi]
    while r < N:
        for i, j, t in pairs:
            if t == 'side chain': continue
            if i == qi and j in search:
                qi = j; chain.append(qi); search = [k for k in search if k not in chain]; break
            elif j == qi and i in search:
                qi = i; chain.append(qi); search = [k for k in search if k not in chain]; break
        r += 1
    return chain


def _s2c_ordered_chain(pairs, N_aa):
    best = []
    for qi in range(N_aa):
        c = _s2c_search_chain(qi, [i for i in range(N_aa) if i != qi], N_aa, pairs)
        if len(c) > len(best): best = c
    chain = best
    if len(chain) == N_aa: return chain
    for _ in range(5):
        if len(chain) > N_aa - 2: break
        ext = []
        for qi in [k for k in range(N_aa) if k not in chain]:
            c = _s2c_search_chain(qi, [i for i in range(N_aa) if i != qi and i not in chain], N_aa, pairs)
            if len(c) > len(ext): ext = c
        chain += ext
    if len(chain) < N_aa: chain += [i for i in range(N_aa) if i not in chain]
    return chain


def _s2c_is_cyclic(pairs, n):
    from collections import defaultdict, deque
    pbs = [(i, j) for i, j, t in pairs if t == 'peptide bond']
    scs = {(min(i, j), max(i, j)) for i, j, t in pairs if t == 'side chain'}
    if n == 2:
        pb_set = {(min(i, j), max(i, j)) for i, j in pbs}
        return bool(pb_set & scs)
    # BFS cycle detection on peptide-bond graph; correctly handles sidechain
    # isopeptide branches (e.g. lipid-modified K→E_g) without false positives
    adj = defaultdict(set)
    for i, j in pbs:
        adj[i].add(j)
        adj[j].add(i)
    visited = set()
    parent = {}
    for start in list(adj):
        if start in visited:
            continue
        queue = deque([start])
        parent[start] = -1
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            for nb in adj[node]:
                if nb not in visited:
                    parent[nb] = node
                    queue.append(nb)
                elif nb != parent.get(node, -1):
                    return True
    return False


def _s2c_strip_n_cap(mol):
    """Remove N-terminal carbonyl-type cap (ac, fmoc, Boc) from backbone N.

    Only strips substituents whose first atom is a carbonyl C (has a =O).
    N-methyl groups (plain CH3) are preserved.
    Returns (stripped_mol, did_strip).
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol

    patt = _C.MolFromSmarts('[N]-[C]-[C](=O)')
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return mol, False

    n_idx, ca_idx = matches[0][0], matches[0][1]
    n_atom = mol.GetAtomWithIdx(n_idx)

    cap_starts = []
    for nb in n_atom.GetNeighbors():
        if nb.GetIdx() == ca_idx:
            continue
        if nb.GetSymbol() == 'C':
            has_carbonyl = any(
                x.GetSymbol() == 'O' and
                mol.GetBondBetweenAtoms(nb.GetIdx(), x.GetIdx()).GetBondTypeAsDouble() == 2.0
                for x in nb.GetNeighbors()
            )
            if has_carbonyl:
                cap_starts.append(nb.GetIdx())

    if not cap_starts:
        return mol, False

    cap_atoms = set()
    queue = list(cap_starts)
    while queue:
        ai = queue.pop()
        if ai in cap_atoms or ai == n_idx:
            continue
        cap_atoms.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            if nb.GetIdx() not in cap_atoms and nb.GetIdx() != n_idx:
                queue.append(nb.GetIdx())

    rw = RWMol(mol)
    for ai in sorted(cap_atoms, reverse=True):
        rw.RemoveAtom(ai)
    try:
        _C.SanitizeMol(rw)
        return rw.GetMol(), True
    except Exception:
        return mol, False


def _s2c_identify_n_cap(mol):
    """Return the abbreviation of the N-terminal cap attached to backbone N.

    Extracts the cap fragment (same atoms that _s2c_strip_n_cap removes), adds a
    dummy N at the attachment point, then matches that small fragment against
    library caps.  Matching against the fragment (not the whole residue) avoids
    false positives from backbone SMARTS patterns overlapping the cap query.
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol, Atom

    # Locate backbone N and the cap atoms (same logic as _s2c_strip_n_cap)
    patt = _C.MolFromSmarts('[N]-[C]-[C](=O)')
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return None
    n_idx, ca_idx = matches[0][0], matches[0][1]
    n_atom = mol.GetAtomWithIdx(n_idx)

    cap_starts = []
    for nb in n_atom.GetNeighbors():
        if nb.GetIdx() == ca_idx:
            continue
        if nb.GetSymbol() == 'C':
            has_carbonyl = any(
                x.GetSymbol() == 'O' and
                mol.GetBondBetweenAtoms(nb.GetIdx(), x.GetIdx()).GetBondTypeAsDouble() == 2.0
                for x in nb.GetNeighbors()
            )
            if has_carbonyl:
                cap_starts.append(nb.GetIdx())
    if not cap_starts:
        return None

    cap_atoms = set()
    queue = list(cap_starts)
    while queue:
        ai = queue.pop()
        if ai in cap_atoms or ai == n_idx:
            continue
        cap_atoms.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            if nb.GetIdx() not in cap_atoms and nb.GetIdx() != n_idx:
                queue.append(nb.GetIdx())

    # Build cap fragment with a placeholder N at the attachment point
    rw = RWMol()
    amap = {}
    for ai in cap_atoms:
        amap[ai] = rw.AddAtom(mol.GetAtomWithIdx(ai))
    n_dummy = rw.AddAtom(Atom('N'))
    for bond in mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in cap_atoms and ei in cap_atoms:
            rw.AddBond(amap[bi], amap[ei], bond.GetBondType())
        elif bi in cap_atoms and ei == n_idx:
            rw.AddBond(amap[bi], n_dummy, bond.GetBondType())
        elif ei in cap_atoms and bi == n_idx:
            rw.AddBond(amap[ei], n_dummy, bond.GetBondType())
    try:
        _C.SanitizeMol(rw)
    except Exception:
        return None
    cap_frag = rw.GetMol()

    # Build library queries: replace backbone_c dummy with N, then match fragment
    _, by_abbr = _load_sdf()
    candidates = []
    for abbr, cap_mol in by_abbr.items():
        p = cap_mol.GetPropsAsDict()
        if p.get('m_type', '') != 'cap':
            continue
        chem_types = p.get('m_chem_types', '')
        if 'backbone_c' not in chem_types or 'backbone_n' in chem_types:
            continue
        r_num = None
        for part in chem_types.split(','):
            part = part.strip()
            if ':' not in part:
                continue
            rn, rtype = part.split(':', 1)
            if rtype.strip() == 'backbone_c':
                try:
                    r_num = int(rn.strip()); break
                except ValueError:
                    pass
        if r_num is None:
            continue
        dummy_idx = next(
            (a.GetIdx() for a in cap_mol.GetAtoms()
             if a.GetAtomicNum() == 0 and a.GetIsotope() == r_num),
            None)
        if dummy_idx is None:
            continue
        rw2 = RWMol(cap_mol)
        rw2.ReplaceAtom(dummy_idx, Atom('N'))
        try:
            _C.SanitizeMol(rw2)
        except Exception:
            continue
        candidates.append((cap_mol.GetNumAtoms(), abbr, rw2.GetMol()))
    candidates.sort(key=lambda x: -x[0])
    for _, abbr, qmol in candidates:
        if cap_frag.HasSubstructMatch(qmol, useChirality=False):
            return abbr
    return None


def _s2c_strip_c_cap(mol):
    """Convert C-terminal cap to carboxylic acid (-COOH) for residue matching.

    Handles am (NH2), ester (OR), and other amide/ester C-caps by finding the
    terminal NCC=O, removing the cap substituent on the carbonyl C, and
    replacing it with OH.  Only touches the BACKBONE terminal C so
    Asn/Gln side-chain amides are safe.
    Returns (converted_mol, did_convert).
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol

    patt = _C.MolFromSmarts('[N]-[C]-[C](=O)')
    matches = mol.GetSubstructMatches(patt)
    for match in matches:
        n_idx, ca_idx, co_idx = match[0], match[1], match[2]
        co_atom = mol.GetAtomWithIdx(co_idx)
        already_acid = False
        cap_root = None
        for nb in co_atom.GetNeighbors():
            ni = nb.GetIdx()
            if ni == ca_idx:
                continue
            bd = mol.GetBondBetweenAtoms(co_idx, ni)
            if nb.GetSymbol() == 'O' and bd.GetBondTypeAsDouble() == 2.0:
                continue  # carbonyl =O, skip
            if nb.GetSymbol() == 'O' and bd.GetBondTypeAsDouble() == 1.0 and nb.GetTotalNumHs() > 0:
                already_acid = True
                break
            cap_root = ni
        if already_acid or cap_root is None:
            continue
        # BFS from cap_root to collect all cap atoms
        cap_atoms: set = set()
        queue = [cap_root]
        while queue:
            ai = queue.pop()
            if ai in cap_atoms or ai == co_idx:
                continue
            cap_atoms.add(ai)
            for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni2 = nb.GetIdx()
                if ni2 not in cap_atoms and ni2 != co_idx:
                    queue.append(ni2)
        rw = RWMol(mol)
        rw.ReplaceAtom(cap_root, _C.Atom('O'))
        for ai in sorted(cap_atoms - {cap_root}, reverse=True):
            rw.RemoveAtom(ai)
        try:
            _C.SanitizeMol(rw)
            return rw.GetMol(), True
        except Exception:
            continue

    return mol, False


# Capped library cache (invalidates when SDF changes)
_s2c_lib_cache: list = []
_s2c_lib_mtime: float = 0.0
_s2c_lib_lock = threading.Lock()

# Raw monomer library (R-atoms preserved, for R-group number detection)
_s2c_raw_lib_cache: dict = {}
_s2c_raw_lib_mtime: float = 0.0

# Full capped library (includes non-backbone monomers, for branch segment matching)
_s2c_full_lib_cache: list = []
_s2c_full_lib_mtime: float = 0.0

# Scaffold crosslink patterns (multi-arm non-amide linkers, e.g. TBMB)
_s2c_scaffold_cache: list | None = None
_s2c_scaffold_mtime: float = 0.0


def _s2c_get_lib():
    """Return (or rebuild) the capped monomer library."""
    from rdkit import Chem as _C
    global _s2c_lib_cache, _s2c_lib_mtime
    mtime = _SDF_PATH.stat().st_mtime
    with _s2c_lib_lock:
        if _s2c_lib_cache and _s2c_lib_mtime == mtime:
            return _s2c_lib_cache
        suppl = _C.SDMolSupplier(str(_SDF_PATH), removeHs=False)
        lib = []
        seen = set()
        for mol in suppl:
            if mol is None: continue
            props = mol.GetPropsAsDict()
            abbr = props.get('m_abbr', '').strip()
            if not abbr or abbr in seen: continue
            ct_str = props.get('m_chem_types', '')
            rg_str = props.get('m_Rgroups', '')
            cts = {}
            for item in ct_str.split(','):
                item = item.strip()
                if ':' in item:
                    k, v = item.split(':', 1)
                    try: cts[int(k)] = v.strip()
                    except ValueError: pass
            bn_r = next((r for r, t in cts.items() if t == 'backbone_n'), None)
            bc_r = next((r for r, t in cts.items() if t == 'backbone_c'), None)
            if bn_r is None or bc_r is None: continue
            bn_atom = next((a for a in mol.GetAtoms()
                            if a.GetAtomicNum() == 0 and a.GetIsotope() == bn_r), None)
            bc_atom = next((a for a in mol.GetAtoms()
                            if a.GetAtomicNum() == 0 and a.GetIsotope() == bc_r), None)
            if bn_atom is None or bc_atom is None: continue
            bn_nbrs = [nb.GetIdx() for nb in bn_atom.GetNeighbors()]
            bc_nbrs = [nb.GetIdx() for nb in bc_atom.GetNeighbors()]
            if not bn_nbrs or not bc_nbrs: continue
            _bb_path = _C.GetShortestPath(mol, bn_nbrs[0], bc_nbrs[0])
            bb_dist = len(_bb_path) - 1 if _bb_path else 0
            smi = _C.MolToSmiles(mol, allHsExplicit=False)
            capped_smi = _s2c_cap_smiles(smi, rg_str)
            try: cm = _C.MolFromSmiles(capped_smi)
            except Exception: cm = None
            if cm is None: continue
            cm = _C.RemoveHs(cm)
            orig_cm = cm                    # preserve stereochemistry
            norm_cm = _s2c_normalize(cm)    # resonance-normalised (strips stereo)
            if norm_cm is None: continue
            from rdkit.Chem import rdMolDescriptors as _RDM
            from rdkit.Chem import rdchem as _RC
            n_rings = _RDM.CalcNumRings(orig_cm)
            has_ez = any(
                b.GetStereo() not in (_RC.BondStereo.STEREONONE, _RC.BondStereo.STEREOANY)
                for b in orig_cm.GetBonds()
                if b.GetBondTypeAsDouble() == 2.0
            )
            lib.append((abbr, orig_cm, norm_cm, norm_cm.GetNumAtoms(), n_rings, has_ez, bb_dist))
            seen.add(abbr)
        lib.sort(key=lambda x: x[3], reverse=True)
        _s2c_lib_cache = lib
        _s2c_lib_mtime = mtime
        return lib


def _s2c_get_raw_lib():
    """Return {abbr: raw_mol} with R-atom dummy atoms (isotope labels) preserved."""
    from rdkit import Chem as _C
    global _s2c_raw_lib_cache, _s2c_raw_lib_mtime
    mtime = _SDF_PATH.stat().st_mtime
    with _s2c_lib_lock:
        if _s2c_raw_lib_cache and _s2c_raw_lib_mtime == mtime:
            return _s2c_raw_lib_cache
        suppl = _C.SDMolSupplier(str(_SDF_PATH), removeHs=False)
        raw = {}
        for mol in suppl:
            if mol is None:
                continue
            props = mol.GetPropsAsDict()
            abbr = props.get('m_abbr', '').strip()
            if abbr and abbr not in raw:
                raw[abbr] = mol
        _s2c_raw_lib_cache = raw
        _s2c_raw_lib_mtime = mtime
        return raw


def _s2c_crosslink_r(abbr, raw_lib):
    """Return the non-backbone R-group number for a crosslink monomer.

    Parses m_chem_types and returns the first R-group that is not backbone_n,
    backbone_c, or backbone_n_mod.  Falls back to 4 if nothing found.
    """
    raw_mol = raw_lib.get(abbr)
    if raw_mol is None:
        return 4
    chem_types = raw_mol.GetPropsAsDict().get('m_chem_types', '')
    _backbone = {'backbone_n', 'backbone_c', 'backbone_n_mod'}
    for part in chem_types.split(','):
        part = part.strip()
        if ':' not in part:
            continue
        r_num, r_type = part.split(':', 1)
        if r_type.strip() not in _backbone:
            try:
                return int(r_num.strip())
            except ValueError:
                pass
    return 4


def _s2c_get_full_lib():
    """Capped library including non-backbone monomers (caps, linkers, etc.)."""
    from rdkit import Chem as _C
    global _s2c_full_lib_cache, _s2c_full_lib_mtime
    mtime = _SDF_PATH.stat().st_mtime
    with _s2c_lib_lock:
        if _s2c_full_lib_cache and _s2c_full_lib_mtime == mtime:
            return _s2c_full_lib_cache
        suppl = _C.SDMolSupplier(str(_SDF_PATH), removeHs=False)
        lib = []
        seen = set()
        for mol in suppl:
            if mol is None:
                continue
            props = mol.GetPropsAsDict()
            abbr = props.get('m_abbr', '').strip()
            if not abbr or abbr in seen:
                continue
            rg_str = props.get('m_Rgroups', '')
            smi = _C.MolToSmiles(mol, allHsExplicit=False)
            capped_smi = _s2c_cap_smiles(smi, rg_str)
            try:
                cm = _C.MolFromSmiles(capped_smi)
            except Exception:
                cm = None
            if cm is None:
                continue
            cm = _C.RemoveHs(cm)
            norm_cm = _s2c_normalize(cm)
            if norm_cm is None:
                continue
            lib.append((abbr, cm, norm_cm, norm_cm.GetNumAtoms()))
            seen.add(abbr)
        lib.sort(key=lambda x: x[3], reverse=True)
        _s2c_full_lib_cache = lib
        _s2c_full_lib_mtime = mtime
        return lib


def _s2c_get_scaffold_patterns():
    """Return scaffold crosslink patterns for all multi-arm non-backbone linkers.

    A scaffold linker is any SDF monomer whose CHUCKLES has ≥ 2 dummy atoms with
    isotope ≥ 4 (non-backbone attachment points) and NO dummy atoms with isotopes
    1 or 2 (backbone N/C connections), making it a pure crosslink scaffold.

    Each entry: (abbr, smarts_mol, {smarts_atom_idx: r_group_num})
    where smarts_atom_idx is the index in the SMARTS mol for each attachment-point
    wildcard (the backbone atoms that bond into the scaffold in the assembled mol).

    Detection strategy: replace [n*] (n ≥ 4) in the CHUCKLES canonical SMILES with
    [*:n] (atom-mapped wildcard) so GetAtomMapNum() recovers the R-group number from
    the SMARTS match, then verify each matched atom sits inside a backbone node.
    """
    import re as _re2
    from rdkit import Chem as _C

    global _s2c_scaffold_cache, _s2c_scaffold_mtime
    mtime = _SDF_PATH.stat().st_mtime
    with _s2c_lib_lock:
        if _s2c_scaffold_cache is not None and _s2c_scaffold_mtime == mtime:
            return _s2c_scaffold_cache

    # Load raw lib OUTSIDE the lock — _s2c_get_raw_lib acquires the same
    # non-reentrant lock and would deadlock if called inside our with-block.
    raw_lib = _s2c_get_raw_lib()
    patterns = []

    for abbr, raw_mol in raw_lib.items():
        if raw_mol is None:
            continue
        atoms = list(raw_mol.GetAtoms())
        # Attachment points: dummy atoms (atomic num 0) with isotope ≥ 4
        attach = [(a.GetIdx(), a.GetIsotope()) for a in atoms
                  if a.GetAtomicNum() == 0 and a.GetIsotope() >= 4]
        if len(attach) < 2:
            continue
        # Exclude backbone linkers (they have R1 or R2 = backbone N/C)
        if any(a.GetAtomicNum() == 0 and a.GetIsotope() in (1, 2) for a in atoms):
            continue

        try:
            chuckles = _C.MolToSmiles(raw_mol, canonical=True)
        except Exception:
            continue

        # [n*] → [*:n] for n ≥ 4; leave backbone R-groups unchanged
        smarts_str = _re2.sub(
            r'\[(\d+)\*\]',
            lambda m: f'[*:{m.group(1)}]' if int(m.group(1)) >= 4 else m.group(0),
            chuckles,
        )
        smarts_mol = _C.MolFromSmarts(smarts_str)
        if smarts_mol is None:
            continue

        r_positions = {
            a.GetIdx(): a.GetAtomMapNum()
            for a in smarts_mol.GetAtoms()
            if a.GetAtomMapNum() >= 4
        }
        if len(r_positions) < 2:
            continue

        global_min_slot = min(iso for _, iso in attach)
        patterns.append((abbr, smarts_mol, r_positions, global_min_slot))

        # For thia_michael_c arms (no leaving group; unreacted form = vinyl CH2=CH-):
        # the full SMARTS CC[*:n] can't match the terminal vinyl because CH2= has no
        # non-H atom for [*:n] to bind.  Generate partial SMARTS for each such arm
        # by removing its [*:n] token, so partial products are detected correctly.
        mol_props = raw_mol.GetPropsAsDict()
        m_chem_str = mol_props.get('m_chem_types', '')
        thia_slots = set()
        for _part in m_chem_str.split(','):
            if ':' in _part:
                _slot_str, _ct = _part.strip().split(':', 1)
                try:
                    if _ct.strip() == 'thia_michael_c':
                        thia_slots.add(int(_slot_str))
                except ValueError:
                    pass

        for _unreacted_slot in thia_slots:
            _partial_str = _re2.sub(
                r'\[\*:' + str(_unreacted_slot) + r'\]', '', smarts_str
            )
            _partial_mol = _C.MolFromSmarts(_partial_str)
            if _partial_mol is None:
                continue
            _partial_r_pos = {
                a.GetIdx(): a.GetAtomMapNum()
                for a in _partial_mol.GetAtoms()
                if a.GetAtomMapNum() >= 4
            }
            if len(_partial_r_pos) < 2:
                continue
            patterns.append((abbr, _partial_mol, _partial_r_pos, global_min_slot))

    with _s2c_lib_lock:
        _s2c_scaffold_cache = patterns
        _s2c_scaffold_mtime = mtime
    return patterns


def _s2c_r_group_at_atom(orig_mol, orig_atom_idx, frag_atoms, abbr, raw_lib):
    """Find the R-group number (1-6) for orig_atom_idx (a junction atom) in monomer abbr.

    Matches the fragment's core atom graph against the raw SDF mol (with dummy R-atoms),
    then returns the isotope label of the dummy atom adjacent to the matched position.
    Returns None if no match found.
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol as _RW

    raw_mol = raw_lib.get(abbr)
    if raw_mol is None:
        return None

    # Collect: for each heavy atom in raw_mol, note adjacent dummy isotopes
    dummy_nbrs = {}
    for a in raw_mol.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetIsotope() > 0:
            for nb in a.GetNeighbors():
                dummy_nbrs.setdefault(nb.GetIdx(), []).append(a.GetIsotope())

    # Build core fragment mol from frag_atoms (bonds within fragment only)
    rw = _RW()
    old_to_new = {}
    for oi in sorted(frag_atoms):
        old_to_new[oi] = rw.AddAtom(_C.Atom(orig_mol.GetAtomWithIdx(oi).GetAtomicNum()))
    for bond in orig_mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in old_to_new and ei in old_to_new:
            rw.AddBond(old_to_new[bi], old_to_new[ei], bond.GetBondType())
    core_frag = rw.GetMol()
    new_junct_idx = old_to_new.get(orig_atom_idx)
    if new_junct_idx is None:
        return None

    # Build stripped raw mol (remove dummy atoms) and track index mapping
    raw_to_stripped = {}
    stripped_to_raw = {}
    si = 0
    for ri in range(raw_mol.GetNumAtoms()):
        if raw_mol.GetAtomWithIdx(ri).GetAtomicNum() != 0:
            raw_to_stripped[ri] = si
            stripped_to_raw[si] = ri
            si += 1
    rw2 = _RW(raw_mol)
    for di in sorted([a.GetIdx() for a in raw_mol.GetAtoms() if a.GetAtomicNum() == 0], reverse=True):
        rw2.RemoveAtom(di)
    stripped = rw2.GetMol()

    # Match core_frag within stripped (find core_frag as sub-graph of stripped)
    matches = stripped.GetSubstructMatches(core_frag, useChirality=False)
    for match in matches:
        if new_junct_idx < len(match):
            stripped_idx = match[new_junct_idx]
            raw_idx = stripped_to_raw.get(stripped_idx)
            if raw_idx is None:
                continue
            isos = dummy_nbrs.get(raw_idx, [])
            if isos:
                return min(isos)

    # Reverse: find stripped within core_frag.
    # Needed when the assembled segment has extra free-cap atoms (e.g. unused R2 α-OH stays
    # when R4 is the chain continuation) — core_frag > stripped so forward match fails.
    rev_matches = core_frag.GetSubstructMatches(stripped, useChirality=False)
    for rev_match in rev_matches:
        for stripped_i, core_i in enumerate(rev_match):
            if core_i == new_junct_idx:
                raw_idx = stripped_to_raw.get(stripped_i)
                if raw_idx is None:
                    continue
                isos = dummy_nbrs.get(raw_idx, [])
                if isos:
                    return min(isos)
    return None


_BRANCH_SEG_SKIP: frozenset = frozenset()  # no deprecated aliases remain


def _s2c_match_branch_segment(seg_frag, full_lib):
    """Match a branch segment mol against the full library; return (abbr, n_matched).

    Tries two directions:
    1. Forward: lib as subgraph of segment (segment >= lib in size).
    2. Reverse: segment as subgraph of lib (lib has one or two extra cap atoms like -OH
       that were consumed when forming the amide bond in the full peptide).

    R-group slot is determined by _s2c_r_group_at_atom from connectivity.
    """
    from rdkit import Chem as _C
    seg_norm = _s2c_normalize(seg_frag)
    if seg_norm is None:
        return None, 0
    best_abbr = None
    best_n = 0
    best_lib_n = float('inf')  # prefer smallest library entry for same n_m (tightest fit)
    seg_n = seg_norm.GetNumAtoms()
    for abbr, cm, norm_cm, n_lib in full_lib:
        if abbr in _BRANCH_SEG_SKIP:
            continue
        n_m = 0
        # Forward: find lib monomer as subgraph of our segment
        if n_lib <= seg_n + 2:
            matches = seg_norm.GetSubstructMatches(norm_cm, useChirality=False)
            if matches:
                n_m = len(matches[0])
        # Reverse: find our segment as subgraph of lib (lib has extra cap -OH/-H atoms).
        # Only try reverse if forward failed — protected variants (e.g. Glu_OAll) share
        # the Glu core but carry extra atoms; tightest-fit tie-break below prefers base E.
        if n_m == 0 and seg_n < n_lib <= seg_n + 4:
            matches = norm_cm.GetSubstructMatches(seg_norm, useChirality=False)
            if matches:
                n_m = len(matches[0])  # = seg_n atoms matched within lib
        # Prefer: (1) higher n_m; (2) smaller n_lib for same n_m (tightest match).
        # This ensures E (n_lib=10) beats Glu_OAll (n_lib=13) when both match 9 atoms.
        if n_m > best_n or (n_m == best_n and n_m > 0 and n_lib < best_lib_n):
            best_n = n_m
            best_abbr = abbr
            best_lib_n = n_lib
            if best_n == seg_n and best_lib_n == seg_n:
                break  # perfect exact-size match — stop searching
    return best_abbr, best_n


def _s2c_walk_branch(mol, glu_atoms, glu_outgoing_n, orphan_atoms, full_lib, raw_lib, anchor_abbr, anchor_n_atom, junction_abbr='E_g'):
    """Walk the branch chain starting from the junction monomer outward.

    Returns list of (abbr, prev_r, cur_r) for each piece in the branch
    (including the junction monomer).
    glu_outgoing_n: junction monomer's alpha-N atom index (outgoing)
    anchor_n_atom: K's epsilon-N atom index (main_at)
    junction_abbr: abbreviation of the junction monomer (detected, not hardcoded)
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol as _RW
    from collections import deque as _dq

    results = []

    # ── Step 1: R-group for K (anchor) → junction monomer connection ───────
    branch_at = None
    for nb in mol.GetAtomWithIdx(anchor_n_atom).GetNeighbors():
        if nb.GetIdx() in glu_atoms:
            branch_at = nb.GetIdx()
            break

    anchor_r = _s2c_r_group_at_atom(mol, anchor_n_atom,
                                     set(range(mol.GetNumAtoms())) - glu_atoms - orphan_atoms,
                                     anchor_abbr, raw_lib)
    glu_in_r = _s2c_r_group_at_atom(mol, branch_at, glu_atoms, junction_abbr, raw_lib) if branch_at is not None else None

    # ── Step 2: Identify junction monomer's outgoing R-group ────────────────
    glu_out_r = _s2c_r_group_at_atom(mol, glu_outgoing_n, glu_atoms, junction_abbr, raw_lib)

    results.append((junction_abbr, anchor_r or 4, glu_in_r or 4))

    if glu_outgoing_n is None or not orphan_atoms:
        return results

    # ── Step 3: Segment orphan atoms at amide N→C=O bonds ───────────────────
    # Walk from the first orphan atom (bonded to glu_outgoing_n) through orphan_atoms
    # Split at each internal amide bond (orphan-N bonded to orphan-C=O)

    # Find orphan entry atom (bonded to glu_outgoing_n)
    entry_orphan = None
    for nb in mol.GetAtomWithIdx(glu_outgoing_n).GetNeighbors():
        if nb.GetIdx() in orphan_atoms:
            entry_orphan = nb.GetIdx()
            break
    if entry_orphan is None:
        return results

    # BFS to find all connected orphan atoms in chain order
    # Build segments by splitting at each amide-N → C=O bond within orphans
    segments = []
    current_seg = set()
    visited = set()
    current_entry = entry_orphan
    prev_junction = glu_outgoing_n  # the non-orphan atom from which we entered

    def is_co(ai):
        atom = mol.GetAtomWithIdx(ai)
        if atom.GetSymbol() != 'C':
            return False
        return any(nb.GetSymbol() == 'O' and
                   mol.GetBondBetweenAtoms(ai, nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                   for nb in atom.GetNeighbors())

    # Walk the linear orphan chain
    queue = _dq([entry_orphan])
    prev_atom = prev_junction
    current_seg_entry = entry_orphan
    junctions = []  # list of (N_atom, CO_atom) internal amide bonds

    while queue:
        ai = queue.popleft()
        if ai in visited:
            continue
        visited.add(ai)
        current_seg.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nb.GetIdx()
            if ni in visited or ni not in orphan_atoms:
                continue
            # Check if this is an amide bond crossing (orphan-N to orphan-C=O or vice versa)
            atom_ai = mol.GetAtomWithIdx(ai)
            atom_ni = mol.GetAtomWithIdx(ni)
            if atom_ai.GetSymbol() == 'N' and is_co(ni):
                junctions.append((ai, ni))
            queue.append(ni)

    # Use junctions to split orphan_atoms into segments
    if not junctions:
        segments = [set(orphan_atoms)]
    else:
        # Rebuild segments by connectivity after removing junction bonds
        remaining = set(orphan_atoms)
        for (n_atom, co_atom) in junctions:
            # Split: everything reachable from entry_orphan without crossing n_atom→co_atom bond
            seg = set()
            q2 = _dq([entry_orphan])
            while q2:
                ai2 = q2.popleft()
                if ai2 in seg:
                    continue
                seg.add(ai2)
                for nb2 in mol.GetAtomWithIdx(ai2).GetNeighbors():
                    ni2 = nb2.GetIdx()
                    if ni2 in seg or ni2 not in remaining:
                        continue
                    if (ai2 == n_atom and ni2 == co_atom):
                        continue  # skip the junction bond
                    q2.append(ni2)
            segments.append(seg)
            entry_orphan = co_atom
            remaining -= seg

        if remaining:
            segments.append(remaining)

    # ── Step 4: Match each segment and find R-groups ─────────────────────────
    prev_out_r = glu_out_r or 1  # E_g's outgoing R (R1 for backbone_n/in_n)
    seg_entry_atoms = []  # entry atom for each segment (bonded to previous piece)

    # Reconstruct entry atom for each segment
    seg_entry = glu_outgoing_n
    for seg_idx, seg in enumerate(segments):
        # Find which atom in seg is bonded to the previous piece's outgoing atom
        seg_in_atom = None
        for ai in seg:
            for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                if nb.GetIdx() == seg_entry or (seg_idx > 0 and nb.GetIdx() in segments[seg_idx - 1]):
                    seg_in_atom = ai
                    break
            if seg_in_atom is not None:
                break
        # Hmm, above is fragile. Better: find atom in seg bonded to atom NOT in seg
        # that is either glu_outgoing_n or in previous segment
        if seg_in_atom is None:
            for ai in seg:
                for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                    ni2 = nb.GetIdx()
                    if ni2 not in seg and ni2 in (set() if seg_idx == 0 else segments[seg_idx - 1] | {glu_outgoing_n}):
                        seg_in_atom = ai
                        break
                if seg_in_atom is not None:
                    break
        seg_entry_atoms.append(seg_in_atom)

    # Find entry atom for first segment (bonded to glu_outgoing_n)
    seg0_in = None
    for ai in segments[0]:
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            if nb.GetIdx() == glu_outgoing_n:
                seg0_in = ai
                break
        if seg0_in is not None:
            break
    if seg0_in is not None:
        seg_entry_atoms[0] = seg0_in

    for seg_idx, seg in enumerate(segments):
        # Build segment mol
        rw3 = _RW()
        s2n = {}
        for ai in sorted(seg):
            s2n[ai] = rw3.AddAtom(_C.Atom(mol.GetAtomWithIdx(ai).GetAtomicNum()))
        for bond in mol.GetBonds():
            bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in s2n and ei in s2n:
                rw3.AddBond(s2n[bi], s2n[ei], bond.GetBondType())
        seg_mol = _C.RemoveHs(rw3.GetMol())

        abbr, n_matched = _s2c_match_branch_segment(seg_mol, full_lib)

        # 1-2 rule: prefer the alternative whose ENTRY slot uses the lowest R-group.
        # E_g (R2=γ-COOH entry) beats E (R4=γ-COOH entry) for γ-linked chains.
        if abbr:
            _seg_in = seg_entry_atoms[seg_idx] if seg_idx < len(seg_entry_atoms) else None
            if _seg_in is not None:
                _cur_in_r = _s2c_r_group_at_atom(mol, _seg_in, seg, abbr, raw_lib)
                if _cur_in_r is not None and _cur_in_r > 2:
                    _seg_norm = _s2c_normalize(seg_mol)
                    _seg_n = _seg_norm.GetNumAtoms() if _seg_norm else 0
                    for _alt_abbr, _alt_cm, _alt_norm, _alt_n in full_lib:
                        if _alt_abbr == abbr or _alt_abbr in _BRANCH_SEG_SKIP:
                            continue
                        _alt_nm = 0
                        if _seg_norm and _alt_n <= _seg_n + 2:
                            _m = _seg_norm.GetSubstructMatches(_alt_norm, useChirality=False)
                            if _m:
                                _alt_nm = len(_m[0])
                        if _alt_nm == 0 and _seg_n < _alt_n <= _seg_n + 4:
                            _m = _alt_norm.GetSubstructMatches(_seg_norm, useChirality=False)
                            if _m:
                                _alt_nm = _seg_n
                        if _alt_nm < n_matched:
                            continue
                        _alt_in_r = _s2c_r_group_at_atom(mol, _seg_in, seg, _alt_abbr, raw_lib)
                        if _alt_in_r is not None and _alt_in_r < _cur_in_r:
                            abbr = _alt_abbr
                            break

        seg_in_atom = seg_entry_atoms[seg_idx] if seg_idx < len(seg_entry_atoms) else None

        # cur_r: R-group on this segment at incoming junction
        cur_r = None
        if seg_in_atom is not None and abbr:
            cur_r = _s2c_r_group_at_atom(mol, seg_in_atom, seg, abbr, raw_lib)

        # Find outgoing junction atom of this segment (bonded to next segment or terminus)
        seg_out_atom = None
        if seg_idx < len(segments) - 1:
            next_seg = segments[seg_idx + 1]
            for ai in seg:
                for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                    if nb.GetIdx() in next_seg:
                        seg_out_atom = ai
                        break
                if seg_out_atom is not None:
                    break

        # Format: (abbr, prev_out_R, cur_R) — AEEA(1,2) means E_g's R1 connects
        # to AEEA's R2; prev_out_R=1 (E_g.in_n), cur_R=2 (AEEA.out_co).
        results.append((abbr or '?', prev_out_r, cur_r or 2))

        if seg_out_atom is not None and abbr:
            prev_out_r = _s2c_r_group_at_atom(mol, seg_out_atom, seg, abbr, raw_lib) or 1
        else:
            prev_out_r = 1

    return results


def _s2c_cip_score(aa_mol, ref_mol, match_tuple):
    """Score CIP stereo agreement for a match.

    Returns (n_agree, n_ref_stereo) where n_agree counts stereocenters with
    matching CIP codes and n_ref_stereo counts defined stereocenters in ref.
    Disagreements score 0 (not -1) because CIP codes are context-dependent:
    the same absolute configuration can produce different R/S labels when the
    molecular graph changes (e.g. isolated fragment vs full peptide chain).
    """
    from rdkit.Chem import AllChem as _AC
    _AC.AssignStereochemistry(aa_mol, cleanIt=True, force=True)
    _AC.AssignStereochemistry(ref_mol, cleanIt=True, force=True)
    n_agree = 0
    n_ref_stereo = 0
    for ref_idx, aa_idx in enumerate(match_tuple):
        ref_cip = ref_mol.GetAtomWithIdx(ref_idx).GetPropsAsDict().get('_CIPCode')
        aa_cip  = aa_mol.GetAtomWithIdx(aa_idx).GetPropsAsDict().get('_CIPCode')
        if ref_cip:
            n_ref_stereo += 1
            if aa_cip and ref_cip == aa_cip:
                n_agree += 1
    return n_agree, n_ref_stereo


def _s2c_match(aa_mol, lib):
    """Match an isolated residue mol against the library.

    Priority: (1) most atoms matched, (2) ring-count agreement (prevents linear
    patterns false-matching cyclic residues like Pip/Aze), (3) CIP stereo score,
    (4) E/Z double-bond stereo score.
    Early break only when all four are perfect (full CIP + E/Z agreement).
    """
    from rdkit.Chem import rdMolDescriptors as _RDM
    from rdkit.Chem import rdchem as _RC
    norm_mol = _s2c_normalize(aa_mol)
    n_q = aa_mol.GetNumAtoms()
    n_rings_aa = _RDM.CalcNumRings(aa_mol)

    # Save E/Z stereo BEFORE _s2c_cip_score (AssignStereochemistry cleanIt=True) corrupts it.
    # Keyed by frozenset({bi, ei}) → stereo enum; used for E/Z scoring in the loop.
    _BS = _RC.BondStereo
    _ez_saved = {}
    for b in aa_mol.GetBonds():
        if b.GetBondTypeAsDouble() == 2.0:
            s = b.GetStereo()
            if s not in (_BS.STEREONONE, _BS.STEREOANY):
                _ez_saved[frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))] = s
    query_has_ez = bool(_ez_saved)

    # D-form tiebreaker: when CIP scores are equal (e.g. flat/achiral input),
    # prefer L-form entries so T beats dT, A beats dA, etc.
    import re as _re2
    def _is_d_form(tok):
        return bool(_re2.match(r'^(D_|[Dd][A-Z])', tok))

    best_abbr = None
    best_atom_n = 0
    best_stereo = -999
    best_n_ref_stereo = 0
    best_ring_match = False
    best_ez = -2
    best_is_d = True   # so the first L-form candidate always beats the initial state

    for abbr, ref_orig, ref_norm, n_ref, n_rings_ref, ref_has_ez, *_ in lib:
        if n_ref > n_q: continue
        match = aa_mol.GetSubstructMatches(ref_orig, useChirality=False)
        if match:
            n_m = len(match[0])
            pairs = [_s2c_cip_score(aa_mol, ref_orig, m) for m in match]
            sc = max(a for a, _ in pairs)
            n_rs = max(r for _, r in pairs)
        else:
            match = norm_mol.GetSubstructMatches(ref_norm, useChirality=False)
            if not match: continue
            n_m = len(match[0])
            sc = 0
            n_rs = 0

        rings_match = (n_rings_ref == n_rings_aa)

        # E/Z scoring using saved pre-CIP state: +1 match, 0 neutral, -1 mismatch
        if query_has_ez and ref_has_ez:
            q_to_ref = {qa: ri for ri, qa in enumerate(match[0])}
            ez_sc = 0
            for b in ref_orig.GetBonds():
                if b.GetBondTypeAsDouble() != 2.0: continue
                ref_s = b.GetStereo()
                if ref_s in (_BS.STEREONONE, _BS.STEREOANY): continue
                ri, rj = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                qi = match[0][ri] if ri < len(match[0]) else None
                qj = match[0][rj] if rj < len(match[0]) else None
                if qi is None or qj is None: continue
                saved = _ez_saved.get(frozenset((qi, qj)))
                if saved is None:
                    ez_sc = 0; break
                ez_sc = 1 if saved == ref_s else -1
        else:
            ez_sc = 0

        is_d = _is_d_form(abbr)
        better = (
            n_m > best_atom_n or
            (n_m == best_atom_n and rings_match and not best_ring_match) or
            (n_m == best_atom_n and rings_match == best_ring_match and sc > best_stereo) or
            (n_m == best_atom_n and rings_match == best_ring_match and sc == best_stereo and n_rs > best_n_ref_stereo) or
            (n_m == best_atom_n and rings_match == best_ring_match and sc == best_stereo and n_rs == best_n_ref_stereo and ez_sc > best_ez) or
            (n_m == best_atom_n and rings_match == best_ring_match and sc == best_stereo and n_rs == best_n_ref_stereo and ez_sc == best_ez and not is_d and best_is_d)
        )
        if better:
            best_atom_n = n_m; best_abbr = abbr; best_stereo = sc
            best_n_ref_stereo = n_rs
            best_ring_match = rings_match; best_ez = ez_sc; best_is_d = is_d
            if n_m == n_q and rings_match:
                if sc == n_rs and (not query_has_ez or ez_sc > 0):
                    break
    return best_abbr, best_atom_n


# ── Coverage library: backbone monomers with tracked in_n / out_co ────────────
from collections import namedtuple as _nt2
_CovEntry = _nt2('_CovEntry', [
    'abbr', 'm_type', 'orig_mol', 'norm_mol', 'smarts_mol',
    'in_n_idx', 'out_co_idx',            # indices into orig_mol / norm_mol
    'smarts_in_n_idx', 'smarts_out_co_idx',  # indices into smarts_mol
    'n_atoms', 'n_rings', 'has_ez',
    'bb_dist',                           # shortest-path bond count N→backbone-C
    # iso SMARTS: carboxyl cap O removed so isopeptide-bonded monomers (E_g→K)
    # match without pulling the backbone N into the match set.  None if no
    # carboxyl R-groups are present.
    'smarts_mol_iso', 'smarts_in_n_idx_iso', 'smarts_out_co_idx_iso',
])
_PlacedNode = _nt2('_PlacedNode', [
    'abbr', 'm_type', 'mol_atoms', 'in_n', 'out_co', 'entry', 'match_tuple',
])

_s2c_cov_lib_cache: list = []
_s2c_cov_lib_mtime: float = 0.0


def _s2c_get_coverage_lib():
    """Backbone monomer library; each entry tracks in_n / out_co atom indices."""
    global _s2c_cov_lib_cache, _s2c_cov_lib_mtime
    mtime = _SDF_PATH.stat().st_mtime
    with _s2c_lib_lock:
        if _s2c_cov_lib_cache and _s2c_cov_lib_mtime == mtime:
            return _s2c_cov_lib_cache
        lib = _s2c_build_cov_lib()
        _s2c_cov_lib_cache = lib
        _s2c_cov_lib_mtime = mtime
        return lib


def _s2c_build_cov_lib():
    """Build coverage library: backbone monomers with atom-map-tracked terminals."""
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol as _RW, rdMolDescriptors as _RDM, rdchem as _RC
    _MAP_IN, _MAP_OUT = 9901, 9902
    suppl = _C.SDMolSupplier(str(_SDF_PATH), removeHs=False)
    lib = []
    seen: set = set()
    for raw_mol in suppl:
        if raw_mol is None:
            continue
        props = raw_mol.GetPropsAsDict()
        abbr = props.get('m_abbr', '').strip()
        if not abbr or abbr in seen:
            continue
        m_type = props.get('m_type', 'aa')
        ct_str = props.get('m_chem_types', '')
        rg_str = props.get('m_Rgroups', '')
        cts: dict = {}
        for item in ct_str.split(','):
            item = item.strip()
            if ':' in item:
                k, v = item.split(':', 1)
                try:
                    cts[int(k)] = v.strip()
                except ValueError:
                    pass
        bn_r = next((r for r, t in cts.items() if t == 'backbone_n'), None)
        bc_r = next((r for r, t in cts.items() if t == 'backbone_c'), None)
        if bn_r is None or bc_r is None:
            continue
        bn_dummy = next((a.GetIdx() for a in raw_mol.GetAtoms()
                         if a.GetAtomicNum() == 0 and a.GetIsotope() == bn_r), None)
        bc_dummy = next((a.GetIdx() for a in raw_mol.GetAtoms()
                         if a.GetAtomicNum() == 0 and a.GetIsotope() == bc_r), None)
        if bn_dummy is None or bc_dummy is None:
            continue
        bn_nbrs = [nb.GetIdx() for nb in raw_mol.GetAtomWithIdx(bn_dummy).GetNeighbors()]
        bc_nbrs = [nb.GetIdx() for nb in raw_mol.GetAtomWithIdx(bc_dummy).GetNeighbors()]
        if not bn_nbrs or not bc_nbrs:
            continue
        in_n_raw = bn_nbrs[0]
        out_co_raw = bc_nbrs[0]
        # Relax H-count SMARTS for atoms adjacent to R-groups that undergo
        # connectivity changes during crosslinking (terminal_alkene → =CH2→=CH-).
        # Restrict to alkene-type R-groups ONLY — other non-backbone R-groups
        # like formamide_c/amide_nh must keep their exact H count or they absorb
        # neighbouring residue fragments (e.g. Lys_For absorbing E_g gamma-CO).
        _MAP_XL = 9904
        _XL_ACTIVE_TYPES = frozenset({'terminal_alkene', 'crosslink_alkene'})
        xl_raw_attach: set = set()
        for _rl_r, _rl_t in cts.items():
            if _rl_t not in _XL_ACTIVE_TYPES:
                continue
            _xl_dum = next((a.GetIdx() for a in raw_mol.GetAtoms()
                           if a.GetAtomicNum() == 0 and a.GetIsotope() == _rl_r), None)
            if _xl_dum is not None:
                for _xl_nb in raw_mol.GetAtomWithIdx(_xl_dum).GetNeighbors():
                    _xl_idx = _xl_nb.GetIdx()
                    if _xl_idx in (in_n_raw, out_co_raw):
                        continue
                    xl_raw_attach.add(_xl_idx)
                    for _xl_bd in raw_mol.GetAtomWithIdx(_xl_idx).GetBonds():
                        if _xl_bd.GetBondTypeAsDouble() == 2.0:
                            _pi = (_xl_bd.GetEndAtomIdx()
                                   if _xl_bd.GetBeginAtomIdx() == _xl_idx
                                   else _xl_bd.GetBeginAtomIdx())
                            if _pi != _xl_dum and _pi not in (in_n_raw, out_co_raw):
                                xl_raw_attach.add(_pi)
        # For carboxyl-type R-groups (isopeptide bond forming), the R-group cap
        # atom (O from [OH] capping) ends up in the SMARTS as [#8] but in the
        # assembled molecule that position is the amide-N of the isopeptide bond.
        # Track attachment atoms with _MAP_CR so we can find and remove the
        # capping O from the SMARTS (same treatment as bc_r_cap_idx for out_co).
        _MAP_CR = 9906
        _CR_ACTIVE_TYPES = frozenset({'carboxyl'})
        cap_remove_raw_attach: set = set()
        for _rl_r, _rl_t in cts.items():
            if _rl_t not in _CR_ACTIVE_TYPES:
                continue
            _cr_dum = next((a.GetIdx() for a in raw_mol.GetAtoms()
                           if a.GetAtomicNum() == 0 and a.GetIsotope() == _rl_r), None)
            if _cr_dum is not None:
                for _cr_nb in raw_mol.GetAtomWithIdx(_cr_dum).GetNeighbors():
                    _cr_idx = _cr_nb.GetIdx()
                    if _cr_idx not in (in_n_raw, out_co_raw):
                        cap_remove_raw_attach.add(_cr_idx)
        rw = _RW(raw_mol)
        rw.GetAtomWithIdx(in_n_raw).SetAtomMapNum(_MAP_IN)
        rw.GetAtomWithIdx(out_co_raw).SetAtomMapNum(_MAP_OUT)
        for _xl_raw in xl_raw_attach:
            rw.GetAtomWithIdx(_xl_raw).SetAtomMapNum(_MAP_XL)
        for _cr_raw in cap_remove_raw_attach:
            rw.GetAtomWithIdx(_cr_raw).SetAtomMapNum(_MAP_CR)
        smi = _C.MolToSmiles(rw.GetMol(), allHsExplicit=False)
        # Cap bn_r with [H] so the query matches in any peptide context.
        # Keep bc_r's SDF default (e.g. [OH] for standard AAs, [H] for _al forms)
        # so that Gly vs Gly_al can be distinguished by bc_r cap presence.
        rg_parts = [rg.strip() for rg in rg_str.split(',')]
        for _ri in range(len(rg_parts)):
            if (_ri + 1) == bn_r:
                rg_parts[_ri] = '[H]'
        cov_rg_str = ', '.join(rg_parts)
        capped_smi = _s2c_cap_smiles(smi, cov_rg_str)
        try:
            cm = _C.MolFromSmiles(capped_smi)
        except Exception:
            continue
        if cm is None:
            continue
        cm = _C.RemoveHs(cm)
        in_n_idx = next((a.GetIdx() for a in cm.GetAtoms() if a.GetAtomMapNum() == _MAP_IN), None)
        out_co_idx = next((a.GetIdx() for a in cm.GetAtoms() if a.GetAtomMapNum() == _MAP_OUT), None)
        if in_n_idx is None or out_co_idx is None:
            continue
        xl_cm_idxs = frozenset(a.GetIdx() for a in cm.GetAtoms()
                               if a.GetAtomMapNum() == _MAP_XL)
        cap_remove_cm_idxs = frozenset(a.GetIdx() for a in cm.GetAtoms()
                                       if a.GetAtomMapNum() == _MAP_CR)
        # Degree-1 single-bond neighbors of carboxyl-attachment atoms are the
        # O atoms added when the carboxyl R-group dummy was replaced by [OH].
        # In the assembled molecule that position is the isopeptide amide-N, so
        # remove these capping O atoms from the SMARTS.
        cap_remove_cap_set: set = set()
        for _cr_cm in cap_remove_cm_idxs:
            for _cr_nb in cm.GetAtomWithIdx(_cr_cm).GetNeighbors():
                _bd = cm.GetBondBetweenAtoms(_cr_cm, _cr_nb.GetIdx())
                if _cr_nb.GetDegree() == 1 and _bd.GetBondTypeAsDouble() != 2.0:
                    cap_remove_cap_set.add(_cr_nb.GetIdx())
        # Build H-count-specific SMARTS while cm still carries _MAP_IN/_MAP_OUT.
        # - in_n → [#7:MAP] (any N; matches free, protonated, or amide N-terminus)
        # - out_co → [#{anum}H{nh}:MAP] using ACTUAL H count from capped mol.
        #   Standard AAs (bc_r=[OH] default) → 0H = distinguishes them from _al forms.
        #   Aldehyde forms (bc_r=[H] default → removed by RemoveHs) → 1H = CHO.
        # - bc_r cap atom (degree-1 neighbor of out_co, not =O) → marked 9903,
        #   then removed from smarts_mol so it doesn't constrain C-terminal context.
        # - All other side-chain atoms → [#{anum}H{nh}] (exact H count) to prevent
        #   neighbouring cap atoms from being absorbed into the match.
        import re as _re2
        _MAP_CAP = 9903
        bc_r_cap_idx = None
        for _nb in cm.GetAtomWithIdx(out_co_idx).GetNeighbors():
            _bd = cm.GetBondBetweenAtoms(out_co_idx, _nb.GetIdx())
            if _bd.GetBondTypeAsDouble() != 2.0 and _nb.GetDegree() == 1:
                bc_r_cap_idx = _nb.GetIdx()
                break
        # Use a per-atom map number to track carboxyl-cap O atoms so we can
        # replace their SMARTS with [#8,#7] (O or N) rather than [#8] alone.
        # This lets E_g match in both free-carboxyl (-OH) and isopeptide (-N)
        # forms while still requiring D/E's carboxyl in free-acid contexts.
        _MAP_CR_CAP = 9907
        rw_s = _RW(cm)
        idx_to_mn: dict = {}
        for _sa in rw_s.GetAtoms():
            _si = _sa.GetIdx()
            if _sa.GetAtomMapNum() in (_MAP_IN, _MAP_OUT):
                idx_to_mn[_si] = _sa.GetAtomMapNum()  # keep 9901/9902 for tracking
            elif _si == bc_r_cap_idx:
                _sa.SetAtomMapNum(_MAP_CAP)
                idx_to_mn[_si] = _MAP_CAP
            elif _si in cap_remove_cap_set:
                _sa.SetAtomMapNum(_MAP_CR_CAP)
                idx_to_mn[_si] = _MAP_CR_CAP
            else:
                _mn = _si + 1  # 1-indexed; won't clash with reserved map nums
                _sa.SetAtomMapNum(_mn)
                idx_to_mn[_si] = _mn
        _smarts_str = _C.MolToSmarts(rw_s.GetMol())
        for _sa in cm.GetAtoms():
            _si = _sa.GetIdx()
            _mn = idx_to_mn[_si]
            _nh = _sa.GetTotalNumHs()
            _ch = _sa.GetFormalCharge()
            _cs = (f'+{_ch}' if _ch > 0 else (str(_ch) if _ch < 0 else ''))
            if _si == in_n_idx:
                _repl = f'[#7:{_MAP_IN}]'   # permissive N
            elif _si == out_co_idx:
                _repl = f'[#{_sa.GetAtomicNum()}H{_nh}{_cs}:{_MAP_OUT}]'  # exact H count
            elif _si == bc_r_cap_idx:
                _repl = f'[*:{_MAP_CAP}]'  # temporary; will be removed below
            elif _si in cap_remove_cap_set:
                # Carboxyl R-group cap atom: keep as [#8] in standard SMARTS
                # (correct for D/E backbone detection).  The iso SMARTS removes
                # this atom entirely so E_g matches in isopeptide form.
                _repl = f'[#{_sa.GetAtomicNum()}{_cs}:{_MAP_CR_CAP}]'
            elif _sa.GetAtomicNum() in (7, 8, 16):
                # Side-chain heteroatoms (N/O/S) may form crosslinks/isopeptide
                # bonds that change their H count; match atom type only.
                _repl = f'[#{_sa.GetAtomicNum()}{_cs}]'
            elif _si in xl_cm_idxs:
                # Atom adjacent to a non-backbone R-group (e.g. terminal alkene
                # for olefin-staple crosslinks).  Crosslinking changes H count
                # (=CH2 → =CH-), so don't pin it.
                _repl = f'[#{_sa.GetAtomicNum()}{_cs}]'
            else:
                _repl = f'[#{_sa.GetAtomicNum()}H{_nh}{_cs}]'
            _smarts_str = _re2.sub(
                r'\[[^\]]*:' + str(_mn) + r'\]', _repl, _smarts_str, count=1)
        smarts_mol = _C.MolFromSmarts(_smarts_str)
        if smarts_mol is None:
            continue
        # Remove bc_r_cap atom ([*:9903]) so the query imposes no constraint on
        # the C-terminal context (amide bond, am cap, free acid all match).
        _cap_idxs_in_sm = sorted(
            [a.GetIdx() for a in smarts_mol.GetAtoms() if a.GetAtomMapNum() == _MAP_CAP],
            reverse=True)
        if _cap_idxs_in_sm:
            rw_sm0 = _RW(smarts_mol)
            for _cap_idx in _cap_idxs_in_sm:
                rw_sm0.RemoveAtom(_cap_idx)
            smarts_mol = rw_sm0.GetMol()
        # Build iso SMARTS by further removing the carboxyl cap atoms (tagged
        # _MAP_CR_CAP).  This lets E_g match in isopeptide form (no -OH at
        # gamma-carboxyl) without pulling K's epsilon-N into the match set.
        _cr_cap_idxs = sorted(
            [a.GetIdx() for a in smarts_mol.GetAtoms()
             if a.GetAtomMapNum() == _MAP_CR_CAP],
            reverse=True)
        if _cr_cap_idxs:
            rw_iso = _RW(smarts_mol)
            for _ci in _cr_cap_idxs:
                rw_iso.RemoveAtom(_ci)
            smarts_mol_iso = rw_iso.GetMol()
            smarts_in_n_idx_iso = next(
                (a.GetIdx() for a in smarts_mol_iso.GetAtoms()
                 if a.GetAtomMapNum() == _MAP_IN), None)
            smarts_out_co_idx_iso = next(
                (a.GetIdx() for a in smarts_mol_iso.GetAtoms()
                 if a.GetAtomMapNum() == _MAP_OUT), None)
            if smarts_in_n_idx_iso is None or smarts_out_co_idx_iso is None:
                smarts_mol_iso = None
                smarts_in_n_idx_iso = None
                smarts_out_co_idx_iso = None
            else:
                rw_iso2 = _RW(smarts_mol_iso)
                for _sa in rw_iso2.GetAtoms():
                    _sa.SetAtomMapNum(0)
                smarts_mol_iso = rw_iso2.GetMol()
        else:
            smarts_mol_iso = None
            smarts_in_n_idx_iso = None
            smarts_out_co_idx_iso = None
        # Strip _MAP_CR_CAP from standard smarts_mol (atoms stay, map num goes).
        for _sa in smarts_mol.GetAtoms():
            if _sa.GetAtomMapNum() == _MAP_CR_CAP:
                _sa.SetAtomMapNum(0)
        # Locate in_n / out_co in smarts_mol by map num, then strip map nums.
        smarts_in_n_idx = next(
            (a.GetIdx() for a in smarts_mol.GetAtoms() if a.GetAtomMapNum() == _MAP_IN), None)
        smarts_out_co_idx = next(
            (a.GetIdx() for a in smarts_mol.GetAtoms() if a.GetAtomMapNum() == _MAP_OUT), None)
        if smarts_in_n_idx is None or smarts_out_co_idx is None:
            continue
        rw_sm = _RW(smarts_mol)
        for _sa in rw_sm.GetAtoms():
            _sa.SetAtomMapNum(0)
        smarts_mol = rw_sm.GetMol()
        rw2 = _RW(cm)
        for a in rw2.GetAtoms():
            if a.GetAtomMapNum() in (_MAP_IN, _MAP_OUT):
                a.SetAtomMapNum(0)
        orig_cm = rw2.GetMol()
        norm_cm = _s2c_normalize(orig_cm)
        if norm_cm is None:
            continue
        n_rings = _RDM.CalcNumRings(orig_cm)
        _BS = _RC.BondStereo
        has_ez = any(
            b.GetStereo() not in (_BS.STEREONONE, _BS.STEREOANY)
            for b in orig_cm.GetBonds() if b.GetBondTypeAsDouble() == 2.0
        )
        _bb_path = _C.GetShortestPath(raw_mol, in_n_raw, out_co_raw)
        _bb_dist = len(_bb_path) - 1 if _bb_path else 0
        lib.append(_CovEntry(
            abbr=abbr, m_type=m_type,
            orig_mol=orig_cm, norm_mol=norm_cm, smarts_mol=smarts_mol,
            in_n_idx=in_n_idx, out_co_idx=out_co_idx,
            smarts_in_n_idx=smarts_in_n_idx, smarts_out_co_idx=smarts_out_co_idx,
            n_atoms=norm_cm.GetNumAtoms(), n_rings=n_rings, has_ez=has_ez,
            bb_dist=_bb_dist,
            smarts_mol_iso=smarts_mol_iso,
            smarts_in_n_idx_iso=smarts_in_n_idx_iso,
            smarts_out_co_idx_iso=smarts_out_co_idx_iso,
        ))
        seen.add(abbr)
    lib.sort(key=lambda e: e.n_atoms, reverse=True)
    return lib


def _s2c_place_monomers(mol, cov_lib):
    """Return all placements of backbone monomers in mol as _PlacedNode objects."""
    placements: list = []
    for entry in cov_lib:
        matches = mol.GetSubstructMatches(entry.smarts_mol, useChirality=False)
        for match in matches:
            placements.append(_PlacedNode(
                abbr=entry.abbr, m_type=entry.m_type,
                mol_atoms=frozenset(match),
                in_n=match[entry.smarts_in_n_idx],
                out_co=match[entry.smarts_out_co_idx],
                entry=entry, match_tuple=match,
            ))
    return placements


def _s2c_build_backbone(mol, placements):
    """Find the longest non-overlapping backbone chain via amide-bond DAG.

    A→B when A.out_co bonds to B.in_n.  Among equal-length chains prefers
    those with the most m_type='aa' members (avoids picking lipid branch chains
    of the same length as the backbone, e.g. C20FA-AEEA-E_g vs K-G-K).
    Returns list of _PlacedNode in N→C order.
    """
    from collections import defaultdict as _dd3
    if not placements:
        return []
    # Deduplicate by atom coverage: L and D amino acids share identical heavy-atom
    # L and D amino acids share identical heavy-atom SMARTS (useChirality=False)
    # so they produce the same placement — _s2c_match (step 4) resolves stereo.
    # HOWEVER beta/gamma backbone variants (E vs E_g, D vs D_b) cover identical
    # atoms but use a DIFFERENT out_co (different carbonyl as R2).  Both must
    # survive deduplication so the chain DAG can select whichever one's R2 is the
    # actual amide-bond exit — that is the (1,2) backbone connection.
    # Dedup key: (mol_atoms, out_co) keeps one entry per distinct backbone exit.
    placements.sort(key=lambda p: (p.m_type != 'aa', -p.entry.n_atoms, p.abbr))
    seen_mol_out: set = set()
    deduped: list = []
    for p in placements:
        key = (p.mol_atoms, p.out_co)
        if key not in seen_mol_out:
            seen_mol_out.add(key)
            deduped.append(p)

    in_n_map = _dd3(list)
    for p in deduped:
        in_n_map[p.in_n].append(p)

    # out_cos from placements that don't share atoms with a given node — these
    # are genuine predecessor exits from OTHER residues.  Excludes the node's own
    # out_co and any alternate-backbone out_co for the same atom set (e.g. E and
    # E_g share mol_atoms; neither should count as the other's predecessor).
    def _external_out_cos(node):
        return {p.out_co for p in deduped if not (p.mol_atoms & node.mol_atoms)}

    def has_predecessor(node):
        ext = _external_out_cos(node)
        for nb in mol.GetAtomWithIdx(node.in_n).GetNeighbors():
            if nb.GetIdx() in ext:
                return True
        return False

    def successors(node, used_atoms):
        result = []
        for nb in mol.GetAtomWithIdx(node.out_co).GetNeighbors():
            for q in in_n_map[nb.GetIdx()]:
                # Check only backbone anchors (in_n, out_co) for overlap — not
                # all mol_atoms — so that RCM-stapled residues sharing the
                # crosslink C=C atoms are not incorrectly excluded.
                if q.in_n not in used_atoms and q.out_co not in used_atoms:
                    result.append(q)
        result.sort(key=lambda q: (q.m_type != 'aa', -q.entry.n_atoms))
        return result

    starts = [p for p in deduped if not has_predecessor(p)]
    is_cyclic_topology = not starts
    if is_cyclic_topology:
        # Pure cyclic peptide: any rotation is valid. Only try one start.
        starts = deduped
        starts.sort(key=lambda p: (p.m_type != 'aa', -p.entry.n_atoms))
        starts = starts[:1]
    else:
        # Mixed topology: linear lipid branches co-exist with a cyclic backbone.
        # Residues on the cyclic backbone all have predecessors → not in starts.
        # Explicitly add them so the DFS can find the (longer) cyclic chain.
        cyclic_cands = [p for p in deduped if has_predecessor(p)]
        starts = starts + cyclic_cands
        starts.sort(key=lambda p: (p.m_type != 'aa', -p.entry.n_atoms))

    best: list = []
    best_score = (-1, -1, -1)

    def chain_score(c):
        # Prefer longer chains, then more AAs, then more total atoms covered
        # (larger atom sets = greedier match, e.g. G over Gly_al).
        return (len(c), sum(1 for nd in c if nd.m_type == 'aa'),
                sum(len(nd.mol_atoms) for nd in c))

    def dfs(node, chain, used):
        nonlocal best, best_score
        sc = chain_score(chain)
        if sc > best_score:
            best[:] = chain[:]
            best_score = sc
        for nxt in successors(node, used)[:3]:  # limit branching factor
            chain.append(nxt)
            dfs(nxt, chain, used | nxt.mol_atoms)
            chain.pop()

    for start in starts:
        dfs(start, [start], start.mol_atoms)

    return best


def _s2c_identify_n_cap_from_atom(mol, n_idx, cap_start_atom, excluded_atoms):
    """Identify N-cap given backbone N atom and first cap atom.

    Extracts the cap fragment (BFS from cap_start_atom, excluding backbone atoms),
    attaches a placeholder N, then matches against library caps.
    Returns cap abbr or None.
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol, Atom
    cap_atoms: set = set()
    queue = [cap_start_atom]
    while queue:
        ai = queue.pop()
        if ai in cap_atoms or ai == n_idx or ai in excluded_atoms:
            continue
        cap_atoms.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nb.GetIdx()
            if ni not in cap_atoms and ni != n_idx and ni not in excluded_atoms:
                queue.append(ni)
    if not cap_atoms:
        return None
    rw = RWMol()
    amap: dict = {}
    for ai in cap_atoms:
        amap[ai] = rw.AddAtom(mol.GetAtomWithIdx(ai))
    n_dummy = rw.AddAtom(Atom('N'))
    for bond in mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in cap_atoms and ei in cap_atoms:
            rw.AddBond(amap[bi], amap[ei], bond.GetBondType())
        elif bi in cap_atoms and ei == n_idx:
            rw.AddBond(amap[bi], n_dummy, bond.GetBondType())
        elif ei in cap_atoms and bi == n_idx:
            rw.AddBond(amap[ei], n_dummy, bond.GetBondType())
    try:
        _C.SanitizeMol(rw)
    except Exception:
        return None
    cap_frag = rw.GetMol()
    _, by_abbr = _load_sdf()
    candidates: list = []
    for abbr, cap_mol in by_abbr.items():
        p = cap_mol.GetPropsAsDict()
        if p.get('m_type', '') != 'cap':
            continue
        chem_types = p.get('m_chem_types', '')
        if 'backbone_c' not in chem_types or 'backbone_n' in chem_types:
            continue
        r_num = None
        for part in chem_types.split(','):
            part = part.strip()
            if ':' not in part:
                continue
            rn, rtype = part.split(':', 1)
            if rtype.strip() == 'backbone_c':
                try:
                    r_num = int(rn.strip()); break
                except ValueError:
                    pass
        if r_num is None:
            continue
        dummy_idx = next((a.GetIdx() for a in cap_mol.GetAtoms()
                          if a.GetAtomicNum() == 0 and a.GetIsotope() == r_num), None)
        if dummy_idx is None:
            continue
        rw2 = RWMol(cap_mol)
        rw2.ReplaceAtom(dummy_idx, Atom('N'))
        try:
            _C.SanitizeMol(rw2)
        except Exception:
            continue
        candidates.append((cap_mol.GetNumAtoms(), abbr, rw2.GetMol()))
    candidates.sort(key=lambda x: -x[0])
    for _, abbr, qmol in candidates:
        if cap_frag.HasSubstructMatch(qmol, useChirality=False):
            return abbr
    return None


def _s2c_identify_c_cap_from_atom(mol, co_idx, cap_start_atom, excluded_atoms):
    """Identify C-cap given backbone carbonyl C atom and first cap atom.

    Extracts the cap fragment (BFS from cap_start_atom, excluding backbone atoms),
    attaches a placeholder C at the attachment point, then matches against library
    C-caps (backbone_n connection type).
    Returns cap abbr or None.
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol, Atom
    cap_atoms: set = set()
    queue = [cap_start_atom]
    while queue:
        ai = queue.pop()
        if ai in cap_atoms or ai == co_idx or ai in excluded_atoms:
            continue
        cap_atoms.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nb.GetIdx()
            if ni not in cap_atoms and ni != co_idx and ni not in excluded_atoms:
                queue.append(ni)
    if not cap_atoms:
        return None
    rw = RWMol()
    amap: dict = {}
    for ai in cap_atoms:
        amap[ai] = rw.AddAtom(mol.GetAtomWithIdx(ai))
    c_dummy = rw.AddAtom(Atom('C'))
    for bond in mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in cap_atoms and ei in cap_atoms:
            rw.AddBond(amap[bi], amap[ei], bond.GetBondType())
        elif bi in cap_atoms and ei == co_idx:
            rw.AddBond(amap[bi], c_dummy, bond.GetBondType())
        elif ei in cap_atoms and bi == co_idx:
            rw.AddBond(amap[ei], c_dummy, bond.GetBondType())
    try:
        _C.SanitizeMol(rw)
    except Exception:
        return None
    cap_frag = rw.GetMol()
    _, by_abbr = _load_sdf()
    candidates: list = []
    for abbr, cap_mol in by_abbr.items():
        p = cap_mol.GetPropsAsDict()
        if p.get('m_type', '') != 'cap':
            continue
        chem_types = p.get('m_chem_types', '')
        if 'backbone_n' not in chem_types or 'backbone_c' in chem_types:
            continue
        r_num = None
        for part in chem_types.split(','):
            part = part.strip()
            if ':' not in part:
                continue
            rn, rtype = part.split(':', 1)
            if rtype.strip() == 'backbone_n':
                try:
                    r_num = int(rn.strip()); break
                except ValueError:
                    pass
        if r_num is None:
            continue
        dummy_idx = next((a.GetIdx() for a in cap_mol.GetAtoms()
                          if a.GetAtomicNum() == 0 and a.GetIsotope() == r_num), None)
        if dummy_idx is None:
            continue
        rw2 = RWMol(cap_mol)
        rw2.ReplaceAtom(dummy_idx, Atom('C'))
        try:
            _C.SanitizeMol(rw2)
        except Exception:
            continue
        candidates.append((cap_mol.GetNumAtoms(), abbr, rw2.GetMol()))
    candidates.sort(key=lambda x: -x[0])
    for _, abbr, qmol in candidates:
        if cap_frag.HasSubstructMatch(qmol, useChirality=False):
            return abbr
    return None


def smiles_to_cabiln_core(smiles: str):
    """Convert a peptide SMILES string to CABILN notation.

    Returns (cabiln_str, match_details) where match_details is a list of
    (abbr, n_matched, n_total_atoms) per residue in chain order.

    Uses library-first backbone detection: matches all backbone monomers against
    the full input mol, builds an amide-bond DAG, and finds the longest chain.
    This handles β-AAs, peptoids, and any future SDF-registered monomer types.
    """
    from rdkit import Chem as _C
    from collections import deque as _dq

    mol = _C.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smiles[:80]}')

    # ── 1. Library-first backbone detection ──────────────────────────────────
    cov_lib = _s2c_get_coverage_lib()
    placements = _s2c_place_monomers(mol, cov_lib)
    # iso_placements: like placements but uses the iso SMARTS (carboxyl cap
    # removed) for monomers that have one.  Used only for branch detection so
    # that isopeptide-bonded monomers (E_g→K) are found without pulling the
    # backbone amide-N into the match set.
    iso_placements: list = []
    for _ie in cov_lib:
        if _ie.smarts_mol_iso is not None:
            _ism, _ini, _ioi = (_ie.smarts_mol_iso,
                                _ie.smarts_in_n_idx_iso,
                                _ie.smarts_out_co_idx_iso)
        else:
            _ism, _ini, _ioi = (_ie.smarts_mol,
                                _ie.smarts_in_n_idx,
                                _ie.smarts_out_co_idx)
        for _imatch in mol.GetSubstructMatches(_ism, useChirality=False):
            iso_placements.append(_PlacedNode(
                abbr=_ie.abbr, m_type=_ie.m_type,
                mol_atoms=frozenset(_imatch),
                in_n=_imatch[_ini],
                out_co=_imatch[_ioi],
                entry=_ie, match_tuple=_imatch,
            ))
    backbone = _s2c_build_backbone(mol, placements)
    if not backbone:
        raise ValueError('Cannot detect peptide backbone (no library monomers matched)')

    n = len(backbone)
    backbone_atoms_all: set = set()
    for node in backbone:
        backbone_atoms_all.update(node.mol_atoms)

    # ── 2. Cyclic detection ───────────────────────────────────────────────────
    cyclic = mol.GetBondBetweenAtoms(backbone[-1].out_co, backbone[0].in_n) is not None

    # ── 3. Cap detection ──────────────────────────────────────────────────────
    n_cap: str | None = None
    c_cap: str | None = None
    if not cyclic:
        # N-cap: carbonyl-bearing C bonded to backbone[0].in_n, outside backbone
        for nb in mol.GetAtomWithIdx(backbone[0].in_n).GetNeighbors():
            nbi = nb.GetIdx()
            if nbi not in backbone_atoms_all and nb.GetSymbol() == 'C':
                has_co = any(
                    x.GetSymbol() == 'O' and
                    mol.GetBondBetweenAtoms(nbi, x.GetIdx()).GetBondTypeAsDouble() == 2.0
                    for x in nb.GetNeighbors()
                )
                if has_co:
                    n_cap = _s2c_identify_n_cap_from_atom(
                        mol, backbone[0].in_n, nbi, backbone_atoms_all
                    ) or 'ac'
                    break
        # C-cap: any non-backbone substituent on out_co (am, NHEt, OEt, OtBu, …)
        for nb in mol.GetAtomWithIdx(backbone[-1].out_co).GetNeighbors():
            nbi = nb.GetIdx()
            if nbi in backbone_atoms_all:
                continue
            # Guard: reject peptide-bond Ns (alpha-N of the next residue whose
            # Cα is itself bonded to a carbonyl C).  A real cap N has no such
            # carbonyl-bearing C neighbour on its non-out_co side.
            if nb.GetSymbol() == 'N':
                is_peptide_n = any(
                    any(x.GetSymbol() == 'O' and
                        mol.GetBondBetweenAtoms(c.GetIdx(), x.GetIdx()).GetBondTypeAsDouble() == 2.0
                        for x in c.GetNeighbors())
                    for c in nb.GetNeighbors()
                    if c.GetIdx() != backbone[-1].out_co and c.GetSymbol() == 'C'
                )
                if is_peptide_n:
                    continue
            identified = _s2c_identify_c_cap_from_atom(
                mol, backbone[-1].out_co, nbi, backbone_atoms_all
            )
            if identified:
                c_cap = identified
                break
            # Fallback: plain am (NH2 with no additional C)
            if nb.GetSymbol() == 'N':
                c_nbrs = [x for x in nb.GetNeighbors() if x.GetSymbol() == 'C']
                if len(c_nbrs) <= 1:
                    c_cap = 'am'
                    break

    # ── 4. Match each backbone residue with stereo-aware matching ────────────
    from rdkit.Chem import AllChem as _AC3
    _AC3.AssignStereochemistry(mol, cleanIt=True, force=True)
    lib = _s2c_get_lib()
    details = []
    abbrs = []
    for pos, node in enumerate(backbone):
        atom_list = list(node.mol_atoms)
        # For free-acid C-terminus (no am cap): include the carboxyl OH in the
        # atom list so isolation gives COOH not CHO (which would match Gly_al).
        if not cyclic and pos == n - 1 and c_cap is None:
            for _nb in mol.GetAtomWithIdx(node.out_co).GetNeighbors():
                _bd = mol.GetBondBetweenAtoms(node.out_co, _nb.GetIdx())
                if (_nb.GetIdx() not in backbone_atoms_all and
                        _nb.GetSymbol() == 'O' and
                        _nb.GetTotalNumHs() > 0 and
                        _bd.GetBondTypeAsDouble() == 1.0):
                    atom_list.append(_nb.GetIdx())
                    break
        aas, _ = _s2c_isolate_residues(mol, [atom_list])
        aa_mol = aas[0]
        # Terminal residues: strip cap artifact introduced by isolate_residues
        # (the cap atoms are outside node.mol_atoms so the isolated mol already
        # has the correct N—H / C=O—OH termini; only strip if isolation added
        # a spurious cap atom from the coverage-lib match itself).
        if not cyclic:
            if pos == 0 and n_cap is not None:
                stripped, did = _s2c_strip_n_cap(aa_mol)
                if did:
                    aa_mol = stripped
            if pos == n - 1 and c_cap is not None:
                stripped, did = _s2c_strip_c_cap(aa_mol)
                if did:
                    aa_mol = stripped
        # Filter lib to entries with the same backbone path length as the
        # backbone-chain placement.  This ensures E_g is preferred over E when
        # the gamma-COOH was the actual chain exit (bb_dist=4 vs 2).  Falls
        # back to the full lib if nothing matches (e.g. novel unregistered dist).
        _node_bb_dist = node.entry.bb_dist
        _lib_filtered = [e for e in lib if e[6] == _node_bb_dist]
        abbr, score = _s2c_match(aa_mol, _lib_filtered or lib)
        if abbr is None:
            raise ValueError(
                f'Unrecognised monomer at backbone position {pos} '
                f'({aa_mol.GetNumAtoms()} heavy atoms); '
                f'register it in the monomer library first'
            )
        details.append((abbr, score, aa_mol.GetNumAtoms()))
        abbrs.append(abbr)

    # ── 5. Branch detection and formatting ────────────────────────────────────
    # Non-backbone placements whose out_co bonds to a backbone atom are branch
    # anchors (e.g. E_g attached to K's epsilon-N via isopeptide bond).
    # Deduplicate branch placements by atom set: L/D pairs and similar monomers
    # (DGlu/E_g/E) cover identical heavy atoms.  When multiple abbrs match
    # the same atoms, prefer standard E/D (whose R4=γ/β-COOH gives E(4,4) etc.)
    # over gamma-form E_g at the isopeptide junction; E_g is preferred only in
    # orphan-segment matching (step 4 of _s2c_walk_branch) via the 1-2 rule.
    # Priority 0 = most preferred; unlisted abbrs get 999.
    _ISOPEP_PRIORITY: dict = {
        'E': 0, 'dE': 1,               # standard Glu: slot annotation e.g. E(4,4)
        'D': 2, 'dD': 3,               # standard Asp
        'DGlu': 4, 'DGlu_g': 5,        # D-Glu variants
        'E_g': 6,                       # gamma-backbone Glu (lower — E preferred at isopeptide junction)
    }
    _bp_by_key: dict = {}  # frozenset(mol_atoms) → best _PlacedNode
    for p in iso_placements:
        if p.mol_atoms & backbone_atoms_all:
            continue  # overlaps backbone — skip
        key = frozenset(p.mol_atoms)
        if key not in _bp_by_key:
            _bp_by_key[key] = p
        else:
            # Keep the placement with higher isopeptide priority
            cur_pri = _ISOPEP_PRIORITY.get(_bp_by_key[key].abbr, 999)
            new_pri = _ISOPEP_PRIORITY.get(p.abbr, 999)
            if new_pri < cur_pri:
                _bp_by_key[key] = p
    branch_placements = list(_bp_by_key.values())
    branch_junctions: dict = {}  # backbone_pos → list of {'main_at', 'glu_node'}

    for bi, bb_node in enumerate(backbone):
        _seen_main_ats: set = set()  # one branch per unique junction atom
        for bp in branch_placements:
            main_at = -1
            # Check ALL atoms of bp for a bond to any backbone atom — not just
            # bp.out_co, because E_g's isopeptide attachment is its gamma-C=O
            # (in mol_atoms) rather than its alpha-C=O (out_co).
            for bpa in bp.mol_atoms:
                for ba in bb_node.mol_atoms:
                    if mol.GetBondBetweenAtoms(bpa, ba) is not None:
                        main_at = ba
                        break
                if main_at != -1:
                    break
            if main_at != -1 and main_at not in _seen_main_ats:
                _seen_main_ats.add(main_at)
                branch_junctions.setdefault(bi, []).append(
                    {'main_at': main_at, 'glu_node': bp}
                )

    if branch_junctions:
        raw_lib = _s2c_get_raw_lib()
        full_lib = _s2c_get_full_lib()
        branch_brackets: dict = {}

        for bi, junction_list in branch_junctions.items():
            anchor_abbr = abbrs[bi]
            bracket_acc = ''
            for jinfo in junction_list:
                main_at = jinfo['main_at']
                glu_node = jinfo['glu_node']
                glu_atoms = set(glu_node.mol_atoms)
                # in_n (alpha-N) is the junction from E_g toward AEEA/C20FA —
                # AEEA's out_co forms an amide bond with E_g's in_n.
                # out_co (alpha-C=O) is a free carboxyl terminus in the assembled
                # molecule when E_g is used as a branch (isopeptide via R4).
                glu_outgoing_n = glu_node.in_n

                # Orphan atoms: reachable from glu_outgoing_n, outside backbone+glu
                orphan_atoms: set = set()
                vis_o = backbone_atoms_all | glu_atoms
                q_o = _dq()
                for nb in mol.GetAtomWithIdx(glu_outgoing_n).GetNeighbors():
                    ni = nb.GetIdx()
                    if ni not in vis_o:
                        q_o.append(ni)
                while q_o:
                    ai = q_o.popleft()
                    if ai in vis_o:
                        continue
                    vis_o.add(ai)
                    orphan_atoms.add(ai)
                    for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                        if nb.GetIdx() not in vis_o:
                            q_o.append(nb.GetIdx())

                chain = _s2c_walk_branch(
                    mol, glu_atoms, glu_outgoing_n, orphan_atoms,
                    full_lib, raw_lib, anchor_abbr, main_at,
                    junction_abbr=glu_node.abbr,
                )
                if chain:
                    parts = '.'.join(f'{a}({p},{c})' for a, p, c in chain)
                    bracket_acc += f'.[{parts}]'
            if bracket_acc:
                branch_brackets[bi] = bracket_acc

        abbrs = [a + branch_brackets.get(i, '') for i, a in enumerate(abbrs)]

    # ── 5.3 Generic sidechain cap detection ──────────────────────────────────
    # After multi-monomer branch detection, atoms not in any backbone residue
    # or branch anchor (e.g. trt on Cys thiol, boc on Lys ε-N) are matched
    # against the full library as single-monomer sidechain caps.
    #
    # _sc_claimed tracks atoms already attributed so we don't double-annotate.
    _sc_claimed: set = set(backbone_atoms_all)
    for _jl in branch_junctions.values():
        for _j in _jl:
            _sc_claimed.update(_j['glu_node'].mol_atoms)

    _sc_unclaimed: set = set(range(mol.GetNumAtoms())) - _sc_claimed
    if _sc_unclaimed:
        from rdkit.Chem import RWMol as _RW2
        _sc_raw = _s2c_get_raw_lib()
        _sc_full = _s2c_get_full_lib()

        for _bi, _node in enumerate(backbone):
            _base_abbr = details[_bi][0]
            # Find all junction points: backbone atoms bonded to unclaimed atoms
            for _ai in list(_node.mol_atoms):
                for _nb in mol.GetAtomWithIdx(_ai).GetNeighbors():
                    _ni = _nb.GetIdx()
                    if _ni not in _sc_unclaimed:
                        continue
                    # BFS to collect the full connected unclaimed cluster
                    _cluster: set = set()
                    _q2 = _dq([_ni])
                    while _q2:
                        _ci = _q2.popleft()
                        if _ci in _sc_claimed or _ci in _cluster:
                            continue
                        _cluster.add(_ci)
                        for _cnb in mol.GetAtomWithIdx(_ci).GetNeighbors():
                            if _cnb.GetIdx() not in _sc_claimed:
                                _q2.append(_cnb.GetIdx())
                    if not _cluster:
                        continue
                    # Build a fragment mol from the cluster
                    _rw_cap = _RW2()
                    _cap_s2n: dict = {}
                    for _cai in sorted(_cluster):
                        _cap_s2n[_cai] = _rw_cap.AddAtom(
                            _C.Atom(mol.GetAtomWithIdx(_cai).GetAtomicNum()))
                    for _bond in mol.GetBonds():
                        _bai, _eai = _bond.GetBeginAtomIdx(), _bond.GetEndAtomIdx()
                        if _bai in _cap_s2n and _eai in _cap_s2n:
                            _rw_cap.AddBond(_cap_s2n[_bai], _cap_s2n[_eai],
                                            _bond.GetBondType())
                    _cap_frag = _C.RemoveHs(_rw_cap.GetMol())
                    if _cap_frag is None or _cap_frag.GetNumAtoms() == 0:
                        continue
                    # Match against library
                    _cap_abbr, _ = _s2c_match_branch_segment(_cap_frag, _sc_full)
                    if not _cap_abbr or _cap_abbr == '?':
                        continue
                    # Skip scaffold monomers — they're handled in step 5.5
                    if _cap_abbr in {a for a, *_ in _s2c_get_scaffold_patterns()}:
                        continue
                    # Determine R-groups at the junction
                    _res_r = _s2c_r_group_at_atom(
                        mol, _ai, _node.mol_atoms, _base_abbr, _sc_raw)
                    _cap_r = _s2c_r_group_at_atom(
                        mol, _ni, _cluster, _cap_abbr, _sc_raw)
                    # Only accept sidechain attachments (R4+).
                    # R1/R2/R3 are backbone N/C connections — those are terminal
                    # caps (ac, am, fmoc, boc) already handled by cap stripping.
                    if _res_r is None or _res_r < 4 or _cap_r is None:
                        continue
                    abbrs[_bi] += f'.{_cap_abbr}({_res_r},{_cap_r})'
                    _sc_claimed.update(_cluster)
                    _sc_unclaimed -= _cluster
                    break  # one cap per junction atom; next junction

    # ── 5.5 Scaffold crosslink detection ─────────────────────────────────────
    # Detect multi-arm non-amide crosslink scaffolds (TBMB, etc.) library-
    # driven: any SDF monomer with ≥2 non-backbone dummy attachment points.
    # For each pattern, [n*] in its CHUCKLES was converted to [*:n] SMARTS;
    # matched atoms at those positions must be backbone atoms.  R-group labels
    # are assigned in backbone chain order (!1, !2, …) for deterministic output.
    scaffold_suffix = ''
    _scaffold_xlinks = 0
    _scaffold_patterns = _s2c_get_scaffold_patterns()

    if _scaffold_patterns:
        _rlib = _s2c_get_raw_lib()
        for _sc_abbr, _sc_smarts, _sc_r_pos, _sc_min_slot in _scaffold_patterns:
            _found = False
            for _match in mol.GetSubstructMatches(_sc_smarts, useChirality=False):
                # Collect R-group positions that map to backbone atoms.
                # Allow partial matches: require ≥2 backbone hits so a scaffold
                # with one unreacted arm is still detected.
                _attach: dict = {}  # smarts_idx → (backbone_pos, r_num)
                for _si, _rnum in _sc_r_pos.items():
                    _mat_idx = _match[_si]
                    _bp = next(
                        (pos for pos, node in enumerate(backbone)
                         if _mat_idx in node.mol_atoms),
                        None,
                    )
                    if _bp is not None:
                        _attach[_si] = (_bp, _rnum)

                _min_arms = max(2, len(_sc_r_pos) - 1)
                if len(_attach) < _min_arms:
                    # Relax the threshold when unreacted arms hit halogen
                    # leaving-group atoms (Br/Cl/I) still present in the SMILES.
                    # A scaffold arm that terminated at a halogen is unambiguously
                    # part of this scaffold — no peptide context produces that.
                    if len(_attach) >= 1:
                        _halogen_nos = {35, 17, 53}  # Br, Cl, I
                        _unmatched = [_match[_si] for _si in _sc_r_pos
                                      if _si not in _attach]
                        if all(mol.GetAtomWithIdx(_ai).GetAtomicNum() in _halogen_nos
                               for _ai in _unmatched):
                            pass  # accept — all unreacted arms are halogens
                        else:
                            continue

                # Sort attachments by backbone chain order.  Renumber R-groups
                # consecutively from the scaffold's minimum slot (e.g. R4,R5,R6
                # for TBMB) so partial matches (2-of-3 arms) don't produce
                # non-consecutive slot numbers like (4,6) when (4,5) is canonical.
                _sorted_bpos = sorted(_attach.values(), key=lambda x: x[0])
                _sorted_r = list(range(_sc_min_slot, _sc_min_slot + len(_sorted_bpos)))
                # Scaffold XL IDs must not reuse !1 when the backbone ring
                # closure already claims it. Start from 2 for cyclic peptides.
                _xl_id = 2 if cyclic else 1
                _tags = []
                for (_bp, _), _scaffold_r in zip(_sorted_bpos, _sorted_r):
                    _res_r = _s2c_crosslink_r(details[_bp][0], _rlib)
                    _tags.append(f'!{_xl_id}')
                    abbrs[_bp] += f'.!{_xl_id}({_res_r},{_scaffold_r})'
                    _xl_id += 1
                scaffold_suffix = f'%{_sc_abbr}.' + '.'.join(_tags)
                # Track the highest XL ID used so sidechain crosslinks start
                # after it (the sidechain counter uses _scaffold_xlinks + 1).
                _scaffold_xlinks = _xl_id - 1
                _found = True
                break

            if _found:
                break  # one scaffold per molecule

    # ── 6. Crosslink detection ────────────────────────────────────────────────
    # Side-chain bonds between pairs of backbone residues that are not the
    # adjacent backbone amide bonds (out_co[i] → in_n[i+1]).
    branch_node_atoms: set = set()
    for jl in branch_junctions.values():
        for jinfo in jl:
            branch_node_atoms.update(jinfo['glu_node'].mol_atoms)

    sc_pairs: list = []
    backbone_amide_pairs = {
        frozenset((backbone[i].out_co, backbone[i + 1].in_n))
        for i in range(n - 1)
    }
    if cyclic:
        backbone_amide_pairs.add(frozenset((backbone[-1].out_co, backbone[0].in_n)))
    for bi in range(n):
        for bj in range(bi + 1, n):
            ni_atoms = backbone[bi].mol_atoms
            nj_atoms = backbone[bj].mol_atoms
            for bond in mol.GetBonds():
                ba, ea = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if (ba in ni_atoms and ea in nj_atoms) or (ba in nj_atoms and ea in ni_atoms):
                    if frozenset((ba, ea)) not in backbone_amide_pairs:
                        sc_pairs.append((bi, bj))
                    break

    if sc_pairs:
        if not branch_junctions:
            raw_lib = _s2c_get_raw_lib()
        # !1 is reserved for the backbone ring-closure marker when cyclic.
        # Sidechain crosslinks must not reuse it, so start from at least 2.
        xlink_ctr = max(1 if cyclic else 0, _scaffold_xlinks) + 1
        for bi, bj in sc_pairs:
            base_i = details[bi][0]
            base_j = details[bj][0]
            r_i = _s2c_crosslink_r(base_i, raw_lib)
            r_j = _s2c_crosslink_r(base_j, raw_lib)
            tag = f'!{xlink_ctr}'
            xlink_ctr += 1
            abbrs[bi] += f'.{tag}({r_i},{r_j})'
            abbrs[bj] += f'.{tag}({r_j},{r_i})'

    # ── 7. Build CABILN string ────────────────────────────────────────────────
    if cyclic:
        cabiln = '!1-' + '-'.join(abbrs) + '-!1'
    else:
        cabiln = '-'.join(abbrs)
        if c_cap:
            cabiln += '-' + c_cap
    if scaffold_suffix:
        cabiln += scaffold_suffix
    if n_cap:
        cabiln = n_cap + '-' + cabiln
    return cabiln, details


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
    _mck = ("monomer", abbr, width, height)
    _mhit = _rc_get(_mck)
    if _mhit is not None:
        return _mhit
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
            _rc_put(_mck, result)
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

def _renumber_xlinks(cabiln: str) -> str:
    """Renumber !X tags so they appear as !1, !2, … in left-to-right reading order."""
    import re as _re2
    seen: dict = {}
    counter = 1
    for m in _re2.finditer(r'!(\d+)', cabiln):
        n = int(m.group(1))
        if n not in seen:
            seen[n] = counter
            counter += 1
    if not seen or seen == {i: i for i in seen}:
        return cabiln  # already in order
    # Two-pass to avoid collision (e.g. !2→!1 while !1 still exists)
    result = cabiln
    for old, new in seen.items():
        result = result.replace(f'!{old}', f'!__T{new}__')
    result = _re2.sub(r'!__T(\d+)__', lambda m: f'!{m.group(1)}', result)
    return result


@app.post("/convert_notation")
async def convert_notation(req: _ConvertReq):
    from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch
    try:
        if req.target == 'bracket':
            return {"result": _renumber_xlinks(cabiln_to_bracket(req.cabiln))}
        elif req.target == 'branch':
            result = cabiln_to_branch(req.cabiln)
            return {"result": _renumber_xlinks(result)}
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
                seq = Sequence(_to_bracket(cabiln))
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
        # Normalize % branch notation to bracket form so all downstream string
        # scanning (bracket spans, crosslink annotation) works uniformly.
        cabiln = _to_bracket(cabiln)

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

            def _annotate_percent_residue(cabiln, res_idx, ann, seq, main_set):
                """Append ann to the right occurrence of a residue inside a % branch."""
                try:
                    abbr = seq.s_monomers[res_idx].get('m_abbr', '')
                except (IndexError, KeyError):
                    return cabiln
                pendant_indices = [i for i in range(len(seq.s_monomers))
                                   if i not in main_set]
                same_abbr = [i for i in pendant_indices
                             if seq.s_monomers[i].get('m_abbr', '') == abbr]
                try:
                    target_occ = same_abbr.index(res_idx)
                except ValueError:
                    return cabiln
                main_part, *branches = cabiln.split('%')
                found = 0
                for bi, branch in enumerate(branches):
                    tokens = branch.split('-')
                    for ti, tok in enumerate(tokens):
                        dot_pos = tok.find('.')
                        base = tok[:dot_pos] if dot_pos != -1 else tok
                        if base == abbr:
                            if found == target_occ:
                                tokens[ti] = tok + ann
                                branches[bi] = '-'.join(tokens)
                                return main_part + '%' + '%'.join(branches)
                            found += 1
                return cabiln

            ann_host   = f'.{_tag}({req.r_host},{req.r_new})'
            ann_target = f'.{_tag}({req.r_new},{req.r_host})'
            # In % branch notation crosslinks use just the tag (R-groups are
            # encoded on the backbone side only, e.g. K.!2(4,4) implies C20FA.!2)
            ann_tag_only = f'.{_tag}'

            result = cabiln
            for _ridx, _ann in [(req.host_residue_idx, ann_host),
                                 (req.target_residue_idx, ann_target)]:
                if _main2 and _ridx in _main2:
                    result = _annotate_main(result, _ridx, _ann)
                elif _seq2 is not None:
                    prev = result
                    result = _annotate_pendant(result, _ridx, _ann, _seq2, _main2)
                    if result == prev and '%' in result:
                        # Residue lives in a % branch — annotate there
                        result = _annotate_percent_residue(
                            result, _ridx, ann_tag_only, _seq2, _main2)
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
            seq = Sequence(_to_bracket(cabiln))
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


class _InsertBackboneReq(BaseModel):
    cabiln: str
    after_idx: int
    new_abbr: str

@app.post("/insert_backbone")
async def insert_backbone(req: _InsertBackboneReq):
    """Insert new_abbr into the main chain immediately after the residue at global index after_idx."""
    def _tok_split(s):
        toks, depth, cur = [], 0, ''
        for ch in s:
            if ch in '[{': depth += 1
            elif ch in ']}': depth -= 1
            if ch == '-' and depth == 0:
                toks.append(cur); cur = ''
            else:
                cur += ch
        if cur: toks.append(cur)
        return toks

    try:
        cabiln = req.cabiln.strip()
        if not cabiln:
            return {"result": req.new_abbr}

        from pyPept.sequence import Sequence
        seq = Sequence(_to_bracket(cabiln))
        chain_ids = seq.s_chains.get('s_monomerIDs', [])
        if not chain_ids:
            return {"error": "No chains found"}

        main_chain = chain_ids[0]
        if req.after_idx not in main_chain:
            return {"error": f"Residue {req.after_idx} is not in the main chain"}
        pos = main_chain.index(req.after_idx)

        # Split branch suffix (% notation) so we only modify the main chain
        pct = cabiln.find('%')
        main_str = cabiln[:pct] if pct >= 0 else cabiln
        suffix   = cabiln[pct:] if pct >= 0 else ''

        toks = _tok_split(main_str)
        toks.insert(pos + 1, req.new_abbr)
        return {"result": '-'.join(toks) + suffix}

    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=500)


class _SmilesToCabilnReq(BaseModel):
    smiles: str
    notation: str = 'percent'  # 'percent' (%) or 'bracket' ([])

@app.post("/smiles_to_cabiln")
async def smiles_to_cabiln_endpoint(req: _SmilesToCabilnReq):
    """Convert a peptide SMILES to CABILN notation using pyPept's monomer library."""
    try:
        from rdkit import Chem as _C
        smi = req.smiles.strip()
        cabiln, details = smiles_to_cabiln_core(smi)
        if req.notation == 'bracket':
            from pyPept.sequence import cabiln_to_bracket
            cabiln = _renumber_xlinks(cabiln_to_bracket(cabiln))
        warning = None
        mol_in = _C.MolFromSmiles(smi)
        if mol_in is not None and _count_defined_stereo(mol_in) == 0:
            warning = "Achiral SMILES detected — chirality inferred from library (L-form assumed). Result may not match intended stereoisomer."
        resp = {
            "cabiln": cabiln,
            "details": [{"abbr": a, "score": s, "total": t} for a, s, t in details],
        }
        if warning:
            resp["warning"] = warning
        return resp
    except Exception as exc:
        return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)


class _ToCabilnReq(BaseModel):
    input: str


@app.post("/to_cabiln")
async def to_cabiln_endpoint(req: _ToCabilnReq):
    """Convert SMILES, BILN, or HELM input to CABILN notation."""
    from rdkit import Chem as _C
    txt = req.input.strip()

    # 1. Try SMILES
    mol = _C.MolFromSmiles(txt)
    if mol is not None:
        try:
            cabiln, details = smiles_to_cabiln_core(txt)
            return {"cabiln": cabiln, "from": "SMILES",
                    "details": [{"abbr": a, "score": s, "total": t} for a, s, t in details]}
        except Exception as exc:
            return JSONResponse({"error": str(exc).split('\n')[0]}, status_code=400)

    # 2. Try HELM
    if 'PEPTIDE' in txt.upper() and '$' in txt:
        try:
            import sys as _sys, pathlib as _pl
            _sys.path.insert(0, str(_pl.Path(__file__).parent.parent / 'src'))
            from pyPept.converter import Converter
            from pyPept.sequence import biln_to_cabiln
            biln = Converter(helm=txt).get_biln()
            cabiln = biln_to_cabiln(biln)
            # Verify it assembles
            from pyPept.molecule import Molecule
            from pyPept.sequence import Sequence
            Molecule(Sequence(cabiln)).get_molecule(fmt='ROMol')
            return {"cabiln": cabiln, "from": "HELM", "details": []}
        except Exception as exc:
            return JSONResponse({"error": f"HELM parse failed: {exc}".split('\n')[0]},
                                status_code=400)

    # 3. Try BILN
    try:
        from pyPept.sequence import biln_to_cabiln, Sequence
        from pyPept.molecule import Molecule
        cabiln = biln_to_cabiln(txt)
        Molecule(Sequence(cabiln)).get_molecule(fmt='ROMol')
        return {"cabiln": cabiln, "from": "BILN", "details": []}
    except Exception as exc:
        pass

    return JSONResponse(
        {"error": "Could not parse input as SMILES, BILN, or HELM"},
        status_code=400)


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
    w, h = max(400, req.width), max(300, req.height)
    _ck = (req.cabiln, w, h, req.seed)
    _hit = _rc_get(_ck)
    if _hit is not None:
        return _hit
    try:
        from rdkit.Chem.Descriptors import ExactMolWt
        from pyPept.molecule import Molecule
        from pyPept.sequence import Sequence

        seq   = Sequence(_to_bracket(req.cabiln))
        mol   = Molecule(seq)
        romol = mol.get_molecule(fmt='ROMol')
        if romol is None:
            return JSONResponse({"error": "Assembly produced no molecule"}, status_code=400)

        svg   = _draw_mol(romol, w, h, seed=req.seed)
        block = _mol_block(romol)

        res_map = mol.get_residue_atom_map()
        residues = [{"idx": i, "abbr": m.get("m_abbr", f"?{i}")}
                    for i, m in enumerate(seq.s_monomers)]

        chain_ids = seq.s_chains.get('s_monomerIDs', [])
        chains = [{"idx": ci, "residues": ids} for ci, ids in enumerate(chain_ids)]

        crosslink_groups = _build_crosslink_groups(seq)
        bracket_groups = _build_bracket_groups(seq, chain_ids, req.cabiln, crosslink_groups)

        result = {"svg": svg, "mol_block": block,
                  "info": f"{romol.GetNumAtoms()} atoms · MW {ExactMolWt(romol):.2f}",
                  "residue_map": {str(k): v for k, v in res_map.items()},
                  "residues": residues,
                  "chains": chains,
                  "bracket_groups": bracket_groups,
                  "crosslink_groups": crosslink_groups,
                  "cabiln_echo": req.cabiln}
        _rc_put(_ck, result)
        return result

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
            seq = Sequence(_to_bracket(txt))
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

        seq        = Sequence(_to_bracket(req.cabiln))
        mol        = Molecule(seq)
        cabiln_mol = mol.get_molecule(fmt='ROMol')
        if cabiln_mol is None:
            return JSONResponse({"error": "CABILN assembly produced no molecule"}, status_code=400)

        smi_canon    = _canon(smiles_mol)
        cabiln_canon = _canon(cabiln_mol)

        warning = None
        if _count_defined_stereo(smiles_mol) == 0:
            # Flat/achiral input — strip stereo from assembled before comparing
            # so connectivity can still be verified meaningfully.
            warning = "Achiral SMILES detected — stereo stripped from CABILN output for comparison. Connectivity only."
            match = _canon_flat(smiles_mol) == _canon_flat(cabiln_mol)
        else:
            # Primary match via InChI (handles tautomers like His imidazole ring).
            # Falls back to canonical SMILES comparison if InChI unavailable.
            try:
                from rdkit.Chem.inchi import MolToInchi
                match = MolToInchi(smiles_mol) == MolToInchi(cabiln_mol)
            except Exception:
                match = smi_canon == cabiln_canon

        resp = {
            "match":            match,
            "smiles_canonical": smi_canon,
            "cabiln_canonical": cabiln_canon,
        }
        if warning:
            resp["warning"] = warning
        return resp

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
    try:
        import sys as _sys, pathlib as _pl
        _sys.path.insert(0, str(_pl.Path(__file__).resolve().parent.parent / 'src'))
        from rdkit.Chem import rdDepictor as _, MolToMolBlock as __
        from rdkit.Chem.Draw import rdMolDraw2D as ___
        from rdkit.Chem.Descriptors import ExactMolWt as ____
        from pyPept.sequence import Sequence as _____
        from pyPept.molecule import Molecule as ______
        print(f"  Imports pre-warmed in {_t.perf_counter()-_t0:.1f}s")
    except Exception as _e:
        print(f"  Import warmup partial: {_e}")
    print(f"  Local     ->  http://localhost:{_port}")
    print("  Tailscale ->  http://100.119.0.78:8732")
    if _port == 8732:
        threading.Thread(target=_open_browser, daemon=True).start()
    uvicorn.run(app, host="0.0.0.0", port=_port, log_level="warning")
