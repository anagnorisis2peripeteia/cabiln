// Monomer discovery, library revisions, filtering and hover previews.
class MonomerLibrary {
  constructor({ getFilters, onUse, onChanged }) {
    this.getFilters = getFilters;
    this.onUse = onUse;
    this.onChanged = onChanged;
    this.loaded = false;
    this.version = null;
    this.monomers = [];
    this.previewCache = {};
    this.previewTimer = null;
    this.reactionPairs = null;
    this.filterActive = false;
    this.panel = document.getElementById('lib-panel');
    this.search = document.getElementById('lib-search');
    this.closeButton = document.getElementById('lib-close');
    this.list = document.getElementById('lib-list');
    this.count = document.getElementById('lib-count');
    this.status = document.getElementById('lib-status');
    this.preview = document.getElementById('lib-preview');
    this.button = document.getElementById('btn-lib');
    this.filterButton = document.getElementById('btn-rxn-filter');
    this.disclosure = createPanel({ panel: this.panel, button: this.button,
      closeButton: this.closeButton, focus: this.search,
      onChange: open => {
        if (open) { this.load(); this.loadReactions(); }
        else this.hidePreview();
      },
    });

    this.search.addEventListener('input', () => this.render());

    this.filterButton.addEventListener('click', async () => {
      if (!Array.isArray(this.reactionPairs)) {
        await this.loadReactions();
        if (!Array.isArray(this.reactionPairs) || !this.getFilters().left) return;
      }
      this.filterActive = !this.filterActive;
      this.filterButton.classList.toggle('active', this.filterActive);
      this.filterButton.setAttribute('aria-pressed', String(this.filterActive));
      this.render();
    });

    window.addEventListener('focus', () => {
      if (this.panel.classList.contains('open')) this.load();
    });

    // Keep one set of listeners while search replaces the rows.
    this.list.addEventListener('click', event => {
      const row = event.target.closest('.lib-row');
      if (row) this.onUse(row.dataset.abbr, !!event.target.closest('.lib-use'));
    });
    this.list.addEventListener('contextmenu', event => {
      const row = event.target.closest('.lib-row');
      if (row && this.getFilters().building) {
        event.preventDefault();
        this.onUse(row.dataset.abbr, true);
      }
    });
    this.list.addEventListener('keydown', event => {
      if (event.target.matches('.lib-row') && ['Enter', ' '].includes(event.key)) {
        event.preventDefault();
        this.onUse(event.target.dataset.abbr);
      }
    });
    for (const type of ['mouseover', 'focusin']) {
      this.list.addEventListener(type, event => {
        const row = event.target.closest('.lib-row');
        if (row && !row.contains(event.relatedTarget)) this.startPreview(row.dataset.abbr, row);
      });
    }
    for (const type of ['mouseout', 'focusout']) {
      this.list.addEventListener(type, event => {
        const row = event.target.closest('.lib-row');
        if (row && !row.contains(event.relatedTarget)) this.hidePreview();
      });
    }
  }

  open(options) { this.disclosure.open(options); }
  close() { this.disclosure.close(); }

  async load() {
    const request = requests.start('library');
    if (!this.loaded) showLoading(this.list, 'Loading monomers…');
    try {
      const res = await fetchCalculation('/monomers', { signal: request.signal });
      const data = await readResponse(res);
      if (!request.current()) return;
      if (!Array.isArray(data)) throw new Error(data.error || 'Invalid monomer library');
      const version = res.headers?.get('X-Library-Version') || null;
      this.status.hidden = true;
      if (this.loaded && version && version === this.version) return;
      if (this.loaded) this.onChanged();
      this.version = version;
      this.monomers = data.map(monomer => ({ ...monomer,
        searchText: [monomer.abbr, monomer.name, monomer.type, monomer.chem_types].join(' ').toLowerCase(),
        attachmentTypes: this.attachmentTypes(monomer.chem_types),
      }));
      this.previewCache = {};
      this.hidePreview();
      this.loaded = true;
      this.render();
    } catch (e) {
      if (!request.current()) return;
      const target = this.loaded ? this.status : this.list;
      target.hidden = false;
      showRetry(target, this.loaded ? 'Could not refresh the library. Previous definitions are still shown.'
        : 'Could not load monomers.', () => this.load());
    } finally {
      request.finish();
    }
  }

  async loadReactions() {
    if (this.reactionPairs !== null) return;
    const request = requests.start('reactions');
    this.status.hidden = true;
    try {
      // Both endpoints share one calculation worker. Let the palette finish so
      // a quick reaction-list request cannot force it into a one-second retry.
      if (requests.has('library') && !await request.waitFor('library')) return;
      const res = await fetchCalculation('/reactions', { signal: request.signal });
      const data = await readResponse(res);
      if (!request.current()) return;
      if (!Array.isArray(data)) throw new Error(data.error || 'Reaction data is unavailable');
      this.reactionPairs = data;
      if (this.filterActive && this.loaded) this.render();
    } catch (error) {
      if (!request.current()) return;
      this.reactionPairs = null;
      showRetry(this.status, 'Reaction filter unavailable. ' + error.message, () => this.loadReactions());
      this.status.hidden = false;
    } finally { request.finish(); }
  }

  attachmentTypes(cts) {
    if (!cts) return [];
    return cts.split(',').map(p => p.includes(':') ? p.split(':')[1].trim() : p.trim()).filter(Boolean);
  }

  render(q = this.search.value.trim().toLowerCase()) {
    const { left: buildLeft, insertBetween: insertBetweenActive, replacements, awaitingSelection } = this.getFilters();
    const swapping = replacements !== null;
    const catalog = swapping ? replacements : this.monomers;
    let filtered = q
      ? catalog.filter(m => m.searchText.includes(q))
      : catalog;

    if (this.filterActive && buildLeft && Array.isArray(this.reactionPairs)) {
      const pairSet = new Set(this.reactionPairs.map(([a, b]) => a + '|' + b));
      let lcts;
      if (buildLeft.selectedSlot !== null) {
        const selRg = buildLeft.rgroups.find(r => r.slot === buildLeft.selectedSlot);
        lcts = selRg && !selRg.used ? [selRg.chem_type] : [];
      } else {
        lcts = buildLeft.rgroups.filter(r => !r.used).map(r => r.chem_type);
      }
      filtered = filtered.filter(m => {
        const mcts = m.attachmentTypes;
        return mcts.some(mct => lcts.some(lct =>
          pairSet.has(lct + '|' + mct) || pairSet.has(mct + '|' + lct)
        ));
      });
    }

    if (insertBetweenActive) {
      filtered = filtered.filter(m => m.backbone_insertable);
    }

    this.count.textContent = swapping
      ? `${filtered.length} / ${catalog.length} replacements with compatible sites`
      : `${filtered.length} / ${this.monomers.length} monomers`;

    if (!filtered.length) {
      this.list.innerHTML = `<div class="placeholder">${swapping && awaitingSelection
        ? 'Select a residue to find replacements' : 'No matches'}</div>`;
      return;
    }

    const rows = filtered.map(m => {
      const badge = m.degenerate
        ? `<span class="lib-badge cap">N/C cap</span>`
        : m.subtype === 'modified' || m.subtype === 'natural'
        ? `<span class="lib-badge aa">${escHtml(m.type)}</span>`
        : m.type === 'cap' && m.subtype === 'protecting'
        ? `<span class="lib-badge protect">cap</span>`
        : `<span class="lib-badge cap">${escHtml(m.type)}</span>`;

      const issues = Array.isArray(m.quality?.issues) ? m.quality.issues : [];
      const qualityText = issues.filter(issue => issue.severity !== 'info')
        .map(issue => issue.message || issue.code).filter(Boolean).join(' · ');
      let lg = m.leaving ? `  LG: ${escHtml(m.leaving)}` : '';
      if (m.degenerate) {
        const parts = [];
        if (m.nterm_abbr) parts.push(`N: ${escHtml(m.nterm_abbr)} (${escHtml(m.nterm_leaving)})`);
        if (m.cterm_abbr) parts.push(`C: ${escHtml(m.cterm_abbr)} (${escHtml(m.cterm_leaving)})`);
        lg = '  ' + parts.join(' | ');
      }
      const label = m.abbr + ': ' + m.name + (qualityText ? '. Library quality: ' + qualityText : '');
      return `<div class="lib-row" data-abbr="${escAttr(m.abbr)}" tabindex="0" aria-label="${escAttr(label)}">
        <div class="lib-abbr">${escHtml(m.abbr)}</div>
        <div class="lib-info">
          <div class="lib-name" title="${escAttr(m.name)}">${escHtml(m.name)}</div>
          <div class="lib-meta">${escHtml(m.chem_types || '')}${lg}</div>
          ${qualityText ? `<div class="lib-quality" title="${escAttr(qualityText)}">${escHtml(qualityText)}</div>` : ''}
        </div>
        ${badge}
        <button type="button" class="lib-use" aria-label="Use ${escAttr(m.abbr)} in builder" title="Choose this monomer in the builder">Use</button>
      </div>`;
    });
    this.list.innerHTML = rows.join('');

  }

  startPreview(abbr, row) {
    clearTimeout(this.previewTimer);
    const request = requests.start('monomer-preview');
    this.previewTimer = setTimeout(async () => {
      try {
        const data = this.previewCache[abbr] || await readResponse(await fetchCalculation(
          `/monomer_svg?abbr=${encodeURIComponent(abbr)}`, { signal: request.signal }
        ));
        if (!request.current()) return;
        if (data.svg) {
          this.previewCache[abbr] = data;
          this.showPreview(data, row);
        }
      } catch (e) { /* Hover previews are optional. */ }
      finally { request.finish(); }
    }, 200);
  }

  showPreview(data, row) {
    let html;
    if (data.degenerate && data.variants) {
      let panels = data.variants.map(v => {
        let pane = `<div class="prev-pane"><span class="prev-label">${v.label}</span>${v.svg}`;
        if (v.reagent) {
          pane += `<span class="prev-rxn">${v.reagent.reaction} (LG: ${v.reagent.reagent_lg})</span>`;
        }
        pane += `</div>`;
        return pane;
      }).join('');
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
    const quality = data.quality || this.monomers.find(item => item.abbr === row.dataset?.abbr)?.quality;
    const issues = Array.isArray(quality?.issues) ? quality.issues : [];
    const notes = issues.filter(issue => issue.severity === 'info');
    const warnings = issues.filter(issue => issue.severity !== 'info');
    if (notes.length) html += `<div class="prev-meta">Library notes: ${notes.map(issue => escHtml(issue.message || issue.code)).join(' · ')}</div>`;
    if (warnings.length) html += `<div class="prev-meta prev-warn">Library quality: ${warnings.map(issue => escHtml(issue.message || issue.code)).join(' · ')}</div>`;
    this.preview.innerHTML = html;
    const hasReagent = !!(data.svg_reagent || (data.variants && data.variants.some(v => v.svg_reagent)));
    this.preview.classList.toggle('has-reagent', hasReagent);
    const rect = row.getBoundingClientRect();
    const previewH = hasReagent ? 240 : 202;
    let top = Math.max(8, rect.top - 40);
    if (top + previewH > window.innerHeight - 8) {
      top = window.innerHeight - 8 - previewH;
    }
    top = Math.max(8, top);
    this.preview.style.left = (rect.right + 8) + 'px';
    this.preview.style.top  = top + 'px';
    this.preview.style.display = 'block';
  }

  hidePreview() {
    clearTimeout(this.previewTimer);
    requests.cancel('monomer-preview');
    this.preview.style.display = 'none';
  }

  get isOpen() { return this.disclosure.isOpen; }

  setDark(dark) { this.preview.classList.toggle('dark', dark); }

  setFilterEnabled(enabled) { this.filterButton.disabled = !enabled; }

  resetFilter() {
    const active = this.filterActive;
    this.filterActive = false;
    this.setFilterEnabled(false);
    this.filterButton.classList.remove('active');
    this.filterButton.setAttribute('aria-pressed', 'false');
    return active;
  }
}
