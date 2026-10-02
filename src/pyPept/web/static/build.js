// Connection and replacement sessions own selection, validation and preview lifetimes.
function createBuildPanel({ library, residueView, getDocument, commitDocument, isStale, onApply }) {
  let connectionValid = false;
  let buildLeft    = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
  let buildRight   = null;  // { abbr, rgroups: [{slot, chem_type, used}], selectedSlot }
  let buildLeftRIdx = null;
  let buildRightRIdx = null;
  let buildReaction = '';
  let insertBetweenActive = false;
  let swapState = null;  // Current source, eligible definitions, reviewed mapping and preview.

  const buildPanel    = document.getElementById('build-panel');
  const btnBuild      = document.getElementById('btn-build');
  const buildClose    = document.getElementById('build-close');
  const buildConnect  = document.getElementById('build-connect');
  const buildSuggestion = document.getElementById('build-suggestion');
  const buildAction = document.getElementById('build-action');
  const swapMapping = document.getElementById('build-swap-mapping');
  const swapSites = document.getElementById('build-swap-sites');
  const buildPreviewButton = document.getElementById('build-preview-button');
  const buildPreviewPanel = document.getElementById('build-preview');
  const buildPreviewStatus = document.getElementById('build-preview-status');
  const buildPreviewReaction = document.getElementById('build-preview-reaction');
  const buildPreviewChanges = document.getElementById('build-preview-changes');
  const buildPreviewFocus = document.getElementById('build-preview-focus');
  const buildPreviewInner = document.getElementById('build-preview-inner');
  const buildPreviewSource = document.getElementById('build-preview-source');
  const buildPreviewNotation = document.getElementById('build-preview-notation');
  const buildStatus   = document.getElementById('build-status');
  const buildHint     = document.getElementById('build-hint');
  const buildLeftAbbr = document.getElementById('build-left-abbr');
  const buildLeftSvg  = document.getElementById('build-left-svg');
  const buildLeftRg   = document.getElementById('build-left-rgroups');
  const buildLeftSite = document.getElementById('build-left-site');
  const buildRightAbbr= document.getElementById('build-right-abbr');
  const buildRightSvg = document.getElementById('build-right-svg');
  const buildRightRg  = document.getElementById('build-right-rgroups');
  const buildRightSite = document.getElementById('build-right-site');
  const buildInsertRow = document.getElementById('build-insert-row');
  const buildInsertInfo = document.getElementById('build-insert-info');
  const buildInsertBtn = document.getElementById('build-insert-btn');
  const buildPreviewViewport = makeZoomable(document.getElementById('build-preview-canvas'), buildPreviewInner);
  buildPreviewFocus.addEventListener('click', () => buildPreviewViewport.focus(
    buildPreviewInner.querySelectorAll('.edit-residue, .edit-connection')));
  const previewPanel = createPanel({ panel: buildPreviewPanel, button: buildPreviewButton,
    closeButton: document.getElementById('build-preview-close'), toggle: false,
    onChange(open) { if (!open) discardPreview(); },
  });

  function useLibraryMonomer(abbr, explicit = false) {
    if (explicit && !panel.isOpen) panel.open();
    if (panel.isOpen && getDocument().notation !== 'cabiln') {
      buildHint.textContent = 'Convert this input to CABILN before building';
    } else if (panel.isOpen && insertBetweenActive) {
      doInsertBetween(abbr);
    } else if (panel.isOpen && !getDocument().text.trim()) {
      commitDocument(abbr, 'cabiln');
      buildHint.textContent = 'First monomer added — select its residue, then choose another monomer';
    } else if (panel.isOpen && isStale()) {
      buildHint.textContent = 'Wait for a valid drawing before choosing attachment sites';
    } else if (panel.isOpen && isSwapMode()) {
      chooseSwapMonomer(abbr);
    } else if (panel.isOpen) {
      loadBuildRight(abbr);
      if (!buildLeft) buildHint.textContent = 'Now select a residue in the sequence';
    } else {
      return false;
    }
    return true;
  }

  // build mode
  function selectBuildResidue(residue) {
    if (!panel.isOpen || isStale()) return;
    const right = !isSwapMode() && buildLeft && buildLeftRIdx !== residue.idx;
    const pending = right
      ? loadBuildRight(residue.abbr, residue.idx)
      : loadBuildLeft(residue.abbr, residue.idx);
    residueView.select(right ? buildLeftRIdx : residue.idx, right ? residue.idx : null);
    return pending;
  }

  const panel = createPanel({ panel: buildPanel, button: btnBuild, closeButton: buildClose,
    focus: buildAction,
    onChange(open) {
      if (open && !library.isOpen) library.open({ focus: false });
      clearBuild();
      if (!open && library.loaded && isSwapMode()) library.render();
    },
  });
  function clearBuild() {
    connectionValid = false;
    clearBuildPreview();
    requests.cancel('build-left', 'build-right', 'bond-check', 'sequence-edit', 'swap-options');
    swapState = null;
    swapMapping.hidden = true;
    swapSites.innerHTML = '';
    buildLeft = null; buildRight = null; buildLeftRIdx = null; buildRightRIdx = null;
    buildReaction = '';
    insertBetweenActive = false;
    buildInsertRow.style.display = 'none';
    buildInsertBtn.textContent = '⊕ Insert Between';
    buildLeftAbbr.textContent = '—';
    document.getElementById('build-left-label').textContent = isSwapMode() ? '1 · Residue to replace' : 'Current residue';
    buildLeftSvg.innerHTML = '<div class="box-placeholder">Select a residue in the drawing or its tile</div>';
    buildLeftRg.innerHTML = '';
    buildLeftSite.hidden = true;
    buildRightAbbr.textContent = '—';
    document.getElementById('build-right-label').textContent = isSwapMode() ? '2 · Eligible replacement' : 'New monomer';
    document.getElementById('build-title').textContent = 'Build';
    document.getElementById('build-swap-tutorial').hidden = !isSwapMode();
    buildPanel.classList.toggle('swapping', isSwapMode());
    buildConnect.textContent = isSwapMode() ? 'Apply swap' : 'Connect';
    buildPreviewButton.textContent = isSwapMode() ? 'Preview swap' : 'Preview';
    const changeReplacement = document.getElementById('build-right-change');
    changeReplacement.textContent = isSwapMode() ? 'Change' : 'Clear';
    changeReplacement.title = isSwapMode() ? 'Choose another replacement' : 'Clear selection';
    changeReplacement.setAttribute('aria-label', isSwapMode() ? 'Change replacement monomer' : 'Clear new monomer selection');
    buildRightSvg.innerHTML = `<div class="box-placeholder">${isSwapMode() ? 'Select a residue to see eligible replacements in Library' : 'Choose Use in the library'}</div>`;
    buildRightRg.innerHTML = '';
    buildRightSite.hidden = true;
    updateControls();
    buildStatus.textContent = isSwapMode() ? 'Select the residue to replace' : 'Choose a residue and a monomer to connect';
    buildStatus.className = 'build-status';
    buildHint.textContent = getDocument().notation !== 'cabiln' ? 'Convert this input to CABILN before building'
      : isSwapMode() ? 'Select the residue to replace; every existing connection will be kept'
      : getDocument().text.trim()
      ? 'Select a residue in the drawing or its tile, then choose Use in the library'
      : 'Choose Use beside a monomer to start a peptide';
    residueView.select();
    const filtered = library.resetFilter();
    if (library.loaded && (filtered || (panel.isOpen && isSwapMode()))) library.render();
  }

  function isSwapMode() { return buildAction.value === 'swap'; }

  buildAction.addEventListener('change', () => {
    const selected = residueView.residue(buildLeftRIdx);
    clearBuild();
    if (library.loaded) library.render();
    if (selected && !isStale()) selectBuildResidue(selected);
  });

  async function loadSwapOptions() {
    const request = requests.start('swap-options', updateControls);
    const source = getDocument().text.trim();
    const index = buildLeftRIdx;
    buildStatus.textContent = 'Finding replacements for all connected sites…';
    try {
      if (requests.has('library', 'reactions') &&
          !await request.waitFor('library', 'reactions')) return;
      const data = await readResponse(await postCalculation('/replacement_options', {
        cabiln: source, residue_idx: index,
      }, request.signal));
      if (!request.current() || source !== getDocument().text.trim() || index !== buildLeftRIdx) return;
      if (data.error) throw new Error(data.error);
      if (data.source_echo !== source || data.residue_idx !== index || !data.context) {
        throw new Error('The replacement list is out of date. Select the residue again.');
      }
      const quality = new Map(library.monomers.map(m => [m.abbr, m.quality]));
      const candidates = data.candidates.map(m => ({ ...m, quality: quality.get(m.abbr),
        searchText: [m.abbr, m.name, m.type, m.chem_types].join(' ').toLowerCase() }));
      swapState = { ...data, candidates, candidate: null, mapping: {}, preview: null };
      buildStatus.textContent = `${candidates.length} eligible replacements in Library`;
      buildHint.textContent = `${data.requirements.length} connections to keep · Choose a replacement, review the mapping, then preview`;
      library.render();
    } catch (error) {
      if (request.current()) {
        buildStatus.textContent = error.message || 'Could not find replacements. Select the residue to retry.';
        buildStatus.className = 'build-status invalid';
      }
    } finally { request.finish(); }
  }

  function chooseSwapMonomer(abbr) {
    const candidate = swapState?.candidates.find(item => item.abbr === abbr);
    if (!candidate) return;
    swapState.candidate = candidate;
    swapState.mapping = { ...candidate.mapping };
    swapState.preview = null;
    renderSwapMapping();
    loadBuildRight(abbr);
  }

  function renderSwapMapping() {
    swapSites.innerHTML = '';
    swapMapping.hidden = !swapState?.candidate;
    if (swapMapping.hidden) return;
    const candidate = swapState.candidate;
    for (const need of swapState.requirements) {
      const label = document.createElement('label');
      label.className = 'swap-site';
      const text = document.createElement('span');
      text.textContent = `${buildLeft.abbr} R${need.slot} →`;
      const select = document.createElement('select');
      select.className = 'hbtn';
      select.dataset.slot = need.slot;
      select.setAttribute('aria-label', `Replacement site for R${need.slot}`);
      for (const slot of candidate.choices[need.slot]) {
        const option = document.createElement('option');
        option.value = slot;
        const site = candidate.sites.find(item => item.slot === slot);
        option.textContent = `R${slot} ${(site?.chem_type || '').replaceAll('_', ' ')}`;
        select.appendChild(option);
      }
      select.value = swapState.mapping[need.slot];
      select.addEventListener('change', () => {
        requests.cancel('sequence-edit');
        swapState.mapping[need.slot] = Number(select.value);
        checkBuildValidity();
        if (buildRight) renderRgroupButtons(buildRightRg, buildRight, 'right');
      });
      const partner = document.createElement('small');
      partner.textContent = need.internal ? `Keeps the internal link to the site mapped from R${need.partner_slot}`
        : `Keeps ${need.partner_abbr} (residue ${need.partner_idx + 1}) R${need.partner_slot}`;
      label.appendChild(text);
      label.appendChild(select);
      label.appendChild(partner);
      swapSites.appendChild(label);
    }
    if (!swapState.requirements.length) swapSites.textContent = 'This monomer has no existing connections.';
  }

  function swapRequest() {
    if (!swapState?.candidate || !buildLeft || !buildRight || isStale() ||
        swapState.source_echo !== getDocument().text.trim() || swapState.residue_idx !== buildLeftRIdx ||
        swapState.candidate.abbr !== buildRight.abbr) return null;
    const slots = Object.values(swapState.mapping);
    if (slots.length !== swapState.requirements.length || new Set(slots).size !== slots.length) return null;
    return { cabiln: swapState.source_echo, residue_idx: swapState.residue_idx,
      new_abbr: swapState.candidate.abbr, slot_map: { ...swapState.mapping }, context: swapState.context };
  }
  function checkAdjacentBackbone() {
    return residueView.insertionAnchor(buildLeftRIdx, buildRightRIdx) !== null;
  }

  function updateInsertBetweenUI() {
    if (panel.isOpen && !isSwapMode() && checkAdjacentBackbone()) {
      const la = (residueView.residue(buildLeftRIdx) || {}).abbr || '?';
      const ra = (residueView.residue(buildRightRIdx) || {}).abbr || '?';
      buildInsertInfo.textContent = `${la} and ${ra} are adjacent on the backbone`;
      buildInsertRow.style.display = 'flex';
    } else {
      buildInsertRow.style.display = 'none';
      if (insertBetweenActive) {
        insertBetweenActive = false;
        buildInsertBtn.textContent = '⊕ Insert Between';
        if (library.loaded) library.render();
      }
    }
  }

  buildInsertBtn.addEventListener('click', () => {
    clearBuildPreview();
    requests.cancel('sequence-edit');
    if (insertBetweenActive) {
      insertBetweenActive = false;
      updateControls();
      buildInsertBtn.textContent = '⊕ Insert Between';
      buildHint.textContent = 'Select a residue, then choose Use beside a library monomer';
      if (library.loaded) library.render();
      return;
    }
    if (!checkAdjacentBackbone()) return;
    insertBetweenActive = true;
    updateControls();
    buildInsertBtn.textContent = '✕ Cancel';
    buildHint.textContent = 'Click a backbone monomer in the library to insert between the selected residues';
    if (!library.isOpen) library.open();
    if (library.loaded) library.render();
  });

  async function doInsertBetween(abbr) {
    clearBuildPreview();
    const after_idx = residueView.insertionAnchor(buildLeftRIdx, buildRightRIdx);
    if (after_idx === null) return;
    const val = getDocument().text.trim();
    const request = requests.start('sequence-edit', updateControls);
    updateControls();
    buildHint.textContent = 'Inserting…';
    try {
      const res = await postCalculation('/insert_backbone', {
        cabiln: val, after_idx, new_abbr: abbr,
      }, request.signal);
      const data = await readResponse(res);
      if (!request.current() || getDocument().text.trim() !== val) return;
      if (data.error) {
        buildHint.textContent = `Could not insert ${abbr}: ${data.error}. Choose another backbone block in Library. Your sequence is unchanged.`;
        return;
      }
      commitDocument(data.result, 'cabiln');
      buildHint.textContent = `${abbr} inserted — select a residue to continue building`;
      if (library.loaded) library.render();
    } catch (e) {
      if (!request.current()) return;
      showRetry(buildHint, `Could not insert ${abbr} between the selected blocks. Your sequence is unchanged.`,
        () => doInsertBetween(abbr));
    } finally {
      request.finish();
    }
  }
  document.getElementById('build-left-change').addEventListener('click', clearBuild);
  document.getElementById('build-right-change').addEventListener('click', () => {
    const residue = isSwapMode() && residueView.residue(buildLeftRIdx);
    if (!residue) { clearBuild(); return; }
    selectBuildResidue(residue);
    if (!library.isOpen) library.open();
  });

  async function loadBuildLeft(abbr, rIdx) {
    clearBuildPreview();
    requests.cancel('bond-check', 'sequence-edit', 'swap-options');
    if (isSwapMode()) {
      requests.cancel('build-right');
      swapState = null;
      swapMapping.hidden = true;
      buildRight = null;
      buildRightRIdx = null;
      buildRightAbbr.textContent = '—';
      buildRightSvg.innerHTML = '<div class="box-placeholder">Choose a replacement in the library</div>';
      buildRightRg.innerHTML = '';
      buildRightSite.hidden = true;
      library.render();
    }
    const request = requests.start('build-left', updateControls);
    const sequence = getDocument().text.trim();
    buildLeftRIdx = rIdx;
    buildLeftAbbr.textContent = abbr;
    if (isSwapMode()) document.getElementById('build-left-label').textContent = `1 · Replace residue ${rIdx + 1}`;
    showLoading(buildLeftSvg, 'Loading attachment sites…');
    buildLeftRg.innerHTML = '';
    buildLeftSite.hidden = true;
    buildLeft = null;
    connectionValid = false;
    updateControls();
    if (library.filterActive && library.loaded) library.render();
    buildStatus.textContent = '';

    try {
      const res = await fetchCalculation(`/monomer_rgroups?abbr=${encodeURIComponent(abbr)}&residue_idx=${rIdx}&cabiln=${encodeURIComponent(sequence)}`, { signal: request.signal });
      const data = await readResponse(res);
      if (!request.current() || getDocument().text.trim() !== sequence) return;
      if (data.error) {
        buildLeftSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
        return;
      }
      buildLeftSvg.innerHTML = data.svg || '';
      buildLeft = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
      renderRgroupButtons(buildLeftRg, buildLeft, 'left');
      buildHint.textContent = buildRight ? 'Choose an attachment site on each side' : 'Choose Use beside a library monomer';
      library.updateFilterStatus();
      if (library.filterActive && library.loaded) library.render();
      updateInsertBetweenUI();
      if (isSwapMode()) loadSwapOptions();
      else checkBuildValidity();
    } catch (e) {
      if (!request.current()) return;
      showRetry(buildLeftSvg, 'Could not load attachment sites.', () => loadBuildLeft(abbr, rIdx));
    } finally {
      request.finish();
    }
  }

  async function loadBuildRight(abbr, rIdx) {
    clearBuildPreview();
    requests.cancel('bond-check', 'sequence-edit');
    const request = requests.start('build-right', updateControls);
    const sequence = getDocument().text.trim();
    buildRightRIdx = rIdx !== undefined ? rIdx : null;
    buildRightAbbr.textContent = abbr;
    document.getElementById('build-right-label').textContent = isSwapMode() ? '2 · Eligible replacement'
      : rIdx !== undefined ? 'Current residue' : 'New monomer';
    showLoading(buildRightSvg, 'Loading attachment sites…');
    buildRightRg.innerHTML = '';
    buildRightSite.hidden = true;
    buildRight = null;
    connectionValid = false;
    updateControls();
    buildStatus.textContent = '';

    const family = !isSwapMode() && rIdx === undefined && library.monomers.find(m => m.abbr === abbr && m.degenerate);
    if (family) {
      buildRightSvg.innerHTML = '<div class="box-placeholder">Choose the attachment form</div>';
      for (const [label, symbol] of [['N-terminal', family.nterm_abbr], ['C-terminal', family.cterm_abbr]]) {
        if (!symbol) continue;
        const button = document.createElement('button');
        button.className = 'rgroup-btn';
        button.textContent = `${label}: ${symbol}`;
        button.addEventListener('click', () => loadBuildRight(symbol));
        buildRightRg.appendChild(button);
      }
      buildHint.textContent = 'Choose a monomer form, then its attachment site';
      request.finish();
      return;
    }

    try {
      let url = `/monomer_rgroups?abbr=${encodeURIComponent(abbr)}`;
      if (rIdx !== undefined) {
        url += `&residue_idx=${rIdx}&cabiln=${encodeURIComponent(sequence)}`;
      }
      const res = await fetchCalculation(url, { signal: request.signal });
      const data = await readResponse(res);
      if (!request.current() || getDocument().text.trim() !== sequence) return;
      if (data.error) {
        buildRightSvg.innerHTML = `<div class="box-placeholder">${escHtml(data.error)}</div>`;
        return;
      }
      buildRightSvg.innerHTML = data.svg || '';
      buildRight = { abbr, rgroups: data.rgroups || [], selectedSlot: null };
      renderRgroupButtons(buildRightRg, buildRight, 'right');
      buildHint.textContent = isSwapMode() ? 'Review each site mapping, then preview the replacement'
        : buildLeft ? 'Choose an attachment site on each side' : 'Select a residue in the sequence';
      updateInsertBetweenUI();
      checkBuildValidity();
    } catch (e) {
      if (!request.current()) return;
      showRetry(buildRightSvg, 'Could not load attachment sites.', () => loadBuildRight(abbr, rIdx));
    } finally {
      request.finish();
    }
  }

  function renderRgroupButtons(container, state, side) {
    container.innerHTML = '';
    state.rgroups.forEach(rg => {
      const swapping = isSwapMode();
      const btn = document.createElement(swapping ? 'span' : 'button');
      btn.className = swapping ? 'site-tag' : 'rgroup-btn';
      if (rg.used) btn.classList.add('used');
      if (!swapping) {
        btn.disabled = !!rg.used;
        btn.setAttribute('aria-pressed', String(state.selectedSlot === rg.slot));
      }
      if (state.selectedSlot === rg.slot) btn.classList.add('selected');
      const mapped = swapping && (rg.used || (side === 'right' && Object.values(swapState?.mapping || {}).includes(rg.slot)));
      btn.classList.toggle('mapped', !!mapped);
      btn.textContent = swapping
        ? `R${rg.slot} · ${mapped ? (side === 'left' ? 'linked' : 'mapped') : 'free'}`
        : `R${rg.slot} ${(rg.chem_type || '').replaceAll('_', ' ')}${rg.used ? ' · used' : ''}`;
      btn.title = `R${rg.slot}: ${rg.chem_type || 'unknown'} · Free-site group: ${rg.leaving || '[H] (implicit)'}${rg.used ? ' — already connected' : ''}`;
      if (!isSwapMode() && !rg.used) {
        btn.addEventListener('click', () => selectRgroup(side, rg.slot));
      }
      container.appendChild(btn);
    });
    const drawing = side === 'left' ? buildLeftSvg : buildRightSvg;
    drawing.querySelectorAll('.site-selected').forEach(path => path.classList.remove('site-selected'));
    const selected = state.rgroups.find(rg => rg.slot === state.selectedSlot);
    const detail = side === 'left' ? buildLeftSite : buildRightSite;
    detail.hidden = !selected;
    if (selected) {
      detail.textContent = `R${selected.slot} · Free-site group: ${selected.leaving || '[H] (implicit)'}`;
      detail.title = 'Group present at this site when it is unconnected. The product preview shows the selected reaction.';
    }
    if (Number.isInteger(selected?.atom_idx)) {
      drawing.querySelectorAll(`.atom-${selected.atom_idx}`).forEach(path => path.classList.add('site-selected'));
    }
    if (isSwapMode()) {
      const occupied = side === 'left' ? state.rgroups.filter(rg => rg.used)
        : state.rgroups.filter(rg => Object.values(swapState?.mapping || {}).includes(rg.slot));
      for (const site of occupied) {
        drawing.querySelectorAll(`.atom-${site.atom_idx}`).forEach(path => path.classList.add('site-selected'));
      }
    }
  }

  function selectRgroup(side, slot) {
    requests.cancel('sequence-edit');
    if (side === 'left' && buildLeft) {
      buildLeft.selectedSlot = buildLeft.selectedSlot === slot ? null : slot;
      renderRgroupButtons(buildLeftRg, buildLeft, 'left');
      if (library.filterActive && library.loaded) library.render();
    } else if (side === 'right' && buildRight) {
      buildRight.selectedSlot = buildRight.selectedSlot === slot ? null : slot;
      renderRgroupButtons(buildRightRg, buildRight, 'right');
    }
    checkBuildValidity();
  }

  function suggestedPair() {
    if (!panel.isOpen || isSwapMode() || isStale() || insertBetweenActive || isPending() ||
        !buildLeft || !buildRight || (buildLeft.selectedSlot && buildRight.selectedSlot)) return null;
    const available = state => state.rgroups.filter(site => !state.selectedSlot || site.slot === state.selectedSlot);
    const pairs = library.compatiblePairs(available(buildLeft), available(buildRight));
    return pairs.length === 1 ? pairs[0] : null;
  }

  buildSuggestion.addEventListener('click', () => {
    const pair = suggestedPair();
    if (!pair) return;
    [buildLeft.selectedSlot, buildRight.selectedSlot] = pair;
    renderRgroupButtons(buildLeftRg, buildLeft, 'left');
    renderRgroupButtons(buildRightRg, buildRight, 'right');
    if (library.filterActive) library.render();
    checkBuildValidity();
  });

  async function checkBuildValidity() {
    clearBuildPreview();
    requests.cancel('bond-check');
    buildReaction = '';
    connectionValid = false;
    updateControls();
    if (isSwapMode()) {
      const valid = !!swapRequest();
      buildStatus.textContent = valid ? 'Preview to validate the complete replacement product'
        : swapState?.candidate && buildRight ? 'Choose a different replacement site for each connection'
        : 'Choose a replacement from the filtered library';
      buildStatus.className = 'build-status';
      updateControls();
      return;
    }
    if (!buildLeft?.selectedSlot || !buildRight?.selectedSlot) {
      buildStatus.textContent = !buildLeft ? 'Select a residue in the sequence'
        : !buildRight ? 'Choose a monomer from the library'
        : !buildLeft.selectedSlot ? 'Choose a site on the current residue'
        : 'Choose a site on the other monomer';
      buildStatus.className = 'build-status';
      return;
    }

    const leftRg = buildLeft.rgroups.find(r => r.slot === buildLeft.selectedSlot);
    const rightRg = buildRight.rgroups.find(r => r.slot === buildRight.selectedSlot);
    if (!leftRg || !rightRg) return;

    const request = requests.start('bond-check', updateControls);

    buildStatus.textContent = 'Checking bond...';
    buildStatus.className = 'build-status';

    try {
      const res = await postCalculation('/validate_bond', {
        chem_type_a: leftRg.chem_type,
        chem_type_b: rightRg.chem_type,
        abbr_a: buildLeft.abbr,
        slot_a: buildLeft.selectedSlot,
        abbr_b: buildRight.abbr,
        slot_b: buildRight.selectedSlot
      }, request.signal);
      const data = await readResponse(res);
      if (!request.current()) return;
      if (!res.ok) throw new Error(data.error || 'Connection check unavailable');
      if (data.valid) {
        buildReaction = data.reaction ? `Reaction: ${data.reaction.replaceAll('_', ' ')}` : 'Bond formation';
        buildStatus.textContent = `Valid: ${data.reaction || 'bond'} (R${buildLeft.selectedSlot}↔R${buildRight.selectedSlot})`;
        buildStatus.className = 'build-status valid';
        connectionValid = true;
      } else {
        buildStatus.textContent = `${buildLeft.abbr} R${buildLeft.selectedSlot} and ${buildRight.abbr} R${buildRight.selectedSlot}: ${data.error || data.reason || 'No supported reaction'}. Choose another free site or monomer.`;
        buildStatus.className = 'build-status invalid';
        connectionValid = false;
      }
    } catch (e) {
      if (!request.current()) return;
      showRetry(buildStatus, `Could not check ${buildLeft.abbr} R${buildLeft.selectedSlot} with ${buildRight.abbr} R${buildRight.selectedSlot}. Your sequence is unchanged.`, checkBuildValidity);
      buildStatus.className = 'build-status invalid';
      connectionValid = false;
    } finally {
      request.finish();
    }
  }

  function isPending() {
    return requests.has('build-left', 'build-right', 'bond-check', 'swap-options', 'build-preview', 'sequence-edit');
  }

  function updateControls() {
    const ready = panel.isOpen && !isStale() && !insertBetweenActive &&
      !isPending() &&
      (isSwapMode() ? !!swapRequest() : connectionValid && !!buildConnection());
    buildPreviewButton.disabled = !ready;
    buildConnect.disabled = !ready || (isSwapMode() && !swapState?.preview);
    const pair = suggestedPair();
    buildSuggestion.hidden = !pair;
    if (pair) {
      buildSuggestion.textContent = `Use ${buildLeft.abbr} R${pair[0]} → ${buildRight.abbr} R${pair[1]}`;
      buildSuggestion.title = 'The only compatible pair for these choices. You can change either site before connecting.';
    }
  }

  function buildConnection() {
    if (isStale() || !buildLeft?.selectedSlot || !buildRight?.selectedSlot) return null;
    return {
      cabiln: getDocument().text.trim(),
      host_residue_idx: buildLeftRIdx ?? 0,
      new_abbr: buildRight.abbr,
      r_host: buildLeft.selectedSlot,
      r_new: buildRight.selectedSlot,
      target_residue_idx: (buildRightRIdx !== null && buildRightRIdx !== buildLeftRIdx)
        ? buildRightRIdx : -1,
    };
  }

  function discardPreview() {
    if (swapState) swapState.preview = null;
    requests.cancel('build-preview');
    updateControls();
    buildPreviewInner.innerHTML = '';
    buildPreviewStatus.textContent = '';
    buildPreviewReaction.textContent = '';
    buildPreviewChanges.textContent = '';
    buildPreviewChanges.hidden = true;
    buildPreviewFocus.disabled = true;
    buildPreviewSource.textContent = '';
    buildPreviewNotation.hidden = true;
  }

  function clearBuildPreview() { previewPanel.close({ focus: false }); }

  function markPreviewChanges(proposal, drawing, swapping) {
    const changes = proposal.change;
    const svg = buildPreviewInner.querySelector('svg');
    if (!changes || !svg) return;
    const atoms = new Set(changes.residues.flatMap(index => drawing.residue_map[index] || []));
    for (const element of svg.querySelectorAll('[class*="atom-"]')) {
      const indices = [...element.classList].filter(cls => /^atom-\d+$/.test(cls))
        .map(cls => Number(cls.slice(5)));
      if (indices.length && indices.every(index => atoms.has(index))) element.classList.add('edit-residue');
    }
    const key = ends => ends.map(end => `${end.residue}:${end.slot}`).sort().join('|');
    const affected = new Set(changes.connections.map(key));
    for (const connection of drawing.connection_bonds || []) {
      if (!affected.has(key(connection.endpoints))) continue;
      for (const index of connection.bonds) {
        svg.querySelectorAll(`.bond-${index}`).forEach(element => element.classList.add('edit-connection'));
      }
    }
    const action = swapping ? `Replace ${buildLeft.abbr} with ${buildRight.abbr}`
      : changes.residues.length ? `Add ${buildRight.abbr}` : 'Join the selected blocks';
    const block = changes.residues.length ? 'Blue glow: changed block. ' : '';
    const connections = changes.connections.length ? `Orange glow: ${swapping ? 'kept' : 'new'} connections.` : 'No existing connections.';
    buildPreviewChanges.textContent = `${action}. ${block}${connections}`;
    buildPreviewChanges.hidden = false;
    buildPreviewFocus.disabled = !svg.querySelector('.edit-residue, .edit-connection');
  }

  buildPreviewButton.addEventListener('click', async () => {
    const swapping = isSwapMode();
    const connection = swapping ? swapRequest() : buildConnection();
    if (!connection || buildPreviewButton.disabled) return;
    clearBuildPreview();
    const request = requests.start('build-preview', () => {
      buildPreviewPanel.setAttribute('aria-busy', 'false');
      updateControls();
    });
    const left = swapping ? `${buildLeft.abbr} (residue ${connection.residue_idx + 1})`
      : `${buildLeft.abbr} (residue ${connection.host_residue_idx + 1}) R${connection.r_host}`;
    const right = swapping ? buildRight.abbr
      : `${buildRight.abbr} (${connection.target_residue_idx < 0 ? 'new' : `residue ${connection.target_residue_idx + 1}`}) R${connection.r_new}`;
    const pair = `${left} ↔ ${right}`;
    updateControls();
    previewPanel.open({ focus: false });
    buildPreviewPanel.setAttribute('aria-busy', 'true');
    buildPreviewStatus.textContent = `${pair} · Preparing preview…`;
    buildPreviewStatus.className = '';
    buildPreviewReaction.textContent = swapping ? Object.entries(connection.slot_map)
      .map(([old, next]) => `R${old} → R${next}`).join(' · ') || 'No existing connections' : buildReaction;
    showLoading(buildPreviewInner, 'Preparing product preview…');
    buildPreviewViewport.reset();
    buildPreviewPanel.scrollIntoView({ block: 'nearest' });
    try {
      // Drawing and verification share one worker. Do not spend preview retries
      // competing with work that is already running for this page.
      if (requests.has('main-render', 'reference-render', 'verify')) {
        buildPreviewStatus.textContent = `${pair} · Waiting for the current drawing or verification…`;
        if (!await request.waitFor('main-render', 'reference-render', 'verify')) return;
      }
      if (!request.current()) return;
      buildPreviewStatus.textContent = `${pair} · Preparing preview…`;
      // These endpoints calculate a candidate; only Connect commits the document.
      const proposal = await readResponse(await postCalculation(
        swapping ? '/replace_monomer' : '/insert_bond', connection, request.signal
      ));
      if (!request.current()) return;
      if (proposal.error) throw new Error(proposal.error);
      const drawing = await readResponse(await postCalculation('/render', {
        cabiln: proposal.result, width: 1000, height: 320,
      }, request.signal));
      if (!request.current()) return;
      if (drawing.error) throw new Error(drawing.error);
      if (swapping && (!CabilnProject.sameContext(connection.context, proposal.context) ||
          !CabilnProject.sameContext(connection.context, drawing.context))) {
        throw new Error('The library changed. Select the residue again to refresh replacements');
      }
      const changedLibrary = getDocument().context?.library_binding && drawing.context?.library_binding
        && !CabilnProject.sameContext(getDocument().context, drawing.context);
      const notice = changedLibrary ? 'Library changed since the current drawing; preview uses current definitions.' : '';
      buildPreviewInner.innerHTML = drawing.svg;
      markPreviewChanges(proposal, drawing, swapping);
      showStructureInfo(buildPreviewStatus, drawing, { lead: pair, notices: [notice] });
      buildPreviewSource.textContent = proposal.result;
      buildPreviewNotation.hidden = false;
      if (swapping) {
        swapState.preview = { result: proposal.result, request: connection };
        buildStatus.textContent = 'Product validated. Apply swap keeps every mapped connection.';
        buildStatus.className = 'build-status valid';
        if (Object.entries(connection.slot_map).some(([old, next]) => Number(old) !== next)) {
          const warning = document.createElement('p');
          warning.className = 'structure-warning';
          warning.textContent = 'Site renumbering can reformat the notation; review Proposed CABILN.';
          buildPreviewStatus.appendChild(warning);
        }
      }
    } catch (error) {
      if (!request.current()) return;
      buildPreviewInner.innerHTML = '';
      showRetry(buildPreviewStatus, `Preview unavailable: ${error.message}. Your sequence is unchanged.`,
        () => buildPreviewButton.click());
      buildPreviewStatus.className = 'error';
    } finally {
      request.finish();
    }
  });

  buildConnect.addEventListener('click', async () => {
    const swapping = isSwapMode();
    const reviewed = swapping ? swapState?.preview : null;
    const connection = swapping ? reviewed?.request : buildConnection();
    if (!connection || buildConnect.disabled) return;
    clearBuildPreview();
    const request = requests.start('sequence-edit', updateControls);

    updateControls();
    buildStatus.textContent = swapping ? 'Checking and applying swap…' : 'Inserting...';
    buildStatus.className = 'build-status';

    try {
      const res = await postCalculation(
        swapping ? '/replace_monomer' : '/insert_bond', connection, request.signal
      );
      const data = await readResponse(res);
      if (!request.current() || getDocument().text.trim() !== connection.cabiln) return;
      if (data.error) {
        buildStatus.textContent = `${data.error}. Your sequence is unchanged. ${swapping
          ? 'Review the replacement and mapping, then preview again.'
          : 'Check the selected blocks and connection points, then try again.'}`;
        buildStatus.className = 'build-status invalid';
        return;
      }
      if (swapping && (data.result !== reviewed.result ||
          !CabilnProject.sameContext(connection.context, data.context))) {
        buildStatus.textContent = 'The replacement changed. Preview it again before applying.';
        buildStatus.className = 'build-status invalid';
        return;
      }
      commitDocument(data.result, 'cabiln');
      onApply?.({ action: swapping ? 'swap' : 'connect', request: connection, result: data.result });
      buildHint.textContent = swapping ? 'Monomer replaced — Undo restores the original peptide'
        : 'Connection added — select a residue to continue building';
    } catch (e) {
      if (!request.current()) return;
      if (swapping) {
        showRetry(buildStatus, 'Could not apply the swap. Your sequence is unchanged. Preview it again to retry.',
          () => buildPreviewButton.click());
      } else {
        showRetry(buildStatus, `Could not connect ${buildLeft.abbr} R${connection.r_host} with ${buildRight.abbr} R${connection.r_new}. Your sequence is unchanged.`,
          () => buildConnect.click());
      }
      buildStatus.className = 'build-status invalid';
    } finally {
      request.finish();
    }
  });

  return {
    refreshChoices: updateControls,
    open: () => panel.open(),
    selectResidue: selectBuildResidue,
    useMonomer: useLibraryMonomer,
    clear: clearBuild,
    get selection() { return [buildLeftRIdx, buildRightRIdx]; },
    get session() {
      return { open: panel.isOpen, action: buildAction.value, busy: isPending(),
        left: buildLeft && { ...buildLeft, index: buildLeftRIdx },
        right: buildRight && { ...buildRight, index: buildRightRIdx },
        ready: !buildPreviewButton.disabled,
        preview: !buildPreviewPanel.hidden ? buildPreviewSource.textContent : null };
    },
    get filters() {
      return { building: panel.isOpen, notation: getDocument().notation, left: buildLeft, insertBetween: insertBetweenActive,
        replacements: panel.isOpen && isSwapMode() ? (swapState?.candidates || []) : null,
        awaitingSelection: !swapState };
    },
  };
}
