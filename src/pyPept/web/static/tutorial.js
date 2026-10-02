function startTutorial({ lesson = 'connect', loadPractice, getExample, getBuildState, isReady, openBuild, openLibrary }) {
  const find = selector => document.querySelector(selector);
  const lessons = {
    connect: { source: 'A-G', name: 'A–G', host: 1, residue: 'glycine', steps: [
      { id: 'load', title: 'Start with two residues', target: '#cabiln-input',
        body: 'Load A–G: alanine joined to glycine. This practice tab leaves your saved draft untouched. Practice edits are not saved automatically.' },
      { id: 'explore', title: 'Find a residue in the drawing',
        body: 'Hover over glycine in the drawing or focus the G tile to light up its atoms. On a touch screen, tap the tile.' },
      { id: 'select', title: 'Select glycine in Build',
        body: 'Open Build, then click or tap glycine in the main drawing or its G tile. Aim at an atom label or bond; drag to pan. Its numbered attachment sites appear in Build.' },
      { id: 'choose', title: 'Choose a monomer to add', target: '#lib-list',
        body: 'Browse or search Library and choose Use beside an amino acid, such as A (alanine) or L (leucine). This exercise uses a monomer with an unused R1 site. Filter can narrow the choices.' },
      { id: 'sites', title: 'Choose the connection', target: '#build-left-rgroups',
        body: 'Choose R2 on glycine, then R1 on your chosen monomer. Used sites cannot be selected. Connect becomes available when the reaction is valid.' },
      { id: 'preview', title: 'Preview the product', target: '#build-preview-button',
        body: 'Choose Preview. Check the proposed molecule, reaction and sequence before applying the connection. The result uses the monomer you chose.' },
      { id: 'apply', title: 'Apply the connection', target: '#build-connect',
        body: 'Choose Connect. The main drawing and residue tiles update to the product.' },
      { id: 'undo', title: 'Undo the edit', target: '#btn-undo',
        body: 'Choose Undo to return to A–G. Redo can restore the connection. Both work with notation conversions and other sequence edits too.' },
      { id: 'finish', title: 'Continue exploring', target: '#btn-project-save',
        body: 'Save project keeps an editable copy on your device. PNG and MOL export the structure. Choose the Swap in Retatrutide lesson above to practise replacing a branched residue, or close this tab to return to your original work.' },
    ] },
    swap: { source: null, name: 'Retatrutide', host: 16, residue: 'K17', steps: [
      { id: 'load', title: 'Try a swap in Retatrutide', target: '#cabiln-input',
        body: 'Load Retatrutide from Examples. You will replace its lipid-bearing lysine while keeping the branch and both backbone connections. This is a practice edit; your saved work stays untouched.' },
      { id: 'select', title: 'Select the branched lysine',
        body: 'Open Build and choose Swap monomer. Select K17: the second of the two adjacent K tiles, immediately before the branch. You can also select its atoms in the drawing. Show control points to the tile.' },
      { id: 'choose', title: 'Choose an eligible replacement', target: '#lib-list',
        body: 'Library now shows replacements that can keep every connection. Search for Orn (L-ornithine) and choose Select. You can try another eligible monomer too; selecting one does not change the peptide.' },
      { id: 'sites', title: 'Review all three connections', target: '#build-swap-mapping',
        body: 'Check the two backbone connections and the AEEA branch. Each current site maps to its own compatible replacement site. The numbers may differ; the menus let you change the mapping.' },
      { id: 'preview', title: 'Preview the replacement', target: '#build-preview-button',
        body: 'Choose Preview swap. Inspect the complete peptide and the site mapping. Expand Proposed CABILN to see the new notation. Apply stays disabled until this product has been validated.' },
      { id: 'apply', title: 'Apply the swap', target: '#build-connect',
        body: 'Choose Apply swap. The selected K becomes your replacement, and the lipid branch remains attached. Changing a mapping or monomer requires a new preview.' },
      { id: 'undo', title: 'Restore Retatrutide', target: '#btn-undo',
        body: 'Choose Undo to restore the original Retatrutide. Redo restores your swap. Back in this guide revisits instructions without changing the molecule.' },
      { id: 'finish', title: 'Try another replacement', target: '#btn-project-save',
        body: 'You can swap other residues, including ones with branches or crosslinks. Eligibility comes from the monomer library and the occupied sites. Save project keeps your edit; close this tab to return to your original work.' },
    ] },
  };
  const panel = find('#tutorial-panel'), button = find('#btn-tutorial');
  const next = find('#tutorial-next'), back = find('#tutorial-back'), load = find('#tutorial-load');
  const status = find('#tutorial-status'), show = find('#tutorial-show'), chooser = find('#tutorial-lesson');
  let index = 0, target = null, frame = null, explored = false, applied = null, reveal = null;
  let loading = false, loadError = '', generation = 0;
  const current = () => lessons[lesson];
  const hostTile = () => `.res-chip[data-residue="${current().host}"]`;
  const observer = new MutationObserver(() => {
    if (frame === null) frame = requestAnimationFrame(() => { frame = null; check(); });
  });
  const disclosure = createPanel({ panel, button, closeButton: find('#tutorial-close'),
    focus: find('#tutorial-title'),
    onChange(open) {
      document.documentElement.classList.toggle('tutorial-open', open);
      if (open) {
        for (const selector of ['#input-bar', '#build-panel', '#main']) {
          observer.observe(find(selector), { subtree: true, childList: true, attributes: true, characterData: true });
        }
        showStep();
      } else {
        observer.disconnect();
        cancelAnimationFrame(frame);
        frame = null;
        target?.classList.remove('tutorial-target');
        target = null;
      }
    },
  });

  function readState() {
    const base = !!current().source && isReady(current().source), build = getBuildState();
    const left = base && build.open && build.action === lesson && build.left?.index === current().host;
    const candidate = left && build.right?.index === null && (lesson === 'swap' ||
      build.right.rgroups.some(site => site.slot === 1 && !site.used));
    const connection = candidate && build.ready && (lesson === 'swap' ||
      (build.left.selectedSlot === 2 && build.right.selectedSlot === 1));
    const product = !!applied && isReady(applied.result);
    explored ||= base && !!find(`${hostTile()}.hover`) && !!find('#render-inner .res-hl');
    return { base, build, left, candidate, connection, preview: connection && !!build.preview,
      product, review: !!applied && (base || product) };
  }

  function recovery(state, id) {
    if (id === 'load' || id === 'finish' || state.review) return null;
    if (!state.base) return { target: '#tutorial-restart', message: `Wait for the drawing, or use Restart to load ${current().name} again.` };
    if (id === 'explore') {
      if (!explored && find('#btn-hl').getAttribute('aria-pressed') === 'false') {
        return { target: '#btn-hl', message: 'Enable Highlight to explore the atoms.', label: 'Enable Highlight', action: () => find('#btn-hl').click() };
      }
      return null;
    }
    if (id === 'undo') return { target: '#build-connect', message: 'Apply the practice edit before trying Undo. Back revisits the earlier steps.' };
    if (!state.build.open) return { target: '#btn-build', message: 'Build is closed. Reopen it to continue.', label: 'Open Build', action: openBuild };
    if (state.build.action !== lesson) return { target: '#build-action', message: `Choose ${lesson === 'swap' ? 'Swap monomer' : 'Connect'} in the Build action menu for this exercise.` };
    if (!state.left) return { target: hostTile(), message: `Select ${current().residue} in the drawing or its tile to continue.` };
    if (id === 'select') return null;
    if (!state.candidate) {
      if (find('#lib-panel').hidden) return { target: '#btn-lib', message: 'Library is closed. Reopen it to choose a monomer.', label: 'Open Library', action: openLibrary };
      return { target: '#lib-list', message: lesson === 'swap' ? 'Choose Select beside an eligible replacement, such as Orn.'
        : 'Browse or search, then choose Use beside a monomer with an unused R1 site, such as A or L.' };
    }
    if (id === 'choose') return null;
    if (!state.connection) return lesson === 'swap'
      ? { target: '#build-swap-mapping', message: 'Give each connection a different replacement site, then wait for the monomer to load.' }
      : { target: state.build.left.selectedSlot === 2 ? '#build-right-rgroups' : '#build-left-rgroups',
        message: 'Choose glycine R2 and your monomer’s R1, then wait for the connection check.' };
    if (id === 'apply' && lesson === 'swap' && !state.preview) return { target: '#build-preview-button', message: 'Preview this mapping before applying the swap.' };
    return null;
  }

  function check() {
    if (!disclosure.isOpen) return;
    const step = current().steps[index], state = readState(), help = recovery(state, step.id);
    const reviewing = state.review && !['undo', 'finish'].includes(step.id);
    const ready = reviewing || ({ load: state.base, explore: explored && state.base, select: state.left,
      choose: state.candidate, sites: state.connection, preview: state.preview, apply: state.product,
      undo: !!applied && state.base, finish: true })[step.id];
    next.disabled = loading || !!loadError || !!help || !ready;
    load.disabled = loading;
    show.disabled = loading;
    status.textContent = loading ? `Loading ${current().name}…` : loadError || help?.message ||
      (reviewing ? `Edit completed${state.base ? ' and undone' : ''}. Back and Next review the instructions; Restart begins again.`
        : next.disabled ? 'Complete this step to continue.' : 'Ready to continue.');
    show.textContent = reviewing ? 'Show drawing' : help?.label || 'Show control';
    reveal = reviewing ? null : help?.action;
    find('#tutorial-body').textContent = reviewing && step.id !== 'load'
      ? `${step.body} You already completed this edit with ${applied.request.new_abbr}; continue to review the next step.` : step.body;
    const element = find(reviewing ? '#render-canvas' : help?.target || step.target || hostTile());
    if (target !== element) {
      target?.classList.remove('tutorial-target');
      target = element;
      target?.classList.add('tutorial-target');
    }
  }

  function showStep() {
    const steps = current().steps;
    chooser.value = lesson;
    find('#tutorial-progress').textContent = `Step ${index + 1} of ${steps.length}`;
    find('#tutorial-title').textContent = steps[index].title;
    back.disabled = index === 0;
    load.hidden = steps[index].id !== 'load';
    load.textContent = `Load ${current().name}`;
    next.textContent = index === steps.length - 1 ? 'Finish' : 'Next';
    check();
  }
  function reset() {
    generation++;
    index = 0;
    explored = false;
    applied = null;
    loading = false;
    loadError = '';
  }
  async function restart() {
    reset();
    const ownGeneration = generation;
    loading = true;
    showStep();
    try {
      const source = lesson === 'swap' ? await getExample('Retatrutide') : 'A-G';
      if (ownGeneration !== generation) return;
      current().source = source;
      loadPractice(source);
    } catch (error) {
      if (ownGeneration === generation) loadError = `Could not load ${current().name}. Choose Load to retry. ${error.message}`;
    } finally {
      if (ownGeneration === generation) { loading = false; check(); }
    }
  }
  next.addEventListener('click', () => {
    check();
    if (next.disabled) return;
    if (index === current().steps.length - 1) disclosure.close({ focus: true });
    else { index++; showStep(); find('#tutorial-title').focus(); }
  });
  back.addEventListener('click', () => { if (index > 0) { index--; showStep(); } });
  chooser.addEventListener('change', () => { lesson = chooser.value; reset(); showStep(); });
  load.addEventListener('click', restart);
  find('#tutorial-restart').addEventListener('click', restart);
  show.addEventListener('click', () => {
    check();
    reveal?.();
    check();
    if (target?.closest('#lib-panel')?.hidden) openLibrary();
    if (target?.closest('#build-panel')?.hidden) openBuild();
    target?.scrollIntoView({ block: 'center', inline: 'nearest' });
    const control = target?.matches('button, input, textarea, select, [tabindex]')
      ? target : target?.querySelector('button:not(:disabled)');
    control?.focus({ preventScroll: true });
  });
  document.addEventListener('click', event => {
    const entry = event.target.closest('.tutorial-start');
    if (!entry) return;
    event.preventDefault();
    if (entry.dataset.lesson && entry.dataset.lesson !== lesson) { lesson = entry.dataset.lesson; reset(); }
    disclosure.open();
    showStep();
  });
  find('#tutorial-launch').hidden = true;
  button.hidden = false;
  disclosure.open({ focus: false });
  return {
    recordApplied(edit) {
      if (edit.action !== lesson || edit.request.cabiln !== current().source) return;
      const matches = lesson === 'swap' ? edit.request.residue_idx === current().host
        : edit.request.host_residue_idx === current().host && edit.request.target_residue_idx === -1 &&
          edit.request.r_host === 2 && edit.request.r_new === 1;
      if (matches) { applied = edit; check(); }
    },
  };
}
