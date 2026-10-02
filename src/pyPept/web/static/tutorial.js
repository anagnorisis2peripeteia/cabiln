function startTutorial({ lesson = 'connect', loadPractice, getExample, getBuildState, isReady, openBuild, openLibrary, guideResidue }) {
  const find = selector => document.querySelector(selector);
  const lessons = {
    connect: { source: 'A-G', name: 'A–G', host: 1, residue: 'G (glycine)', steps: [
      { id: 'load', title: 'Build a chain of three blocks', target: '#cabiln-input',
        body: 'A peptide is a chain of building blocks called amino acids. Click Load A–G to start with two blocks, A and G. You will add a third. This practice tab leaves your saved work untouched; practice edits are not saved automatically.' },
      { id: 'explore', title: 'Match a tile to the drawing',
        body: 'The coloured tiles above the drawing each represent one building block. Hover over the G tile, or tap it on a touch screen. Its part of the drawing lights up. G is short for glycine; you do not need to recognise its chemical structure.' },
      { id: 'select', title: 'Choose the marked G tile',
        body: 'Open Build, then click the G tile marked “choose” above the drawing. G appears under Current residue: “residue” means a building block already in your chain. You can also click a line or letter inside the outlined part of the main drawing.' },
      { id: 'choose', title: 'Choose a block to add', target: '#lib-list',
        body: 'Library lists building blocks, called monomers in the app. Find A (alanine) by browsing or searching, then click Use on its row. A appears under New monomer. L (leucine) works for this exercise too.' },
      { id: 'sites', title: 'Choose the connection', target: '#build-left-rgroups',
        body: 'R1, R2 and R3 label connection points on a block. Under Current residue, click R2. Under New monomer, click R1, or use the suggested pair if offered. This joins the end of G to the start of your new block. Crossed-out points are already in use.' },
      { id: 'preview', title: 'Check the chain before adding', target: '#build-preview-button',
        body: 'Click Preview. Blue glow marks {replacement}, the new block after G. Orange glow marks its new link. Under Proposed CABILN, check for A-G-{replacement}: CABILN is the text form of the drawing. Your working chain stays A-G until you click Connect.' },
      { id: 'apply', title: 'Add the new block', target: '#build-connect',
        body: 'Click Connect. The tiles above the main drawing should now read A, G, {replacement}. The new block is part of your chain.' },
      { id: 'undo', title: 'Undo the edit', target: '#btn-undo',
        body: 'Click Undo beside the sequence box. The third tile disappears and A–G returns. Redo adds it again. Back in this guide only revisits instructions; it does not undo edits.' },
      { id: 'finish', title: 'Continue exploring', target: '#btn-project-save',
        body: 'You have selected, added and removed a building block. Save project downloads an editable copy; PNG downloads a picture. Choose Swap in Retatrutide above to practise replacing a block in a larger chain, or close this tab to return to your original work.' },
    ] },
    swap: { source: null, name: 'Retatrutide', host: 16, residue: 'K at position 17', steps: [
      { id: 'load', title: 'Try a swap in Retatrutide', target: '#cabiln-input',
        body: 'Click Load Retatrutide to open an example peptide: a chain of building blocks. You will replace one block while keeping everything attached to it. The next step marks the block for you. This practice tab leaves your saved work untouched.' },
      { id: 'select', title: 'Choose the marked K tile',
        body: 'Open Build. Set Action to Swap monomer: “monomer” means a building block. Then click the K tile marked “choose” above the drawing. This is block 17, called lysine. Its part of the drawing is also outlined; you can click a line or letter inside it.' },
      { id: 'choose', title: 'Choose a replacement block', target: '#lib-list',
        body: 'Search Library for Orn and click Select on its row. Orn is the code for ornithine, another building block. Only replacements that can keep all the current connections are listed. Your chain stays unchanged while you choose.' },
      { id: 'sites', title: 'Review all three connections', target: '#build-swap-mapping',
        body: 'The three rows keep the links to the previous block K, the next block A, and the side arm labelled AEEA. R numbers label connection points. For Orn, leave the suggested R1 → R1, R2 → R2 and R4 → R4 choices, then click Next.' },
      { id: 'preview', title: 'Preview the replacement', target: '#build-preview-button',
        body: 'Click Preview swap, then Focus change to see the replacement in blue and its three retained links in orange. Reset view shows the whole chain and side arm again. The app checks all three connections before enabling Apply swap. Your working chain stays unchanged until you apply.' },
      { id: 'apply', title: 'Apply the swap', target: '#build-connect',
        body: 'Click Apply swap. The marked K tile becomes {replacement}, and the long side arm stays attached. If you change the replacement or any connection point, preview again before applying.' },
      { id: 'undo', title: 'Restore Retatrutide', target: '#btn-undo',
        body: 'Click Undo beside the sequence box. K returns in place of {replacement}. Redo restores your swap. Back in this guide only revisits instructions; it does not undo edits.' },
      { id: 'finish', title: 'Try another replacement', target: '#btn-project-save',
        body: 'You have replaced a block without disconnecting its neighbours or side arm. Try another block, or choose Build a peptide above to practise adding one. Save project downloads an editable copy. Close this practice tab to return to your original work.' },
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
        guideResidue(null);
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
    if (state.build.busy) return { target: '#build-status', message: 'Build is loading or checking your choices. Please wait.' };
    if (!state.left) return { target: hostTile(), message: `Click the marked ${current().residue} tile above the drawing.` };
    if (id === 'select') return null;
    if (!state.candidate) {
      if (find('#lib-panel').hidden) return { target: '#btn-lib', message: 'Library is closed. Reopen it to choose a building block.', label: 'Open Library', action: openLibrary };
      return { target: '#lib-list', message: lesson === 'swap' ? 'Search for Orn, then click Select on its row.'
        : 'Find A (alanine) or L (leucine) in Library, then click Use on its row.' };
    }
    if (id === 'choose') return null;
    if (!state.connection) return lesson === 'swap'
      ? { target: '#build-swap-mapping', message: 'Each row needs a different connection point. For Orn, use R1, R2 and R4 in order.' }
      : { target: state.build.left.selectedSlot === 2 ? '#build-right-rgroups' : '#build-left-rgroups',
        message: state.build.left.selectedSlot === 2 ? 'Click R1 under New monomer, then wait for the connection check.'
          : 'Click R2 under Current residue (G).' };
    if (id === 'apply' && lesson === 'swap' && !state.preview) return { target: '#build-preview-button', message: 'Click Preview swap before applying this replacement.' };
    return null;
  }

  function check() {
    if (!disclosure.isOpen) return;
    const step = current().steps[index], state = readState(), help = recovery(state, step.id);
    const reviewing = state.review && !['undo', 'finish'].includes(step.id);
    find('#tutorial-title').textContent = `${reviewing ? 'Review: ' : ''}${step.title}`;
    guideResidue(state.base && !state.review && (step.id === 'select' ||
      (!['load', 'explore', 'undo', 'finish'].includes(step.id) && !state.left)) ? current().host : null);
    const ready = reviewing || ({ load: state.base, explore: explored && state.base, select: state.left,
      choose: state.candidate, sites: state.connection, preview: state.preview, apply: state.product,
      undo: !!applied && state.base, finish: true })[step.id];
    next.disabled = loading || !!loadError || !!help || !ready;
    load.disabled = loading;
    show.disabled = loading;
    const replacement = (state.review ? applied.request.new_abbr : state.build.right?.abbr) || 'your chosen block';
    const confirmation = { load: `${current().name} is loaded. Click Next.`,
      explore: 'G lights up in the drawing. Click Next.',
      select: `${current().residue} is selected in Build. Click Next.`,
      choose: `${replacement} is selected. Click Next.`,
      sites: 'The connection choices are ready. Click Next.',
      preview: 'Preview ready. Your chain is unchanged. Click Next.',
      apply: `${replacement} is now in your chain. Click Next.`,
      undo: `${current().name} is restored. Click Next.` };
    status.textContent = loading ? `Loading ${current().name}…` : loadError || help?.message ||
      (reviewing ? `Edit completed${state.base ? ' and undone' : ''}. Back and Next review the instructions; Restart begins again.`
        : next.disabled ? 'Complete this step to continue.' : confirmation[step.id] || 'Ready to finish.');
    show.textContent = reviewing ? 'Show drawing' : help?.label || 'Show control';
    reveal = reviewing ? null : help?.action;
    const body = step.body.replaceAll('{replacement}', replacement);
    find('#tutorial-body').textContent = reviewing && step.id !== 'load'
      ? `${body} You already completed this edit with ${applied.request.new_abbr}; continue to review the next step.` : body;
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
