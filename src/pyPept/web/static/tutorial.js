function startTutorial({ lesson = 'connect', loadPractice, getExample, getBuildState, isReady, openBuild, openLibrary, guideResidue }) {
  const monomers = createMonomerTutorial(isReady);
  const find = selector => monomers.element(selector);
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
        body: 'Click Preview below the two blocks. Blue glow marks {replacement}, your new block; orange marks its link to G. Open Reaction and notation to see A-G-{replacement}, the text form of the drawing. Your chain stays A-G until you click Connect.' },
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
  lessons.monomer = { ...lessons.connect, steps: monomers.steps(lessons.connect.steps) };
  const panel = find('#tutorial-panel'), button = find('#btn-tutorial');
  const next = find('#tutorial-next'), back = find('#tutorial-back'), load = find('#tutorial-load');
  const status = find('#tutorial-status'), show = find('#tutorial-show'), chooser = find('#tutorial-lesson');
  const cue = createTutorialCue(panel, show);
  let index = 0, target = null, frame = null, explored = false, applied = null, reveal = null;
  let loading = false, loadError = '', generation = 0, clearanceFrame = null;
  const current = () => lessons[lesson];
  const buildAction = () => lesson === 'swap' ? 'swap' : 'connect';
  const suggestedMonomer = () => lesson === 'monomer' ? monomers.abbr : lesson === 'swap' ? 'Orn' : 'A';
  const hostTile = () => `.res-chip[data-residue="${current().host}"]`;
  function scheduleCheck() {
    if (frame === null) frame = requestAnimationFrame(() => { frame = null; check(); });
  }
  const observer = new MutationObserver(scheduleCheck);
  const formObserver = new MutationObserver(scheduleCheck);
  const sizeObserver = new ResizeObserver(() => {
    const height = `${Math.ceil(panel.getBoundingClientRect().height) + 24}px`;
    cancelAnimationFrame(clearanceFrame);
    clearanceFrame = requestAnimationFrame(() => {
      clearanceFrame = null;
      document.documentElement.style.setProperty('--tutorial-clearance', height);
    });
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
        sizeObserver.observe(panel);
      } else {
        observer.disconnect();
        formObserver.disconnect();
        sizeObserver.disconnect();
        cancelAnimationFrame(clearanceFrame);
        clearanceFrame = null;
        document.documentElement.style.removeProperty('--tutorial-clearance');
        cancelAnimationFrame(frame);
        frame = null;
        target?.classList.remove('tutorial-target');
        target = null;
        cue.clear();
        guideResidue(null);
      }
      monomers.moveGuide(panel, find('#tutorial-cue'));
    },
  });

  function readState() {
    const base = !!current().source && isReady(current().source), build = getBuildState();
    const left = base && build.open && build.action === buildAction() && build.left?.index === current().host;
    const candidate = left && build.right?.index === null && (lesson === 'swap' ||
      build.right.rgroups.some(site => site.slot === 1 && !site.used)) &&
      (lesson !== 'monomer' || build.right.abbr === monomers.abbr);
    const connection = candidate && build.ready && (lesson === 'swap' ||
      (build.left.selectedSlot === 2 && build.right.selectedSlot === 1));
    const product = !!applied && isReady(applied.result);
    explored ||= base && !!find(`${hostTile()}.hover`) && !!find('#render-inner .res-hl');
    return { base, build, left, candidate, connection, preview: connection && !!build.preview,
      product, review: !!applied && (base || product) };
  }

  function recovery(state, id) {
    if (id === 'load' || id === 'finish' || state.review) return null;
    if (!state.base) {
      if (find('#render-pane').getAttribute('aria-busy') === 'true') {
        return { target: '#render-progress', message: 'Updating the drawing. You can continue when it is ready.' };
      }
      return { target: '#tutorial-restart', message: `Use Restart to load ${current().name} again.` };
    }
    if (id === 'explore') {
      if (!explored && find('#btn-hl').getAttribute('aria-pressed') === 'false') {
        return { target: '#btn-hl', message: 'Enable Highlight to explore the atoms.', label: 'Enable Highlight', action: () => find('#btn-hl').click() };
      }
      return null;
    }
    if (id === 'undo') return { target: '#build-connect', message: 'Apply the practice edit before trying Undo. Back revisits the earlier steps.' };
    if (!state.build.open) return { target: '#btn-build', message: 'Build is closed. Reopen it to continue.', label: 'Open Build', action: openBuild };
    if (state.build.action !== buildAction()) return { target: '#build-action', message: `Choose ${lesson === 'swap' ? 'Swap monomer' : 'Connect'} in the Build action menu for this exercise.` };
    if (state.build.busy) return { target: '#build-status', message: 'Build is loading or checking your choices. Please wait.' };
    if (!state.left) return { target: hostTile(), message: `Click the marked ${current().residue} tile above the drawing.` };
    if (id === 'select') return null;
    if (!state.candidate) {
      if (find('#lib-panel').hidden) return { target: '#btn-lib', message: 'Library is closed. Reopen it to choose a building block.', label: 'Open Library', action: openLibrary };
      return { target: '#lib-list', message: lesson === 'swap' ? 'Search for Orn, then click Select on its row.'
        : lesson === 'monomer' ? `Find ${monomers.abbr} in Library, then click Use on its row.`
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

  function actionFor(step, help) {
    let selector = help?.target || step.target || hostTile();
    let label = {
      load: `Click Load ${current().name}`,
      explore: 'Hover over G, or tap it',
      select: `Click the marked ${lesson === 'swap' ? 'K' : 'G'} tile`,
      choose: lesson === 'swap' ? 'Click Select beside Orn' : 'Click Use beside A',
      sites: 'Look over the three connections',
      preview: lesson === 'swap' ? 'Click Preview swap' : 'Click Preview',
      apply: lesson === 'swap' ? 'Click Apply swap' : 'Click Connect',
      undo: 'Click Undo', finish: 'You did it! Click Finish',
    }[step.id];
    if (step.id === 'load') selector = '#tutorial-load';
    if (selector === '#lib-list') {
      const abbr = suggestedMonomer();
      const choice = `.lib-row[data-abbr="${CSS.escape(abbr)}"] .lib-use`;
      selector = find(choice) ? choice : '#lib-search';
      label = selector === '#lib-search' ? `Search for ${abbr}`
        : `Click ${lesson === 'swap' ? 'Select' : 'Use'} beside ${abbr}`;
    }
    if (selector === '#build-left-rgroups' || selector === '#build-right-rgroups') {
      const slot = selector === '#build-left-rgroups' ? 2 : 1;
      label = `Click R${slot} under ${slot === 2 ? 'Current residue' : 'New monomer'}`;
      selector += ` button[title^="R${slot}:"]`;
    }
    label = ({ '#btn-build': 'Click Build', '#btn-lib': 'Click Library',
      '#btn-hl': 'Turn Highlight on', '#tutorial-restart': 'Click Restart to try again',
      '#build-action': `Choose ${lesson === 'swap' ? 'Swap monomer' : 'Connect'}`,
      '#build-preview-button': lesson === 'swap' ? 'Click Preview swap' : 'Click Preview',
      '#build-connect': lesson === 'swap' ? 'Click Apply swap' : 'Click Connect',
      '#build-status': 'Wait for the connection check', '#render-progress': 'Wait for the drawing',
    })[selector] || label;
    if (selector === hostTile() && step.id !== 'explore') label = `Click the marked ${lesson === 'swap' ? 'K' : 'G'} tile`;
    return { element: find(selector), label };
  }

  function check(revealTarget = false) {
    if (!disclosure.isOpen) return;
    const step = current().steps[index], state = readState();
    const detail = lesson === 'monomer' ? monomers.inspect(step.id, state) : null;
    const help = detail ? null : recovery(state, step.id);
    const reviewing = !detail && state.review && !['undo', 'finish'].includes(step.id);
    find('#tutorial-title').textContent = `${reviewing ? 'Review: ' : ''}${step.title}`;
    guideResidue(state.base && !state.review && (step.id === 'select' ||
      (['choose', 'sites', 'preview', 'apply'].includes(step.id) && !state.left)) ? current().host : null);
    const ready = detail ? detail.ready : reviewing || ({ load: state.base, explore: explored && state.base, select: state.left,
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
    status.textContent = loading ? `Loading ${current().name}…` : loadError || detail?.message || help?.message ||
      (reviewing ? `Edit completed${state.base ? ' and undone' : ''}. Back and Next review the instructions; Restart begins again.`
        : next.disabled ? 'Complete this step to continue.' : confirmation[step.id] || 'Ready to finish.');
    const previewReady = step.id === 'preview' && state.preview;
    const showLabel = reviewing ? 'Show drawing' : previewReady ? 'Show preview' : help?.label || 'Show control';
    // Replacing the text during a press can cancel WebKit's click event.
    if (show.textContent !== showLabel) show.textContent = showLabel;
    reveal = reviewing ? null : help?.action;
    const body = step.body.replaceAll('{replacement}', replacement).replaceAll('{custom}', monomers.abbr);
    find('#tutorial-body').textContent = reviewing && step.id !== 'load'
      ? `${body} You already completed this edit with ${applied.request.new_abbr}; continue to review the next step.` : body;
    const element = find(reviewing ? '#render-canvas' : detail?.target || help?.target || (previewReady ? '#build-preview-canvas' : step.target) || hostTile());
    const action = detail ? { element, label: detail.label } : actionFor(step, help), complete = !next.disabled;
    const reviewConnections = lesson === 'swap' && step.id === 'sites' && complete && !reviewing;
    const reviewAttachments = lesson === 'monomer' && step.id === 'attachments' && complete && !!monomers.registrationDocument;
    const reviewDetails = reviewConnections || reviewAttachments;
    const focusTarget = complete ? element : action.element || element;
    if (target !== focusTarget) {
      target?.classList.remove('tutorial-target');
      target = focusTarget;
      target?.classList.add('tutorial-target');
      if (reviewAttachments) target?.scrollIntoView({ block: 'center', inline: 'nearest' });
    }
    panel.classList.toggle('step-complete', complete);
    find('#tutorial-instruction').textContent = reviewConnections ? 'Check the three links, then click Next.'
      : reviewAttachments ? 'Find R1 and R2, then click Next.' : complete
      ? (step.id === 'finish' ? 'You did it! Click Finish.' : 'Click Next to keep going.') : action.label;
    find('#tutorial-meter-fill').style.transform = `scaleX(${(index + Number(complete)) / current().steps.length})`;
    monomers.moveGuide(panel, find('#tutorial-cue'));
    if (loading || state.build.busy || (!state.base && find('#render-pane').getAttribute('aria-busy') === 'true')) cue.clear();
    else cue.point(reviewDetails ? element : complete ? next : action.element,
      reviewConnections ? 'These three links stay attached' : reviewAttachments ? 'Your block’s connection points'
        : complete ? (step.id === 'finish' ? 'Finish!' : 'Next step →') : action.label,
      { complete, reveal: revealTarget });
  }

  function showStep() {
    const steps = current().steps;
    chooser.value = lesson;
    find('#tutorial-progress').textContent = `Step ${index + 1} of ${steps.length}`;
    back.disabled = index === 0;
    load.hidden = steps[index].id !== 'load';
    load.textContent = `Load ${current().name}`;
    next.textContent = index === steps.length - 1 ? 'Finish' : 'Next';
    check(true);
  }
  function reset() {
    generation++;
    index = 0;
    explored = false;
    applied = null;
    loading = false;
    loadError = '';
    monomers.reset();
    monomers.closeForm();
    lessons.monomer.steps = monomers.steps(lessons.connect.steps);
    const url = new URL(window.location.href);
    url.searchParams.set('tutorial', lesson === 'connect' ? '1' : lesson);
    url.searchParams.delete('reopen');
    url.hash = '';
    window.history.replaceState(null, '', url);
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
      ? target : target?.querySelector('button:not(:disabled), input:not(:disabled)');
    if (control && control.ownerDocument !== document) control.ownerDocument.defaultView.focus();
    control?.focus({ preventScroll: true });
  });
  function registrationChanged() {
    monomers.recordRegistration();
    const document = monomers.registrationDocument;
    formObserver.disconnect();
    if (document?.body && disclosure.isOpen) {
      formObserver.observe(document.body, { subtree: true, childList: true, attributes: true, characterData: true });
    }
    check();
  }
  window.addEventListener('cabiln-library-changed', () => { monomers.recordRegistration(); check(); });
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
    registrationChanged,
    get lesson() { return lesson; },
    recordProject(kind, project) {
      if (lesson === 'monomer') { monomers.recordProject(kind, project); check(); }
    },
    recordApplied(edit) {
      if (edit.action !== buildAction() || edit.request.cabiln !== current().source) return;
      const matches = lesson === 'swap' ? edit.request.residue_idx === current().host
        : edit.request.host_residue_idx === current().host && edit.request.target_residue_idx === -1 &&
          edit.request.r_host === 2 && edit.request.r_new === 1;
      if (matches && (lesson !== 'monomer' || edit.request.new_abbr === monomers.abbr)) { applied = edit; check(); }
    },
  };
}

function createTutorialCue(panel, show) {
  const cue = document.getElementById('tutorial-cue');
  const ring = document.getElementById('tutorial-cue-ring');
  const pointer = document.getElementById('tutorial-cue-pointer');
  const label = document.getElementById('tutorial-cue-label');
  const arrow = pointer.querySelector('svg'), path = arrow.querySelector('path');
  const reduced = matchMedia('(prefers-reduced-motion: reduce)');
  let target = null, caption = '', frame = null, animatedTarget = null;
  let targetDocument = null;

  function stopMotion() {
    for (const animation of cue.getAnimations({ subtree: true })) animation.cancel();
  }
  function schedule() {
    if (target && frame === null) frame = requestAnimationFrame(() => { frame = null; position(); });
  }
  const resize = new ResizeObserver(schedule);

  function screenBox(element) {
    const box = element.getBoundingClientRect();
    const host = element.ownerDocument.defaultView?.frameElement;
    const offset = host?.getBoundingClientRect();
    const x = offset ? offset.left + host.clientLeft : 0;
    const y = offset ? offset.top + host.clientTop : 0;
    return { left: box.left + x, top: box.top + y, right: box.right + x,
      bottom: box.bottom + y, width: box.width, height: box.height };
  }

  function visibleBox(element) {
    if (!element?.isConnected || !element.getClientRects().length) return null;
    const box = screenBox(element), host = element.ownerDocument.defaultView?.frameElement;
    let left = Math.max(8, box.left), top = Math.max(8, box.top);
    let right = Math.min(innerWidth - 8, box.right), bottom = Math.min(innerHeight - 8, box.bottom);
    for (let parent = element.parentElement; parent; parent = parent.parentElement) {
      const style = getComputedStyle(parent), bounds = screenBox(parent);
      if (/(auto|scroll|hidden|clip)/.test(style.overflowX)) {
        left = Math.max(left, bounds.left); right = Math.min(right, bounds.right);
      }
      if (/(auto|scroll|hidden|clip)/.test(style.overflowY)) {
        top = Math.max(top, bounds.top); bottom = Math.min(bottom, bounds.bottom);
      }
    }
    if (host) {
      const bounds = host.getBoundingClientRect();
      left = Math.max(left, bounds.left); right = Math.min(right, bounds.right);
      top = Math.max(top, bounds.top); bottom = Math.min(bottom, bounds.bottom);
    }
    if (!panel.contains(element)) {
      const guide = panel.getBoundingClientRect();
      if (innerWidth > 960) right = Math.min(right, guide.left - 8);
      else bottom = Math.min(bottom, guide.top - 8);
    }
    if (right - left < Math.min(24, box.width * .75) || bottom - top < Math.min(20, box.height * .75)) return null;
    return { left, top, right, bottom };
  }

  function position() {
    if (!target || document.hidden) { cue.hidden = true; stopMotion(); return; }
    let element = target, box = visibleBox(element);
    if (!box) { element = show; box = visibleBox(element); }
    if (!box) { cue.hidden = true; stopMotion(); return; }
    label.textContent = element === target ? caption : 'Show me where';
    cue.hidden = false;
    const guide = panel.getBoundingClientRect(), inside = panel.contains(element);
    cue.classList.toggle('internal', inside);
    const area = { left: inside ? guide.left + 8 : 8, top: inside ? guide.top + 8 : 8,
      right: inside ? guide.right - 8 : innerWidth > 960 ? guide.left - 8 : innerWidth - 8,
      bottom: inside ? guide.bottom - 8 : innerWidth > 960 ? innerHeight - 8 : guide.top - 8 };
    const width = pointer.offsetWidth, height = pointer.offsetHeight;
    const below = box.top - height - 10 < area.top;
    const x = Math.max(area.left, Math.min((box.left + box.right - width) / 2, area.right - width));
    const y = below ? Math.min(box.bottom + 10, area.bottom - height) : box.top - height - 10;
    pointer.style.left = `${x}px`;
    pointer.style.top = `${Math.max(area.top, y)}px`;
    pointer.classList.toggle('below', below);
    arrow.style.marginLeft = `${Math.max(0, Math.min(width - 24, (box.left + box.right) / 2 - x - 12))}px`;
    path.setAttribute('d', below ? 'M12 30V2M5 9l7-7 7 7' : 'M12 2v28M5 23l7 7 7-7');
    Object.assign(ring.style, { left: `${box.left - 4}px`, top: `${box.top - 4}px`,
      width: `${box.right - box.left + 8}px`, height: `${box.bottom - box.top + 8}px` });
    if (animatedTarget !== element) {
      stopMotion();
      animatedTarget = element;
      if (!reduced.matches) {
        // Three nudges draw attention, then leave a steady arrow to follow.
        const timing = { duration: 1100, iterations: 3, easing: 'ease-in-out' };
        ring.animate([{ opacity: .35 }, { opacity: 1 }, { opacity: .35 }], timing);
        if (!inside) arrow.animate([{ transform: 'translateY(0)' }, { transform: `translateY(${below ? -4 : 4}px)` }, { transform: 'translateY(0)' }], timing);
      }
    }
  }
  window.addEventListener('resize', schedule);
  document.addEventListener('scroll', schedule, true);
  document.addEventListener('visibilitychange', () => {
    if (document.hidden) { cue.hidden = true; stopMotion(); }
    else schedule();
  });
  reduced.addEventListener('change', () => { stopMotion(); schedule(); });
  return {
    point(element, text, { complete, reveal }) {
      if (target !== element) {
        resize.disconnect();
        targetDocument?.removeEventListener('scroll', schedule, true);
        target = element;
        targetDocument = target?.ownerDocument !== document ? target?.ownerDocument : null;
        targetDocument?.addEventListener('scroll', schedule, true);
        animatedTarget = null;
        for (const item of [target, target?.ownerDocument.defaultView?.frameElement, panel,
          document.getElementById('input-bar'), document.getElementById('build-panel'), document.getElementById('main')]) {
          if (item) resize.observe(item);
        }
      }
      caption = text;
      cue.classList.toggle('complete', complete);
      if (!target) { cue.hidden = true; stopMotion(); return; }
      if (reveal) target.scrollIntoView({ block: 'nearest', inline: 'nearest' });
      schedule();
    },
    clear() {
      target = animatedTarget = null;
      cancelAnimationFrame(frame);
      frame = null;
      resize.disconnect();
      targetDocument?.removeEventListener('scroll', schedule, true);
      targetDocument = null;
      cue.hidden = true;
      stopMotion();
    },
  };
}
