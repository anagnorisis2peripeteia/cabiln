function startTutorial({ loadPractice, isReady, openBuild, openLibrary }) {
  const find = selector => document.querySelector(selector);
  const text = selector => find(selector)?.textContent;
  const site = (side, slot) => [...find(`#build-${side}-rgroups`).querySelectorAll('button')]
    .find(button => button.textContent.startsWith(`R${slot} `));
  const selected = (side, slot) => site(side, slot)?.classList.contains('selected');
  const steps = [
    { title: 'Start with two residues', target: '#cabiln-input',
      body: 'Load A–G: alanine joined to glycine. This practice tab leaves your saved draft untouched. Practice edits are not saved automatically.' },
    { title: 'Find a residue in the drawing', target: '.res-chip[data-residue="1"]',
      body: 'Hover over glycine in the drawing or focus the G tile to light up its atoms. On a touch screen, tap the tile.' },
    { title: 'Select glycine in Build', target: '.res-chip[data-residue="1"]',
      body: 'Open Build, then click or tap glycine in the main drawing or its G tile. Aim at an atom label or bond; drag to pan. Its numbered attachment sites appear in Build.' },
    { title: 'Choose a monomer to add', target: '#lib-list',
      body: 'Browse or search Library and choose Use beside an amino acid, such as A (alanine) or L (leucine). This exercise uses a monomer with an unused R1 site. Filter can narrow the choices.' },
    { title: 'Choose the connection', target: '#build-left-rgroups',
      body: 'Choose R2 on glycine, then R1 on your chosen monomer. Used sites cannot be selected. Connect becomes available when the reaction is valid.' },
    { title: 'Preview the product', target: '#build-preview-button',
      body: 'Choose Preview. Check the proposed molecule, reaction and sequence before applying the connection. The result uses the monomer you chose.' },
    { title: 'Apply the connection', target: '#build-connect',
      body: 'Choose Connect. The main drawing and residue tiles update to the product you previewed.' },
    { title: 'Undo the edit', target: '#btn-undo',
      body: 'Choose Undo to return to A–G. Redo can restore the connection. Both work with notation conversions and other sequence edits too.' },
    { title: 'Continue exploring', target: '#btn-project-save',
      body: 'Save project keeps an editable copy on your device. PNG and MOL export the structure. Try Examples, Swap monomer in Build, or Verify against a reference. Close this practice tab to return to your original work.' },
  ];
  const panel = find('#tutorial-panel');
  const button = find('#btn-tutorial');
  const next = find('#tutorial-next');
  const back = find('#tutorial-back');
  const load = find('#tutorial-load');
  const status = find('#tutorial-status');
  const show = find('#tutorial-show');
  let index = 0, target = null, frame = null, explored = false, proposal = null, reveal = null;
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
    const base = isReady('A-G');
    const left = base && text('#build-left-abbr') === 'G' && site('left', 2);
    const candidate = text('#build-right-label') === 'New monomer' && site('right', 1)?.disabled === false;
    const abbr = text('#build-right-abbr');
    const connection = left && candidate && selected('left', 2) && selected('right', 1) && !find('#build-connect').disabled;
    const preview = connection && !find('#build-preview').hidden && find('#build-preview-inner svg') && text('#build-preview-source');
    if (candidate && proposal?.abbr !== abbr) proposal = null;
    if (preview) proposal = { abbr, source: preview };
    explored ||= base && !!find('.res-chip[data-residue="1"].hover') && !!find('#render-inner .res-hl');
    return { base, left, candidate, connection, preview, applied: !!proposal && isReady(proposal.source) };
  }

  function recovery(state) {
    if (index === 0 || index > 6 || state.applied) return null;
    if (!state.base) return { target: '#tutorial-restart', message: 'Wait for the drawing, or use Restart to load A–G again.' };
    if (index === 1) {
      if (!explored && find('#btn-hl').getAttribute('aria-pressed') === 'false') {
        return { target: '#btn-hl', message: 'Enable Highlight to explore the atoms.', label: 'Enable Highlight', action: () => find('#btn-hl').click() };
      }
      return null;
    }
    if (find('#build-panel').hidden) return { target: '#btn-build', message: 'Build is closed. Reopen it to continue.', label: 'Open Build', action: openBuild };
    if (find('#build-action').value !== 'connect') return { target: '#build-action', message: 'Choose Connect in the Build action menu for this exercise.' };
    if (!state.left) return { target: '.res-chip[data-residue="1"]', message: 'Select glycine in the drawing or its G tile to continue.' };
    if (index >= 3 && !state.candidate) {
      if (find('#lib-panel').hidden) return { target: '#btn-lib', message: 'Library is closed. Reopen it to choose a monomer.', label: 'Open Library', action: openLibrary };
      return { target: '#lib-list', message: 'Browse or search, then choose Use beside a monomer with an unused R1 site, such as A or L.' };
    }
    if (index >= 4 && !state.connection) return {
      target: selected('left', 2) ? '#build-right-rgroups' : '#build-left-rgroups',
      message: 'Choose glycine R2 and your monomer’s R1, then wait for the connection check.',
    };
    if (index === 6 && !proposal) return { target: '#build-preview-button', message: 'The selection changed. Preview this product before connecting.' };
    return null;
  }

  function check() {
    const state = readState(), help = recovery(state);
    const ready = (state.applied && index >= 2 && index <= 6) ||
      [state.base, explored && state.base, state.left, state.candidate,
        state.connection, state.preview, state.applied, state.base, true][index];
    next.disabled = !!help || !ready;
    status.textContent = help?.message || (next.disabled ? 'Complete this step to continue.' : 'Ready to continue.');
    show.textContent = help?.label || 'Show control';
    reveal = help?.action;
    const element = find(help?.target || steps[index].target);
    if (target !== element) {
      target?.classList.remove('tutorial-target');
      target = element;
      target?.classList.add('tutorial-target');
    }
  }

  function showStep() {
    find('#tutorial-progress').textContent = `Step ${index + 1} of ${steps.length}`;
    find('#tutorial-title').textContent = steps[index].title;
    find('#tutorial-body').textContent = steps[index].body;
    back.disabled = index === 0;
    load.hidden = index !== 0;
    next.textContent = index === steps.length - 1 ? 'Finish' : 'Next';
    check();
  }
  next.addEventListener('click', () => {
    check();
    if (next.disabled) return;
    if (index === steps.length - 1) disclosure.close({ focus: true });
    else { index++; showStep(); find('#tutorial-title').focus(); }
  });
  back.addEventListener('click', () => { if (index > 0) { index--; showStep(); } });
  function restart() {
    index = 0;
    explored = false;
    proposal = null;
    loadPractice();
    showStep();
  }
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
    if (event.target.closest('.tutorial-start')) { event.preventDefault(); disclosure.open(); }
  });
  find('#tutorial-launch').hidden = true;
  button.hidden = false;
  disclosure.open({ focus: false });
}
