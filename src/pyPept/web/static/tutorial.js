function startTutorial({ loadPractice, isReady }) {
  const find = selector => document.querySelector(selector);
  const text = selector => find(selector)?.textContent;
  const selected = (side, slot) => find(`#build-${side}-rgroups .selected`)?.textContent.startsWith(`R${slot} `);
  const steps = [
    { title: 'Start with two residues', target: '#cabiln-input',
      body: 'Load A–G: alanine joined to glycine. This practice tab leaves your saved draft untouched. Practice edits are not saved automatically.',
      ready: () => isReady('A-G') },
    { title: 'Find a residue in the drawing', target: '.res-chip[data-residue="1"]',
      body: 'Hover over or focus the G tile. Its atoms light up in the drawing. On a touch screen, tap the tile.',
      ready: () => isReady('A-G') && find('.res-chip[data-residue="1"].hover') && find('#render-inner .res-hl') },
    { title: 'Select glycine in Build',
      target: () => find('#build-panel').hidden ? '#btn-build' : '.res-chip[data-residue="1"]',
      body: 'Open Build, then select the G tile. The selected residue and its numbered attachment sites appear in Build.',
      ready: () => isReady('A-G') && !find('#build-panel').hidden && text('#build-left-abbr') === 'G' && find('#build-left-rgroups button') },
    { title: 'Choose a monomer to add', target: '#lib-search',
      body: 'In Library, search Alanine and choose Use beside A (L-Alanine). Filter can narrow the library to compatible attachment sites.',
      ready: () => text('#build-left-abbr') === 'G' && text('#build-right-abbr') === 'A' && find('#build-right-rgroups button') },
    { title: 'Choose the connection', target: '#build-left-rgroups',
      body: 'Choose R2 on glycine, then R1 on alanine. Used sites cannot be selected. Connect becomes available when the reaction is valid.',
      ready: () => isReady('A-G') && selected('left', 2) && selected('right', 1) && !find('#build-connect').disabled },
    { title: 'Preview the product', target: '#build-preview-button',
      body: 'Choose Preview. Check the proposed molecule, reaction and A–G–A sequence before applying the connection.',
      ready: () => !find('#build-preview').hidden && find('#build-preview-inner svg') && text('#build-preview-source') === 'A-G-A' },
    { title: 'Apply the connection', target: '#build-connect',
      body: 'Choose Connect. The main drawing and residue tiles update to A–G–A.',
      ready: () => isReady('A-G-A') },
    { title: 'Undo the edit', target: '#btn-undo',
      body: 'Choose Undo to return to A–G. Redo can restore the connection. Both work with notation conversions and other sequence edits too.',
      ready: () => isReady('A-G') },
    { title: 'Continue exploring', target: '#btn-project-save',
      body: 'Save project keeps an editable copy on your device. PNG and MOL export the structure. Try Examples, Swap monomer in Build, or Verify against a reference. Close this practice tab to return to your original work.',
      ready: () => true },
  ];
  const panel = find('#tutorial-panel');
  const button = find('#btn-tutorial');
  const next = find('#tutorial-next');
  const back = find('#tutorial-back');
  const load = find('#tutorial-load');
  const status = find('#tutorial-status');
  const completed = new Set();
  let index = 0, target = null, frame = null;
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

  function check() {
    const step = steps[index];
    if (step.ready()) completed.add(index);
    next.disabled = !completed.has(index);
    status.textContent = next.disabled ? 'Complete this step to continue.' : 'Ready to continue.';
    const selector = typeof step.target === 'function' ? step.target() : step.target;
    const element = find(selector);
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
    if (next.disabled) return;
    if (index === steps.length - 1) disclosure.close({ focus: true });
    else { index++; showStep(); find('#tutorial-title').focus(); }
  });
  back.addEventListener('click', () => { if (index > 0) { index--; showStep(); } });
  load.addEventListener('click', loadPractice);
  find('#tutorial-restart').addEventListener('click', () => {
    index = 0;
    completed.clear();
    loadPractice();
    showStep();
  });
  find('#tutorial-show').addEventListener('click', () => {
    target?.scrollIntoView({ block: 'center', inline: 'nearest' });
    const control = target?.matches('button, input, textarea, [tabindex]')
      ? target : target?.querySelector('button:not(:disabled)');
    control?.focus({ preventScroll: true });
  });
  document.addEventListener('click', event => {
    if (event.target.closest('.tutorial-start')) { event.preventDefault(); disclosure.open(); }
  });
  button.hidden = false;
  disclosure.open({ focus: false });
}
