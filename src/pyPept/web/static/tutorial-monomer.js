function createMonomerTutorial(isReady) {
  const find = selector => document.querySelector(selector);
  const dialog = find('#registration-dialog'), frame = find('#registration-frame');
  const example = find('#tutorial-example'), freshTab = find('#tutorial-open-tab');
  const params = new URLSearchParams(location.search);
  let restoring = params.get('tutorial') === 'monomer' && params.get('reopen') === '1';
  let abbr = restoring ? new URLSearchParams(location.hash.slice(1)).get('monomer') || 'MyBlock' : 'MyBlock';
  if (!restoring) {
    try { abbr = sessionStorage.getItem('cabiln.tutorial.monomer') || abbr; }
    catch (_) { /* The lesson also works without browser storage. */ }
  }
  let pending = null, definition = null, saved = false, opened = false;
  if (!restoring) definition = CabilnLibrary.monomers.find(item => item.abbr === abbr) || null;
  const registration = () => dialog.open ? frame.contentWindow?.CabilnRegistration : null;
  const product = () => `A-G-${abbr}`;
  const result = (ready, target, label, message) => ({ ready, target, label, message });

  const reopenStep = { id: 'reopen', title: 'Restore your custom block', target: '#tutorial-open-tab',
    body: 'Open a fresh practice tab, then use Open project to choose peptide.cabiln.json from your downloads. The file contains your custom block and the chain. Other users’ libraries stay unchanged.' };
  const finishStep = { id: 'finish', title: 'Your project carries its building blocks', target: '#btn-project-save',
    body: 'You restored the chain and its custom block from a file. Additions belong to this tab; refreshing usually keeps them, but closing it can lose them. Save project keeps their definitions for later. Local installation is a separate option outside this lesson.' };

  example.addEventListener('click', () => registration()?.useExample({
    smiles: 'N[C@@H](CCCS)C(=O)O', abbr, name: 'Practice amino acid',
  }));

  function read() {
    const form = registration()?.state;
    if (form?.registering && form.destination === 'session') pending = form.submitted;
    const added = form?.registered || pending;
    if (added) {
      const found = CabilnLibrary.monomers.find(item => item.abbr === added.abbr &&
        Object.keys(added).every(key => JSON.stringify(item[key]) === JSON.stringify(added[key])));
      if (found && (!definition || definition.abbr !== found.abbr || definition.chuckles !== found.chuckles)) {
        definition = found;
        abbr = found.abbr;
        try { sessionStorage.setItem('cabiln.tutorial.monomer', abbr); }
        catch (_) { /* Saving a project remains available. */ }
      }
    }
    return { form, registered: !!definition };
  }

  function inspect(id, state) {
    const { form, registered } = read();
    example.hidden = !['molecule', 'attachments', 'register'].includes(id) || registered || !form;
    freshTab.hidden = id !== 'reopen' || restoring || !saved;
    freshTab.href = `/?tutorial=monomer&reopen=1#monomer=${encodeURIComponent(abbr)}`;
    if (['molecule', 'attachments', 'register'].includes(id)) {
      if (registered) return result(true, dialog.open ? '#registration-close' : '#btn-lib', '',
        `${abbr} is in this tab’s library. Click Next.`);
      if (!dialog.open) return result(false, '#register-link', 'Click Add monomer',
        'Open Library, then Add monomer. Your chain stays here while you create its new block.');
      if (!form) return result(false, '#registration-frame', 'Wait for the form', 'Loading Add monomer…');
      if (!form.smiles) return result(false, '#tutorial-example', 'Click Use practice molecule',
        'This fills a molecule and suggested name. You do not need to write SMILES.');
      if (id === 'molecule') return result(true, 'registration:#smiles-in', '', 'The molecule is entered. Click Next to detect its connection points.');
      if (form.choices) return result(false, 'registration:#attachment-choice-list', 'Choose an attachment option', 'Choose which numbered attachment sites this block should use.');
      if (!form.payload) return result(false, 'registration:#btn-preview', 'Click Preview & detect R-groups',
        form.previewing ? 'Detecting connection points. Please wait.' : 'Preview the current molecule. If options appear, choose the intended attachment sites.');
      if (id === 'attachments') return result(true, 'registration:#detected-display', '',
        'Find R1 and R2 beside the drawing. These are the connection points we will use. Click Next.');
      if (!form.payload.abbr) return result(false, 'registration:#abbr-in', 'Enter an abbreviation', 'Give your block a short, unique label.');
      if (!form.payload.name) return result(false, 'registration:#name-in', 'Enter a name', 'Give your block a readable name.');
      return result(false, 'registration:#btn-register', 'Click Add to this tab',
        form.registering ? 'Adding your block. Please wait.' : 'Review the abbreviation and name, then add the block. Save project will keep its definition.');
    }
    if (id === 'return') {
      if (!registered) return result(false, '#register-link', 'Click Add monomer', 'Add your practice block before returning to Build. Back revisits the form instructions.');
      return result(!dialog.open, 'registration:#registration-done', 'Click Show in Library',
        dialog.open ? 'Return to Library to use your new block.' : `${abbr} is ready to use in Build. Click Next.`);
    }
    if (id === 'save') {
      if (!state.product) return result(false, '#tutorial-back', 'Go Back to finish building', 'Connect your custom block to A–G before saving.');
      return result(!!saved, '#btn-project-save', 'Click Save project', saved
        ? 'The download includes your custom block. Keep the file, then click Next.'
        : 'Save the project to your device. A picture or MOL export does not contain the editable monomer definitions.');
    }
    if (id === 'reopen') {
      if (opened && isReady(product())) return result(true, '#render-canvas', '', `${abbr} and the three-block chain are restored. Click Next.`);
      if (!restoring && !saved) return result(false, '#btn-project-save', 'Click Save project', 'Save your completed practice chain before opening a fresh tab.');
      return result(false, restoring ? '#btn-project-open' : '#tutorial-open-tab',
        restoring ? 'Click Open project' : 'Click Open a fresh tab', restoring
          ? `Choose the downloaded project containing A-G-${abbr}. This fresh tab starts with its own library.`
          : 'Continue in the new tab. Your downloaded project will restore the custom block there.');
    }
    if (id === 'finish') return result(!!opened && isReady(product()), '#render-canvas', '', 'The project and its custom monomer are restored.');
    return null;
  }

  return {
    get abbr() { return abbr; },
    steps(connectionSteps) {
      if (restoring) return [reopenStep, finishStep];
      return [
        { id: 'load', title: 'Create a block of your own', target: '#tutorial-load',
          body: 'We will add a practice molecule to your tab, connect it to A–G, then save and reopen the project. Click Load A–G. Practice edits are not saved automatically.' },
        { id: 'molecule', title: 'Enter a practice molecule', target: '#register-link',
          body: 'Library → Add monomer opens the form. Click Use practice molecule to fill an example and its name. SMILES is a text description of the molecule.' },
        { id: 'attachments', title: 'See where the block can connect', target: 'registration:#detected-display',
          body: 'Click Preview & detect R-groups. The drawing and list show numbered connection points. R1 and R2 join this block to a peptide. You do not need to recognise the structure for this exercise.' },
        { id: 'register', title: 'Add the block to this tab', target: 'registration:#btn-register',
          body: 'Review the abbreviation and name, then click Add to this tab. You can change the suggested name. The abbreviation must be unique. This lesson always uses temporary storage, including on a local installation.' },
        { id: 'return', title: 'Find your block in Library', target: 'registration:#registration-done',
          body: 'Click Show in Library. The search shows your new block, while the sequence remains A–G. Closing the form also keeps a block that was successfully added.' },
        ...connectionSteps.filter(step => ['select', 'choose', 'sites', 'preview', 'apply'].includes(step.id))
          .map(step => step.id === 'choose' ? { ...step, body: 'Find {custom} in Library and click Use. Your custom block appears under New monomer.' } : step),
        { id: 'save', title: 'Keep the block with your project', target: '#btn-project-save',
          body: 'Click Save project. Keep the downloaded peptide.cabiln.json file. It includes the sequence and custom monomer definitions. Temporary additions are not installed in the shared library.' },
        reopenStep, finishStep,
      ];
    },
    inspect,
    recordRegistration: read,
    reset() { restoring = false; saved = opened = false; example.hidden = freshTab.hidden = true; },
    recordProject(kind, project) {
      const included = project.monomers?.find(item => item.abbr === abbr);
      if (project.document.text !== product() || !included || (definition && included.chuckles !== definition.chuckles)) return;
      if (kind === 'save') saved = true;
      else { opened = true; definition = included; }
    },
    element(selector) {
      return selector?.startsWith('registration:') ? frame.contentDocument?.querySelector(selector.slice(13)) : find(selector);
    },
    get registrationDocument() { return dialog.open ? frame.contentDocument : null; },
    moveGuide(panel, cue) {
      const inside = !panel.hidden && dialog.open;
      dialog.classList.toggle('with-tutorial', inside);
      const parent = inside ? dialog : document.body;
      if (panel.parentElement !== parent) parent.append(panel, cue);
    },
    closeForm() { if (dialog.open) dialog.close(); },
  };
}
