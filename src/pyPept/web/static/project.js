// Portable project data is kept separate from rendering and request lifetimes.
const CabilnProject = (() => {
  const FORMAT = 'cabiln-project';
  const VERSION = 1;
  const MAX_FILE_BYTES = 2 * 1024 * 1024;
  const MODES = ['cabiln', 'smiles', 'biln', 'helm'];
  const clone = value => value == null ? null : JSON.parse(JSON.stringify(value));

  function object(value) {
    return value && typeof value === 'object' && !Array.isArray(value);
  }

  function quality(value) {
    if (value == null) return null;
    if (!object(value) || !['complete', 'partial', 'unresolved'].includes(value.recognition_status)) {
      throw new Error('The project contains invalid import-quality information.');
    }
    if (typeof value.search_complete !== 'boolean' || typeof value.inferred_stereo !== 'boolean') {
      throw new Error('The project contains invalid import-quality flags.');
    }
    const assignments = value.assignments || [];
    const index = item => Number.isSafeInteger(item) && item >= 0;
    if (!Array.isArray(assignments) || assignments.length > 20000 || assignments.some(item =>
      !object(item) || typeof item.symbol !== 'string' || typeof item.recognized !== 'boolean' ||
      !index(item.residue_index) || !Array.isArray(item.source_atoms) || !item.source_atoms.every(index) ||
      !Array.isArray(item.attachments) || item.attachments.some(pair =>
        !Array.isArray(pair) || pair.length !== 2 || !index(pair[0]) || pair[0] < 1 || !index(pair[1])))) {
      throw new Error('The project contains invalid monomer assignments.');
    }
    const synthetic = value.synthetic_components || [];
    const warnings = value.warnings || [];
    if (!Array.isArray(synthetic) || !synthetic.every(index) || !Array.isArray(warnings) ||
        warnings.some(item => typeof item !== 'string')) {
      throw new Error('The project contains invalid import-quality details.');
    }
    return {
      recognition_status: value.recognition_status,
      search_complete: value.search_complete,
      inferred_stereo: value.inferred_stereo,
      assignments: clone(assignments), synthetic_components: clone(synthetic), warnings: [...warnings],
    };
  }

  function document(value, fallbackNotation = 'cabiln') {
    if (!object(value) || typeof value.text !== 'string') {
      throw new Error('The project is missing its editable sequence.');
    }
    const notation = value.notation || fallbackNotation;
    if (!MODES.includes(notation)) throw new Error('The project uses an unsupported input notation.');
    if (value.text.length > 20000) throw new Error('A project sequence exceeds the 20,000 character limit.');
    return {
      text: value.text, notation,
      warning: typeof value.warning === 'string' ? value.warning : '',
      quality: quality(value.quality),
      canonical: object(value.canonical) ? clone(value.canonical) : null,
      context: object(value.context) ? clone(value.context) : null,
    };
  }

  function reference(value, fallbackContext = null) {
    if (typeof value === 'string') value = { text: value };
    const result = { text: typeof value?.text === 'string' ? value.text : '', original: null };
    if (result.text.length > 20000) throw new Error('The project reference exceeds 20,000 characters.');
    if (object(value?.original)) {
      const original = value.original;
      if (!['text', 'mol'].includes(original.kind) || typeof original.content !== 'string') {
        throw new Error('The project contains an invalid original reference.');
      }
      result.original = { kind: original.kind, content: original.content,
        name: typeof original.name === 'string' ? original.name : '' };
      if (original.content.length > 1000000 || result.original.name.length > 256) {
        throw new Error('The original reference or its filename exceeds the project limit.');
      }
    }
    if (value?.context != null && !object(value.context)) {
      throw new Error('The project contains an invalid reference binding.');
    }
    result.context = clone(value?.context || ((result.text || result.original) ? fallbackContext : null));
    return result;
  }

  function drafts(value) {
    const result = {};
    if (value != null && (!object(value) || Object.keys(value).some(mode => !MODES.includes(mode)))) {
      throw new Error('The project contains unsupported notation drafts.');
    }
    for (const mode of MODES) {
      let item = value?.[mode];
      if (typeof item === 'string') item = { text: item };
      if (item != null) {
        result[mode] = document(item, mode);
        if (result[mode].notation !== mode) throw new Error('A project draft has the wrong notation.');
      }
    }
    return result;
  }

  function read(value, requireContext = true) {
    if (!object(value) || value.format !== FORMAT || value.version !== VERSION) {
      throw new Error('Choose a CABILN project (.cabiln.json) saved by this version of the app.');
    }
    if (requireContext && !object(value.context)) {
      throw new Error('This project has no library binding. Open its JSON to recover the source text, then save a new project.');
    }
    return { format: FORMAT, version: VERSION,
      document: document(value.document), drafts: drafts(value.drafts),
      reference: reference(value.reference, value.context), context: clone(value.context),
      saved_at: typeof value.saved_at === 'string' ? value.saved_at : '' };
  }

  function fromConversion(value) {
    return { quality: value.recognition_status ? quality(value) : null,
      canonical: clone(value.canonical), context: clone(value.context) };
  }

  function sorted(value) {
    if (!object(value)) return value;
    return Object.fromEntries(Object.keys(value).sort().map(key => [key, sorted(value[key])]));
  }

  function sameBinding(left, right) {
    return !!left?.library_binding && !!right?.library_binding &&
      JSON.stringify(sorted(left.library_binding)) === JSON.stringify(sorted(right.library_binding));
  }

  function sameContext(left, right) {
    return sameBinding(left, right) && left.project_version === right.project_version &&
      JSON.stringify(sorted(left.canonical)) === JSON.stringify(sorted(right.canonical));
  }

  function qualitySummary(value) {
    if (!value) return '';
    const label = { complete: 'Complete recognition', partial: 'Partial recognition',
      unresolved: 'Unresolved recognition' }[value.recognition_status];
    const parts = ['Import quality: ' + label];
    const assignments = value.assignments || [];
    if (assignments.length) {
      parts.push(`${assignments.filter(item => item.recognized === true).length}/${assignments.length} library matches`);
    }
    if (value.search_complete === false) parts.push('search limited');
    if (value.inferred_stereo === true) parts.push('stereochemistry inferred');
    return parts.join(' · ');
  }

  return { FORMAT, VERSION, MAX_FILE_BYTES, clone, document, drafts, reference,
    read, quality, fromConversion, sameBinding, sameContext, qualitySummary };
})();
