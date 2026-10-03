// A tab owns its definitions; the server token only addresses an expiring copy.
const CabilnLibrary = (() => {
  try {
    if (window.parent && window.parent !== window && window.parent.CabilnLibrary) {
      return window.parent.CabilnLibrary;
    }
  } catch (_) { /* An unrelated embedding page does not own this library. */ }
  const KEY = 'cabiln.monomers.v1';
  let state = { token: '', monomers: [] };
  let revision = 0;
  let recovery = null;
  let retained = true;
  try {
    const saved = JSON.parse(window.sessionStorage.getItem(KEY));
    if (saved && Array.isArray(saved.monomers) && saved.monomers.length <= 64 &&
        typeof saved.token === 'string') state = saved;
  } catch (_) { retained = false; }

  const clone = value => JSON.parse(JSON.stringify(value));
  const signature = value => JSON.stringify(value, (key, item) =>
    item && typeof item === 'object' && !Array.isArray(item)
      ? Object.fromEntries(Object.keys(item).sort().map(name => [name, item[name]])) : item);

  function store(next) {
    state = { token: next.token, monomers: clone(next.monomers) };
    try { window.sessionStorage.setItem(KEY, JSON.stringify(state)); retained = true; }
    catch (_) { retained = false; }
  }

  async function create(monomers, signal) {
    const response = await fetchCalculation('/session_library', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ monomers }), signal, library: null,
    });
    const data = await readResponse(response);
    if (data.error || !data.token || !Array.isArray(data.monomers)) {
      throw new Error(data.error || 'Could not prepare the temporary library.');
    }
    return data;
  }

  async function stage(additions = [], signal) {
    if (!Array.isArray(additions)) throw new Error('Project monomers must be a list.');
    const expectedRevision = revision;
    const merged = new Map(state.monomers.map(item => [item.abbr, item]));
    for (const item of additions) {
      if (!item || typeof item.abbr !== 'string') throw new Error('A project monomer has no abbreviation.');
      const previous = merged.get(item.abbr);
      if (previous && signature(previous) !== signature(item)) {
        throw new Error(`${item.abbr} already has another definition in this tab. Open this project in a new tab.`);
      }
      merged.set(item.abbr, item);
    }
    if (merged.size > 64) throw new Error('A tab can hold up to 64 custom monomers.');
    const monomers = [...merged.values()];
    const candidate = signature(monomers) === signature(state.monomers)
      ? clone(state) : await create(monomers, signal);
    return { ...candidate, expectedRevision };
  }

  function commit(candidate, notify = true) {
    if (candidate.expectedRevision !== revision) {
      throw new Error('The tab library changed during this operation. Try again.');
    }
    const changed = signature(candidate.monomers) !== signature(state.monomers);
    if (changed) revision++;
    store(candidate);
    if (changed && notify) window.dispatchEvent(new Event('cabiln-library-changed'));
  }

  async function add(monomer, signal) {
    if (state.monomers.some(item => item.abbr === monomer.abbr)) {
      throw new Error(`${monomer.abbr} is already in this tab’s library.`);
    }
    const candidate = await stage([monomer], signal);
    commit(candidate);
    return candidate;
  }

  async function request(url, { library, ...options } = {}) {
    let selected = library === undefined ? state : library;
    const send = () => fetch(url, {
      ...options, headers: { ...options.headers,
        ...(selected?.token ? { 'X-Cabiln-Library': selected.token } : {}) },
    });
    if (!selected?.monomers.length) return send();
    let response = selected.token ? await send() : null;
    if (response && response.headers?.get('X-Cabiln-Library-Expired') !== '1') return response;
    await response?.body?.cancel();
    if (selected === state) {
      const original = state;
      if (!recovery) {
        recovery = create(original.monomers).then(restored => {
          if (state === original) store(restored);
          return restored;
        }).finally(() => { recovery = null; });
      }
      selected = await recovery;
    } else Object.assign(selected, await create(selected.monomers, options.signal));
    return send();
  }

  return { request, stage, commit, add,
    get monomers() { return clone(state.monomers); },
    get temporary() { return state.monomers.length > 0; },
    get retained() { return retained; },
  };
})();
window.CabilnLibrary = CabilnLibrary;
