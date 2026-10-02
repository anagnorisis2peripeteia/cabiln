// The shipped scripts run in a minimal DOM with manually completed fetches.
// Aborted requests can still resolve, matching responses already queued by a browser.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

class Element {
  constructor() {
    this.value = '';
    this._html = '';
    this._text = '';
    this.style = {};
    this.dataset = {};
    this.children = [];
    this.listeners = new Map();
    this.attributes = new Map();
    this.disabled = false;
    this.hidden = true;
    const classes = new Set();
    this.classList = {
      add: name => classes.add(name),
      remove: name => classes.delete(name),
      contains: name => classes.has(name),
      toggle(name, enabled = !classes.has(name)) {
        if (enabled) classes.add(name);
        else classes.delete(name);
      },
    };
  }
  get innerHTML() { return this._html; }
  set innerHTML(html) {
    this._html = html;
    this._text = html.replace(/<[^>]*>/g, '').replace(/&amp;/g, '&').replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"');
    this.children = [];
  }
  get textContent() { return this._text + this.children.map(child => child.textContent).join(''); }
  set textContent(text) { this._text = text; this._html = ''; this.children = []; }
  contains(element) { return this === element || this.children.some(child => child.contains(element)); }
  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, []);
    this.listeners.get(type).push(listener);
  }
  dispatchEvent(event) {
    event.target ??= this;
    event.preventDefault ??= () => { event.defaultPrevented = true; };
    return Promise.all((this.listeners.get(event.type) || []).map(fn => fn(event)));
  }
  querySelectorAll() { return []; }
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
  getAttribute(name) { return this.attributes.get(name) ?? null; }
  querySelector() { return null; }
  appendChild(element) { this.children.push(element); }
  click() { return this.dispatchEvent({ type: 'click' }); }
  focus() {}
  scrollIntoView() {}
}

function page(script, { registration = false, storedDraft = null, now = Date.now } = {}) {
  const elements = new Map();
  const element = id => {
    if (!elements.has(id)) {
      const node = new Element();
      node.id = id;
      elements.set(id, node);
    }
    return elements.get(id);
  };
  element('notation-select').value = 'cabiln';
  element('type-in').value = 'aa';
  element('subtype-in').value = 'modified';
  for (const id of ['btn-register', 'btn-png', 'btn-mol']) element(id).disabled = true;
  const requests = [];
  const timers = new Map();
  let timerId = 0;
  const downloads = [];
  const storage = new Map(storedDraft ? [['cabiln.draft.v1', JSON.stringify(storedDraft)]] : []);
  const window = new Element();
  window.location = { search: '' };
  window.localStorage = {
    getItem: key => storage.get(key) || null,
    setItem: (key, value) => storage.set(key, value),
    removeItem: key => storage.delete(key),
  };
  const context = vm.createContext({
    console, AbortController, AbortSignal, Blob, URLSearchParams,
    Date: class extends Date { static now() { return now(); } },
    URL: { createObjectURL(blob) { downloads.push(blob); return 'blob:project'; }, revokeObjectURL() {} },
    Event: class { constructor(type) { this.type = type; } },
    document: {
      getElementById: element,
      querySelectorAll: () => [],
      querySelector: () => null,
      createElement: () => new Element(),
    },
    window,
    setTimeout(callback) { timers.set(++timerId, callback); return timerId; },
    clearTimeout(id) { timers.delete(id); },
    setInterval() {},
    fetch(url, options = {}) {
      if (url === '/capabilities') {
        return Promise.resolve({ ok: true, json: async () => ({ registration }) });
      }
      return new Promise((resolve, reject) => {
        requests.push({
          url, options, reject,
          resolve(data, ok = true, status = ok ? 200 : 400, headers = {}) {
            resolve({ ok, status, headers: { get: key => headers[key] || null }, json: async () => data });
          },
        });
      });
    },
  });
  const filename = path.join(__dirname, '../src/pyPept/web/static', script);
  const html = fs.readFileSync(path.join(path.dirname(filename), script === 'builder.js' ? 'index.html' : 'register.html'), 'utf8');
  for (const [, name] of html.matchAll(/<script src="\/static\/([^"]+)"[^>]*><\/script>/g)) {
    const dependency = path.join(path.dirname(filename), name);
    vm.runInContext(fs.readFileSync(dependency, 'utf8'), context, { filename: dependency });
  }
  return {
    element, requests, downloads, storage,
    run: source => vm.runInContext(source, context),
    async input(id, value) {
      element(id).value = value;
      await element(id).dispatchEvent({ type: 'input' });
    },
    async timers() {
      const pending = [...timers.values()];
      timers.clear();
      for (const callback of pending) callback();
      await Promise.resolve();
    },
  };
}

const flush = () => new Promise(setImmediate);

async function openBuilder(ui, source = 'G') {
  ui.run(`replaceDocument({text: ${JSON.stringify(source)}, notation: 'cabiln'})`);
  if (!ui.run('build.filters.building')) {
    await ui.element('btn-build').click();
    ui.requests.findLast(request => request.url === '/monomers').resolve([], true, 200, { 'X-Library-Version': 'old' });
    await flush();
    ui.requests.findLast(request => request.url === '/reactions').resolve([]);
    await flush();
  }
  ui.requests.length = 0;
}

async function loadSelection(ui, side, abbr, index, rgroups = []) {
  const pending = side === 'left' || index !== null
    ? ui.run(`build.selectResidue({abbr:${JSON.stringify(abbr)}, idx:${index}})`)
    : ui.run(`build.useMonomer(${JSON.stringify(abbr)})`);
  const request = ui.requests.at(-1);
  assertRequest(request, '/monomer_rgroups');
  request.resolve({ svg: `<svg>${abbr}</svg>`, rgroups });
  await pending;
  await flush();
}

function assertRequest(request, prefix) {
  if (!request?.url.startsWith(prefix)) throw new Error(`Expected ${prefix}, got ${request?.url}`);
}

async function chooseSite(ui, side, slot) {
  const button = ui.element(`build-${side}-rgroups`).children.find(child => child.textContent.startsWith(`R${slot} `));
  if (!button) throw new Error(`Missing ${side} R${slot}`);
  await button.click();
}

async function prepareConnection(ui, { source = 'G', left = 'G', right = 'A', leftSlot = 2, rightSlot = 1, index = 0 } = {}) {
  await openBuilder(ui, source);
  await loadSelection(ui, 'left', left, index, [{slot: leftSlot, chem_type: 'backbone_c'}]);
  await loadSelection(ui, 'right', right, null, [{slot: rightSlot, chem_type: 'backbone_n'}]);
  await chooseSite(ui, 'left', leftSlot);
  // Each contract starts at the validation request, after setup has settled.
  ui.requests.length = 0;
  await chooseSite(ui, 'right', rightSlot);
  assertRequest(ui.requests[0], '/validate_bond');
}

module.exports = { Element, page, flush, openBuilder, loadSelection, chooseSite, prepareConnection };
