// The shipped scripts run in a minimal DOM with manually completed fetches.
// Aborted requests can still resolve, matching responses already queued by a browser.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

class Element {
  constructor() {
    this.value = '';
    this.innerHTML = '';
    this.textContent = '';
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
  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, []);
    this.listeners.get(type).push(listener);
  }
  dispatchEvent(event) {
    event.target ??= this;
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
    if (!elements.has(id)) elements.set(id, new Element());
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
  window.localStorage = {
    getItem: key => storage.get(key) || null,
    setItem: (key, value) => storage.set(key, value),
    removeItem: key => storage.delete(key),
  };
  const context = vm.createContext({
    console, AbortController, Blob,
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
  const dependencies = script === 'builder.js'
    ? ['requests.js', 'project.js', 'document.js'] : ['requests.js'];
  for (const name of dependencies) {
    const dependency = path.join(path.dirname(filename), name);
    vm.runInContext(fs.readFileSync(dependency, 'utf8'), context, { filename: dependency });
  }
  vm.runInContext(fs.readFileSync(filename, 'utf8'), context, { filename });
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

module.exports = { Element, page };
