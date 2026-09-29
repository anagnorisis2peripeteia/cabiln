// Editable source, its evidence, and per-notation history have one owner.
// The last successful drawing and in-flight requests have separate lifetimes.
class CabilnDocument {
  #present = { text: '', notation: 'cabiln', warning: '', quality: null,
    canonical: null, context: null };
  #past = [];
  #future = [];
  #drafts = {};
  #typedAt = 0;

  get present() { return this.#present; }
  get past() { return this.#past; }
  get future() { return this.#future; }
  get drafts() { return this.#drafts; }

  // Server normalization and binding updates amend the current revision.
  // They must not insert an Undo step or discard Redo history.
  replace(patch, drafts = this.#drafts) {
    this.#present = { ...this.#present, ...patch };
    this.#drafts = { ...drafts, [this.#present.notation]: this.#present };
  }

  commit(next, typing = false) {
    if (Object.keys(this.#present).every(key => next[key] === this.#present[key])) {
      return false;
    }
    const now = Date.now();
    if (!typing || !this.#typedAt || now - this.#typedAt > 750 ||
        next.notation !== this.#present.notation) {
      this.#past.push(this.#present);
      if (this.#past.length > 100) this.#past.shift();
    }
    this.#future = [];
    this.#typedAt = typing ? now : 0;
    this.replace(next);
    return true;
  }

  travel(direction) {
    const from = direction === 'undo' ? this.#past : this.#future;
    const to = direction === 'undo' ? this.#future : this.#past;
    if (!from.length) return false;
    to.push(this.#present);
    this.#typedAt = 0;
    this.replace(from.pop());
    return true;
  }

  restoreDrafts(drafts) {
    this.#drafts = { ...drafts };
  }
}
