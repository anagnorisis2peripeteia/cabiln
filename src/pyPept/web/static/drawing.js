// Keep the accepted molecule, depiction and atom order together, even while a
// later edit is pending. Only a ready CABILN drawing can supply exports or reuse.
class DrawingState {
  #displayed = null;
  #ready = false;
  #stale = false;
  #seed = 0;

  get hasDrawing() { return this.#displayed !== null; }
  get stale() { return this.#stale; }
  get seed() { return this.#seed; }
  get snapshot() {
    return this.#ready && this.#displayed?.notation === 'cabiln' ? this.#displayed : null;
  }
  get cabiln() { return this.snapshot?.source || ''; }
  get svg() { return this.snapshot?.data.svg || ''; }
  get molBlock() { return this.snapshot?.data.mol_block || ''; }

  matches(source, notation) {
    return this.#displayed?.source === source && this.#displayed?.notation === notation;
  }
  begin() { this.#ready = false; }
  invalidate() {
    this.#ready = false;
    this.#stale = true;
    this.#seed = 0;
  }
  clear() {
    this.invalidate();
    this.#displayed = null;
    this.#stale = false;
  }
  accept(data, source, notation, view = {}) {
    this.#displayed = { data, source, notation, ...view };
    this.#seed = view.seed || 0;
    this.#stale = false;
    this.#ready = true;
  }
  nextLayout() {
    this.#seed = this.#seed === 0 ? 2 : this.#seed + 1;
    return this.#seed;
  }
}
