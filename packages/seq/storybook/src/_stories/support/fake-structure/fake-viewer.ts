/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MVSData } from '@molstar/mvs-builder';
import type { ComponentExpressionT } from '@molstar/mvs-builder/tree/mvs/param-types';
import { State } from '@molstar/core/reactive';
import {
  type Cell,
  type ChainTable,
  hasCell,
  normalizeSpans,
  observedSpans,
  type Row,
  rowsOf,
  selectorsOfSpans,
  type Span,
  spanOfCell,
  spansOfSelector,
  toggleSpans,
  withoutAuthChain,
} from './chain-table.js';
import { BaseColor, colorNodesOf, stateIssues, structureState } from './mvs-state.js';

// The fake structure viewer (structure-syncing, The fake structure viewer): a DOM-free model of a structure drawn as
// rows of residue cells, with the grammar of a viewer of its own and none of Seq*: hover; click selects; Cmd/Ctrl+click
// toggles; Esc clears; double-click focuses; "select chain" and "within 5 Å" make selections no sequence made. What it
// is told and tells is MolViewSpec's words (selectors); its state is MolViewSpec, validated before it is taken. Its
// messages to the host are delayed by `latency` (a stale hover: the host's pointer may have moved on). Views are thin
// (`fake-view.ts`); tests drive the model.

/** The camera frame: the cells it is on and the time it takes to get there. */
export interface Camera {
  readonly spans: readonly Span[];
  readonly durationMs: number;
}

export interface FakeStructureOptions {
  /** Milliseconds a message to the host is delayed (0, the default: at once). */
  readonly latency?: number;
  /** How long its own focus takes. */
  readonly focusMs?: number;
}

/** What the viewer tells its host: selectors of the cells in question. */
export type Selectors = readonly ComponentExpressionT[];

/** Cells within `RadiusCells` of a selected cell on its row are "within 5 Å" (the fake has no coordinates). */
export const RadiusCells = 3;

export class FakeStructureViewer {
  private _rows: readonly Row[];
  private _table: ChainTable;
  /** Counts the table's edits (`dropChain`); `onEdit` tells a host that holds the model. */
  readonly edits = new State(0);
  private readonly _editFns = new Set<() => void>();
  readonly entry: string;
  readonly latency: number;
  /** Its own pointer's hover. */
  readonly hover = new State<readonly Span[]>([]);
  /** What the host lit (apart from the hover: a host's `highlight(null)` never clears the pointer's). */
  readonly highlight = new State<readonly Span[]>([]);
  readonly selection = new State<readonly Span[]>([]);
  /** The camera frame; null: it has not moved. */
  readonly camera = new State<Camera | null>(null);
  /** The MolViewSpec state: its `color` nodes are what the cells are drawn with. */
  readonly mvs: State<MVSData>;
  private _colors = new Map<number, string>();
  private _base: string = BaseColor;
  private readonly _hoverFns = new Set<(s: Selectors, at: number) => void>();
  private readonly _selectFns = new Set<(s: Selectors, at: number) => void>();
  private readonly _focusFns = new Set<(s: Selectors, at: number) => void>();
  private readonly _timers = new Set<ReturnType<typeof setTimeout>>();
  private readonly _abort = new AbortController();
  private readonly _focusMs: number;
  /** Messages sent to the host, by kind, in order (what a story prints and a twin reads). */
  readonly sent: { readonly kind: 'hover' | 'select' | 'focus'; readonly selectors: Selectors }[] = [];

  constructor(table: ChainTable, options: FakeStructureOptions = {}) {
    this._table = table;
    this._rows = rowsOf(table);
    this.entry = table.entry;
    this.latency = options.latency ?? 0;
    this._focusMs = options.focusMs ?? 400;
    this.mvs = new State<MVSData>(structureState(table.entry, []));
  }

  /** The rows now: a chain instance each, as the table (edited or not) has them. */
  get rows(): readonly Row[] {
    return this._rows;
  }

  /**
   * What the hovered cells are called, for a card beside them (a first user looks for the names: `A 12 GLY`, a range
   * `A 12 GLY to 15 ALA`, the unobserved ones said so): the first span's, undefined with no hover.
   */
  hoverName(): string | undefined {
    const span = this.hover.value[0];
    if (!span) return undefined;
    const { residues, asym } = this._rows[span.row];
    const name = (i: number) => {
      const r = residues[i];
      return `${r.auth}${r.ins} ${r.comp}${r.observed ? '' : ' (not modeled)'}`;
    };
    return span.end - span.start > 1
      ? `${asym} ${name(span.start)} to ${name(span.end - 1)}`
      : `${asym} ${name(span.start)}`;
  }

  /** The chain table now (the host's adapter reads its model from it). */
  get table(): ChainTable {
    return this._table;
  }

  // ---- its own edits, told to a host that holds its model ----

  /** Tells `fn` after each edit of the table (`dropChain`); returns the unsubscribe. */
  onEdit(fn: () => void): () => void {
    this._editFns.add(fn);
    return () => void this._editFns.delete(fn);
  }

  /**
   * "Drop chain": the chain whose author chain id is `auth` goes, with what is bound to it (its ligands). The cells
   * the pointer, the host and the selection held are gone with it, the colors are applied again to the rows left.
   * Returns false for a chain it does not hold.
   */
  dropChain(auth: string): boolean {
    const next = withoutAuthChain(this._table, auth);
    if (next === this._table) return false;
    this._table = next;
    this._rows = rowsOf(next);
    this.hover.set([]);
    this.highlight.set([]);
    this.selection.set([]);
    this.camera.set(null);
    // the state's colors are by selector: read again over the rows that are left
    this.setState(this.mvs.value);
    this.edits.set(this.edits.value + 1);
    for (const f of [...this._editFns]) f();
    return true;
  }

  /** The author chain id of the chain the selection or the pointer is on, else the last polymer chain's. */
  chainToDrop(): string | undefined {
    const picked = this.selection.value.length ? this.selection.value : this.hover.value;
    const row = picked.length ? this._rows[picked[0].row] : [...this._rows].reverse().find((r) => r.polymer);
    return row?.auth;
  }

  /** Aborted when the viewer goes: what a host's link ends with. */
  get signal(): AbortSignal {
    return this._abort.signal;
  }

  // ---- the host's commands: MolViewSpec's words ----

  /** Takes a MolViewSpec state when it is valid (its issues are returned and nothing changes otherwise). */
  setState(mvs: MVSData): string[] | undefined {
    const issues = stateIssues(mvs);
    if (issues) return issues;
    const colors = new Map<number, string>();
    let base: string = BaseColor;
    for (const n of colorNodesOf(mvs)) {
      if (!n.selector) {
        base = n.color;
        colors.clear();
        continue;
      }
      for (const s of n.selector.flatMap((sel) => spansOfSelector(this.rows, sel)))
        for (let i = s.start; i < s.end; i++) colors.set(key(s.row, i), n.color);
    }
    this._colors = colors;
    this._base = base;
    this.mvs.set(mvs);
    return undefined;
  }

  /** The color a cell is drawn with. */
  colorOf(row: number, index: number): string {
    return this._colors.get(key(row, index)) ?? this._base;
  }

  /** The host's highlight, selection and focus (null: cleared); unobserved cells are never taken. */
  setHighlight(selectors: Selectors | null) {
    this.highlight.set(this.spansOf(selectors));
  }
  setSelection(selectors: Selectors | null) {
    this.selection.set(this.spansOf(selectors));
  }
  focusOn(selectors: Selectors, durationMs: number) {
    const spans = this.spansOf(selectors);
    if (spans.length) this.camera.set({ spans, durationMs });
  }

  private spansOf(selectors: Selectors | null): Span[] {
    return selectors ? normalizeSpans(selectors.flatMap((s) => spansOfSelector(this.rows, s))) : [];
  }

  // ---- its own events, told to the host ----

  onHover(fn: (s: Selectors, at: number) => void): () => void {
    this._hoverFns.add(fn);
    return () => void this._hoverFns.delete(fn);
  }
  onSelect(fn: (s: Selectors, at: number) => void): () => void {
    this._selectFns.add(fn);
    return () => void this._selectFns.delete(fn);
  }
  onFocus(fn: (s: Selectors, at: number) => void): () => void {
    this._focusFns.add(fn);
    return () => void this._focusFns.delete(fn);
  }

  private tell(
    kind: 'hover' | 'select' | 'focus',
    spans: readonly Span[],
    fns: Set<(s: Selectors, at: number) => void>,
  ) {
    const selectors = selectorsOfSpans(this.rows, spans);
    // when it happened: a late message says so (a host drops what its pointer has overtaken)
    const at = performance.now();
    const send = () => {
      if (this._abort.signal.aborted) return;
      this.sent.push({ kind, selectors });
      for (const f of [...fns]) f(selectors, at);
    };
    if (this.latency <= 0) return send();
    const t = setTimeout(() => {
      this._timers.delete(t);
      send();
    }, this.latency);
    this._timers.add(t);
  }

  // ---- its grammar ----

  /** The pointer over a cell (null: left): an unobserved cell is a gap and is nothing. */
  pointerOver(cell: Cell | null) {
    const spans = cell && this.rows[cell.row].residues[cell.index].observed ? [spanOfCell(cell)] : [];
    if (sameSpans(this.hover.value, spans)) return;
    this.hover.set(spans);
    this.tell('hover', spans, this._hoverFns);
  }

  /** A click selects the cell; with Cmd/Ctrl it is toggled in the selection. */
  click(cell: Cell, toggle = false) {
    if (!this.rows[cell.row].residues[cell.index].observed) return;
    const one = [spanOfCell(cell)];
    const next = toggle ? toggleSpans(this.rows, this.selection.value, one) : one;
    this.setSpans(next);
  }

  /** Esc: the selection is cleared. */
  clear() {
    if (this.selection.value.length) this.setSpans([]);
  }

  /** A double-click focuses the cell: the camera moves over `focusMs`. */
  doubleClick(cell: Cell) {
    if (!this.rows[cell.row].residues[cell.index].observed) return;
    const spans = [spanOfCell(cell)];
    this.camera.set({ spans, durationMs: this._focusMs });
    this.tell('focus', spans, this._focusFns);
  }

  /** "Select chain": every observed residue of the rows of the selection's chains (with none: the hovered cell's, else the first). */
  selectChain() {
    const picked = this.selection.value.length ? this.selection.value : this.hover.value;
    const rows = new Set(picked.map((s) => s.row));
    const chains = [...(rows.size ? rows : [0])].map((r) => this.rows[r]);
    const same = this.rows.filter((r) => chains.some((c) => c.asym === r.asym && c.operator === r.operator));
    this.setSpans(
      observedSpans(
        this.rows,
        same.map((r) => ({ row: r.index, start: 0, end: r.residues.length })),
      ),
    );
  }

  /**
   * "Within 5 Å": the selection grown by `RadiusCells` residues along each row, and by the ligands of its chain in its
   * copy (the fake has no coordinates, so a neighbor is a neighbor along the chain).
   */
  within() {
    const sel = this.selection.value;
    if (!sel.length) return;
    const grown: Span[] = sel.map((s) => ({
      row: s.row,
      start: Math.max(0, s.start - RadiusCells),
      end: Math.min(this.rows[s.row].residues.length, s.end + RadiusCells),
    }));
    const copies = new Set(sel.map((s) => `${this.rows[s.row].auth} ${this.rows[s.row].operator}`));
    for (const r of this.rows)
      if (!r.polymer && copies.has(`${r.auth} ${r.operator}`))
        grown.push({ row: r.index, start: 0, end: r.residues.length });
    this.setSpans(observedSpans(this.rows, normalizeSpans(grown)));
  }

  /** Whether the cell is selected. */
  isSelected(cell: Cell): boolean {
    return hasCell(this.selection.value, cell);
  }

  private setSpans(spans: readonly Span[]) {
    const next = normalizeSpans(spans);
    if (sameSpans(this.selection.value, next)) return;
    this.selection.set(next);
    this.tell('select', next, this._selectFns);
  }

  dispose() {
    this._abort.abort();
    for (const t of this._timers) clearTimeout(t);
    this._timers.clear();
  }
}

const key = (row: number, index: number) => row * 1_000_000 + index;

function sameSpans(a: readonly Span[], b: readonly Span[]): boolean {
  return a.length === b.length && a.every((s, i) => s.row === b[i].row && s.start === b[i].start && s.end === b[i].end);
}
