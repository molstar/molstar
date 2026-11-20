/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The query forms' model over the static model and an overlay (_spec/README.md "Overlay", The query model): the
// spec's systems and annotations, then the overlay's (item systems, dynamic annotations), and `(node :key k)` naming a
// node the store holds. An address of a dynamic node resolves while the node is held and waits while it is not. With
// the store's tiles, `(row :key k)` also names a row of an alignment read from a source once its rows are held, and
// waits, asking for them, while they are not.

import type { Atom, Node } from '@molstar/core/repr/tree';
import type { CellsTable, DataStore, Want } from '@molstar/seqx-core/data';
import type { QueryAnnotation, QueryModel, QuerySystem } from '@molstar/seqx-core/query';
import { pairsQueryAnnotation, queryAnnotation, type StaticModel } from './model.js';
import { atomsOf, componentsOf, moleculeOf, unitsOf } from './molecule-query.js';
import { nodeAddress, type Overlay } from './overlay.js';

/** What a row form reads of an alignment's rows read from a source: the store's loaded tiles, and how to ask for them. */
export interface RowsSource {
  readonly data: Pick<DataStore, 'coverage' | 'keysOf' | 'get'> & Partial<Pick<DataStore, 'tooMany'>>;
  /** Asks for the rows' tiles a form named (whoever resolves asks: a row form needs no layer to have read them). */
  readonly want?: (want: Want) => void;
}

const models = new WeakMap<StaticModel, WeakMap<Overlay, Map<RowsSource | undefined, QueryModel>>>();

/**
 * The query model of the static model with an overlay's dynamic nodes; one per model, overlay and rows source (kept).
 * With `rows`, a row form reads the rows of an alignment from a source in the tiles the store holds (`storeRows`),
 * pending while they are not.
 */
export function overlayQueryModel(model: StaticModel, overlay: Overlay, rows?: RowsSource): QueryModel {
  let byOverlay = models.get(model);
  if (!byOverlay) models.set(model, (byOverlay = new WeakMap()));
  let byReader = byOverlay.get(overlay);
  if (!byReader) byOverlay.set(overlay, (byReader = new Map()));
  let q: QueryModel | undefined = byReader.get(rows);
  if (q) return q;
  const systems = new Map<string, QuerySystem | undefined>();
  const annotations = new Map<string, QueryAnnotation | undefined>();
  q = {
    system(address) {
      const s = model.system(address);
      if (s) return s;
      if (systems.has(address)) return systems.get(address);
      const e = overlay.system(address);
      // a dynamic unit's numberings are its own inline labels; its molecule's number it where it has none of an id
      const own = overlay.molecular.numberings?.get(address);
      const parent = e?.parent && q!.system(e.parent.system);
      const out: QuerySystem | undefined = e
        ? {
            address,
            length: e.length ?? model.config.openQueryLength,
            cyclic: e.cyclic,
            numbering: (id) => own?.get(id),
            numberingIds: () => (own ? [...own.keys()] : []),
            ...(parent ? { parent: { system: parent, start: e.parent!.start } } : {}),
          }
        : undefined;
      systems.set(address, out);
      return out;
    },
    annotation(address) {
      const a = model.annotation(address);
      if (a) return a;
      if (annotations.has(address)) return annotations.get(address);
      const e = overlay.annotation(address);
      const out =
        e?.shape === 'ranges'
          ? queryAnnotation(e, model.systems.get(e.target) ?? overlay.system(e.target), model.config)
          : e?.shape === 'pairs'
            ? pairsQueryAnnotation(e, (s) => model.systems.get(s) ?? overlay.system(s), model.config)
            : undefined;
      annotations.set(address, out);
      return out;
    },
    components: (s) =>
      model.components?.(s) ??
      componentsOf(overlay.molecular, overlay.system(s), (a) => model.systems.get(a) ?? overlay.system(a)),
    compoundClass: (code) => model.compoundClass?.(code) ?? overlay.molecular.compounds?.get(code)?.class,
    annotationsOn: (s) => [
      ...(model.annotationsOn?.(s) ?? []),
      ...overlay.nodes.flatMap((n) => {
        const a = overlay.annotation(n.address);
        return a?.shape === 'ranges' && a.target === s ? [a.address] : [];
      }),
    ],
    units: (c) => model.units?.(c) ?? unitsOf(overlay.molecular, c),
    molecule: (a) => model.molecule?.(a) ?? moleculeOf(overlay.molecular, a),
    atoms: (s) => model.atoms?.(s) ?? atomsOf(overlay.molecular, s),
    node(owner: string, key: Atom): Node | undefined {
      // the store's keys are strings: a number never names one (`1` is not `"1"`)
      return typeof key === 'string' ? overlay.node(nodeAddress(owner, key))?.node : undefined;
    },
    chain: () => undefined,
    row: (at, key) => {
      const own = model.row?.(at, key);
      if (own !== undefined || !rows) return own;
      const held = storeRows(model, rows, at);
      return held?.status === 'held' ? rowOf(held.tables, String(key)) : undefined;
    },
    rowsPending: (at) => (rows ? storeRows(model, rows, at)?.status === 'pending' : false),
  };
  byReader.set(rows, q);
  return q;
}

type HeldRows = { readonly status: 'pending' } | { readonly status: 'held'; readonly tables: readonly CellsTable[] };

const Pending: HeldRows = { status: 'pending' };

/**
 * The tiles of an alignment's rows read from a source (`rows`, its rows system), asked for whole; pending while any is
 * not loaded (or the source not attached); undefined when the rows are not a source's (EC-35: rows are positions).
 */
function storeRows(model: StaticModel, src: RowsSource, rows: string): HeldRows | undefined {
  for (const a of model.annotations.values()) {
    if (a.shape !== 'cells' || a.rows !== rows || a.access?.kind !== 'source') continue;
    const { source, system } = a.access;
    const n = src.data.coverage(source)?.systems.find((s) => s.name === system)?.length;
    if (n === undefined) return Pending;
    const want = { source, system, start: 0, end: n, level: 1 };
    src.want?.(want);
    const keys = src.data.keysOf(want);
    // none listed past the store's bound on a want's tiles (never asked): pending, as before they are held
    if (!keys.length && src.data.tooMany?.(want)) return Pending;
    const tables: CellsTable[] = [];
    for (const k of keys) {
      const e = src.data.get(k);
      if (e?.status !== 'loaded') return Pending;
      const t = e.tables.cells;
      if (t) tables.push(t);
    }
    return { status: 'held', tables };
  }
  return undefined;
}

/** Each tile's rows by label (the first row of a label), made once a table. */
const RowIndex = new WeakMap<CellsTable, ReadonlyMap<string, number>>();

/** The first row of a label in the tiles (each indexed once), undefined when none has it. */
function rowOf(tables: readonly CellsTable[], label: string): number | undefined {
  let best: number | undefined;
  for (const t of tables) {
    if (t.form !== 'codes' || !t.label) continue;
    let index = RowIndex.get(t);
    if (!index) {
      const m = new Map<string, number>();
      const n = Math.floor(t.codes.length / Math.max(1, t.columns));
      for (let i = 0; i < n; i++) {
        const l = t.label.get(i);
        if (l !== null && !m.has(l)) m.set(l, t.start + i);
      }
      RowIndex.set(t, (index = m));
    }
    const r = index.get(label);
    if (r !== undefined && (best === undefined || r < best)) best = r;
  }
  return best;
}
