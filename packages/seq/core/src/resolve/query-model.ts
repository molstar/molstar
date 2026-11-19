/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The query forms' model over the static model and an overlay (_spec/README.md "Overlay", The query model): the
// spec's systems and annotations, then the overlay's (item systems, dynamic annotations), and `(node :key k)` naming a
// node the store holds. An address of a dynamic node resolves while the node is held and waits while it is not. With
// the store's tiles, `(row :key k)` also names a row of an alignment read from a source, once its rows are held.

import type { Atom, Node } from '@molstar/core/repr/tree';
import type { DataStore } from '@molstar/seqx-core/data';
import type { QueryAnnotation, QueryModel, QuerySystem } from '@molstar/seqx-core/query';
import { queryAnnotation, type StaticModel } from './model.js';
import { nodeAddress, type Overlay } from './overlay.js';

/** What a row form reads of a source alignment's rows: the store's loaded tiles (a `DataStore`). */
export type RowsReader = Pick<DataStore, 'coverage' | 'keysOf' | 'get'>;

const models = new WeakMap<StaticModel, WeakMap<Overlay, Map<RowsReader | undefined, QueryModel>>>();

/**
 * The query model of the static model with an overlay's dynamic nodes; one per model, overlay and reader (kept). With
 * `rows`, a row form reads the rows of an alignment from a source in the tiles it holds (`storeRow`).
 */
export function overlayQueryModel(model: StaticModel, overlay: Overlay, rows?: RowsReader): QueryModel {
  let byOverlay = models.get(model);
  if (!byOverlay) models.set(model, (byOverlay = new WeakMap()));
  let byReader = byOverlay.get(overlay);
  if (!byReader) byOverlay.set(overlay, (byReader = new Map()));
  let q = byReader.get(rows);
  if (q) return q;
  const systems = new Map<string, QuerySystem | undefined>();
  const annotations = new Map<string, QueryAnnotation | undefined>();
  q = {
    system(address) {
      const s = model.system(address);
      if (s) return s;
      if (systems.has(address)) return systems.get(address);
      const e = overlay.system(address);
      const out: QuerySystem | undefined = e
        ? {
            address,
            length: e.length ?? model.config.openQueryLength,
            cyclic: e.cyclic,
            numbering: () => undefined,
            numberingIds: () => [],
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
          : undefined;
      annotations.set(address, out);
      return out;
    },
    node(owner: string, key: Atom): Node | undefined {
      // the store's keys are strings: a number never names one (`1` is not `"1"`)
      return typeof key === 'string' ? overlay.node(nodeAddress(owner, key))?.node : undefined;
    },
    chain: () => undefined,
    row: (at, key) => model.row?.(at, key) ?? (rows ? storeRow(model, rows, at, key) : undefined),
  };
  byReader.set(rows, q);
  return q;
}

/**
 * A row of an alignment read from a source by its key (its label), in the tiles of its rows the store has loaded: the
 * first row of the label; undefined while its rows are not held (EC-35: rows are positions).
 */
function storeRow(model: StaticModel, reader: RowsReader, rows: string, key: Atom): number | undefined {
  for (const a of model.annotations.values()) {
    if (a.shape !== 'cells' || a.rows !== rows || a.access?.kind !== 'source') continue;
    const { source, system } = a.access;
    const n = reader.coverage(source)?.systems.find((s) => s.name === system)?.length;
    if (n === undefined) return undefined;
    const label = String(key);
    let best: number | undefined;
    for (const k of reader.keysOf({ source, system, start: 0, end: n, level: 1 })) {
      const e = reader.get(k);
      const t = e?.status === 'loaded' ? e.tables.cells : undefined;
      if (t?.form !== 'codes' || !t.label) continue;
      const m = Math.floor(t.codes.length / Math.max(1, t.columns));
      for (let i = 0; i < m; i++)
        if (t.label.get(i) === label && (best === undefined || t.start + i < best)) {
          best = t.start + i;
          break;
        }
    }
    return best;
  }
  return undefined;
}
