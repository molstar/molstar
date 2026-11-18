/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The query forms' model over the static model and an overlay (_spec/README.md "Overlay", The query model): the
// spec's systems and annotations, then the overlay's (item systems, dynamic annotations), and `(node :key k)` naming a
// node the store holds. An address of a dynamic node resolves while the node is held and waits while it is not.

import type { Atom, Node } from '@molstar/core/repr/tree';
import type { QueryAnnotation, QueryModel, QuerySystem } from '@molstar/seqx-core/query';
import { queryAnnotation, type StaticModel } from './model.js';
import { nodeAddress, type Overlay } from './overlay.js';

const models = new WeakMap<StaticModel, WeakMap<Overlay, QueryModel>>();

/** The query model of the static model with an overlay's dynamic nodes; one per model and overlay (kept). */
export function overlayQueryModel(model: StaticModel, overlay: Overlay): QueryModel {
  let byOverlay = models.get(model);
  if (!byOverlay) models.set(model, (byOverlay = new WeakMap()));
  let q = byOverlay.get(overlay);
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
  };
  byOverlay.set(overlay, q);
  return q;
}
