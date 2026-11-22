/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover and selection through a map read by region (_spec/README.md "Render", A liftover: projection): the viewer's
// `RegionProjector` over its store, its work in `lift.ts`, loaded with the genome chunk (a liftover is over assemblies,
// so a model holding one waits for it; before it lands nothing is carried). Its answers are kept for the store's
// version, so a warm hover across the map reads them (nothing allocated).

import type { RegionProjector } from '@molstar/seqx-core/view';
import { Lift } from './layer-data.js';
import type { LiftCache, RegionReader } from './lift.js';

/** A projector through the maps read by region over `data` (a store: the tiles it holds now, its version). */
export function regionProjector(
  data: RegionReader & { readonly version: { readonly value: number } },
): RegionProjector {
  let cache: LiftCache | undefined;
  return {
    version: data.version,
    reads: (model) => (Lift.module?.regionMaps(model.spec).length ?? 0) > 0,
    project(model, from, ranges, to, out) {
      const m = Lift.module;
      if (!m) return false;
      cache ??= m.newLiftCache();
      return m.projectLifted(cache, data.version.value, model, data, from, ranges, to, out);
    },
    places(model, view, system, start, end) {
      const m = Lift.module;
      if (!m) return undefined;
      cache ??= m.newLiftCache();
      return m.placesLifted(cache, data.version.value, model, data, view, system, start, end);
    },
  };
}
