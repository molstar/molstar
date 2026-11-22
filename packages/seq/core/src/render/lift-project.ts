/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover and selection through a map read by region (_spec/README.md "Render", A liftover: projection): the viewer's
// `RegionProjector` over its store, its work in `lift.ts`, loaded with the genome chunk (a liftover is over assemblies,
// so a model holding one waits for it; before it lands nothing is carried).

import { withGenome } from '@molstar/seqx-core/resolve';
import type { RegionProjector } from '@molstar/seqx-core/view';
import type { RegionReader } from './lift.js';

const Lift: { project?: typeof import('./lift.js').projectLifted } = {};
withGenome(() =>
  import('./lift.js').then((m) => {
    Lift.project = m.projectLifted;
  }),
);

/** A projector through the maps read by region over `data` (a store: the tiles it holds now, its version). */
export function regionProjector(
  data: RegionReader & { readonly version: { readonly value: number } },
): RegionProjector {
  return {
    get version() {
      return data.version.value;
    },
    project: (model, from, ranges, to, out) => Lift.project?.(model, data, from, ranges, to, out) ?? false,
  };
}
