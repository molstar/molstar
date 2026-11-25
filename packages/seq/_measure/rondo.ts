/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Size entry, the genome ring extension: what composing it adds to a viewer (its own code: the ribbons, the hues,
// `rondoRing`, the viewer's commands and binding; the core's extension entry and the viewer it imports are every
// viewer's already, so `sizes.mjs` leaves them out of this row). Bundled by `sizes.mjs`, never run.

export { rondo, rondoRing } from '../ext-rondo/src/index.ts';
export { rondoViewer } from '../ext-rondo/src/viewer.ts';
