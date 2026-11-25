/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Size entry, the BED family's runtime: what composing it adds to a viewer (its own code and io's BGZF and tabix; the protocol, issues and tree it
// imports are every viewer's already, so `sizes.mjs` leaves them out of this row). Bundled by `sizes.mjs`, never run.

export { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
