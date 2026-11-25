/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Size entry, the mmCIF runtime: what composing it adds to a viewer before any source opens (its own code and Seq io's
// converter; the protocol, issues and tree it imports are every viewer's already, so `sizes.mjs` leaves them out of
// this row). `@molstar/io`'s parser is a chunk loaded when the first mmCIF source opens, listed apart by `sizes.mjs`.
// Bundled by `sizes.mjs`, never run.

export { MmcifRuntime } from '@molstar/seqx-core/runtimes/mmcif';
