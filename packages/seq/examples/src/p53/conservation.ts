/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// p53's conservation signal over the gene's region (`RegionSystem`, `RegionLength`): a bigWig read by region and zoom
// level through the `bigwig` runtime (`@molstar/seqx-core/runtimes/bigwig`, opt in). Generated from a seed by
// `seqx-core`'s test-only bigWig writer, which compares the committed file with what it writes; never recorded.

/** The committed bigWig, as a spec names it. */
export const ConservationFile = 'fixtures/p53/gene-region-conservation.bw';
/** Its format: the application composes the `bigwig` runtime to read it. */
export const ConservationFormat = 'bigwig';
/** Its zoom levels' reductions (the levels it serves besides 1). */
export const ConservationZooms = [16, 256] as const;
/** The bases at its end without data (past the alignment the signal comes from). */
export const ConservationUncovered = 200;
