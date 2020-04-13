/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { DSSPContext, DSSPType, BridgeType, Bridge } from './common';

/**
 * Two nonoverlapping stretches of three residues each, i - 1, i, i + 1 and j - 1, j, j + 1,
 * form either a parallel or antiparallel bridge, depending on which of
 * two basic patterns is matched. We assign a bridge between residues i and j
 * if there are two H bonds characteristic of P-structure; in particular,
 *
 * Parallel Bridge(i, j) =:
 *      [Hbond(i - 1, j) and Hbond(j, i + 1)] or
 *      [Hbond(j - 1, i) and Hbond(i, j + 1)]
 *
 * Antiparallel Bridge(i, j) =:
 *      [Hbond(i, j) and Hbond(j, i)] or
 *      [Hbond(i - 1, j + 1) and Hbond(j - 1, i + l)]
 *
 * Type: B
 */
export function assignBridges(ctx: DSSPContext) {
    const { proteinInfo, hbonds, flags, bridges } = ctx;

    const { offset, b } = hbonds;
    let i: number, j: number;

    for (let k = 0, kl = proteinInfo.residueIndices.length; k < kl; ++k) {
        for (let t = offset[k], _t = offset[k + 1]; t < _t; t++) {
            const l = b[t];
            if (k > l) continue;

            // Parallel Bridge(i, j) =: [Hbond(i - 1, j) and Hbond(j, i + 1)]
            i = k + 1; // k is i - 1
            j = l;
            if (i !== j && hbonds.getDirectedEdgeIndex(j, i + 1) !== -1) {
                flags[i] |= DSSPType.Flag.B;
                flags[j] |= DSSPType.Flag.B;
                // TODO move to constructor, actually omit object all together
                bridges[bridges.length] = new Bridge(i, j, BridgeType.PARALLEL);
            }

            // Parallel Bridge(i, j) =: [Hbond(j - 1, i) and Hbond(i, j + 1)]
            i = k;
            j = l - 1; // l is j + 1
            if (i !== j && hbonds.getDirectedEdgeIndex(j - 1, i) !== -1) {
                flags[i] |= DSSPType.Flag.B;
                flags[j] |= DSSPType.Flag.B;
                bridges[bridges.length] = new Bridge(j, i, BridgeType.PARALLEL);
            }

            // Antiparallel Bridge(i, j) =: [Hbond(i, j) and Hbond(j, i)]
            i = k;
            j = l;
            if (i !== j && hbonds.getDirectedEdgeIndex(j, i) !== -1) {
                flags[i] |= DSSPType.Flag.B;
                flags[j] |= DSSPType.Flag.B;
                bridges[bridges.length] = new Bridge(j, i, BridgeType.ANTI_PARALLEL);
            }

            // Antiparallel Bridge(i, j) =: [Hbond(i - 1, j + 1) and Hbond(j - 1, i + l)]
            i = k + 1;
            j = l - 1;
            if (i !== j && hbonds.getDirectedEdgeIndex(j - 1, i + 1) !== -1) {
                flags[i] |= DSSPType.Flag.B;
                flags[j] |= DSSPType.Flag.B;
                bridges[bridges.length] = new Bridge(j, i, BridgeType.ANTI_PARALLEL);
            }
        }
    }

    bridges.sort((a, b) => a.partner1 > b.partner1 ? 1 : a.partner1 < b.partner1 ? -1 : 0);
}