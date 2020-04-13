/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { DSSPContext, DSSPType } from './common';

/**
 * A minimal helix is defined by two consecutive n-turns.
 * For example, a 4-helix, of minimal length 4 from residues i to i + 3,
 * requires 4-turns at residues i - 1 and i,
 *
 *      3-helix(i,i + 2)=: [3-turn(i - 1) and 3-turn(i)]
 *      4-helix(i,i + 3)=: [4-turn(i - 1) and 4-turn(i)]
 *      5-helix(i,i + 4)=: [5-turn(i - 1) and 5-turn(i)]
 *
 * Type: G (n=3), H (n=4), I (n=5)
 */
export function assignHelices(ctx: DSSPContext) {
    const { proteinInfo, flags } = ctx;
    const residueCount = proteinInfo.residueIndices.length;

    const turnFlag = [DSSPType.Flag.T3S, DSSPType.Flag.T4S, DSSPType.Flag.T5S, DSSPType.Flag.T3, DSSPType.Flag.T4, DSSPType.Flag.T5];
    const helixFlag = [0, 0, 0, DSSPType.Flag.G, DSSPType.Flag.H, DSSPType.Flag.I];

    const helixCheckOrder = ctx.params.oldOrdering ? [4, 3, 5] : [3, 4, 5];
    for (let ni = 0; ni < helixCheckOrder.length; ni++) {
        const n = helixCheckOrder[ni];

        for (let i = 1, il = residueCount - n; i < il; i++) {
            const fI = DSSPType.create(flags[i]);
            const fI1 = DSSPType.create(flags[i - 1]);
            const fI2 = DSSPType.create(flags[i + 1]);

            // TODO rework to elegant solution which will not break instantly
            if (ctx.params.oldOrdering) {
                if ((n === 3 && (DSSPType.is(fI, DSSPType.Flag.H) || DSSPType.is(fI2, DSSPType.Flag.H)) || // for 3-10 yield to alpha helix
                    (n === 5 && ((DSSPType.is(fI, DSSPType.Flag.H) || DSSPType.is(fI, DSSPType.Flag.G)) || (DSSPType.is(fI2, DSSPType.Flag.H) || DSSPType.is(fI2, DSSPType.Flag.G)))))) { // for pi yield to all other helices
                    continue;
                }
            } else {
                if ((n === 4 && (DSSPType.is(fI, DSSPType.Flag.G) || DSSPType.is(fI2, DSSPType.Flag.G)) || // for alpha helix yield to 3-10
                    (n === 5 && ((DSSPType.is(fI, DSSPType.Flag.H) || DSSPType.is(fI, DSSPType.Flag.G)) || (DSSPType.is(fI2, DSSPType.Flag.H) || DSSPType.is(fI2, DSSPType.Flag.G)))))) { // for pi yield to all other helices
                    continue;
                }
            }

            if (DSSPType.is(fI, turnFlag[n]) && DSSPType.is(fI, turnFlag[n - 3]) && // check fI for turn start of proper type
                DSSPType.is(fI1, turnFlag[n]) && DSSPType.is(fI1, turnFlag[n - 3])) { // check fI1 accordingly
                if (ctx.params.oldDefinition) {
                    for (let k = 0; k < n; k++) {
                        flags[i + k] |= helixFlag[n];
                    }
                } else {
                    for (let k = -1; k <= n; k++) {
                        flags[i + k] |= helixFlag[n];
                    }
                }
            }
        }
    }
}