/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { DSSPContext, DSSPType } from './common';


/**
 * The basic turn pattern is a single H bond of type (i, i + n).
 * We assign an n-turn at residue i if there is an H bond from CO(i) to NH(i + n),
 * i.e., “n-turn(i)=: Hbond(i, i + n), n = 3, 4, 5.”
 *
 * Type: T
 */
export function assignTurns(ctx: DSSPContext) {
    const { proteinInfo, hbonds, flags } = ctx;

    const turnFlag = [DSSPType.Flag.T3S, DSSPType.Flag.T4S, DSSPType.Flag.T5S, DSSPType.Flag.T3, DSSPType.Flag.T4, DSSPType.Flag.T5];

    for (let idx = 0; idx < 3; idx++) {
        for (let i = 0, il = proteinInfo.residueIndices.length - 1; i < il; ++i) {

            // check if hbond exists
            if (hbonds.getDirectedEdgeIndex(i, i + idx + 3) !== -1) {
                flags[i] |= turnFlag[idx + 3] | turnFlag[idx];
                if (ctx.params.oldDefinition) {
                    for (let k = 1; k < idx + 3; ++k) {
                        flags[i + k] |= turnFlag[idx + 3] | DSSPType.Flag.T;
                    }
                } else {
                    for (let k = 0; k <= idx + 3; ++k) {
                        flags[i + k] |= turnFlag[idx + 3] | DSSPType.Flag.T;
                    }
                }
            }
        }
    }
}