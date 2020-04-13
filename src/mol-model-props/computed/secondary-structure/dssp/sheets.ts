/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { DSSPContext, DSSPType, BridgeType } from './common';

function isHelixType(f: DSSPType) {
    return DSSPType.is(f, DSSPType.Flag.G) || DSSPType.is(f, DSSPType.Flag.H) || DSSPType.is(f, DSSPType.Flag.I);
}

/**
 * sheet=: set of one or more ladders connected by shared residues
 *
 * Type: E
 */
export function assignSheets(ctx: DSSPContext) {
    const { ladders, flags } = ctx;
    for (let ladderIndex = 0; ladderIndex < ladders.length; ladderIndex++) {
        const ladder = ladders[ladderIndex];
        for (let lcount = ladder.firstStart; lcount <= ladder.firstEnd; lcount++) {
            const diff = ladder.firstStart - lcount;
            const l2count = ladder.secondStart - diff;

            if (ladder.firstStart !== ladder.firstEnd) {
                flags[lcount] |= DSSPType.Flag.E;
                flags[l2count] |= DSSPType.Flag.E;
            } else {
                if (!isHelixType(flags[lcount]) && DSSPType.is(flags[lcount], DSSPType.Flag.E)) {
                    flags[lcount] |= DSSPType.Flag.B;
                }
                if (!isHelixType(flags[l2count]) && DSSPType.is(flags[l2count], DSSPType.Flag.E)) {
                    flags[l2count] |= DSSPType.Flag.B;
                }
            }
        }

        if (ladder.nextLadder === 0) continue;

        const conladder = ladders[ladder.nextLadder];
        for (let lcount = ladder.firstStart; lcount <= conladder.firstEnd; lcount++) {
            flags[lcount] |= DSSPType.Flag.E;
        }
        if (ladder.type === BridgeType.PARALLEL) {
            for (let lcount = ladder.secondStart; lcount <= conladder.secondEnd; lcount++) {
                flags[lcount] |= DSSPType.Flag.E;
            }
        } else {
            for (let lcount = conladder.secondEnd; lcount <= ladder.secondStart; lcount++) {
                flags[lcount] |= DSSPType.Flag.E;
            }
        }
    }
}