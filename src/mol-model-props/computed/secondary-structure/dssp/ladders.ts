/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { DSSPContext, Ladder, BridgeType, Bridge } from './common';

function shouldExtendLadder(ladder: Ladder, bridge: Bridge): boolean {
    // in order to extend ladders, same type must be present
    if (bridge.type !== ladder.type) return false;

    // only extend if residue 1 is sequence neighbor with regard to ladder
    if (bridge.partner1 !== ladder.firstEnd + 1) return false;

    if (bridge.type === BridgeType.PARALLEL) {
        if (bridge.partner2 === ladder.secondEnd + 1) {
            return true;
        }
    } else {
        if (bridge.partner2 === ladder.secondStart - 1) {
            return true;
        }
    }

    return false;
}

/**
 * For beta structures, we define: a bulge-linked ladder consists of two ladders or bridges of the same type
 * connected by at most one extra residue of one strand and at most four extra residues  on the other strand,
 * all residues in bulge-linked ladders are marked E, including any extra residues.
 */
function resemblesBulge(ladder1: Ladder, ladder2: Ladder) {
    if (!(ladder1.type === ladder2.type && ladder2.firstStart - ladder1.firstEnd < 6 &&
        ladder1.firstStart < ladder2.firstStart && ladder2.nextLadder === 0)) return false;

    if (ladder1.type === BridgeType.PARALLEL) {
        return bulgeCriterion2(ladder1, ladder2);
    } else {
        return bulgeCriterion2(ladder2, ladder1);
    }
}

function bulgeCriterion2(ladder1: Ladder, ladder2: Ladder) {
    return ladder2.secondStart - ladder1.secondEnd > 0 && ((ladder2.secondStart - ladder1.secondEnd < 6 &&
        ladder2.firstStart - ladder1.firstEnd < 3) || ladder2.secondStart - ladder1.secondEnd < 3);
}

/**
 * ladder=: set of one or more consecutive bridges of identical type
 *
 * Type: E
 */
export function assignLadders(ctx: DSSPContext) {
    const { bridges, ladders } = ctx;

    // create ladders
    for (let bridgeIndex = 0; bridgeIndex < bridges.length; bridgeIndex++) {
        const bridge = bridges[bridgeIndex];
        let found = false;
        for (let ladderIndex = 0; ladderIndex < ladders.length; ladderIndex++) {
            const ladder = ladders[ladderIndex];
            if (shouldExtendLadder(ladder, bridge)) {
                found = true;
                ladder.firstEnd++;
                if (bridge.type === BridgeType.PARALLEL) {
                    ladder.secondEnd++;
                } else {
                    ladder.secondStart--;
                }
            }
        }

        // no suitable assignment: create new ladder with single bridge as content
        if (!found) {
            ladders[ladders.length] = {
                previousLadder: 0,
                nextLadder: 0,
                firstStart: bridge.partner1,
                firstEnd: bridge.partner1,
                secondStart: bridge.partner2,
                secondEnd: bridge.partner2,
                type: bridge.type
            };
        }
    }

    // connect ladders
    for (let ladderIndex1 = 0; ladderIndex1 < ladders.length; ladderIndex1++) {
        const ladder1 = ladders[ladderIndex1];
        for (let ladderIndex2 = ladderIndex1; ladderIndex2 < ladders.length; ladderIndex2++) {
            const ladder2 = ladders[ladderIndex2];
            if (resemblesBulge(ladder1, ladder2)) {
                ladder1.nextLadder = ladderIndex2;
                ladder2.previousLadder = ladderIndex1;
            }
        }
    }
}