/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';
import { createRangeArray, makeBuckets } from 'mol-data/util';
import { Column } from 'mol-data/db';
import { RuntimeContext } from 'mol-task';

export async function sortAtomSite(ctx: RuntimeContext, atom_site: mmCIF_Database['atom_site'], start: number, end: number) {
    const indices = createRangeArray(start, end - 1);

    const { label_entity_id, label_asym_id, label_seq_id } = atom_site;
    const entityBuckets = makeBuckets(indices, label_entity_id.value, false);
    if (ctx.shouldUpdate) await ctx.update();
    for (let ei = 0, _eI = entityBuckets.length - 1; ei < _eI; ei++) {
        const chainBuckets = makeBuckets(indices, label_asym_id.value, false, entityBuckets[ei], entityBuckets[ei + 1]);
        for (let cI = 0, _cI = chainBuckets.length - 1; cI < _cI; cI++) {
            const aI = chainBuckets[cI];
            // are we in HETATM territory?
            if (label_seq_id.valueKind(aI) !== Column.ValueKind.Present) continue;

            makeBuckets(indices, label_seq_id.value, true, aI, chainBuckets[cI + 1]);
            if (ctx.shouldUpdate) await ctx.update();
        }
        if (ctx.shouldUpdate) await ctx.update();
    }

    return indices;
}