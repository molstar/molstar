/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ShrakeRupleyContext } from './common';
import { isPolymer, MaxAsa, DefaultMaxAsa } from '../../../../mol-model/structure/model/types';

export async function normalizeAccessibleSurfaceArea(ctx: ShrakeRupleyContext) {
    const updateChunk = ctx.updateChunk / 10;
    const residueCount = ctx.relativeAccessibleSurfaceArea.length;
    for (let i = 0; i < residueCount; i += updateChunk) {
        computeRange(ctx, i, Math.min(i + updateChunk, residueCount));
    }
}

function computeRange(ctx: ShrakeRupleyContext, begin: number, end: number) {
    const { accessibleSurfaceArea, relativeAccessibleSurfaceArea, structure } = ctx;
    const { residues, derived } = structure.model.atomicHierarchy;

    for (let i = begin; i < end; ++i) {
        // skip entities not part of a polymer chain
        if (!ctx.nonPolymer) {
            if (!isPolymer(derived.residue.moleculeType[i])) continue;
        }

        const maxAsa = (MaxAsa as any)[residues.label_comp_id.value(i)];
        const rasa = accessibleSurfaceArea[i] / (maxAsa === undefined ? DefaultMaxAsa : maxAsa);
        relativeAccessibleSurfaceArea[i] = rasa;
    }
}