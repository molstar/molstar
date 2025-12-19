/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Diego del Alamo <diego.delalamo@gmail.com>
 *
 * Structure-level TM-align wrapper
 */

import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';
import { TMAlign } from '../../../../mol-math/linear-algebra/3d/tm-align';
import { StructureElement } from '../element';
import { getPositionTable } from './superposition';

export { tmAlign, tmAlignMultiple };

export type TMAlignResult = TMAlign.Result;

/**
 * Perform TM-align on two structure element loci.
 * Aligns structure B onto structure A (A is the reference).
 *
 * @param a Reference structure loci (will not be transformed)
 * @param b Mobile structure loci (transformation returned)
 * @returns TM-align result with transformation, scores, and alignment
 */
function tmAlign(a: StructureElement.Loci, b: StructureElement.Loci): TMAlignResult {
    const lenA = StructureElement.Loci.size(a);
    const lenB = StructureElement.Loci.size(b);

    if (lenA === 0 || lenB === 0) {
        return {
            bTransform: Mat4.identity(),
            tmScoreA: 0,
            tmScoreB: 0,
            rmsd: 0,
            alignedLength: 0,
            sequenceIdentity: 0,
            alignmentA: [],
            alignmentB: []
        };
    }

    const posA = getPositionTable(a, lenA);
    const posB = getPositionTable(b, lenB);

    return TMAlign.compute({ a: posA, b: posB });
}

/**
 * Perform TM-align on multiple structure element loci.
 * The first structure is used as the reference; all others are aligned to it.
 *
 * @param xs Array of structure element loci (first is reference)
 * @returns Array of TM-align results (length = xs.length - 1)
 */
function tmAlignMultiple(xs: StructureElement.Loci[]): TMAlignResult[] {
    const results: TMAlignResult[] = [];
    if (xs.length < 2) return results;

    const refLoci = xs[0];
    const lenRef = StructureElement.Loci.size(refLoci);
    const posRef = getPositionTable(refLoci, lenRef);

    for (let i = 1; i < xs.length; i++) {
        const mobileLoci = xs[i];
        const lenMobile = StructureElement.Loci.size(mobileLoci);
        const posMobile = getPositionTable(mobileLoci, lenMobile);

        results.push(TMAlign.compute({ a: posRef, b: posMobile }));
    }

    return results;
}
