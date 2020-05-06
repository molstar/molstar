/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MinimizeRmsd } from '../../../../mol-math/linear-algebra/3d/minimize-rmsd';
import StructureElement from '../element';
import { OrderedSet } from '../../../../mol-data/int';
import { AlignSequences } from '../../../sequence/alignment/sequence';
import StructureProperties from '../properties';

export function superpose(xs: StructureElement.Loci[]): MinimizeRmsd.Result[] {
    const ret: MinimizeRmsd.Result[] = [];
    if (xs.length <= 0) return ret;

    const n = getMinSize(xs);
    const input: MinimizeRmsd.Input = {
        a: getPositionTable(xs[0], n),
        b: getPositionTable(xs[1], n)
    };

    ret[0] = MinimizeRmsd.compute(input);
    for (let i = 2; i < xs.length; i++) {
        input.b = getPositionTable(xs[i], n);
        input.centerB = void 0;
        ret.push(MinimizeRmsd.compute(input));
    }

    return ret;
}

type AlignAndSuperposeResult = MinimizeRmsd.Result & { alignmentScore: number };
const reProtein = /(polypeptide|cyclic-pseudo-peptide)/i;

export function alignAndSuperpose(xs: StructureElement.Loci[]): AlignAndSuperposeResult[] {
    const ret: AlignAndSuperposeResult[] = [];
    if (xs.length <= 0) return ret;

    const l = StructureElement.Loci.getFirstLocation(xs[0])!;
    const subtype = StructureProperties.entity.subtype(l);
    const substMatrix = subtype.match(reProtein) ? 'blosum62' : 'default';

    for (let i = 1; i < xs.length; i++) {
        const { a, b, score } = AlignSequences.compute({
            a: xs[0].elements[0],
            b: xs[i].elements[0],
        }, { substMatrix });

        const lociA = StructureElement.Loci(xs[0].structure, [a]);
        const lociB = StructureElement.Loci(xs[i].structure, [b]);
        const n = OrderedSet.size(a.indices);

        ret.push({
            ...MinimizeRmsd.compute({
                a: getPositionTable(lociA, n),
                b: getPositionTable(lociB, n)
            }),
            alignmentScore: score
        });
    }

    return ret;
}

function getPositionTable(xs: StructureElement.Loci, n: number): MinimizeRmsd.Positions {
    const ret = MinimizeRmsd.Positions.empty(n);
    let o = 0;
    for (const u of xs.elements) {
        const { unit, indices } = u;
        const { elements } = unit;
        const { x, y, z } = unit.conformation;
        for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
            const e = elements[OrderedSet.getAt(indices, i)];
            ret.x[o] = x(e);
            ret.y[o] = y(e);
            ret.z[o] = z(e);
            o++;
            if (o >= n) break;
        }
        if (o >= n) break;
    }
    return ret;
}

function getMinSize(xs: StructureElement.Loci[]) {
    if (xs.length === 0) return 0;
    let s = StructureElement.Loci.size(xs[0]);
    for (let i = 1; i < xs.length; i++) {
        const t = StructureElement.Loci.size(xs[i]);
        if (t < s) s = t;
    }
    return s;
}