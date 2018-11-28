/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure'
import { Link } from './structure/structure/unit/links'
import { Shape } from './shape';
import { Sphere3D } from 'mol-math/geometry';
import { CentroidHelper } from 'mol-math/geometry/centroid-helper';
import { Vec3 } from 'mol-math/linear-algebra';
import { OrderedSet } from 'mol-data/int';
import { Structure } from './structure/structure';

/** A Loci that includes every loci */
export const EveryLoci = { kind: 'every-loci' as 'every-loci' }
export type EveryLoci = typeof EveryLoci
export function isEveryLoci(x: any): x is EveryLoci {
    return !!x && x.kind === 'every-loci';
}

/** A Loci that is empty */
export const EmptyLoci = { kind: 'empty-loci' as 'empty-loci' }
export type EmptyLoci = typeof EmptyLoci
export function isEmptyLoci(x: any): x is EmptyLoci {
    return !!x && x.kind === 'empty-loci';
}

export interface DataLoci {
    readonly kind: 'data-loci',
    readonly data: any,
    readonly tag: string
    readonly indices: OrderedSet<number>
}
export function isDataLoci(x: any): x is DataLoci {
    return !!x && x.kind === 'data-loci';
}
export function areDataLociEqual(a: DataLoci, b: DataLoci) {
    return a.data === b.data && a.tag === b.tag && OrderedSet.areEqual(a.indices, b.indices)
}
export function createDataLoci(data: any, tag: string, indices: OrderedSet<number>): DataLoci {
    return { kind: 'data-loci', data, tag, indices }
}

export function areLociEqual(lociA: Loci, lociB: Loci) {
    if (isEveryLoci(lociA) && isEveryLoci(lociB)) return true
    if (isEmptyLoci(lociA) && isEmptyLoci(lociB)) return true
    if (isDataLoci(lociA) && isDataLoci(lociB)) {
        return areDataLociEqual(lociA, lociB)
    }
    if (Structure.isLoci(lociA) && Structure.isLoci(lociB)) {
        return Structure.areLociEqual(lociA, lociB)
    }
    if (StructureElement.isLoci(lociA) && StructureElement.isLoci(lociB)) {
        return StructureElement.areLociEqual(lociA, lociB)
    }
    if (Link.isLoci(lociA) && Link.isLoci(lociB)) {
        return Link.areLociEqual(lociA, lociB)
    }
    if (Shape.isLoci(lociA) && Shape.isLoci(lociB)) {
        return Shape.areLociEqual(lociA, lociB)
    }
    return false
}


export { Loci }

type Loci = StructureElement.Loci | Structure.Loci | Link.Loci | EveryLoci | EmptyLoci | DataLoci | Shape.Loci

namespace Loci {

    const sphereHelper = new CentroidHelper(), tempPos = Vec3.zero();

    export function getBoundingSphere(loci: Loci): Sphere3D | undefined {
        if (loci.kind === 'every-loci' || loci.kind === 'empty-loci') return void 0;

        sphereHelper.reset();
        if (loci.kind === 'structure-loci') {
            return Sphere3D.clone(loci.structure.boundary.sphere)
        } else if (loci.kind === 'element-loci') {
            for (const e of loci.elements) {
                const { indices } = e;
                const pos = e.unit.conformation.position;
                const { elements } = e.unit;
                for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                    pos(elements[OrderedSet.getAt(indices, i)], tempPos);
                    sphereHelper.includeStep(tempPos);
                }
            }
            sphereHelper.finishedIncludeStep();
            for (const e of loci.elements) {
                const { indices } = e;
                const pos = e.unit.conformation.position;
                const { elements } = e.unit;
                for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                    pos(elements[OrderedSet.getAt(indices, i)], tempPos);
                    sphereHelper.radiusStep(tempPos);
                }
            }
        } else if (loci.kind === 'link-loci') {
            for (const e of loci.links) {
                let pos = e.aUnit.conformation.position;
                pos(e.aUnit.elements[e.aIndex], tempPos);
                sphereHelper.includeStep(tempPos);
                pos = e.bUnit.conformation.position;
                pos(e.bUnit.elements[e.bIndex], tempPos);
                sphereHelper.includeStep(tempPos);
            }
            sphereHelper.finishedIncludeStep();
            for (const e of loci.links) {
                let pos = e.aUnit.conformation.position;
                pos(e.aUnit.elements[e.aIndex], tempPos);
                sphereHelper.radiusStep(tempPos);
                pos = e.bUnit.conformation.position;
                pos(e.bUnit.elements[e.bIndex], tempPos);
                sphereHelper.radiusStep(tempPos);
            }
        } else if (loci.kind === 'group-loci') {
            // TODO
            return void 0;
        }

        return Sphere3D.create(Vec3.clone(sphereHelper.center), Math.sqrt(sphereHelper.radiusSq));
    }
}