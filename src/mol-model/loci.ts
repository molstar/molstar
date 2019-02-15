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

/** A generic data loci */
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

export { Loci }

type Loci = StructureElement.Loci | Structure.Loci | Link.Loci | EveryLoci | EmptyLoci | DataLoci | Shape.Loci

namespace Loci {
    export function areEqual(lociA: Loci, lociB: Loci) {
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

    const sphereHelper = new CentroidHelper(), tempPos = Vec3.zero();

    export function getBoundingSphere(loci: Loci, boundingSphere?: Sphere3D): Sphere3D | undefined {
        if (loci.kind === 'every-loci' || loci.kind === 'empty-loci') return void 0;

        if (!boundingSphere) boundingSphere = Sphere3D.zero()
        sphereHelper.reset();

        if (loci.kind === 'structure-loci') {
            return Sphere3D.copy(boundingSphere, loci.structure.boundary.sphere)
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
                e.aUnit.conformation.position(e.aUnit.elements[e.aIndex], tempPos);
                sphereHelper.includeStep(tempPos);
                e.bUnit.conformation.position(e.bUnit.elements[e.bIndex], tempPos);
                sphereHelper.includeStep(tempPos);
            }
            sphereHelper.finishedIncludeStep();
            for (const e of loci.links) {
                e.aUnit.conformation.position(e.aUnit.elements[e.aIndex], tempPos);
                sphereHelper.radiusStep(tempPos);
                e.aUnit.conformation.position(e.bUnit.elements[e.bIndex], tempPos);
                sphereHelper.radiusStep(tempPos);
            }
        } else if (loci.kind === 'group-loci') {
            // TODO
            return void 0;
        } else if (loci.kind === 'data-loci') {
            // TODO maybe add loci.getBoundingSphere()???
            return void 0;
        }

        Vec3.copy(boundingSphere.center, sphereHelper.center)
        boundingSphere.radius = Math.sqrt(sphereHelper.radiusSq)
        return boundingSphere
    }

    const tmpSphere3D = Sphere3D.zero()
    export function getCenter(loci: Loci, center?: Vec3): Vec3 | undefined {
        const boundingSphere = getBoundingSphere(loci, tmpSphere3D)
        return boundingSphere ? Vec3.copy(center || Vec3.zero(), boundingSphere.center) : undefined
    }
}