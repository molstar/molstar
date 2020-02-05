/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure'
import { Bond } from './structure/structure/unit/bonds'
import { Shape, ShapeGroup } from './shape';
import { Sphere3D } from '../mol-math/geometry';
import { CentroidHelper } from '../mol-math/geometry/centroid-helper';
import { Vec3 } from '../mol-math/linear-algebra';
import { OrderedSet } from '../mol-data/int';
import { Structure } from './structure/structure';
import { PrincipalAxes } from '../mol-math/linear-algebra/matrix/principal-axes';
import { ParamDefinition } from '../mol-util/param-definition';
import { Interactions } from '../mol-model-props/computed/interactions/interactions';
import { Features } from '../mol-model-props/computed/interactions/features';

/** A Loci that includes every loci */
export const EveryLoci = { kind: 'every-loci' as 'every-loci' }
export type EveryLoci = typeof EveryLoci
export function isEveryLoci(x?: Loci): x is EveryLoci {
    return !!x && x.kind === 'every-loci';
}

/** A Loci that is empty */
export const EmptyLoci = { kind: 'empty-loci' as 'empty-loci' }
export type EmptyLoci = typeof EmptyLoci
export function isEmptyLoci(x?: Loci): x is EmptyLoci {
    return !!x && x.kind === 'empty-loci';
}

/** A generic data loci */
export interface DataLoci<T = unknown> {
    readonly kind: 'data-loci',
    readonly data: T,
    readonly tag: string
    readonly indices: OrderedSet<number>
}
export function isDataLoci(x?: Loci): x is DataLoci {
    return !!x && x.kind === 'data-loci';
}
export function areDataLociEqual(a: DataLoci, b: DataLoci) {
    return a.data === b.data && a.tag === b.tag && OrderedSet.areEqual(a.indices, b.indices)
}
export function isDataLociEmpty(loci: DataLoci) {
    return OrderedSet.size(loci.indices) === 0 ? true : false
}
export function createDataLoci<T = unknown>(data: T, tag: string, indices: OrderedSet<number>): DataLoci<T> {
    return { kind: 'data-loci', data, tag, indices }
}

export { Loci }

type Loci = StructureElement.Loci | Structure.Loci | Bond.Loci | Interactions.Loci | EveryLoci | EmptyLoci | DataLoci | Shape.Loci | ShapeGroup.Loci

namespace Loci {
    interface FiniteArray<T, L extends number = number> extends ReadonlyArray<T> { length: L };
    export interface Bundle<L extends number> { loci: FiniteArray<Loci, L> }

    export function areEqual(lociA: Loci, lociB: Loci) {
        if (isEveryLoci(lociA) && isEveryLoci(lociB)) return true
        if (isEmptyLoci(lociA) && isEmptyLoci(lociB)) return true
        if (isDataLoci(lociA) && isDataLoci(lociB)) {
            return areDataLociEqual(lociA, lociB)
        }
        if (Structure.isLoci(lociA) && Structure.isLoci(lociB)) {
            return Structure.areLociEqual(lociA, lociB)
        }
        if (StructureElement.Loci.is(lociA) && StructureElement.Loci.is(lociB)) {
            return StructureElement.Loci.areEqual(lociA, lociB)
        }
        if (Bond.isLoci(lociA) && Bond.isLoci(lociB)) {
            return Bond.areLociEqual(lociA, lociB)
        }
        if (Interactions.isLoci(lociA) && Interactions.isLoci(lociB)) {
            return Interactions.areLociEqual(lociA, lociB)
        }
        if (Shape.isLoci(lociA) && Shape.isLoci(lociB)) {
            return Shape.areLociEqual(lociA, lociB)
        }
        if (ShapeGroup.isLoci(lociA) && ShapeGroup.isLoci(lociB)) {
            return ShapeGroup.areLociEqual(lociA, lociB)
        }
        return false
    }

    export function isEvery(loci?: Loci): loci is EveryLoci {
        return !!loci && loci.kind === 'every-loci';
    }

    export function isEmpty(loci: Loci): loci is EmptyLoci {
        if (isEveryLoci(loci)) return false
        if (isEmptyLoci(loci)) return true
        if (isDataLoci(loci)) return isDataLociEmpty(loci)
        if (Structure.isLoci(loci)) return Structure.isLociEmpty(loci)
        if (StructureElement.Loci.is(loci)) return StructureElement.Loci.isEmpty(loci)
        if (Bond.isLoci(loci)) return Bond.isLociEmpty(loci)
        if (Interactions.isLoci(loci)) return Interactions.isLociEmpty(loci)
        if (Shape.isLoci(loci)) return Shape.isLociEmpty(loci)
        if (ShapeGroup.isLoci(loci)) return ShapeGroup.isLociEmpty(loci)
        return false
    }

    export function remap<T>(loci: Loci, data: T) {
        if (data instanceof Structure) {
            if (StructureElement.Loci.is(loci)) {
                loci = StructureElement.Loci.remap(loci, data)
            } else if (Structure.isLoci(loci)) {
                loci = Structure.remapLoci(loci, data)
            } else if (Bond.isLoci(loci)) {
                loci = Bond.remapLoci(loci, data)
            } else if (Interactions.isLoci(loci)) {
                // TODO might be too expensive
                // loci = Interactions.remapLoci(loci, data)
            }
        }
        return loci
    }

    const sphereHelper = new CentroidHelper(), tmpPos = Vec3.zero();

    export function getBoundingSphere(loci: Loci, boundingSphere?: Sphere3D): Sphere3D | undefined {
        if (loci.kind === 'every-loci' || loci.kind === 'empty-loci') return void 0;

        if (!boundingSphere) boundingSphere = Sphere3D()
        sphereHelper.reset();

        if (loci.kind === 'structure-loci') {
            return Sphere3D.copy(boundingSphere, loci.structure.boundary.sphere)
        } else if (loci.kind === 'element-loci') {
            return Sphere3D.copy(boundingSphere, StructureElement.Loci.getBoundary(loci).sphere);
        } else if (loci.kind === 'bond-loci') {
            for (const e of loci.bonds) {
                e.aUnit.conformation.position(e.aUnit.elements[e.aIndex], tmpPos);
                sphereHelper.includeStep(tmpPos);
                e.bUnit.conformation.position(e.bUnit.elements[e.bIndex], tmpPos);
                sphereHelper.includeStep(tmpPos);
            }
            sphereHelper.finishedIncludeStep();
            for (const e of loci.bonds) {
                e.aUnit.conformation.position(e.aUnit.elements[e.aIndex], tmpPos);
                sphereHelper.radiusStep(tmpPos);
                e.aUnit.conformation.position(e.bUnit.elements[e.bIndex], tmpPos);
                sphereHelper.radiusStep(tmpPos);
            }
        } else if (loci.kind === 'interaction-loci') {
            const { unitsFeatures } = loci.interactions
            for (const e of loci.contacts) {
                Features.setPosition(tmpPos, e.unitA, e.indexA, unitsFeatures.get(e.unitA.id))
                sphereHelper.includeStep(tmpPos)
                Features.setPosition(tmpPos, e.unitB, e.indexB, unitsFeatures.get(e.unitB.id))
                sphereHelper.includeStep(tmpPos);
            }
            sphereHelper.finishedIncludeStep();
            for (const e of loci.contacts) {
                Features.setPosition(tmpPos, e.unitA, e.indexA, unitsFeatures.get(e.unitA.id))
                sphereHelper.radiusStep(tmpPos)
                Features.setPosition(tmpPos, e.unitB, e.indexB, unitsFeatures.get(e.unitB.id))
                sphereHelper.radiusStep(tmpPos);
            }
        } else if (loci.kind === 'shape-loci') {
            return Sphere3D.copy(boundingSphere, loci.shape.geometry.boundingSphere)
        } else if (loci.kind === 'group-loci') {
            return ShapeGroup.getBoundingSphere(loci, boundingSphere)
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

    export function getPrincipalAxes(loci: Loci): PrincipalAxes | undefined {
        if (loci.kind === 'every-loci' || loci.kind === 'empty-loci') return void 0;

        if (loci.kind === 'structure-loci') {
            return StructureElement.Loci.getPrincipalAxes(Structure.toStructureElementLoci(loci.structure))
        } else if (loci.kind === 'element-loci') {
            return StructureElement.Loci.getPrincipalAxes(loci)
        } else if (loci.kind === 'bond-loci') {
            // TODO
            return void 0;
        } else if (loci.kind === 'interaction-loci') {
            // TODO
            return void 0;
        } else if (loci.kind === 'shape-loci') {
            // TODO
            return void 0;
        } else if (loci.kind === 'group-loci') {
            // TODO
            return void 0;
        } else if (loci.kind === 'data-loci') {
            // TODO maybe add loci.getPrincipalAxes()???
            return void 0;
        }
    }

    //

    const Granularity = {
        'element': (loci: Loci) => loci,
        'residue': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeResidues(loci, true)
                : loci
        },
        'chain': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeChains(loci)
                : loci
        },
        'elementInstances': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToAllInstances(loci)
                : loci
        },
        'residueInstances': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToAllInstances(StructureElement.Loci.extendToWholeResidues(loci, true))
                : loci
        },
        'chainInstances': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToAllInstances(StructureElement.Loci.extendToWholeChains(loci))
                : loci
        },
        'entity': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeEntities(loci)
                : loci
        },
        'model': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeModels(loci)
                : loci
        },
        'structure': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? Structure.toStructureElementLoci(loci.structure)
                : loci
        },
        'shape': (loci: Loci) => {
            return ShapeGroup.isLoci(loci)
                ? Shape.Loci(loci.shape)
                : loci
        },
    }
    export type Granularity = keyof typeof Granularity
    export const GranularityOptions = ParamDefinition.objectToOptions(Granularity);

    export function applyGranularity(loci: Loci, granularity: Granularity) {
        return Granularity[granularity](loci)
    }

    /**
     * Converts structure related loci to StructureElement.Loci and applies
     * granularity if given
    */
    export function normalize(loci: Loci, granularity?: Granularity) {
        if (granularity !== 'element' && Bond.isLoci(loci)) {
            // convert Bond.Loci to a StructureElement.Loci so granularity can be applied
            loci = Bond.toStructureElementLoci(loci)
        }
        if (Structure.isLoci(loci)) {
            // convert to StructureElement.Loci
            loci = Structure.toStructureElementLoci(loci.structure)
        }
        if (StructureElement.Loci.is(loci)) {
            // ensure the root structure is used
            loci = StructureElement.Loci.remap(loci, loci.structure.root)
        }
        if (granularity) {
            // needs to be applied AFTER remapping to root
            loci = applyGranularity(loci, granularity)
        }
        return loci
    }
}