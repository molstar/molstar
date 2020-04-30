/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure';
import { Bond } from './structure/structure/unit/bonds';
import { Shape, ShapeGroup } from './shape';
import { Sphere3D } from '../mol-math/geometry';
import { Vec3 } from '../mol-math/linear-algebra';
import { Structure } from './structure/structure';
import { PrincipalAxes } from '../mol-math/linear-algebra/matrix/principal-axes';
import { ParamDefinition } from '../mol-util/param-definition';
import { shallowEqual } from '../mol-util';
import { FiniteArray } from '../mol-util/type-helpers';
import { BoundaryHelper } from '../mol-math/geometry/boundary-helper';
import { stringToWords } from '../mol-util/string';
import { Volume } from './volume/volume';

/** A Loci that includes every loci */
export const EveryLoci = { kind: 'every-loci' as 'every-loci' };
export type EveryLoci = typeof EveryLoci
export function isEveryLoci(x?: Loci): x is EveryLoci {
    return !!x && x.kind === 'every-loci';
}

/** A Loci that is empty */
export const EmptyLoci = { kind: 'empty-loci' as 'empty-loci' };
export type EmptyLoci = typeof EmptyLoci
export function isEmptyLoci(x?: Loci): x is EmptyLoci {
    return !!x && x.kind === 'empty-loci';
}

/** A generic data loci */
export interface DataLoci<T = unknown, E = unknown> {
    readonly kind: 'data-loci',
    readonly tag: string
    readonly data: T,
    readonly elements: ReadonlyArray<E>

    getBoundingSphere(boundingSphere: Sphere3D): Sphere3D
    getLabel(): string
}
export function isDataLoci(x?: Loci): x is DataLoci {
    return !!x && x.kind === 'data-loci';
}
export function areDataLociEqual(a: DataLoci, b: DataLoci) {
    // use shallowEqual to allow simple data objects that are contructed on-the-fly
    if (!shallowEqual(a.data, b.data) || a.tag !== b.tag) return false;
    if (a.elements.length !== b.elements.length) return false;
    for (let i = 0, il = a.elements.length; i < il; ++i) {
        if (!shallowEqual(a.elements[i], b.elements[i])) return false;
    }
    return true;
}
export function isDataLociEmpty(loci: DataLoci) {
    return loci.elements.length === 0 ? true : false;
}
export function DataLoci<T = unknown, E = unknown>(tag: string, data: T, elements: ReadonlyArray<E>, getBoundingSphere: DataLoci<T, E>['getBoundingSphere'], getLabel: DataLoci<T, E>['getLabel']): DataLoci<T, E> {
    return { kind: 'data-loci', tag, data, elements, getBoundingSphere, getLabel };
}

export { Loci };

type Loci = StructureElement.Loci | Structure.Loci | Bond.Loci | EveryLoci | EmptyLoci | DataLoci | Shape.Loci | ShapeGroup.Loci | Volume.Loci | Volume.Isosurface.Loci | Volume.Cell.Loci

namespace Loci {
    export interface Bundle<L extends number> { loci: FiniteArray<Loci, L> }

    const boundaryHelper = new BoundaryHelper('98');
    export function getBundleBoundingSphere(bundle: Bundle<any>): Sphere3D {
        const spheres = bundle.loci.map(l => getBoundingSphere(l)).filter(s => !!s) as Sphere3D[];
        boundaryHelper.reset();
        for (const s of spheres) boundaryHelper.includePositionRadius(s.center, s.radius);
        boundaryHelper.finishedIncludeStep();
        for (const s of spheres) boundaryHelper.radiusPositionRadius(s.center, s.radius);
        return boundaryHelper.getSphere();
    }

    export function areEqual(lociA: Loci, lociB: Loci) {
        if (isEveryLoci(lociA) && isEveryLoci(lociB)) return true;
        if (isEmptyLoci(lociA) && isEmptyLoci(lociB)) return true;
        if (isDataLoci(lociA) && isDataLoci(lociB)) {
            return areDataLociEqual(lociA, lociB);
        }
        if (Structure.isLoci(lociA) && Structure.isLoci(lociB)) {
            return Structure.areLociEqual(lociA, lociB);
        }
        if (StructureElement.Loci.is(lociA) && StructureElement.Loci.is(lociB)) {
            return StructureElement.Loci.areEqual(lociA, lociB);
        }
        if (Bond.isLoci(lociA) && Bond.isLoci(lociB)) {
            return Bond.areLociEqual(lociA, lociB);
        }
        if (Shape.isLoci(lociA) && Shape.isLoci(lociB)) {
            return Shape.areLociEqual(lociA, lociB);
        }
        if (ShapeGroup.isLoci(lociA) && ShapeGroup.isLoci(lociB)) {
            return ShapeGroup.areLociEqual(lociA, lociB);
        }
        if (Volume.isLoci(lociA) && Volume.isLoci(lociB)) {
            return Volume.areLociEqual(lociA, lociB);
        }
        if (Volume.Isosurface.isLoci(lociA) && Volume.Isosurface.isLoci(lociB)) {
            return Volume.Isosurface.areLociEqual(lociA, lociB);
        }
        if (Volume.Cell.isLoci(lociA) && Volume.Cell.isLoci(lociB)) {
            return Volume.Cell.areLociEqual(lociA, lociB);
        }
        return false;
    }

    export function isEvery(loci?: Loci): loci is EveryLoci {
        return !!loci && loci.kind === 'every-loci';
    }

    export function isEmpty(loci: Loci): loci is EmptyLoci {
        if (isEveryLoci(loci)) return false;
        if (isEmptyLoci(loci)) return true;
        if (isDataLoci(loci)) return isDataLociEmpty(loci);
        if (Structure.isLoci(loci)) return Structure.isLociEmpty(loci);
        if (StructureElement.Loci.is(loci)) return StructureElement.Loci.isEmpty(loci);
        if (Bond.isLoci(loci)) return Bond.isLociEmpty(loci);
        if (Shape.isLoci(loci)) return Shape.isLociEmpty(loci);
        if (ShapeGroup.isLoci(loci)) return ShapeGroup.isLociEmpty(loci);
        if (Volume.isLoci(loci)) return Volume.isLociEmpty(loci);
        if (Volume.Isosurface.isLoci(loci)) return Volume.Isosurface.isLociEmpty(loci);
        if (Volume.Cell.isLoci(loci)) return Volume.Cell.isLociEmpty(loci);
        return false;
    }

    export function remap<T>(loci: Loci, data: T) {
        if (data instanceof Structure) {
            if (StructureElement.Loci.is(loci)) {
                loci = StructureElement.Loci.remap(loci, data);
            } else if (Structure.isLoci(loci)) {
                loci = Structure.remapLoci(loci, data);
            } else if (Bond.isLoci(loci)) {
                loci = Bond.remapLoci(loci, data);
            }
        }
        return loci;
    }

    export function getBoundingSphere(loci: Loci, boundingSphere?: Sphere3D): Sphere3D | undefined {
        if (loci.kind === 'every-loci' || loci.kind === 'empty-loci') return void 0;

        if (!boundingSphere) boundingSphere = Sphere3D();

        if (loci.kind === 'structure-loci') {
            return Sphere3D.copy(boundingSphere, loci.structure.boundary.sphere);
        } else if (loci.kind === 'element-loci') {
            return Sphere3D.copy(boundingSphere, StructureElement.Loci.getBoundary(loci).sphere);
        } else if (loci.kind === 'bond-loci') {
            return Bond.getBoundingSphere(loci, boundingSphere);
        } else if (loci.kind === 'shape-loci') {
            return Sphere3D.copy(boundingSphere, loci.shape.geometry.boundingSphere);
        } else if (loci.kind === 'group-loci') {
            return ShapeGroup.getBoundingSphere(loci, boundingSphere);
        } else if (loci.kind === 'data-loci') {
            return loci.getBoundingSphere(boundingSphere);
        } else if (loci.kind === 'volume-loci') {
            return Volume.getBoundingSphere(loci.volume, boundingSphere);
        } else if (loci.kind === 'isosurface-loci') {
            return Volume.Isosurface.getBoundingSphere(loci.volume, loci.isoValue, boundingSphere);
        } else if (loci.kind === 'cell-loci') {
            return Volume.Cell.getBoundingSphere(loci.volume, loci.indices, boundingSphere);
        }
    }

    const tmpSphere3D = Sphere3D.zero();
    export function getCenter(loci: Loci, center?: Vec3): Vec3 | undefined {
        const boundingSphere = getBoundingSphere(loci, tmpSphere3D);
        return boundingSphere ? Vec3.copy(center || Vec3(), boundingSphere.center) : undefined;
    }

    export function getPrincipalAxes(loci: Loci): PrincipalAxes | undefined {
        if (loci.kind === 'every-loci' || loci.kind === 'empty-loci') return void 0;

        if (loci.kind === 'structure-loci') {
            return StructureElement.Loci.getPrincipalAxes(Structure.toStructureElementLoci(loci.structure));
        } else if (loci.kind === 'element-loci') {
            return StructureElement.Loci.getPrincipalAxes(loci);
        } else if (loci.kind === 'bond-loci') {
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
        } else if (loci.kind === 'volume-loci') {
            // TODO
            return void 0;
        } else if (loci.kind === 'isosurface-loci') {
            // TODO
            return void 0;
        } else if (loci.kind === 'cell-loci') {
            // TODO
            return void 0;
        }
    }

    //

    const Granularity = {
        'element': (loci: Loci) => loci,
        'residue': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeResidues(loci, true)
                : loci;
        },
        'chain': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeChains(loci)
                : loci;
        },
        'entity': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeEntities(loci)
                : loci;
        },
        'model': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToWholeModels(loci)
                : loci;
        },
        'structure': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? Structure.toStructureElementLoci(loci.structure)
                : ShapeGroup.isLoci(loci)
                    ? Shape.Loci(loci.shape)
                    : loci;
        },
        'elementInstances': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToAllInstances(loci)
                : loci;
        },
        'residueInstances': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToAllInstances(StructureElement.Loci.extendToWholeResidues(loci, true))
                : loci;
        },
        'chainInstances': (loci: Loci) => {
            return StructureElement.Loci.is(loci)
                ? StructureElement.Loci.extendToAllInstances(StructureElement.Loci.extendToWholeChains(loci))
                : loci;
        },
    };
    export type Granularity = keyof typeof Granularity
    export const GranularityOptions = ParamDefinition.objectToOptions(Granularity, k => {
        switch (k) {
            case 'element': return'Atom/Coarse Element';
            case 'elementInstances': return ['Atom/Coarse Element Instances', 'With Symmetry'];
            case 'structure': return'Structure/Shape';
            default: return k.indexOf('Instances')
                ? [stringToWords(k), 'With Symmetry'] : stringToWords(k);
        }
    });

    export function applyGranularity(loci: Loci, granularity: Granularity) {
        return Granularity[granularity](loci);
    }

    /**
     * Converts structure related loci to StructureElement.Loci and applies
     * granularity if given
     */
    export function normalize(loci: Loci, granularity?: Granularity) {
        if (granularity !== 'element' && Bond.isLoci(loci)) {
            // convert Bond.Loci to a StructureElement.Loci so granularity can be applied
            loci = Bond.toStructureElementLoci(loci);
        }
        if (Structure.isLoci(loci)) {
            // convert to StructureElement.Loci
            loci = Structure.toStructureElementLoci(loci.structure);
        }
        if (StructureElement.Loci.is(loci)) {
            // ensure the root structure is used
            loci = StructureElement.Loci.remap(loci, loci.structure.root);
        }
        if (granularity) {
            // needs to be applied AFTER remapping to root
            loci = applyGranularity(loci, granularity);
        }
        return loci;
    }
}