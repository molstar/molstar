/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, SortedArray } from 'mol-data/int'
import Unit from './unit'
import { ElementIndex } from '../model';
import { ResidueIndex, ChainIndex } from '../model/indexing';
import Structure from './structure';
import { Boundary } from './util/boundary';
import { BoundaryHelper } from 'mol-math/geometry/boundary-helper';
import { Vec3 } from 'mol-math/linear-algebra';

interface StructureElement<U = Unit> {
    readonly kind: 'element-location',
    unit: U,
    /** Index into element (atomic/coarse) properties of unit.model */
    element: ElementIndex
}

namespace StructureElement {
    export function create(unit?: Unit, element?: ElementIndex): StructureElement {
        return { kind: 'element-location', unit: unit!, element: element || (0 as ElementIndex) };
    }

    // TODO: when nominal types are available, make this indexed by UnitIndex
    export type Set = SortedArray<ElementIndex>

    /** Index into Unit.elements */
    export type UnitIndex = { readonly '@type': 'structure-element-index' } & number

    export interface Property<T> { (location: StructureElement): T }
    export interface Predicate extends Property<boolean> { }

    export function property<T>(p: Property<T>) { return p; }

    function _wrongUnitKind(kind: string) { throw new Error(`Property only available for ${kind} models.`); }
    export function atomicProperty<T>(p: (location: StructureElement<Unit.Atomic>) => T) {
        return property(l => Unit.isAtomic(l.unit) ? p(l as StructureElement<Unit.Atomic>) : _wrongUnitKind('atomic') );
    }

    export function coarseProperty<T>(p: (location: StructureElement<Unit.Spheres | Unit.Gaussians>) => T) {
        return property(l => Unit.isCoarse(l.unit) ? p(l as StructureElement<Unit.Spheres | Unit.Gaussians>) : _wrongUnitKind('coarse') );
    }

    /** Represents multiple element index locations */
    export interface Loci {
        readonly kind: 'element-loci',
        readonly structure: Structure,
        /** Access i-th element as unit.elements[indices[i]] */
        readonly elements: ReadonlyArray<{
            unit: Unit,
            /**
             * Indices into the unit.elements array.
             * Can use OrderedSet.forEach to iterate (or OrderedSet.size + OrderedSet.getAt)
             */
            indices: OrderedSet<UnitIndex>
        }>
    }

    export function Loci(structure: Structure, elements: ArrayLike<{ unit: Unit, indices: OrderedSet<UnitIndex> }>): Loci {
        return { kind: 'element-loci', structure, elements: elements as Loci['elements'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'element-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.elements.length !== b.elements.length) return false
        for (let i = 0, il = a.elements.length; i < il; ++i) {
            const elementA = a.elements[i]
            const elementB = b.elements[i]
            if (elementA.unit.id !== elementB.unit.id) return false
            if (!OrderedSet.areEqual(elementA.indices, elementB.indices)) return false
        }
        return true
    }

    export function isLocation(x: any): x is StructureElement {
        return !!x && x.kind === 'element-location';
    }

    export function residueIndex(e: StructureElement) {
        if (Unit.isAtomic(e.unit)) {
            return e.unit.residueIndex[e.element];
        } else {
            // TODO: throw error instead?
            return -1 as ResidueIndex;
        }
    }

    export function chainIndex(e: StructureElement) {
        if (Unit.isAtomic(e.unit)) {
            return e.unit.chainIndex[e.element];
        } else {
            // TODO: throw error instead?
            return -1 as ChainIndex;
        }
    }

    export function entityIndex(l: StructureElement) {
        switch (l.unit.kind) {
            case Unit.Kind.Atomic:
                return l.unit.model.atomicHierarchy.index.getEntityFromChain(l.unit.chainIndex[l.element])
            case Unit.Kind.Spheres:
                return l.unit.model.coarseHierarchy.spheres.entityKey[l.element]
            case Unit.Kind.Gaussians:
                return l.unit.model.coarseHierarchy.gaussians.entityKey[l.element]
        }
    }

    export namespace Loci {
        export function all(structure: Structure): Loci {
            return Loci(structure, structure.units.map(unit => ({
                unit,
                indices: OrderedSet.ofRange<UnitIndex>(0 as UnitIndex, unit.elements.length as UnitIndex)
            })));
        }

        export function remap(loci: Loci, structure: Structure): Loci {
            return Loci(structure, loci.elements.map(e => ({
                unit: structure.unitMap.get(e.unit.id)!,
                indices: e.indices
            })));
        }

        export function union(xs: Loci, ys: Loci): Loci {
            if (xs.elements.length > ys.elements.length) return union(ys, xs);
            if (xs.elements.length === 0) return ys;

            const map = new Map<number, OrderedSet<UnitIndex>>();

            for (const e of xs.elements) map.set(e.unit.id, e.indices);

            const elements: Loci['elements'][0][] = [];
            for (const e of ys.elements) {
                if (map.has(e.unit.id)) {
                    elements[elements.length] = { unit: e.unit, indices: OrderedSet.union(map.get(e.unit.id)!, e.indices) };
                } else {
                    elements[elements.length] = e;
                }
            }

            return Loci(xs.structure, elements);
        }

        export function subtract(xs: Loci, ys: Loci): Loci {
            const map = new Map<number, OrderedSet<UnitIndex>>();
            for (const e of ys.elements) map.set(e.unit.id, e.indices);

            const elements: Loci['elements'][0][] = [];
            for (const e of xs.elements) {
                if (map.has(e.unit.id)) {
                    const indices = OrderedSet.subtract(e.indices, map.get(e.unit.id)!);
                    if (OrderedSet.size(indices) === 0) continue;
                    elements[elements.length] = { unit: e.unit, indices };
                } else {
                    elements[elements.length] = e;
                }
            }

            return Loci(xs.structure, elements);
        }

        export function areIntersecting(xs: Loci, ys: Loci): boolean {
            if (xs.elements.length > ys.elements.length) return areIntersecting(ys, xs);
            if (xs.elements.length === 0) return ys.elements.length === 0;

            const map = new Map<number, OrderedSet<UnitIndex>>();

            for (const e of xs.elements) map.set(e.unit.id, e.indices);
            for (const e of ys.elements) {
                if (!map.has(e.unit.id)) continue;
                if (OrderedSet.areIntersecting(map.get(e.unit.id)!, e.indices)) return true;
            }

            return false;
        }

        export function extendToWholeResidues(loci: Loci): Loci {
            const elements: Loci['elements'][0][] = [];

            for (const lociElement of loci.elements) {
                if (lociElement.unit.kind !== Unit.Kind.Atomic) elements[elements.length] = lociElement;

                const unitElements = lociElement.unit.elements;
                const h = lociElement.unit.model.atomicHierarchy;

                const { index: residueIndex, offsets: residueOffsets } = h.residueAtomSegments;

                const newIndices: UnitIndex[] = [];
                const indices = lociElement.indices, len = OrderedSet.size(indices);
                let i = 0;
                while (i < len) {
                    const rI = residueIndex[unitElements[OrderedSet.getAt(indices, i)]];
                    while (i < len && residueIndex[unitElements[OrderedSet.getAt(indices, i)]] === rI) {
                        i++;
                    }

                    for (let j = residueOffsets[rI], _j = residueOffsets[rI + 1]; j < _j; j++) {
                        const idx = OrderedSet.indexOf(unitElements, j);
                        if (idx >= 0) newIndices[newIndices.length] = idx as UnitIndex;
                    }
                }

                elements[elements.length] = { unit: lociElement.unit, indices: SortedArray.ofSortedArray(newIndices) };
            }

            return Loci(loci.structure, elements);
        }

        const boundaryHelper = new BoundaryHelper(), tempPos = Vec3.zero();
        export function getBoundary(loci: Loci): Boundary {
            boundaryHelper.reset(0);

            for (const e of loci.elements) {
                const { indices } = e;
                const pos = e.unit.conformation.position, r = e.unit.conformation.r;
                const { elements } = e.unit;
                for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                    const eI = elements[OrderedSet.getAt(indices, i)];
                    pos(eI, tempPos);
                    boundaryHelper.boundaryStep(tempPos, r(eI));
                }
            }
            boundaryHelper.finishBoundaryStep();
            for (const e of loci.elements) {
                const { indices } = e;
                const pos = e.unit.conformation.position, r = e.unit.conformation.r;
                const { elements } = e.unit;
                for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                    const eI = elements[OrderedSet.getAt(indices, i)];
                    pos(eI, tempPos);
                    boundaryHelper.extendStep(tempPos, r(eI));
                }
            }

            return { box: boundaryHelper.getBox(), sphere: boundaryHelper.getSphere() };
        }
    }
}

export default StructureElement