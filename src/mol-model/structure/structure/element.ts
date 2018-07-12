/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, SortedArray } from 'mol-data/int'
import Unit from './unit'
import { ElementIndex } from '../model';

interface StructureElement<U = Unit> {
    unit: U,
    /** Index into element (atomic/coarse) properties of unit.model */
    element: ElementIndex
}

namespace StructureElement {
    export function create(unit?: Unit, element?: ElementIndex): StructureElement { return { unit: unit as any, element: element || (0 as ElementIndex) }; }

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

    export function Loci(elements: ArrayLike<{ unit: Unit, indices: OrderedSet<UnitIndex> }>): Loci {
        return { kind: 'element-loci', elements: elements as Loci['elements'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'element-loci';
    }
}

export default StructureElement