/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, SortedArray } from 'mol-data/int'
import Unit from './unit'
import { ElementIndex } from '../model';

/** Element index in Model */
// type Element = { readonly '@type': 'element' } & number

namespace Element {
    export type Set = SortedArray<ElementIndex>

    /** Index into Unit.elements */
    export type Index = { readonly '@type': 'element-index' } & number

    /** All the information required to access element properties */
    export interface Location<U = Unit> {
        unit: U,
        /** Index into element (atomic/coarse) properties of unit.model */
        element: ElementIndex
    }
    export function Location(unit?: Unit, element?: ElementIndex): Location { return { unit: unit as any, element: element || (0 as ElementIndex) }; }
    export interface Property<T> { (location: Location): T }
    export interface Predicate extends Property<boolean> { }

    export function property<T>(p: Property<T>) { return p; }

    function _wrongUnitKind(kind: string) { throw new Error(`Property only available for ${kind} models.`); }
    export function atomicProperty<T>(p: (location: Location<Unit.Atomic>) => T) {
        return property(l => Unit.isAtomic(l.unit) ? p(l as Location<Unit.Atomic>) : _wrongUnitKind('atomic') );
    }

    export function coarseProperty<T>(p: (location: Location<Unit.Spheres | Unit.Gaussians>) => T) {
        return property(l => Unit.isCoarse(l.unit) ? p(l as Location<Unit.Spheres | Unit.Gaussians>) : _wrongUnitKind('coarse') );
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
            indices: OrderedSet<Index>
        }>
    }

    export function Loci(elements: ArrayLike<{ unit: Unit, indices: OrderedSet<Index> }>): Loci {
        return { kind: 'element-loci', elements: elements as Loci['elements'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'element-loci';
    }
}

export default Element