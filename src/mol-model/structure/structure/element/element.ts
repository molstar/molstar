/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SortedArray } from '../../../../mol-data/int';
import { ElementIndex, ResidueIndex, ChainIndex } from '../../model';
import Unit from '../unit';
import { Location } from './location';
import StructureProperties from '../properties';

// TODO: when nominal types are available, make this indexed by UnitIndex
export type Set = SortedArray<ElementIndex>

/** Index into Unit.elements */
export type UnitIndex = { readonly '@type': 'unit-element-index' } & number

export interface Property<T> { (location: Location): T }
export interface Predicate extends Property<boolean> { }

export function property<T>(p: Property<T>) { return p; }

function _wrongUnitKind(kind: string) { throw new Error(`Property only available for ${kind} models.`); }
export function atomicProperty<T>(p: (location: Location<Unit.Atomic>) => T) {
    return property(l => Unit.isAtomic(l.unit) ? p(l as Location<Unit.Atomic>) : _wrongUnitKind('atomic'));
}

export function coarseProperty<T>(p: (location: Location<Unit.Spheres | Unit.Gaussians>) => T) {
    return property(l => Unit.isCoarse(l.unit) ? p(l as Location<Unit.Spheres | Unit.Gaussians>) : _wrongUnitKind('coarse'));
}

export function residueIndex(e: Location) {
    if (Unit.isAtomic(e.unit)) {
        return e.unit.residueIndex[e.element];
    } else {
        // TODO: throw error instead?
        return -1 as ResidueIndex;
    }
}

export function chainIndex(e: Location) {
    if (Unit.isAtomic(e.unit)) {
        return e.unit.chainIndex[e.element];
    } else {
        // TODO: throw error instead?
        return -1 as ChainIndex;
    }
}

export function entityIndex(l: Location) {
    return StructureProperties.entity.key(l);
}