/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from '../structure'
import Atom from './atom'
import Unit from './unit'
import Model from '../model'
import Column from '../../mol-base/collections/column'

/** Atom pointer */
interface Property<T> { (location: Property.Location): T }

namespace Property {
    export interface Predicate extends Property<boolean> { }

    /** All the information required to access atom properties */
    export interface Location {
        //structure: Structure,
        unit: Unit,
        // hierarchy: Model['hierarchy'],
        // conformation: Model['conformation'],
        atomIndex: number,
        // residueIndex: number,
        // chainIndex: number
    }

    export function createLocation(): Location {
        return {
            //structure: structure!,
            unit: {} as any,
            //hierarchy: (!unit ? void 0 : unit.model.hierarchy)!,
            //conformation: (!unit ? void 0 : unit.model.conformation)!,
            atomIndex: 0
            //residueIndex: 0,
            //chainIndex: 0
        };
    }

    export function update(l: Location, structure: Structure, atom: Atom) {
        // const u = structure.units[Atom.unit(atom)];
        // const i = Atom.index(atom);
        // l.structure = structure;
        // l.unit = u;
        // l.atomIndex = i;
        // l.residueIndex = u.residueIndex[i];
        // l.chainIndex = u.chainIndex[i];
        throw 'not implemented'
    }

    export function naive<T>(p: Property<T>) { return p; }

    export function cachedAtomColumn<T>(col: (model: Model) => Column<T>): Property<T> {
        let lastUnit: Unit | undefined = void 0;
        let cached: ((row: number) => T) | undefined = void 0;
        return function (location) {
            if (location.unit === lastUnit && !!cached) return cached(location.atomIndex);
            lastUnit = location.unit;
            cached = col(lastUnit.model).value;
            return cached(location.atomIndex);
        };
    }

    // TODO: cached versions of other properties.
}

export default Property