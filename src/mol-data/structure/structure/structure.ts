/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model, Format } from '../model'
import Unit from './unit'
import Operator from './operator'
import AtomSet from './atom/set'
import Atom from './atom'
import { OrderedSet, Iterator } from 'mol-base/collections/integer'

interface Structure extends Readonly<{
    units: { readonly [id: number]: Unit },
    atoms: AtomSet
}> { }

namespace Structure {
    export const Empty = { units: {}, atoms: AtomSet.Empty };

    export function create(units: Structure['units'], atoms: AtomSet): Structure {
        return { units, atoms };
    }

    export function ofData(format: Format) {
        const models = Model.create(format);
        return models.map(ofModel);
    }

    export function ofModel(model: Model): Structure {
        const chains = model.hierarchy.chainSegments;
        const builder = Builder();

        for (let c = 0; c < chains.count; c++) {
            const unit = Unit.create(model, Operator.Identity);
            builder.addUnit(unit);
            builder.addAtoms(unit.id, OrderedSet.ofBounds(chains.segments[c], chains.segments[c + 1]));
        }

        return builder.getStructure();
    }

    export interface Builder {
        addUnit(unit: Unit): void,
        addAtoms(unitId: number, atoms: OrderedSet): void,
        getStructure(): Structure,
        readonly atomCount: number
    }

    class BuilderImpl implements Builder {
        private units = Object.create(null);
        private atoms = Object.create(null);
        atomCount = 0;

        addUnit(unit: Unit) { this.units[unit.id] = unit; }
        addAtoms(unitId: number, atoms: OrderedSet) { this.atoms[unitId] = atoms; this.atomCount += OrderedSet.size(atoms); }
        getStructure(): Structure { return this.atomCount > 0 ? Structure.create(this.units, AtomSet.create(this.atoms)) : Empty; }
    }

    export function Builder(): Builder { return new BuilderImpl(); }

    /** Transient = location gets overwritten when move() is called. */
    export function atomLocationsTransient(s: Structure): Iterator<Atom.Location> {
        const l = Atom.Location();
        const update = Atom.updateLocation;
        return Iterator.map(AtomSet.atoms(s.atoms), a => update(s, l, a));
    }

    // TODO: "lift" atom set operators?
    // TODO: "diff"
}

export default Structure