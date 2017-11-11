/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, Iterator } from 'mol-data/int'
import { UniqueArray } from 'mol-data/util'
import SymmetryOperator from 'mol-math/geometry/symmetry-operator'
import { Model, Format } from '../model'
import Unit from './unit'
import AtomSet from './atom/set'
import AtomGroup from './atom/group'
import Atom from './atom'


interface Structure extends Readonly<{
    units: Unit[],
    atoms: AtomSet
}> { }

namespace Structure {
    export const Empty: Structure = { units: [], atoms: AtomSet.Empty };

    export function create(units: Unit[], atoms: AtomSet): Structure {
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
            const group = AtomGroup.createNew(OrderedSet.ofBounds(chains.segments[c], chains.segments[c + 1]));
            const unit = Unit.create(model, SymmetryOperator.Default, group);
            builder.add(unit, OrderedSet.ofBounds(chains.segments[c], chains.segments[c + 1]));
        }

        return builder.getStructure();
    }

    export interface Builder {
        add(unit: Unit, atoms: OrderedSet): void,
        addUnit(unit: Unit): void,
        setAtoms(unitId: number, atoms: OrderedSet): void,
        getStructure(): Structure,
        readonly atomCount: number
    }

    class BuilderImpl implements Builder {
        private _unitId = 0;
        private units: Unit[] = [];
        private atoms = AtomSet.Generator();
        atomCount = 0;

        add(unit: Unit, atoms: OrderedSet) { const id = this.addUnit(unit); this.setAtoms(id, atoms); }
        addUnit(unit: Unit) { const id = this._unitId++; this.units[id] = unit; return id; }
        setAtoms(unitId: number, atoms: OrderedSet) { this.atoms.add(unitId, atoms); this.atomCount += OrderedSet.size(atoms); }
        getStructure(): Structure { return this.atomCount > 0 ? Structure.create(this.units, this.atoms.getSet()) : Empty; }
    }

    export function Builder(): Builder { return new BuilderImpl(); }

    /** Transient = location gets overwritten when move() is called. */
    export function atomLocationsTransient(s: Structure): Iterator<Atom.Location> {
        const l = Atom.Location();
        const update = Atom.updateLocation;
        return Iterator.map(AtomSet.atoms(s.atoms), a => update(s, l, a));
    }

    export function getModels(s: Structure) {
        const { units, atoms } = s;
        const arr = UniqueArray.create<Model['id'], Model>();
        const ids = AtomSet.unitIds(atoms);
        for (let i = 0; i < ids.length; i++) {
            const u = units[ids[i]];
            UniqueArray.add(arr, u.model.id, u.model);
        }
        return arr.array;
    }

    // TODO: "lift" atom set operators?
    // TODO: "diff"
}

export default Structure