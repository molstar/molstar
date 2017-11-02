/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from '../model'
import Unit from './unit'
import Operator from './operator'
import AtomSet from './atom/set'
import { OrderedSet } from 'mol-base/collections/integer'

interface Structure extends Readonly<{
    units: { readonly [id: number]: Unit },
    atoms: AtomSet
}> { }

namespace Structure {
    export const Empty = { units: {}, atoms: AtomSet.Empty };

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
        getStructure(): Structure { return this.atomCount > 0 ? { units: this.units, atoms: AtomSet.create(this.atoms) } : Empty; }
    }

    export function Builder(): Builder { return new BuilderImpl(); }


    // TODO: "lift" atom set operators?
    // TODO: "diff"
}

export default Structure