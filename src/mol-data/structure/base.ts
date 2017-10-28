/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Model from '../model'
import Structure from '../structure'
import Unit from './unit'
import Operator from './operator'
import AtomSet from './atom-set'
import OrderedSet from '../../mol-base/collections/integer/ordered-set'

class Builder {
    private units = Object.create(null);
    private atoms = Object.create(null);

    addUnit(unit: Unit) { this.units[unit.id] = unit; }
    addAtoms(unitId: number, atoms: OrderedSet) { this.atoms[unitId] = atoms; }

    getStructure(): Structure { return { units: this.units, atoms: AtomSet.create(this.atoms) } }
}

export const Empty: Structure = { units: {}, atoms: AtomSet.Empty };

export function ofModel(model: Model): Structure {
    const chains = model.segments.chains;
    const builder = new Builder();

    for (let c = 0; c < chains.count; c++) {
        const unit = Unit.create(model, Operator.Identity);
        builder.addUnit(unit);
        builder.addAtoms(unit.id, OrderedSet.ofBounds(chains.segments[c], chains.segments[c + 1]));
    }

    return builder.getStructure();
}