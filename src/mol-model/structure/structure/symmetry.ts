/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from './structure'
import AtomSet from './atom/set'
import Unit from './unit'
import { Selection } from '../query'
import { ModelSymmetry } from '../model'

namespace Symmetry {
    export const  buildAssembly = buildAssemblyImpl;
}

export default Symmetry;

function buildAssemblyImpl(structure: Structure, name: string) {
    const models = Structure.getModels(structure);
    if (models.length !== 1) throw new Error('Can only build assemblies from structures based on 1 model.');

    const assembly = ModelSymmetry.findAssembly(models[0], name);
    if (!assembly) throw new Error(`Assembly '${name}' is not defined.`);

    const assembler = Structure.Builder();

    let unitId = 0;
    for (const g of assembly.operatorGroups) {
        const selection = g.selector(structure);
        if (Selection.structureCount(selection) === 0) continue;
        const { units, atoms } = Selection.union(selection);

        const unitIds = AtomSet.unitIds(atoms);

        for (const oper of g.operators) {
            for (let uI = 0, _uI = unitIds.length; uI < _uI; uI++) {
                const unit = units.get(unitIds[uI]);
                assembler.add(Unit.create(unitId++, unit.model, oper, unit.naturalGroup), AtomSet.unitGetByIndex(atoms, uI));
            }
        }
    }

    return assembler.getStructure();
}