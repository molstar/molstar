/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from './structure'
import { Selection } from '../query'
import { ModelSymmetry } from '../model'
import { Task } from 'mol-task';

namespace Symmetry {
    export const buildAssembly = buildAssemblyImpl;
}

export default Symmetry;

function buildAssemblyImpl(structure: Structure, name: string) {
    return Task.create('Build Symmetry', async ctx => {
        const models = Structure.getModels(structure);
        if (models.length !== 1) throw new Error('Can only build assemblies from structures based on 1 model.');

        const assembly = ModelSymmetry.findAssembly(models[0], name);
        if (!assembly) throw new Error(`Assembly '${name}' is not defined.`);

        const assembler = Structure.Builder();

        for (const g of assembly.operatorGroups) {
            const selection = await ctx.runChild(g.selector(structure));
            if (Selection.structureCount(selection) === 0) {
                continue;
            }
            const { units } = Selection.unionStructure(selection);

            for (const oper of g.operators) {
                for (const unit of units) {
                    assembler.addWithOperator(unit, oper);
                }
            }
        }

        return assembler.getStructure();
    });
}