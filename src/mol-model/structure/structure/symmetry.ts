/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from './structure'
import { Selection } from '../query'
import { ModelSymmetry } from '../model'
import { Task } from 'mol-task';
import { SortedArray } from 'mol-data/int';
import Unit from './unit';
import { EquivalenceClasses } from 'mol-data/util';

namespace StructureSymmetry {
    // Units that have the same elements but differ with operator only.
    export type UnitGroup = { readonly elements: SortedArray, readonly units: ReadonlyArray<Unit> }
    export type TransformGroups = ReadonlyArray<UnitGroup>

    export function buildAssembly(structure: Structure, name: string) {
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

    export function getTransformGroups(s: Structure): StructureSymmetry.TransformGroups {
        // group everything by the "invariantId"
        const invariantGroups = EquivalenceClasses<number, Unit>(u => u.invariantId, (a, b) => a.invariantId === b.invariantId && a.model.id === b.model.id);
        for (const u of s.units) invariantGroups.add(u.id, u);

        const ret: UnitGroup[] = [];
        // group everything by the "element array"
        for (const group of invariantGroups.groups) {
            const setGrouping = EquivalenceClasses<number, Unit>(u => SortedArray.hashCode(u.elements), (a, b) => SortedArray.areEqual(a.elements, b.elements));

            for (const id of group) {
                const unit = s.unitMap.get(id);
                setGrouping.add(unit.id, unit);
            }

            for (const eqUnits of setGrouping.groups) {
                const first = s.unitMap.get(eqUnits[0]);
                ret.push({ elements: first.elements, units: eqUnits.map(id => s.unitMap.get(id)) });
            }
        }

        return ret;
    }
}

export default StructureSymmetry;