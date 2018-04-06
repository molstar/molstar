/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementGroup, ElementSet, Structure, Unit } from 'mol-model/structure';
import { RenderObject } from 'mol-gl/renderer';
import { EquivalenceClasses } from 'mol-data/util';
import { OrderedSet } from 'mol-data/int'
import { Task } from 'mol-task'

export interface RepresentationProps {

}

export interface UnitRepresentation {
    create: (unit: Unit, elementGroup: ElementGroup, props?: Partial<RepresentationProps>) => Task<RenderObject[]>,
    update: (props: RepresentationProps) => boolean,
}

// export interface StructureRepresentation {
//     create: (structure: Structure, props: RepresentationProps) => boolean,
//     update: (props: RepresentationProps) => boolean
// }

export class StructureRepresentation {
    constructor(private repr: UnitRepresentation) {

    }
    create(structure: Structure) {
        return Task.create('S. repr.', async ctx => {

            const { elements, units } = structure;
            const uniqueGroups = EquivalenceClasses<number, ElementGroup>(
                ElementGroup.hashCode,
                (a, b) => units[a.id].model.id === units[b.id].model.id && OrderedSet.areEqual(a.elements, b.elements));

            for (let i = 0, _i = ElementSet.unitCount(elements); i < _i; i++) {
                const group = ElementSet.unitGetByIndex(elements, i);
                uniqueGroups.add(group.id, group);

            }

            const u = this.repr.create(units[0], 0 as any, 0 as any);
            await ctx.update('Building units...');
            await ctx.runChild(u);

            //ctx.runChild(...)
            //uniqueGroups.groups

            return true
        });
    }
    update(elements: ElementSet, props: RepresentationProps) {
        // TODO check model.id, conformation.id, unit.id, elementGroup(.hashCode/.areEqual)
        return false
    }
}