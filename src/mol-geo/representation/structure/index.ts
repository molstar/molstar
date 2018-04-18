/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementGroup, ElementSet, Structure, Unit } from 'mol-model/structure';
import { EquivalenceClasses } from 'mol-data/util';
import { OrderedSet } from 'mol-data/int'
import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/scene';

export interface RepresentationProps {

}

export interface UnitsRepresentation<Props> {
    renderObjects: ReadonlyArray<RenderObject>,
    create: (units: ReadonlyArray<Unit>, elementGroup: ElementGroup, props: Props) => Task<void>,
    update: (props: RepresentationProps) => boolean,
}

export interface StructureRepresentation<Props> {
    renderObjects: ReadonlyArray<RenderObject>,
    create: (structure: Structure, props?: Props) => Task<void>,
    update: (elements: ElementSet, props: Props) => boolean
}

export function StructureRepresentation<Props>(reprCtor: () => UnitsRepresentation<Props>): StructureRepresentation<Props> {
    const renderObjects: RenderObject[] = []
    const unitReprs: UnitsRepresentation<Props>[] = []

    return {
        renderObjects,
        create(structure: Structure, props: Props = {} as Props) {
            return Task.create('StructureRepresentation', async ctx => {
                const { elements, units } = structure;
                const uniqueGroups = EquivalenceClasses<number, ElementGroup>(
                    ElementGroup.hashCode,
                    (a, b) => {
                        console.log(units, a.id, b.id)
                        return units[a.id].model.id === units[b.id].model.id && OrderedSet.areEqual(a.elements, b.elements)
                    }
                );

                for (let i = 0, _i = ElementSet.unitCount(elements); i < _i; i++) {
                    const group = ElementSet.unitGetByIndex(elements, i);
                    uniqueGroups.add(i, group);
                    console.log(i, group)
                }

                for (let i = 0, _i = uniqueGroups.groups.length; i < _i; i++) {
                    const groupUnits: Unit[] = []
                    const group = uniqueGroups.groups[i]
                    // console.log('group', i)
                    for (let j = 0, _j = group.length; j < _j; j++) {
                        groupUnits.push(units[group[j]])
                        // console.log(units[group[j]].operator.matrix)
                    }
                    const elementGroup = ElementSet.unitGetByIndex(elements, group[0])
                    const repr = reprCtor()
                    unitReprs.push(repr)
                    await ctx.update({ message: 'Building units...', current: i, max: _i });
                    await ctx.runChild(repr.create(groupUnits, elementGroup, props));
                    renderObjects.push(...repr.renderObjects)
                }

                // console.log(ElementSet.unitCount(elements))
                // console.log(uniqueGroups)

                // console.log({ elements, units })

                // const repr = reprCtor()
                // unitReprs.push(repr)
                // await ctx.update('Building units...');
                // await ctx.runChild(repr.create(units, elements, props));
                // renderObjects.push(...repr.renderObjects)
            });
        },
        update(elements: ElementSet, props: RepresentationProps) {
            // TODO check model.id, conformation.id, unit.id, elementGroup(.hashCode/.areEqual)
            return false
        }
    }
}