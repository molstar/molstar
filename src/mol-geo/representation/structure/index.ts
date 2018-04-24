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
// import { Mat4, EPSILON } from 'mol-math/linear-algebra';

export interface RepresentationProps {

}

export interface UnitsRepresentation<Props = {}> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (units: ReadonlyArray<Unit>, elementGroup: ElementGroup, props: Props) => Task<void>
    update: (props: RepresentationProps) => Task<boolean>
}

export interface StructureRepresentation<Props = {}> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (structure: Structure, props?: Props) => Task<void>
    update: (props: Props) => Task<void>
}

interface GroupRepresentation<T> {
    repr: UnitsRepresentation<T>
    units: Unit[]
    elementGroup: ElementGroup
}

export function StructureRepresentation<Props>(reprCtor: () => UnitsRepresentation<Props>): StructureRepresentation<Props> {
    const renderObjects: RenderObject[] = []
    const groupReprs: GroupRepresentation<Props>[] = []

    return {
        renderObjects,
        create(structure: Structure, props: Props = {} as Props) {
            return Task.create('StructureRepresentation.create', async ctx => {
                const { elements, units } = structure;
                const uniqueGroups = EquivalenceClasses<number, { unit: Unit, group: ElementGroup }>(
                    ({ unit, group }) => ElementGroup.hashCode(group),
                    (a, b) => a.unit.model.id === b.unit.model.id && OrderedSet.areEqual(a.group.elements, b.group.elements)
                );

                // const uniqueTransformations = EquivalenceClasses<number, { unit: Unit, group: ElementGroup }>(
                //     ({ unit, group }) => unit.operator.matrix.join(','),
                //     (a, b) => Mat4.areEqual(a.unit.operator.matrix, b.unit.operator.matrix, EPSILON.Value)
                // );

                const unitIndices = ElementSet.unitIndices(elements);
                for (let i = 0, _i = unitIndices.length; i < _i; i++) {
                    const unitIndex = unitIndices[i];
                    const group = ElementSet.groupFromUnitIndex(elements, unitIndex);
                    const unit = units[unitIndex]
                    uniqueGroups.add(unitIndex, { unit, group });
                    // uniqueTransformations.add(unitIndex, { unit, group });
                }

                // console.log({ uniqueGroups, uniqueTransformations })

                for (let i = 0, il = uniqueGroups.groups.length; i < il; i++) {
                    const groupUnits: Unit[] = []
                    const group = uniqueGroups.groups[i]
                    // console.log('group', i)
                    for (let j = 0, jl = group.length; j < jl; j++) {
                        groupUnits.push(units[group[j]])
                    }
                    const elementGroup = ElementSet.groupFromUnitIndex(elements, group[0])
                    const repr = reprCtor()
                    groupReprs.push({ repr, units: groupUnits, elementGroup })
                    await ctx.update({ message: 'Building units...', current: i, max: il });
                    await ctx.runChild(repr.create(groupUnits, elementGroup, props));
                    renderObjects.push(...repr.renderObjects)
                }
            });
        },
        update(props: Props) {
            return Task.create('StructureRepresentation.update', async ctx => {
                // TODO check model.id, conformation.id, unit.id, elementGroup(.hashCode/.areEqual)
                renderObjects.length = 0 // clear
                for (let i = 0, il = groupReprs.length; i < il; ++i) {
                    const groupRepr = groupReprs[i]
                    const { repr, units, elementGroup } = groupRepr
                    await ctx.update({ message: 'Updating units...', current: i, max: il });
                    if (!await ctx.runChild(repr.update(props))) {
                        await ctx.runChild(repr.create(units, elementGroup, props))
                    }
                    renderObjects.push(...repr.renderObjects)
                }
            })
        }
    }
}