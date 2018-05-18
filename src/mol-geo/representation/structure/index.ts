/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSymmetry, Unit } from 'mol-model/structure';
import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/scene';
import { Representation, RepresentationProps } from '..';
// import { Mat4, EPSILON } from 'mol-math/linear-algebra';


export interface UnitsRepresentation<P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (group: Unit.SymmetryGroup, props: P) => Task<void>
    update: (props: P) => Task<boolean>
}

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (structure: Structure, props?: P) => Task<void>
    update: (props: P) => Task<void>
}

interface GroupRepresentation<T> {
    repr: UnitsRepresentation<T>
    group: Unit.SymmetryGroup
}

export function StructureRepresentation<P>(reprCtor: () => UnitsRepresentation<P>): StructureRepresentation<P> {
    const renderObjects: RenderObject[] = []
    const groupReprs: GroupRepresentation<P>[] = []

    return {
        renderObjects,
        create(structure: Structure, props: P = {} as P) {
            return Task.create('StructureRepresentation.create', async ctx => {
                const groups = StructureSymmetry.getTransformGroups(structure);
                for (let i = 0; i < groups.length; i++) {
                    const group = groups[i];
                    const repr = reprCtor()
                    groupReprs.push({ repr, group })
                    await repr.create(group, props).runAsChild(ctx, { message: 'Building structure unit representations...', current: i, max: groups.length });
                    renderObjects.push(...repr.renderObjects)
                }
            });
        },
        update(props: P) {
            return Task.create('StructureRepresentation.update', async ctx => {
                // TODO check model.id, conformation.id, unit.id, elementGroup(.hashCode/.areEqual)
                renderObjects.length = 0 // clear
                for (let i = 0, il = groupReprs.length; i < il; ++i) {
                    const groupRepr = groupReprs[i]
                    const { repr, group } = groupRepr
                    const state = { message: 'Updating structure unit representations...', current: i, max: il };
                    if (!await repr.update(props).runAsChild(ctx, state)) {
                        await repr.create(group, props).runAsChild(ctx, state)
                    }
                    renderObjects.push(...repr.renderObjects)
                }
            })
        }
    }
}