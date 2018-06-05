/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSymmetry, Unit, Element, Queries } from 'mol-model/structure';
import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { Representation, RepresentationProps } from '..';
import { ColorTheme } from '../../theme';
import { PickingId, PickingInfo } from '../../util/picking';

export interface UnitsRepresentation<P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (group: Unit.SymmetryGroup, props: P) => Task<void>
    update: (props: P) => Task<boolean>
    getLocation: (pickingId: PickingId) => Element.Location | null
}

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (structure: Structure, props?: P) => Task<void>
    update: (props: P) => Task<void>
    getLocation: (pickingId: PickingId) => Element.Location | null
    getLabel: (pickingId: PickingId) => PickingInfo | null
}

interface GroupRepresentation<T> {
    repr: UnitsRepresentation<T>
    group: Unit.SymmetryGroup
}

function label(loc: Element.Location) {
    const model = loc.unit.model.label
    const instance = loc.unit.conformation.operator.name
    let element = ''

    if (Unit.isAtomic(loc.unit)) {
        const asym_id = Queries.props.chain.auth_asym_id(loc)
        const seq_id = Queries.props.residue.auth_seq_id(loc)
        const comp_id = Queries.props.residue.auth_comp_id(loc)
        const atom_id = Queries.props.atom.auth_atom_id(loc)
        element = `[${comp_id}]${seq_id}:${asym_id}.${atom_id}`
    } else if (Unit.isCoarse(loc.unit)) {
        const asym_id = Queries.props.coarse.asym_id(loc)
        const seq_id_begin = Queries.props.coarse.seq_id_begin(loc)
        const seq_id_end = Queries.props.coarse.seq_id_end(loc)
        if (seq_id_begin === seq_id_end) {
            const entityKey = Queries.props.coarse.entityKey(loc)
            const seq = loc.unit.model.sequence.byEntityKey[entityKey]
            const comp_id = seq.compId.value(seq_id_begin)
            element = `[${comp_id}]${seq_id_begin}:${asym_id}`
        } else {
            element = `${seq_id_begin}-${seq_id_end}:${asym_id}`
        }
    } else {
        element = 'unknown'
    }

    return { label: `${model} ${instance} ${element}` }
}

export const DefaultStructureProps = {
    colorTheme: { name: 'instance-index' } as ColorTheme,
    alpha: 1,
    visible: true,
    doubleSided: false,
    depthMask: true,
    hoverSelection: { objectId: -1, instanceId: -1, elementId: -1 } as PickingId
}
export type StructureProps = Partial<typeof DefaultStructureProps>

export function StructureRepresentation<P extends StructureProps>(reprCtor: () => UnitsRepresentation<P>): StructureRepresentation<P> {
    const renderObjects: RenderObject[] = []
    const groupReprs: GroupRepresentation<P>[] = []
    // let currentProps: typeof DefaultStructureProps

    function getLocation(pickingId: PickingId) {
        for (let i = 0, il = groupReprs.length; i < il; ++i) {
            const loc = groupReprs[i].repr.getLocation(pickingId)
            if (loc) return loc
        }
        return null
    }

    return {
        renderObjects,
        create(structure: Structure, props: P = {} as P) {
            // currentProps = Object.assign({}, DefaultStructureProps, props)

            return Task.create('StructureRepresentation.create', async ctx => {
                // const { query } = currentProps
                // const qs = await query(structure).runAsChild(ctx)
                // const subStructure = Selection.unionStructure(qs)

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
                renderObjects.length = 0 // clear

                for (let i = 0, il = groupReprs.length; i < il; ++i) {
                    const groupRepr = groupReprs[i]
                    const { repr, group } = groupRepr
                    const state = { message: 'Updating structure unit representations...', current: i, max: il };
                    if (!await repr.update(props).runAsChild(ctx, state)) {
                        console.log('update failed, need to rebuild')
                        await repr.create(group, props).runAsChild(ctx, state)
                    }
                    renderObjects.push(...repr.renderObjects)
                }
            })
        },
        getLocation,
        getLabel(pickingId: PickingId) {
            const loc = getLocation(pickingId)
            return loc ? label(loc) : null
        }
    }
}