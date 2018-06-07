/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSymmetry, Unit } from 'mol-model/structure';
import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { Representation, RepresentationProps } from '..';
import { ColorTheme } from '../../theme';
import { PickingId } from '../../util/picking';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';

export interface UnitsRepresentation<P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (group: Unit.SymmetryGroup, props: P) => Task<void>
    update: (props: P) => Task<boolean>
    getLoci: (pickingId: PickingId) => Loci | null
    mark: (loci: Loci, action: MarkerAction) => void
}

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

interface GroupRepresentation<T> {
    repr: UnitsRepresentation<T>
    group: Unit.SymmetryGroup
}

export const DefaultStructureProps = {
    colorTheme: { name: 'instance-index' } as ColorTheme,
    alpha: 1,
    visible: true,
    doubleSided: false,
    depthMask: true
}
export type StructureProps = Partial<typeof DefaultStructureProps>

export function StructureRepresentation<P extends StructureProps>(reprCtor: () => UnitsRepresentation<P>): StructureRepresentation<P> {
    const renderObjects: RenderObject[] = []
    const groupReprs: GroupRepresentation<P>[] = []
    // let currentProps: typeof DefaultStructureProps

    function getLoci(pickingId: PickingId) {
        for (let i = 0, il = groupReprs.length; i < il; ++i) {
            const loc = groupReprs[i].repr.getLoci(pickingId)
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
        getLoci,
        mark(loci: Loci, action: MarkerAction) {
            for (let i = 0, il = groupReprs.length; i < il; ++i) {
                groupReprs[i].repr.mark(loci, action)
            }
        }
    }
}