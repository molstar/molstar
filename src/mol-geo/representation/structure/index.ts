/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSymmetry, Unit } from 'mol-model/structure';
import { Task, RuntimeContext } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { Representation, RepresentationProps } from '..';
import { ColorTheme } from '../../theme';
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';

export interface UnitsVisual<P> {
    readonly renderObjects: ReadonlyArray<RenderObject>
    create: (ctx: RuntimeContext, group: Unit.SymmetryGroup, props: P) => Promise<void>
    update: (ctx: RuntimeContext, props: P) => Promise<boolean>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => void
}

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

interface GroupVisual<T> {
    readonly visual: UnitsVisual<T>
    readonly group: Unit.SymmetryGroup
}

export const DefaultStructureProps = {
    colorTheme: { name: 'instance-index' } as ColorTheme,
    alpha: 1,
    visible: true,
    doubleSided: false,
    depthMask: true,
    useFog: true,
}
export type StructureProps = Partial<typeof DefaultStructureProps>

export function StructureRepresentation<P extends StructureProps>(visualCtor: () => UnitsVisual<P>): StructureRepresentation<P> {
    const renderObjects: RenderObject[] = []
    const groupVisuals: GroupVisual<P>[] = []
    // let currentProps: typeof DefaultStructureProps

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
                    if(ctx.shouldUpdate) await ctx.update({
                        message: 'Building structure unit visuals...',
                        current: i,
                        max: groups.length
                    })
                    const group = groups[i];
                    const visual = visualCtor()
                    groupVisuals.push({ visual, group })
                    await visual.create(ctx, group, props)
                    renderObjects.push(...visual.renderObjects)
                }
            });
        },
        update(props: P) {
            return Task.create('StructureRepresentation.update', async ctx => {
                renderObjects.length = 0 // clear

                for (let i = 0, il = groupVisuals.length; i < il; ++i) {
                    if(ctx.shouldUpdate) await ctx.update({
                        message: 'Updating structure unit visuals...',
                        current: i,
                        max: il
                    })
                    const groupVisual = groupVisuals[i]
                    const { visual, group } = groupVisual

                    if (!await visual.update(ctx, props)) {
                        console.log('update failed, need to rebuild')
                        visual.create(ctx, group, props)
                    }
                    renderObjects.push(...visual.renderObjects)
                }
            })
        },
        getLoci(pickingId: PickingId) {
            for (let i = 0, il = groupVisuals.length; i < il; ++i) {
                const loc = groupVisuals[i].visual.getLoci(pickingId)
                if (!isEmptyLoci(loc)) return loc
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            for (let i = 0, il = groupVisuals.length; i < il; ++i) {
                groupVisuals[i].visual.mark(loci, action)
            }
        }
    }
}