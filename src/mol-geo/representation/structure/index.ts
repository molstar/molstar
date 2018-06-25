/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSymmetry, Unit } from 'mol-model/structure';
import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { Representation, RepresentationProps, Visual } from '..';
import { ColorTheme, SizeTheme } from '../../theme';
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<Unit.SymmetryGroup, P> { }
export interface  StructureVisual<P extends RepresentationProps = {}> extends Visual<Structure, P> { }

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export const DefaultStructureProps = {
    colorTheme: { name: 'instance-index' } as ColorTheme,
    sizeTheme: { name: 'physical' } as SizeTheme,
    alpha: 1,
    visible: true,
    doubleSided: false,
    depthMask: true,
    useFog: true,
}
export type StructureProps = Partial<typeof DefaultStructureProps>

export function StructureRepresentation<P extends StructureProps>(unitsVisualCtor: () => UnitsVisual<P>, structureVisualCtor?: () => StructureVisual<P>): StructureRepresentation<P> {
    let unitsVisuals = new Map<number, { group: Unit.SymmetryGroup, visual: UnitsVisual<P> }>()
    let structureVisual: StructureVisual<P> | undefined

    let _props: Required<P>
    let _structure: Structure
    let _groups: ReadonlyArray<Unit.SymmetryGroup>

    function create(structure: Structure, props: P = {} as P) {
        _props = Object.assign({}, DefaultStructureProps, _props, props)

        return Task.create('Creating StructureRepresentation', async ctx => {
            if (!_structure) {
                _groups = structure.symmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visual = unitsVisualCtor()
                    await visual.create(ctx, group, _props)
                    unitsVisuals.set(group.hashCode, { visual, group })
                }

                if (structureVisualCtor) {
                    structureVisual = structureVisualCtor()
                    await structureVisual.create(ctx, structure, _props)
                }
            } else {
                if (_structure.hashCode === structure.hashCode) {
                    await update(_props)
                } else {
                    _groups = structure.symmetryGroups;
                    const newGroups: Unit.SymmetryGroup[] = []
                    const oldUnitsVisuals = unitsVisuals
                    unitsVisuals = new Map()
                    for (let i = 0; i < _groups.length; i++) {
                        const group = _groups[i];
                        const visualGroup = oldUnitsVisuals.get(group.hashCode)
                        if (visualGroup) {
                            const { visual, group } = visualGroup
                            if (!await visual.update(ctx, _props)) {
                                await visual.create(ctx, group, _props)
                            }
                            oldUnitsVisuals.delete(group.hashCode)
                        } else {
                            newGroups.push(group)
                            const visual = unitsVisualCtor()
                            await visual.create(ctx, group, _props)
                            unitsVisuals.set(group.hashCode, { visual, group })
                        }
                    }
    
                    // for new groups, reuse leftover visuals
                    const unusedVisuals: UnitsVisual<P>[] = []
                    oldUnitsVisuals.forEach(({ visual }) => unusedVisuals.push(visual))
                    newGroups.forEach(async group => {
                        const visual = unusedVisuals.pop() || unitsVisualCtor()
                        await visual.create(ctx, group, _props)
                        unitsVisuals.set(group.hashCode, { visual, group })
                    })
                    unusedVisuals.forEach(visual => visual.destroy())

                    if (structureVisual) {
                        if (!await structureVisual.update(ctx, _props)) {
                            await structureVisual.create(ctx, _structure, _props)
                        }
                    }
                }
            }
            _structure = structure
        });
    }

    function update(props: P) {
        return Task.create('Updating StructureRepresentation', async ctx => {
            _props = Object.assign({}, DefaultStructureProps, _props, props)

            unitsVisuals.forEach(async ({ visual, group }) => {
                if (!await visual.update(ctx, _props)) {
                    await visual.create(ctx, group, _props)
                }
            })

            if (structureVisual) {
                if (!await structureVisual.update(ctx, _props)) {
                    await structureVisual.create(ctx, _structure, _props)
                }
            }
        })
    }

    function getLoci(pickingId: PickingId) {
        let loci: Loci = EmptyLoci
        unitsVisuals.forEach(({ visual }) => {
            const _loci = visual.getLoci(pickingId)
            if (!isEmptyLoci(_loci)) loci = _loci
        })
        return loci
    }

    function mark(loci: Loci, action: MarkerAction) {
        unitsVisuals.forEach(({ visual }) => visual.mark(loci, action))
    }

    function destroy() {
        unitsVisuals.forEach(({ visual }) => visual.destroy())
        unitsVisuals.clear()
    }

    return {
        get renderObjects() {
            const renderObjects: RenderObject[] = []
            unitsVisuals.forEach(({ visual }) => renderObjects.push(...visual.renderObjects))
            return renderObjects
        },
        create,
        update,
        getLoci,
        mark,
        destroy
    }
}