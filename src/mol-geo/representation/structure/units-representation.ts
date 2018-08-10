/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Unit } from 'mol-model/structure';
import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { Representation, RepresentationProps, Visual } from '..';
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';
import { getQualityProps } from '../util';
import { DefaultStructureProps, StructureProps } from '.';

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<Unit.SymmetryGroup, P> { }
export interface  StructureVisual<P extends RepresentationProps = {}> extends Visual<Structure, P> { }

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export function UnitsRepresentation<P extends StructureProps>(visualCtor: () => UnitsVisual<P>): StructureRepresentation<P> {
    let visuals = new Map<number, { group: Unit.SymmetryGroup, visual: UnitsVisual<P> }>()

    let _props: Required<P>
    let _structure: Structure
    let _groups: ReadonlyArray<Unit.SymmetryGroup>

    function create(structure: Structure, props: P = {} as P) {
        _props = Object.assign({}, DefaultStructureProps, _props, props, getQualityProps(props, structure))
        _props.colorTheme!.structure = structure

        return Task.create('Creating StructureRepresentation', async ctx => {
            if (!_structure) {
                _groups = structure.unitSymmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visual = visualCtor()
                    await visual.create(ctx, group, _props)
                    visuals.set(group.hashCode, { visual, group })
                }
            } else {
                if (_structure.hashCode === structure.hashCode) {
                    await update(_props)
                } else {
                    _groups = structure.unitSymmetryGroups;
                    const newGroups: Unit.SymmetryGroup[] = []
                    const oldUnitsVisuals = visuals
                    visuals = new Map()
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
                            const visual = visualCtor()
                            await visual.create(ctx, group, _props)
                            visuals.set(group.hashCode, { visual, group })
                        }
                    }

                    // for new groups, reuse leftover visuals
                    const unusedVisuals: UnitsVisual<P>[] = []
                    oldUnitsVisuals.forEach(({ visual }) => unusedVisuals.push(visual))
                    newGroups.forEach(async group => {
                        const visual = unusedVisuals.pop() || visualCtor()
                        await visual.create(ctx, group, _props)
                        visuals.set(group.hashCode, { visual, group })
                    })
                    unusedVisuals.forEach(visual => visual.destroy())
                }
            }
            _structure = structure
        });
    }

    function update(props: P) {
        return Task.create('Updating StructureRepresentation', async ctx => {
            _props = Object.assign({}, DefaultStructureProps, _props, props, getQualityProps(props, _structure))
            _props.colorTheme!.structure = _structure

            visuals.forEach(async ({ visual, group }) => {
                if (!await visual.update(ctx, _props)) {
                    await visual.create(ctx, group, _props)
                }
            })
        })
    }

    function getLoci(pickingId: PickingId) {
        let loci: Loci = EmptyLoci
        visuals.forEach(({ visual }) => {
            const _loci = visual.getLoci(pickingId)
            if (!isEmptyLoci(_loci)) loci = _loci
        })
        return loci
    }

    function mark(loci: Loci, action: MarkerAction) {
        visuals.forEach(({ visual }) => visual.mark(loci, action))
    }

    function destroy() {
        visuals.forEach(({ visual }) => visual.destroy())
        visuals.clear()
    }

    return {
        get renderObjects() {
            const renderObjects: RenderObject[] = []
            visuals.forEach(({ visual }) => renderObjects.push(visual.renderObject))
            return renderObjects
        },
        get props() {
            return _props
        },
        create,
        update,
        getLoci,
        mark,
        destroy
    }
}