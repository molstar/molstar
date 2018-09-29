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
import { PickingId } from '../../geometry/picking';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../geometry/marker-data';
import { StructureProps, StructureParams } from '.';
import { StructureGroup } from './units-visual';

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> { }
export interface  StructureVisual<P extends RepresentationProps = {}> extends Visual<Structure, P> { }

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export function UnitsRepresentation<P extends StructureProps>(label: string, visualCtor: () => UnitsVisual<P>): StructureRepresentation<P> {
    let visuals = new Map<number, { group: Unit.SymmetryGroup, visual: UnitsVisual<P> }>()

    let _props: P
    let _structure: Structure
    let _groups: ReadonlyArray<Unit.SymmetryGroup>

    function createOrUpdate(props: Partial<P> = {}, structure?: Structure) {
        _props = Object.assign({}, _props, props)

        return Task.create('Creating or updating UnitsRepresentation', async ctx => {
            if (!_structure && !structure) {
                throw new Error('missing structure')
            } else if (structure && !_structure) {
                // console.log('initial structure')
                // First call with a structure, create visuals for each group.
                _groups = structure.unitSymmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visual = visualCtor()
                    await visual.createOrUpdate(ctx, _props, { group, structure })
                    visuals.set(group.hashCode, { visual, group })
                }
            } else if (structure && _structure.hashCode !== structure.hashCode) {
                // console.log('_structure.hashCode !== structure.hashCode')
                // Tries to re-use existing visuals for the groups of the new structure.
                // Creates additional visuals if needed, destroys left-over visuals.
                _groups = structure.unitSymmetryGroups;
                // const newGroups: Unit.SymmetryGroup[] = []
                const oldVisuals = visuals
                visuals = new Map()
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visualGroup = oldVisuals.get(group.hashCode)
                    if (visualGroup) {
                        const { visual } = visualGroup
                        await visual.createOrUpdate(ctx, _props, { group, structure })
                        visuals.set(group.hashCode, { visual, group })
                        oldVisuals.delete(group.hashCode)
                    } else {
                        // newGroups.push(group)
                        const visual = visualCtor()
                        await visual.createOrUpdate(ctx, _props, { group, structure })
                        visuals.set(group.hashCode, { visual, group })
                    }
                }
                oldVisuals.forEach(({ visual }) => visual.destroy())

                // TODO review logic
                // For new groups, re-use left-over visuals
                // const unusedVisuals: UnitsVisual<P>[] = []
                // oldVisuals.forEach(({ visual }) => unusedVisuals.push(visual))
                // newGroups.forEach(async group => {
                //     const visual = unusedVisuals.pop() || visualCtor()
                //     await visual.createOrUpdate(ctx, _props, group)
                //     visuals.set(group.hashCode, { visual, group })
                // })
                // unusedVisuals.forEach(visual => visual.destroy())
            } else if (structure && _structure.hashCode === structure.hashCode) {
                // console.log('_structure.hashCode === structure.hashCode')
                // Expects that for structures with the same hashCode,
                // the unitSymmetryGroups are the same as well.
                // Re-uses existing visuals for the groups of the new structure.
                _groups = structure.unitSymmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visualGroup = visuals.get(group.hashCode)
                    if (visualGroup) {
                        await visualGroup.visual.createOrUpdate(ctx, _props, { group, structure })
                        visualGroup.group = group
                    } else {
                        throw new Error(`expected to find visual for hashCode ${group.hashCode}`)
                    }
                }
            } else {
                // console.log('no new structure')
                // No new structure given, just update all visuals with new props.
                const visualsList: [ UnitsVisual<P>, Unit.SymmetryGroup ][] = [] // TODO avoid allocation
                visuals.forEach(({ visual, group }) => visualsList.push([ visual, group ]))
                for (let i = 0, il = visualsList.length; i < il; ++i) {
                    const [ visual, group ] = visualsList[i]
                    await visual.createOrUpdate(ctx, _props, { group, structure: _structure })
                }
            }
            if (structure) _structure = structure
        });
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
        let changed = false
        visuals.forEach(({ visual }) => {
            changed = visual.mark(loci, action) || changed
        })
        return changed
    }

    function destroy() {
        visuals.forEach(({ visual }) => visual.destroy())
        visuals.clear()
    }

    return {
        label,
        params: StructureParams, // TODO
        get renderObjects() {
            const renderObjects: RenderObject[] = []
            visuals.forEach(({ visual }) => {
                if (visual.renderObject) renderObjects.push(visual.renderObject)
            })
            return renderObjects
        },
        get props() {
            return _props
        },
        createOrUpdate,
        getLoci,
        mark,
        destroy
    }
}