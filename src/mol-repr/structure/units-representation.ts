/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Unit } from 'mol-model/structure';
import { Task } from 'mol-task'
import { GraphicsRenderObject, getNextMaterialId } from 'mol-gl/render-object';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Visual } from '../visual';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { StructureGroup } from './units-visual';
import { StructureRepresentation, StructureParams, StructureRepresentationState, StructureRepresentationStateBuilder } from './representation';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { Theme, createEmptyTheme } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { UnitKind, UnitKindOptions } from './visual/util/common';
import { Subject } from 'rxjs';

export const UnitsParams = {
    ...StructureParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic', 'spheres'], UnitKindOptions),
}
export type UnitsParams = typeof UnitsParams

export interface UnitsVisual<P extends UnitsParams> extends Visual<StructureGroup, P> { }

export function UnitsRepresentation<P extends UnitsParams>(label: string, ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, P>, visualCtor: (materialId: number) => UnitsVisual<P>): StructureRepresentation<P> {
    let version = 0
    const updated = new Subject<number>()
    const materialId = getNextMaterialId()
    const renderObjects: GraphicsRenderObject[] = []
    const _state = StructureRepresentationStateBuilder.create()
    let visuals = new Map<number, { group: Unit.SymmetryGroup, visual: UnitsVisual<P> }>()

    let _structure: Structure
    let _groups: ReadonlyArray<Unit.SymmetryGroup>
    let _params: P
    let _props: PD.Values<P>
    let _theme = createEmptyTheme()

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, structure?: Structure) {
        if (structure && structure !== _structure) {
            _params = getParams(ctx, structure)
            if (!_props) _props = PD.getDefaultValues(_params)
        }
        _props = Object.assign({}, _props, props)

        return Task.create('Creating or updating UnitsRepresentation', async runtime => {
            if (!_structure && !structure) {
                throw new Error('missing structure')
            } else if (structure && !_structure) {
                // console.log(label, 'initial structure')
                // First call with a structure, create visuals for each group.
                _groups = structure.unitSymmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visual = visualCtor(materialId)
                    const promise = visual.createOrUpdate({ webgl: ctx.webgl, runtime }, _theme, _props, { group, structure })
                    if (promise) await promise
                    visuals.set(group.hashCode, { visual, group })
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length })
                }
            } else if (structure && !Structure.areEquivalent(structure, _structure)) {
                // console.log(label, 'structure not equivalent')
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
                        // console.log(label, 'found visualGroup to reuse')
                        // console.log('old', visualGroup.group)
                        // console.log('new', group)
                        const { visual } = visualGroup
                        const promise = visual.createOrUpdate({ webgl: ctx.webgl, runtime }, _theme, _props, { group, structure })
                        if (promise) await promise
                        visuals.set(group.hashCode, { visual, group })
                        oldVisuals.delete(group.hashCode)
                    } else {
                        // console.log(label, 'not found visualGroup to reuse, creating new')
                        // newGroups.push(group)
                        const visual = visualCtor(materialId)
                        const promise = visual.createOrUpdate({ webgl: ctx.webgl, runtime }, _theme, _props, { group, structure })
                        if (promise) await promise
                        visuals.set(group.hashCode, { visual, group })
                    }
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length })
                }
                oldVisuals.forEach(({ visual }) => {
                    // console.log(label, 'removed unused visual')
                    visual.destroy()
                })

                // TODO review logic
                // For new groups, re-use left-over visuals
                // const unusedVisuals: UnitsVisual<P>[] = []
                // oldVisuals.forEach(({ visual }) => unusedVisuals.push(visual))
                // newGroups.forEach(async group => {
                //     const visual = unusedVisuals.pop() || visualCtor()
                //     await visual.createOrUpdate({ ...ctx, runtime }, _props, group)
                //     visuals.set(group.hashCode, { visual, group })
                // })
                // unusedVisuals.forEach(visual => visual.destroy())
            } else if (structure && structure !== _structure && Structure.areEquivalent(structure, _structure)) {
                // console.log(label, 'structures equivalent but not identical')
                // Expects that for structures with the same hashCode,
                // the unitSymmetryGroups are the same as well.
                // Re-uses existing visuals for the groups of the new structure.
                _groups = structure.unitSymmetryGroups;
                // console.log('new', structure.unitSymmetryGroups)
                // console.log('old', _structure.unitSymmetryGroups)
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visualGroup = visuals.get(group.hashCode)
                    if (visualGroup) {
                        const promise = visualGroup.visual.createOrUpdate({ webgl: ctx.webgl, runtime }, _theme, _props, { group, structure })
                        if (promise) await promise
                        visualGroup.group = group

                    } else {
                        throw new Error(`expected to find visual for hashCode ${group.hashCode}`)
                    }
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length })
                }
            } else {
                // console.log(label, 'no new structure')
                // No new structure given, just update all visuals with new props.
                const visualsList: [ UnitsVisual<P>, Unit.SymmetryGroup ][] = [] // TODO avoid allocation
                visuals.forEach(({ visual, group }) => visualsList.push([ visual, group ]))
                for (let i = 0, il = visualsList.length; i < il; ++i) {
                    const [ visual ] = visualsList[i]
                    const promise = visual.createOrUpdate({ webgl: ctx.webgl, runtime }, _theme, _props)
                    if (promise) await promise
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: il })
                }
            }
            // update list of renderObjects
            renderObjects.length = 0
            visuals.forEach(({ visual }) => {
                if (visual.renderObject) renderObjects.push(visual.renderObject)
            })
            // set new structure
            if (structure) _structure = structure
            // increment version
            updated.next(version++)
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

    function setState(state: Partial<StructureRepresentationState>) {
        const { visible, alphaFactor, pickable, overpaint, transform, unitTransforms } = state
        if (visible !== undefined) visuals.forEach(({ visual }) => visual.setVisibility(visible))
        if (alphaFactor !== undefined) visuals.forEach(({ visual }) => visual.setAlphaFactor(alphaFactor))
        if (pickable !== undefined) visuals.forEach(({ visual }) => visual.setPickable(pickable))
        if (overpaint !== undefined) visuals.forEach(({ visual }) => visual.setOverpaint(overpaint))
        if (transform !== undefined) visuals.forEach(({ visual }) => visual.setTransform(transform))
        if (unitTransforms !== undefined) {
            visuals.forEach(({ visual, group }) => {
                if (unitTransforms) {
                    // console.log(group.hashCode, unitTransforms.getSymmetryGroupTransforms(group))
                    visual.setTransform(undefined, unitTransforms.getSymmetryGroupTransforms(group))
                } else {
                    visual.setTransform(undefined, null)
                }
            })
        }

        StructureRepresentationStateBuilder.update(_state, state)
    }

    function setTheme(theme: Theme) {
        _theme = theme
    }

    function destroy() {
        visuals.forEach(({ visual }) => visual.destroy())
        visuals.clear()
    }

    return {
        label,
        get groupCount() {
            let groupCount = 0
            visuals.forEach(({ visual }) => {
                if (visual.renderObject) groupCount += visual.groupCount
            })
            return groupCount
        },
        get props() { return _props },
        get params() { return _params },
        get state() { return _state },
        get theme() { return _theme },
        renderObjects,
        updated,
        createOrUpdate,
        setState,
        setTheme,
        getLoci,
        mark,
        destroy
    }
}