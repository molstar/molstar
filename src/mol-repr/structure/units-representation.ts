/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Unit } from 'mol-model/structure';
import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { Visual, RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { StructureGroup } from './units-visual';
import { StructureRepresentation, StructureParams } from './representation';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { Theme, ThemeProps, createTheme } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { UnitKind, UnitKindOptions } from './visual/util/common';
import { BehaviorSubject } from 'rxjs';

export const UnitsParams = {
    ...StructureParams,
    unitKinds: PD.MultiSelect<UnitKind>('Unit Kind', '', ['atomic', 'spheres'], UnitKindOptions),
}
export type UnitsParams = typeof UnitsParams

export interface UnitsVisual<P extends UnitsParams> extends Visual<StructureGroup, P> { }

export function UnitsRepresentation<P extends UnitsParams>(label: string, getParams: RepresentationParamsGetter<Structure, P>, visualCtor: () => UnitsVisual<P>): StructureRepresentation<P> {
    const updated = new BehaviorSubject(0)
    let visuals = new Map<number, { group: Unit.SymmetryGroup, visual: UnitsVisual<P> }>()

    let _structure: Structure
    let _groups: ReadonlyArray<Unit.SymmetryGroup>
    let _params: P
    let _props: PD.Values<P>
    let _theme: Theme

    function createOrUpdate(ctx: RepresentationContext, props: Partial<PD.Values<P>> = {}, themeProps: ThemeProps = {}, structure?: Structure) {
        if (structure && structure !== _structure) {
            _params = getParams(ctx, structure)
            if (!_props) _props = PD.getDefaultValues(_params)
        }
        _props = Object.assign({}, _props, props)
        _theme = createTheme(ctx, { structure: structure || _structure }, props, themeProps, _theme)

        return Task.create('Creating or updating UnitsRepresentation', async runtime => {
            if (!_structure && !structure) {
                throw new Error('missing structure')
            } else if (structure && !_structure) {
                // console.log('initial structure')
                // First call with a structure, create visuals for each group.
                _groups = structure.unitSymmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visual = visualCtor()
                    await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, { group, structure })
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
                        await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, { group, structure })
                        visuals.set(group.hashCode, { visual, group })
                        oldVisuals.delete(group.hashCode)
                    } else {
                        // newGroups.push(group)
                        const visual = visualCtor()
                        await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, { group, structure })
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
                //     await visual.createOrUpdate({ ...ctx, runtime }, _props, group)
                //     visuals.set(group.hashCode, { visual, group })
                // })
                // unusedVisuals.forEach(visual => visual.destroy())
            } else if (structure && structure !== _structure && _structure.hashCode === structure.hashCode) {
                // console.log('_structure.hashCode === structure.hashCode')
                // Expects that for structures with the same hashCode,
                // the unitSymmetryGroups are the same as well.
                // Re-uses existing visuals for the groups of the new structure.
                _groups = structure.unitSymmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visualGroup = visuals.get(group.hashCode)
                    if (visualGroup) {
                        await visualGroup.visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, { group, structure })
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
                    const [ visual ] = visualsList[i]
                    await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props)
                }
            }
            if (structure) _structure = structure
            updated.next(updated.getValue() + 1)
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

    function setVisibility(value: boolean) {
        visuals.forEach(({ visual }) => {
            visual.setVisibility(value)
        })
    }

    function destroy() {
        visuals.forEach(({ visual }) => visual.destroy())
        visuals.clear()
    }

    return {
        label,
        get renderObjects() {
            const renderObjects: RenderObject[] = []
            visuals.forEach(({ visual }) => {
                if (visual.renderObject) renderObjects.push(visual.renderObject)
            })
            return renderObjects
        },
        get props() { return _props },
        get params() { return _params },
        get updated() { return updated },
        createOrUpdate,
        getLoci,
        mark,
        setVisibility,
        destroy
    }
}