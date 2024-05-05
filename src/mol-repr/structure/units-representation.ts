/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationState } from './representation';
import { Visual } from '../visual';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Structure, Unit, StructureElement, Bond } from '../../mol-model/structure';
import { Subject } from 'rxjs';
import { getNextMaterialId, GraphicsRenderObject } from '../../mol-gl/render-object';
import { Theme } from '../../mol-theme/theme';
import { Task } from '../../mol-task';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci, EmptyLoci, isEmptyLoci, isEveryLoci, isDataLoci, EveryLoci } from '../../mol-model/loci';
import { MarkerAction, MarkerActions, applyMarkerAction } from '../../mol-util/marker-action';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { Mat4, EPSILON } from '../../mol-math/linear-algebra';
import { Interval } from '../../mol-data/int';
import { StructureParams } from './params';
import { Clipping } from '../../mol-theme/clipping';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { StructureGroup } from './visual/util/common';
import { Substance } from '../../mol-theme/substance';
import { LocationCallback } from '../util';
import { Emissive } from '../../mol-theme/emissive';

export interface UnitsVisual<P extends StructureParams> extends Visual<StructureGroup, P> { }

export function UnitsRepresentation<P extends StructureParams>(label: string, ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, P>, visualCtor: (materialId: number, structure: Structure, props: PD.Values<P>, webgl?: WebGLContext) => UnitsVisual<P>): StructureRepresentation<P> {
    let version = 0;
    const { webgl } = ctx;
    const updated = new Subject<number>();
    const materialId = getNextMaterialId();
    const renderObjects: GraphicsRenderObject[] = [];
    const geometryState = new Representation.GeometryState();
    const _state = StructureRepresentationStateBuilder.create();
    let visuals = new Map<number, { group: Unit.SymmetryGroup, visual: UnitsVisual<P> }>();

    let _structure: Structure;
    let _groups: ReadonlyArray<Unit.SymmetryGroup>;
    let _params: P;
    let _props: PD.Values<P>;
    let _theme = Theme.createEmpty();

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, structure?: Structure) {
        if (structure && structure !== _structure) {
            _params = getParams(ctx, structure);
            if (!_props) _props = PD.getDefaultValues(_params);
        }
        _props = Object.assign({}, _props, props);

        return Task.create('Creating or updating UnitsRepresentation', async runtime => {
            if (!_structure && !structure) {
                throw new Error('missing structure');
            } else if (structure && !_structure) {
                // console.log(label, 'initial structure');
                // First call with a structure, create visuals for each group.
                _groups = structure.unitSymmetryGroups;
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visual = visualCtor(materialId, structure, _props, webgl);
                    const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, { group, structure });
                    if (promise) await promise;
                    setVisualState(visual, group, _state); // current state for new visual
                    visuals.set(group.hashCode, { visual, group });
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length });
                }
            } else if (structure && (!Structure.areUnitIdsAndIndicesEqual(structure, _structure) || structure.child !== _structure.child)) {
                // console.log(label, 'structures not equivalent');
                // Tries to re-use existing visuals for the groups of the new structure.
                // Creates additional visuals if needed, destroys left-over visuals.
                _groups = structure.unitSymmetryGroups;
                // const newGroups: Unit.SymmetryGroup[] = []
                const oldVisuals = visuals;
                visuals = new Map();
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visualGroup = oldVisuals.get(group.hashCode);
                    if (visualGroup) {
                        // console.log(label, 'found visualGroup to reuse');
                        // console.log('old', visualGroup.group)
                        // console.log('new', group)
                        let { visual } = visualGroup;
                        if (visual.mustRecreate?.({ group, structure }, _props, webgl)) {
                            visual.destroy();
                            visual = visualCtor(materialId, structure, _props, webgl);
                            const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, { group, structure });
                            if (promise) await promise;
                            setVisualState(visual, group, _state); // current state for new visual
                        } else {
                            const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, { group, structure });
                            if (promise) await promise;
                        }
                        visuals.set(group.hashCode, { visual, group });
                        oldVisuals.delete(group.hashCode);

                        // Remove highlight
                        // TODO: remove selection too??
                        if (visual.renderObject) {
                            const arr = visual.renderObject.values.tMarker.ref.value.array;
                            applyMarkerAction(arr, Interval.ofBounds(0, arr.length), MarkerAction.RemoveHighlight);
                        }
                    } else {
                        // console.log(label, 'did not find visualGroup to reuse, creating new');
                        // newGroups.push(group)
                        const visual = visualCtor(materialId, structure, _props, webgl);
                        const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, { group, structure });
                        if (promise) await promise;
                        setVisualState(visual, group, _state); // current state for new visual
                        visuals.set(group.hashCode, { visual, group });
                    }
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length });
                }
                oldVisuals.forEach(({ visual }) => {
                    // console.log(label, 'removed unused visual');
                    visual.destroy();
                });
            } else if (structure && structure !== _structure && Structure.areUnitIdsAndIndicesEqual(structure, _structure)) {
                // console.log(label, 'structures equivalent but not identical');
                // Expects that for structures with the same hashCode,
                // the unitSymmetryGroups are the same as well.
                // Re-uses existing visuals for the groups of the new structure.
                _groups = structure.unitSymmetryGroups;
                // console.log('new', structure.unitSymmetryGroups)
                // console.log('old', _structure.unitSymmetryGroups)
                for (let i = 0; i < _groups.length; i++) {
                    const group = _groups[i];
                    const visualGroup = visuals.get(group.hashCode);
                    if (visualGroup) {
                        let { visual } = visualGroup;
                        if (visual.mustRecreate?.({ group, structure }, _props, ctx.webgl)) {
                            visual.destroy();
                            visual = visualCtor(materialId, structure, _props, ctx.webgl);
                            visualGroup.visual = visual;
                            const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, { group, structure });
                            if (promise) await promise;
                            setVisualState(visual, group, _state); // current state for new visual
                        } else {
                            const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, { group, structure });
                            if (promise) await promise;
                        }
                        visualGroup.group = group;
                    } else {
                        throw new Error(`expected to find visual for hashCode ${group.hashCode}`);
                    }
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length });
                }
            } else {
                // console.log(label, 'no new structure');
                // No new structure given, just update all visuals with new props.
                const visualsList: { group: Unit.SymmetryGroup, visual: UnitsVisual<P> }[] = []; // TODO avoid allocation
                visuals.forEach(vg => visualsList.push(vg));
                for (let i = 0, il = visualsList.length; i < il; ++i) {
                    let { visual, group } = visualsList[i];
                    if (visual.mustRecreate?.({ group, structure: _structure }, _props, ctx.webgl)) {
                        visual.destroy();
                        visual = visualCtor(materialId, _structure, _props, webgl);
                        visualsList[i].visual = visual;
                        const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, { group, structure: _structure });
                        if (promise) await promise;
                        setVisualState(visual, group, _state); // current state for new visual
                    } else {
                        const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props);
                        if (promise) await promise;
                    }
                    if (runtime.shouldUpdate) await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: il });
                }
            }
            // update list of renderObjects
            renderObjects.length = 0;
            visuals.forEach(({ visual }) => {
                if (visual.renderObject) {
                    renderObjects.push(visual.renderObject);
                    geometryState.add(visual.renderObject.id, visual.geometryVersion);
                }
            });
            geometryState.snapshot();
            // set new structure
            if (structure) _structure = structure;
            // increment version
            updated.next(version++);
        });
    }

    function getLoci(pickingId: PickingId) {
        let loci: Loci = EmptyLoci;
        visuals.forEach(({ visual }) => {
            const _loci = visual.getLoci(pickingId);
            if (!isEmptyLoci(_loci)) loci = _loci;
        });
        return loci;
    }

    function eachLocation(cb: LocationCallback) {
        visuals.forEach(({ visual }) => {
            visual.eachLocation(cb);
        });
    }

    function getAllLoci() {
        return [Structure.Loci(_structure.child ?? _structure)];
    }

    function mark(loci: Loci, action: MarkerAction) {
        if (!_structure) return false;
        if (!MarkerActions.is(_state.markerActions, action)) return false;
        if (Structure.isLoci(loci) || StructureElement.Loci.is(loci) || Bond.isLoci(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, _structure)) return false;
            // Remap `loci` from equivalent structure to the current `_structure`
            loci = Loci.remap(loci, _structure);
            if (Structure.isLoci(loci) || (StructureElement.Loci.is(loci) && StructureElement.Loci.isWholeStructure(loci))) {
                // Change to `EveryLoci` to allow for downstream optimizations
                loci = EveryLoci;
            }
        } else if (!isEveryLoci(loci) && !isDataLoci(loci)) {
            return false;
        }
        if (Loci.isEmpty(loci)) return false;

        let changed = false;
        visuals.forEach(({ visual }) => {
            changed = visual.mark(loci, action) || changed;
        });
        return changed;
    }

    function setVisualState(visual: UnitsVisual<P>, group: Unit.SymmetryGroup, state: Partial<StructureRepresentationState>) {
        const { visible, alphaFactor, pickable, overpaint, transparency, emissive, substance, clipping, themeStrength, transform, unitTransforms } = state;

        if (visible !== undefined) visual.setVisibility(visible);
        if (alphaFactor !== undefined) visual.setAlphaFactor(alphaFactor);
        if (pickable !== undefined) visual.setPickable(pickable);
        if (overpaint !== undefined) visual.setOverpaint(overpaint, webgl);
        if (transparency !== undefined) visual.setTransparency(transparency, webgl);
        if (emissive !== undefined) visual.setEmissive(emissive, webgl);
        if (substance !== undefined) visual.setSubstance(substance, webgl);
        if (clipping !== undefined) visual.setClipping(clipping);
        if (themeStrength !== undefined) visual.setThemeStrength(themeStrength);
        if (transform !== undefined) {
            if (transform !== _state.transform || !Mat4.areEqual(transform, _state.transform, EPSILON)) {
                visual.setTransform(transform);
            }
        }
        if (unitTransforms !== undefined) {
            if (unitTransforms) {
                // console.log(group.hashCode, unitTransforms.getSymmetryGroupTransforms(group))
                visual.setTransform(undefined, unitTransforms.getSymmetryGroupTransforms(group));
            } else if (unitTransforms !== _state.unitTransforms) {
                visual.setTransform(undefined, null);
            }
        }
    }

    function setState(state: Partial<StructureRepresentationState>) {
        const { visible, alphaFactor, pickable, overpaint, transparency, emissive, substance, clipping, themeStrength, transform, unitTransforms, syncManually, markerActions } = state;
        const newState: Partial<StructureRepresentationState> = {};

        if (visible !== undefined) newState.visible = visible;
        if (alphaFactor !== undefined) newState.alphaFactor = alphaFactor;
        if (pickable !== undefined) newState.pickable = pickable;
        if (overpaint !== undefined && _structure) {
            newState.overpaint = Overpaint.remap(overpaint, _structure);
        }
        if (transparency !== undefined && _structure) {
            newState.transparency = Transparency.remap(transparency, _structure);
        }
        if (emissive !== undefined && _structure) {
            newState.emissive = Emissive.remap(emissive, _structure);
        }
        if (substance !== undefined && _structure) {
            newState.substance = Substance.remap(substance, _structure);
        }
        if (clipping !== undefined && _structure) {
            newState.clipping = Clipping.remap(clipping, _structure);
        }
        if (themeStrength !== undefined) newState.themeStrength = themeStrength;
        if (transform !== undefined && !Mat4.areEqual(transform, _state.transform, EPSILON)) {
            newState.transform = transform;
        }
        if (unitTransforms !== _state.unitTransforms || unitTransforms?.version !== _state.unitTransformsVersion) {
            newState.unitTransforms = unitTransforms;
            _state.unitTransformsVersion = unitTransforms ? unitTransforms?.version : -1;
        }
        if (syncManually !== undefined) newState.syncManually = syncManually;
        if (markerActions !== undefined) newState.markerActions = markerActions;

        visuals.forEach(({ visual, group }) => setVisualState(visual, group, newState));

        StructureRepresentationStateBuilder.update(_state, newState);
    }

    function setTheme(theme: Theme) {
        _theme = theme;
    }

    function destroy() {
        visuals.forEach(({ visual }) => visual.destroy());
        visuals.clear();
    }

    return {
        label,
        get groupCount() {
            let groupCount = 0;
            visuals.forEach(({ visual }) => {
                if (visual.renderObject) groupCount += visual.groupCount;
            });
            return groupCount;
        },
        get geometryVersion() { return geometryState.version; },
        get props() { return _props; },
        get params() { return _params; },
        get state() { return _state; },
        get theme() { return _theme; },
        renderObjects,
        updated,
        createOrUpdate,
        setState,
        setTheme,
        getLoci,
        getAllLoci,
        eachLocation,
        mark,
        destroy
    };
}