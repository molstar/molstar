/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ComplexVisual, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationState } from './representation';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Structure, StructureElement, Bond } from '../../mol-model/structure';
import { Subject } from 'rxjs';
import { getNextMaterialId, GraphicsRenderObject } from '../../mol-gl/render-object';
import { Theme } from '../../mol-theme/theme';
import { Task } from '../../mol-task';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci, isEveryLoci, isDataLoci, EveryLoci } from '../../mol-model/loci';
import { MarkerAction, MarkerActions } from '../../mol-util/marker-action';
import { Overpaint } from '../../mol-theme/overpaint';
import { StructureParams } from './params';
import { Clipping } from '../../mol-theme/clipping';
import { Transparency } from '../../mol-theme/transparency';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Substance } from '../../mol-theme/substance';
import { LocationCallback } from '../util';
import { Emissive } from '../../mol-theme/emissive';

export function ComplexRepresentation<P extends StructureParams>(label: string, ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, P>, visualCtor: (materialId: number, structure: Structure, props: PD.Values<P>, webgl?: WebGLContext) => ComplexVisual<P>): StructureRepresentation<P> {
    let version = 0;
    const { webgl } = ctx;
    const updated = new Subject<number>();
    const geometryState = new Representation.GeometryState();
    const materialId = getNextMaterialId();
    const renderObjects: GraphicsRenderObject[] = [];
    const _state = StructureRepresentationStateBuilder.create();
    let visual: ComplexVisual<P> | undefined;

    let _structure: Structure;
    let _params: P;
    let _props: PD.Values<P>;
    let _theme = Theme.createEmpty();

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, structure?: Structure) {
        if (structure && structure !== _structure) {
            _params = getParams(ctx, structure);
            _structure = structure;
            if (!_props) _props = PD.getDefaultValues(_params);
        }
        _props = Object.assign({}, _props, props);

        return Task.create('Creating or updating ComplexRepresentation', async runtime => {
            let newVisual = false;
            if (!visual) {
                visual = visualCtor(materialId, _structure, _props, webgl);
                newVisual = true;
            } else if (visual.mustRecreate?.(_structure, _props, webgl)) {
                visual.destroy();
                visual = visualCtor(materialId, _structure, _props, webgl);
                newVisual = true;
            }
            const promise = visual.createOrUpdate({ webgl, runtime }, _theme, _props, structure);
            if (promise) await promise;
            if (newVisual) setState(_state); // current state for new visual
            // update list of renderObjects
            renderObjects.length = 0;
            if (visual && visual.renderObject) {
                renderObjects.push(visual.renderObject);
                geometryState.add(visual.renderObject.id, visual.geometryVersion);
            }
            geometryState.snapshot();
            // increment version
            version += 1;
            updated.next(version);
        });
    }

    function getLoci(pickingId: PickingId) {
        return visual ? visual.getLoci(pickingId) : EmptyLoci;
    }

    function getAllLoci() {
        return [Structure.Loci(_structure.child ?? _structure)];
    }

    function eachLocation(cb: LocationCallback) {
        visual?.eachLocation(cb);
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
        return visual ? visual.mark(loci, action) : false;
    }

    function setState(state: Partial<StructureRepresentationState>) {
        StructureRepresentationStateBuilder.update(_state, state);

        if (state.visible !== undefined && visual) {
            // hide visual when _unitTransforms is set and not the identity
            visual.setVisibility(state.visible && (_state.unitTransforms === null || _state.unitTransforms.isIdentity));
        }
        if (state.alphaFactor !== undefined && visual) visual.setAlphaFactor(state.alphaFactor);
        if (state.pickable !== undefined && visual) visual.setPickable(state.pickable);
        if (state.overpaint !== undefined && visual) {
            // Remap loci from equivalent structure to the current structure
            const remappedOverpaint = Overpaint.remap(state.overpaint, _structure);
            visual.setOverpaint(remappedOverpaint, webgl);
        }
        if (state.transparency !== undefined && visual) {
            // Remap loci from equivalent structure to the current structure
            const remappedTransparency = Transparency.remap(state.transparency, _structure);
            visual.setTransparency(remappedTransparency, webgl);
        }
        if (state.emissive !== undefined && visual) {
            // Remap loci from equivalent structure to the current structure
            const remappedEmissive = Emissive.remap(state.emissive, _structure);
            visual.setEmissive(remappedEmissive, webgl);
        }
        if (state.substance !== undefined && visual) {
            // Remap loci from equivalent structure to the current structure
            const remappedSubstance = Substance.remap(state.substance, _structure);
            visual.setSubstance(remappedSubstance, webgl);
        }
        if (state.clipping !== undefined && visual) {
            // Remap loci from equivalent structure to the current structure
            const remappedClipping = Clipping.remap(state.clipping, _structure);
            visual.setClipping(remappedClipping);
        }
        if (state.themeStrength !== undefined && visual) visual.setThemeStrength(state.themeStrength);
        if (state.transform !== undefined && visual) visual.setTransform(state.transform);
        if (state.unitTransforms !== undefined && visual) {
            // Since ComplexVisuals always renders geometries between units, the application
            // of `unitTransforms` does not make sense. When given here and not the identity,
            // it is ignored and sets the visual's visibility to `false`.
            visual.setVisibility(_state.visible && (state.unitTransforms === null || state.unitTransforms.isIdentity));
        }
    }

    function setTheme(theme: Theme) {
        _theme = theme;
    }

    function destroy() {
        if (visual) visual.destroy();
    }

    return {
        label,
        get groupCount() {
            return visual ? visual.groupCount : 0;
        },
        get props() { return _props; },
        get params() { return _params; },
        get state() { return _state; },
        get theme() { return _theme; },
        get geometryVersion() { return geometryState.version; },
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