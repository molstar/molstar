/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { shallowMerge2 } from '@molstar/core/util/object';
import {
  type ComplexVisual,
  type StructureRepresentation,
  StructureRepresentationStateBuilder,
  type StructureRepresentationState,
} from './representation.js';
import { Representation, type RepresentationContext, type RepresentationParamsGetter } from '../representation.js';
import { Structure, StructureElement, Bond } from '@molstar/model/model/structure';
import { Subject } from 'rxjs';
import { getNextMaterialId, type GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Theme } from '@molstar/graphics/theme/theme';
import { Task } from '@molstar/core/task';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { EmptyLoci, Loci, isEveryLoci, isDataLoci, EveryLoci } from '@molstar/model/model/loci';
import { type MarkerAction, MarkerActions } from '@molstar/core/util/marker-action';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import type { StructureParams } from './params.js';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { Transparency } from '@molstar/graphics/theme/transparency';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Substance } from '@molstar/graphics/theme/substance';
import type { LocationCallback } from '../util.js';
import { Emissive } from '@molstar/graphics/theme/emissive';
import { Wiggle } from '@molstar/graphics/theme/wiggle';
import type { SignatureContext } from '../signatures.js';
import { createVisual, getFactorySignatures, type VisualFactory } from '../visual-factory.js';

export function ComplexRepresentation<P extends StructureParams>(
  label: string,
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, P>,
  visualCtor: VisualFactory<Structure, P>,
): StructureRepresentation<P> {
  let version = 0;
  let { graphics } = ctx;
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

  /** updates one after another, see `Representation.createUpdateQueue` */
  const serialize = Representation.createUpdateQueue();
  function createOrUpdate(props: Partial<PD.Values<P>> = {}, structure?: Structure) {
    if (structure && structure !== _structure) {
      _params = getParams(ctx, structure);
      _structure = structure;
      if (!_props) _props = PD.getDefaultValues(_params);
    }
    // the same object while no prop changes: visuals skip what depends on the props alone
    _props = _props ? shallowMerge2(_props, props) : Object.assign({}, _props, props);

    const task = Task.create('Creating or updating ComplexRepresentation', async (runtime) => {
      let newVisual = false;
      if (!visual) {
        visual = createVisual(visualCtor, materialId, _structure, _props, graphics);
        newVisual = true;
      } else if (visual.mustRecreate?.(_structure, _props, graphics)) {
        visual.destroy();
        visual = createVisual(visualCtor, materialId, _structure, _props, graphics);
        newVisual = true;
      }
      const promise = visual.createOrUpdate({ graphics, runtime }, _theme, _props, structure);
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
    return serialize(task);
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
    if (state.pickable !== undefined && visual) visual.setPickable(state.pickable);
    if (state.overpaint !== undefined && visual) {
      // Remap loci from equivalent structure to the current structure
      const remappedOverpaint = Overpaint.remap(state.overpaint, _structure);
      visual.setOverpaint(remappedOverpaint, graphics);
    }
    if (state.transparency !== undefined && visual) {
      // Remap loci from equivalent structure to the current structure
      const remappedTransparency = Transparency.remap(state.transparency, _structure);
      visual.setTransparency(remappedTransparency, graphics);
    }
    if (state.emissive !== undefined && visual) {
      // Remap loci from equivalent structure to the current structure
      const remappedEmissive = Emissive.remap(state.emissive, _structure);
      visual.setEmissive(remappedEmissive, graphics);
    }
    if (state.substance !== undefined && visual) {
      // Remap loci from equivalent structure to the current structure
      const remappedSubstance = Substance.remap(state.substance, _structure);
      visual.setSubstance(remappedSubstance, graphics);
    }
    if (state.clipping !== undefined && visual) {
      // Remap loci from equivalent structure to the current structure
      const remappedClipping = Clipping.remap(state.clipping, _structure);
      visual.setClipping(remappedClipping);
    }
    if (state.wiggle !== undefined && visual) {
      // Remap loci from equivalent structure to the current structure
      const remappedWiggle = Wiggle.remap(state.wiggle, _structure);
      visual.setWiggle(remappedWiggle, graphics);
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
    get props() {
      return _props;
    },
    get params() {
      return _params;
    },
    get state() {
      return _state;
    },
    get theme() {
      return _theme;
    },
    get geometryVersion() {
      return geometryState.version;
    },
    renderObjects,
    updated,
    createOrUpdate,
    setState,
    setTheme,
    getLoci,
    getAllLoci,
    eachLocation,
    mark,
    setGraphics: (value: GraphicsContext) => {
      graphics = value;
    },
    getSignatures: (props: PD.Values<P>, ctx: SignatureContext) => getFactorySignatures(visualCtor, props, ctx),
    destroy,
  };
}
