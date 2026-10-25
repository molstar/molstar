/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { shallowMerge2 } from '@molstar/core/util/object';
import {
  type StructureRepresentation,
  StructureRepresentationStateBuilder,
  type StructureRepresentationState,
} from './representation.js';
import type { Visual } from '../visual.js';
import type { VisualContext } from '../visual.js';
import { isPromiseLike } from '@molstar/core/util/type-helpers';
import { Representation, type RepresentationContext, type RepresentationParamsGetter } from '../representation.js';
import { Structure, Unit, StructureElement, Bond, StructureSymmetry } from '@molstar/model/model/structure';
import { Subject } from 'rxjs';
import { getNextMaterialId, type GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Theme } from '@molstar/graphics/theme/theme';
import { type RuntimeContext, Task } from '@molstar/core/task';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { Loci, EmptyLoci, isEveryLoci, isDataLoci, EveryLoci } from '@molstar/model/model/loci';
import { MarkerAction, MarkerActions, applyMarkerAction } from '@molstar/core/util/marker-action';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { Transparency } from '@molstar/graphics/theme/transparency';
import { Mat4, EPSILON } from '@molstar/core/math/linear-algebra';
import { Interval } from '@molstar/core/data/int';
import type { StructureParams } from './params.js';
import { Clipping } from '@molstar/graphics/theme/clipping';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { StructureGroup } from './visual/util/common.js';
import { isWholeStructureLoci } from './visual/util/loci-index.js';
import { Substance } from '@molstar/graphics/theme/substance';
import type { LocationCallback } from '../util.js';
import { Emissive } from '@molstar/graphics/theme/emissive';
import { Wiggle } from '@molstar/graphics/theme/wiggle';
import { HashMap } from '@molstar/core/util/map';
import { hash2 } from '@molstar/core/data/util';
import type { GpuCullingValue, RenderMergedValue } from '@molstar/graphics/geo/geometry/base';
import { canCullOnGpu } from '@molstar/graphics/gl/renderable/merged';
import { mergeRenderObjectsByType, pickMergedLoci, type MergedEntry } from '../merged.js';
import { type SignatureContext, withMergedSignatures } from '../signatures.js';
import { createVisual, getFactorySignatures, type VisualFactory } from '../visual-factory.js';

function createVisualsMap<P extends StructureParams>() {
  return new HashMap<Unit.SymmetryGroup, { group: Unit.SymmetryGroup; visual: UnitsVisual<P> }>(
    (group) => hash2(group.hashCode, group.transformHash),
    Unit.SymmetryGroup.areInvariantElementsEqual,
  );
}

export interface UnitsVisual<P extends StructureParams> extends Visual<StructureGroup, P> {}

export function UnitsRepresentation<P extends StructureParams>(
  label: string,
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, P>,
  /** the visual of a symmetry group, which it may choose by the group (`selectVisual`) */
  visualCtor: VisualFactory<StructureGroup, P>,
): StructureRepresentation<P> {
  let version = 0;
  let { graphics } = ctx;
  const updated = new Subject<number>();
  const materialId = getNextMaterialId();
  const renderObjects: GraphicsRenderObject[] = [];
  const geometryState = new Representation.GeometryState();
  const _state = StructureRepresentationStateBuilder.create();
  let visuals = createVisualsMap<P>();

  let mergedByType = new Map<string, MergedEntry<UnitsVisual<P>>>();
  /** visuals whose render object is a member of a merged one: their layers are not smoothed (perf-audit 07 P1, N2) */
  const mergedVisuals = new Set<UnitsVisual<P>>();
  /** whether the layers of a visual were last set with smoothing, to set them again when its merging changes */
  const layersSmoothed = new WeakMap<UnitsVisual<P>, boolean>();
  /** visuals by the id of their render object, to resolve a pick without asking every visual */
  const visualsByObjectId = new Map<number, UnitsVisual<P>>();
  /**
   * The entries of `visuals` by the ids of their units, to mark only the visuals a small loci has (perf-audit 07 P3);
   * reset whenever an entry is added, an entry's visual replaced in place stays found
   */
  let visualsByUnitId: Map<number, { visual: UnitsVisual<P> }> | undefined;
  let visualCount = 0;

  function getVisualsByUnitId() {
    if (visualsByUnitId) return visualsByUnitId;
    visualsByUnitId = new Map();
    visualCount = 0;
    visuals.forEach((entry) => {
      ++visualCount;
      for (const u of entry.group.units) visualsByUnitId!.set(u.id, entry);
    });
    return visualsByUnitId;
  }

  let _structure: Structure;
  let _groups: ReadonlyArray<Unit.SymmetryGroup>;
  let _params: P;
  let _props: PD.Values<P>;
  let _theme = Theme.createEmpty();

  /** the visuals' context of the update running, one object for all of its visuals (a new one with new `graphics`) */
  let visualContext: VisualContext | undefined;
  const getVisualContext = (runtime: RuntimeContext) =>
    visualContext?.runtime === runtime && visualContext.graphics === graphics
      ? visualContext
      : (visualContext = { graphics, runtime });

  /**
   * `visual` updated, or a visual of `structureGroup` created when there is none or it must be recreated
   * (`mustRecreate`); `withData` false updates `visual` by the props and theme only, its group unchanged. A promise
   * only where the visual's update is one: a recolor of many units updates them all synchronously, without a promise
   * each (callers await only a promise)
   */
  function createOrUpdateVisual(
    runtime: RuntimeContext,
    visual: UnitsVisual<P> | undefined,
    structureGroup: StructureGroup,
    withData = true,
  ): UnitsVisual<P> | Promise<UnitsVisual<P>> {
    const ctx = getVisualContext(runtime);
    if (visual && !visual.mustRecreate?.(structureGroup, _props, graphics)) {
      const promise = visual.createOrUpdate(ctx, _theme, _props, withData ? structureGroup : undefined);
      return promise ? promise.then(() => visual) : visual;
    }
    visual?.destroy();
    const created = createVisual(visualCtor, materialId, structureGroup, _props, graphics);
    const finish = () => {
      setVisualState(created, structureGroup.group, _state, true); // current state for new visual
      return created;
    };
    const promise = created.createOrUpdate(ctx, _theme, _props, structureGroup);
    return promise ? promise.then(finish) : finish();
  }

  /** updates one after another, see `Representation.createUpdateQueue` */
  const serialize = Representation.createUpdateQueue();
  function createOrUpdate(props: Partial<PD.Values<P>> = {}, structure?: Structure) {
    if (structure && structure !== _structure) {
      _params = getParams(ctx, structure);
      if (!_props) _props = PD.getDefaultValues(_params);
    }
    // the same object while no prop changes: visuals skip what depends on the props alone
    _props = _props ? shallowMerge2(_props, props) : Object.assign({}, _props, props);

    const task = Task.create('Creating or updating UnitsRepresentation', async (runtime) => {
      if (!_structure && !structure) {
        throw new Error('missing structure');
      } else if (structure && !_structure) {
        // console.log(label, 'initial structure');
        // First call with a structure, create visuals for each group.
        _groups = StructureSymmetry.computeChildAwareTransformGroups(structure);
        for (let i = 0; i < _groups.length; i++) {
          const group = _groups[i];
          const created = createOrUpdateVisual(runtime, undefined, { group, structure });
          const visual = isPromiseLike(created) ? await created : created;
          visuals.set(group, { visual, group });
          visualsByUnitId = undefined;
          if (runtime.shouldUpdate)
            await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length });
        }
      } else if (
        structure &&
        (!Structure.areUnitIdsAndIndicesEqual(structure, _structure) || structure.child !== _structure.child)
      ) {
        // console.log(label, 'structures not equivalent');
        // Tries to re-use existing visuals for the groups of the new structure.
        // Creates additional visuals if needed, destroys left-over visuals.
        _groups = StructureSymmetry.computeChildAwareTransformGroups(structure);
        // const newGroups: Unit.SymmetryGroup[] = []
        const oldVisuals = visuals;
        visuals = createVisualsMap<P>();
        visualsByUnitId = undefined;
        for (let i = 0; i < _groups.length; i++) {
          const group = _groups[i];
          const visualGroup = oldVisuals.get(group);
          const updated = createOrUpdateVisual(runtime, visualGroup?.visual, { group, structure });
          const visual = isPromiseLike(updated) ? await updated : updated;
          visuals.set(group, { visual, group });
          visualsByUnitId = undefined;
          if (visualGroup) {
            oldVisuals.delete(group);

            // Remove highlight
            // TODO: remove selection too??
            if (visual.renderObject) {
              const { array, count } = visual.renderObject.values.tMarker.ref.value;
              applyMarkerAction(array, Interval.ofBounds(0, count), MarkerAction.RemoveHighlight);
            }
          }
          if (runtime.shouldUpdate)
            await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length });
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
        _groups = StructureSymmetry.computeChildAwareTransformGroups(structure);
        // console.log('new', _groups)
        // console.log('old', _structure.unitSymmetryGroups)
        for (let i = 0; i < _groups.length; i++) {
          const group = _groups[i];
          const visualGroup = visuals.get(group);
          if (visualGroup) {
            const updated = createOrUpdateVisual(runtime, visualGroup.visual, { group, structure });
            visualGroup.visual = isPromiseLike(updated) ? await updated : updated;
            visualGroup.group = group;
          } else {
            throw new Error(`expected to find visual for hashCode ${group.hashCode}`);
          }
          if (runtime.shouldUpdate)
            await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: _groups.length });
        }
      } else {
        // console.log(label, 'no new structure');
        // No new structure given, just update all visuals with new props.
        const visualsList: { group: Unit.SymmetryGroup; visual: UnitsVisual<P> }[] = []; // TODO avoid allocation
        visuals.forEach((vg) => visualsList.push(vg));
        for (let i = 0, il = visualsList.length; i < il; ++i) {
          const entry = visualsList[i];
          const structureGroup = { group: entry.group, structure: _structure };
          const updated = createOrUpdateVisual(runtime, entry.visual, structureGroup, false);
          entry.visual = isPromiseLike(updated) ? await updated : updated;
          if (runtime.shouldUpdate)
            await runtime.update({ message: 'Creating or updating UnitsVisual', current: i, max: il });
        }
      }
      // update list of renderObjects
      renderObjects.length = 0;
      const entries: { visual: UnitsVisual<P>; renderObject: GraphicsRenderObject }[] = [];
      visuals.forEach(({ visual }) => {
        if (visual.renderObject) {
          entries.push({ visual, renderObject: visual.renderObject });
          geometryState.add(visual.renderObject.id, visual.geometryVersion);
        }
      });
      updateRenderObjectList(entries);
      geometryState.snapshot();
      // set new structure
      if (structure) _structure = structure;
      // increment version
      updated.next(version++);
    });
    return serialize(task);
  }

  function updateRenderObjectList(entries: { visual: UnitsVisual<P>; renderObject: GraphicsRenderObject }[]) {
    // from BaseGeometry.Params, part of every units visual's props
    const renderMerged = (_props as { renderMerged?: RenderMergedValue }).renderMerged ?? false;
    const gpuCulling = (_props as { gpuCulling?: GpuCullingValue }).gpuCulling ?? false;
    mergedByType = mergeRenderObjectsByType(
      entries,
      renderMerged,
      mergedByType,
      renderObjects,
      gpuCulling,
      canCullOnGpu(ctx.graphics),
    );
    mergedVisuals.clear();
    mergedByType.forEach(({ visuals }) => {
      for (const v of visuals) mergedVisuals.add(v);
    });
    visualsByObjectId.clear();
    const { overpaint, transparency, emissive, substance } = _state;
    const hasLayers =
      overpaint.layers.length > 0 ||
      transparency.layers.length > 0 ||
      emissive.layers.length > 0 ||
      substance.layers.length > 0;
    for (const { visual, renderObject } of entries) {
      visualsByObjectId.set(renderObject.id, visual);
      // merged or unmerged now: the layers of new visuals, and those set under the previous merging
      if (hasLayers && (!layersSmoothed.has(visual) || layersSmoothed.get(visual) !== !mergedVisuals.has(visual))) {
        setLayers(visual, _state);
      }
    }
  }

  function getMergedLoci(pickingId: PickingId): Loci | undefined {
    return pickMergedLoci(mergedByType, pickingId, (visual, index, localInstanceId, object) =>
      visual.getLoci({ ...pickingId, objectId: object.members[index].id, instanceId: localInstanceId }),
    );
  }

  function getLoci(pickingId: PickingId) {
    const mergedLoci = getMergedLoci(pickingId);
    if (mergedLoci !== undefined) return mergedLoci;

    const visual = visualsByObjectId.get(pickingId.objectId);
    return visual ? visual.getLoci(pickingId) : EmptyLoci;
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
      if (Structure.isLoci(loci) || (StructureElement.Loci.is(loci) && isWholeStructureLoci(loci))) {
        // Change to `EveryLoci` to allow for downstream optimizations
        loci = EveryLoci;
      }
    } else if (!isEveryLoci(loci) && !isDataLoci(loci)) {
      return false;
    }
    if (Loci.isEmpty(loci)) return false;

    let changed = false;
    const marked = getMarkedVisuals(loci);
    if (marked) {
      for (const visual of marked) changed = visual.mark(loci, action) || changed;
    } else {
      visuals.forEach(({ visual }) => {
        changed = visual.mark(loci, action) || changed;
      });
    }
    return changed;
  }

  /**
   * The visuals of the units in `loci` when it has fewer units than there are visuals, `undefined` for all: the others
   * have none of it (their `lociApply` keys by unit), and a hover then costs the same with 10 or 10,000 visuals
   */
  function getMarkedVisuals(loci: Loci): Set<UnitsVisual<P>> | undefined {
    const byUnitId = getVisualsByUnitId();
    if (StructureElement.Loci.is(loci)) {
      if (loci.elements.length >= visualCount) return;
      const marked = new Set<UnitsVisual<P>>();
      for (const e of loci.elements) {
        const entry = byUnitId.get(e.unit.id);
        if (entry) marked.add(entry.visual);
      }
      return marked;
    }
    if (Bond.isLoci(loci)) {
      if (loci.bonds.length >= visualCount) return;
      const marked = new Set<UnitsVisual<P>>();
      for (const b of loci.bonds) {
        const a = byUnitId.get(b.aUnit.id);
        if (a) marked.add(a.visual);
        const c = byUnitId.get(b.bUnit.id);
        if (c) marked.add(c.visual);
      }
      return marked;
    }
  }

  /** the layers that can be smoothed: not on members of a merged object, which keeps one layer type for all */
  function setLayers(visual: UnitsVisual<P>, state: Partial<StructureRepresentationState>) {
    const { overpaint, transparency, emissive, substance } = state;
    if (!overpaint && !transparency && !emissive && !substance) return;
    const smooth = !mergedVisuals.has(visual);
    if (overpaint !== undefined) visual.setOverpaint(overpaint, graphics, smooth);
    if (transparency !== undefined) visual.setTransparency(transparency, graphics, smooth);
    if (emissive !== undefined) visual.setEmissive(emissive, graphics, smooth);
    if (substance !== undefined) visual.setSubstance(substance, graphics, smooth);
    layersSmoothed.set(visual, smooth);
  }

  /** `isNew`: a visual just made, at identity, with `state` the representation's current `_state` */
  function setVisualState(
    visual: UnitsVisual<P>,
    group: Unit.SymmetryGroup,
    state: Partial<StructureRepresentationState>,
    isNew = false,
  ) {
    const { visible, pickable, clipping, wiggle, themeStrength, transform, unitTransforms } = state;

    if (visible !== undefined) visual.setVisibility(visible);
    if (pickable !== undefined) visual.setPickable(pickable);
    // a new visual's layers wait for `updateRenderObjectList`: smoothed before the merging is decided, they would keep
    // the group from merging (perf-audit-2 07/P16)
    if (!isNew) setLayers(visual, state);
    if (clipping !== undefined) visual.setClipping(clipping);
    if (wiggle !== undefined) visual.setWiggle(wiggle, graphics);
    if (themeStrength !== undefined) visual.setThemeStrength(themeStrength);
    if (transform !== undefined) {
      // a new visual's transform is the identity, not `_state.transform` (perf-audit-2 07/B2)
      const changed = isNew
        ? !Mat4.isIdentity(transform)
        : transform !== _state.transform || !Mat4.areEqual(transform, _state.transform, EPSILON);
      if (changed) visual.setTransform(transform);
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
    const {
      visible,
      pickable,
      overpaint,
      transparency,
      emissive,
      substance,
      clipping,
      wiggle,
      themeStrength,
      transform,
      unitTransforms,
      syncManually,
      markerActions,
    } = state;
    const newState: Partial<StructureRepresentationState> = {};

    if (visible !== undefined) newState.visible = visible;
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
    if (wiggle !== undefined && _structure) {
      newState.wiggle = Wiggle.remap(wiggle, _structure);
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
    mergedByType.clear();
    mergedVisuals.clear();
    visualsByObjectId.clear();
    visualsByUnitId = undefined;
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
    get geometryVersion() {
      return geometryState.version;
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
    getSignatures: (props: PD.Values<P>, ctx: SignatureContext) => {
      // render objects of the same type are merged by `renderMerged` (groups of one stay unmerged unless culled on the
      // GPU); 'auto' merges above `AutoRenderMergedThreshold` units, taken as unmerged: such objects compile while they
      // build
      const p = props as { renderMerged?: RenderMergedValue; gpuCulling?: GpuCullingValue };
      return withMergedSignatures(getFactorySignatures(visualCtor, props, ctx), p);
    },
    destroy,
  };
}
