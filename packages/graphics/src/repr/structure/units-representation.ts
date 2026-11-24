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
import { Visual } from '../visual.js';
import type { VisualContext } from '../visual.js';
import { isPromiseLike } from '@molstar/core/util/type-helpers';
import { Representation, type RepresentationContext, type RepresentationParamsGetter } from '../representation.js';
import { Structure, Unit, StructureElement, Bond, StructureSymmetry } from '@molstar/model/model/structure';
import { Subject } from 'rxjs';
import { getNextMaterialId, type GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Theme } from '@molstar/graphics/theme/theme';
import { type RuntimeContext, Task } from '@molstar/core/task';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { getLodPixelLevelCount } from '@molstar/graphics/geo/geometry/lod-levels';
import { Loci, EmptyLoci, isEveryLoci, isDataLoci, EveryLoci } from '@molstar/model/model/loci';
import { MarkerAction, MarkerActions } from '@molstar/core/util/marker-action';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { Transparency } from '@molstar/graphics/theme/transparency';
import { Mat4, EPSILON } from '@molstar/core/math/linear-algebra';
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
import { type GpuCullingValue, type RenderMergedValue, resolveRenderMerged } from '@molstar/graphics/geo/geometry/base';
import {
  canCullOnGpu,
  canMergeValues,
  isMergeableType,
  type MergeableValues,
} from '@molstar/graphics/gl/renderable/merged';
import { mergeRenderObjectsByType, pickMergedLoci, type MergedEntry } from '../merged.js';
import { type SignatureContext, withMergedSignatures } from '../signatures.js';
import {
  choosesByProps,
  createVisual,
  getFactorySignatures,
  getKindFactory,
  getVisualKind,
  type VisualFactory,
} from '../visual-factory.js';
import type { UnitsVisual as UnitsBuilderVisual, UnitsVisualGeometryBuilder } from './units-visual.js';
import { UnitsMergedVisual } from './units-merged-visual.js';

function createVisualsMap<P extends StructureParams>() {
  // `kind`: the visual's kind (`getVisualKind`)
  return new HashMap<Unit.SymmetryGroup, { group: Unit.SymmetryGroup; visual: UnitsVisual<P>; kind: string }>(
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

  /**
   * The visuals over the groups of a kind where they merge (`UnitsMergedVisual`; a kind the factory's choice of visual
   * for a group, `getVisualKind`: one for a factory that chooses by the props alone), per-group visuals for the others
   * (a kind of a type that does not merge, e.g. GPU surfaces' texture meshes, or one found not to). A kind found not to
   * merge (`unmergeableAt`) is tried again only where its per-group visuals' own values merge after an update
   * (`mergesAgain`), not by a geometry made to find out, and not at the props, theme and structure it failed at
   */
  const mergedByKind = new Map<string, UnitsMergedVisual<P>>();
  const unmergeableAt = new Map<string, { props: PD.Values<P>; theme: Theme; structure: Structure }>();
  const failedAt = (kind: string, structure: Structure) => {
    const at = unmergeableAt.get(kind);
    return !!at && at.props === _props && at.theme === _theme && at.structure === structure;
  };
  /**
   * each kind's groups of the last update, the same array while they are the same groups (a merged visual's update of
   * the same groups takes its fast paths)
   */
  let groupsByKind = new Map<string, ReadonlyArray<Unit.SymmetryGroup>>();

  /** the groups of `groups` by their kind (`getVisualKind`), each kind's the last update's array where the same */
  function getGroupsByKind(structure: Structure, groups: ReadonlyArray<Unit.SymmetryGroup>) {
    const byKind = new Map<string, Unit.SymmetryGroup[]>();
    if (choosesByProps(visualCtor) && groups.length > 0) {
      // one kind for all: the groups themselves
      byKind.set(getVisualKind(visualCtor, { group: groups[0], structure }, _props, graphics), groups as any);
    } else {
      for (const group of groups) {
        const kind = getVisualKind(visualCtor, { group, structure }, _props, graphics);
        let list = byKind.get(kind);
        if (!list) byKind.set(kind, (list = []));
        list.push(group);
      }
    }
    const out = new Map<string, ReadonlyArray<Unit.SymmetryGroup>>();
    byKind.forEach((list, kind) => {
      const before = groupsByKind.get(kind);
      const same = !!before && before.length === list.length && before.every((g, i) => g === list[i]);
      out.set(kind, same ? before! : list);
    });
    groupsByKind = out;
    return out;
  }

  function hasVisualsOfKind(kind: string) {
    let has = false;
    visuals.forEach((entry) => {
      if (entry.kind === kind) has = true;
    });
    return has;
  }

  /**
   * The builder of a merged visual over `groups`, of `kind`: where `renderMerged` merges them (two or more) and the
   * kind's visual is of a builder of a type that merges
   */
  function getMergedBuilder(
    kind: string,
    structure: Structure,
    groups: ReadonlyArray<Unit.SymmetryGroup>,
  ): UnitsVisualGeometryBuilder<P, any> | undefined {
    if (groups.length < 2) return;
    const renderMerged = (_props as { renderMerged?: RenderMergedValue }).renderMerged ?? false;
    if (!resolveRenderMerged(renderMerged, groups.length)) return;
    // found not to merge: the per-group visuals decide after their update (`mergesAgain`)
    if (unmergeableAt.has(kind) && (hasVisualsOfKind(kind) || failedAt(kind, structure))) return;
    const probe = getKindFactory(visualCtor, kind)(-1) as UnitsBuilderVisual<P>;
    const builder = probe.builder as UnitsVisualGeometryBuilder<P, any> | undefined;
    if (!builder || !isMergeableType(builder.geometryUtils.createEmpty().kind)) return;
    return builder;
  }

  /**
   * The builder of a merged visual over `groups` of `kind` where its per-group visuals, made since a merged one did not
   * merge, now do by their own values (`canMergeValues` of each, as the merged visual checks its parts): undefined
   * where they do not, or where the merged one failed at these props and structure
   */
  function mergesAgain(kind: string, structure: Structure, groups: ReadonlyArray<Unit.SymmetryGroup>) {
    if (!unmergeableAt.has(kind) || failedAt(kind, structure)) return;
    if (groups.length < 2) return;
    const renderMerged = (_props as { renderMerged?: RenderMergedValue }).renderMerged ?? false;
    if (!resolveRenderMerged(renderMerged, groups.length)) return;
    let builder: UnitsVisualGeometryBuilder<P, any> | undefined;
    let type = '';
    let merges = true;
    visuals.forEach(({ visual, kind: k }) => {
      if (!merges || k !== kind) return;
      const b = (visual as UnitsBuilderVisual<P>).builder as UnitsVisualGeometryBuilder<P, any> | undefined;
      if (!b) {
        merges = false;
        return;
      }
      if (!builder) {
        builder = b;
        type = b.geometryUtils.createEmpty().kind;
      } else if (b.createGeometry !== builder.createGeometry || b.geometryUtils !== builder.geometryUtils) {
        merges = false;
        return;
      }
      // a group without an object (nothing drawn) none to check; the layers not: smoothed while not merged
      // (`setLayers`), set again unsmoothed by the merged visual
      const ro = visual.renderObject;
      if (ro && !canMergeValues(type, [ro.values as unknown as MergeableValues], false)) merges = false;
    });
    return merges ? builder : undefined;
  }

  /**
   * The merged visual of `kind` updated to `structure` and its `groups`; false where they do not merge (the per-group
   * visuals instead)
   */
  async function createOrUpdateMerged(
    runtime: RuntimeContext,
    kind: string,
    builder: UnitsVisualGeometryBuilder<P, any>,
    structure: Structure,
    groups: ReadonlyArray<Unit.SymmetryGroup>,
  ) {
    let mergedVisual = mergedByKind.get(kind);
    if (
      mergedVisual &&
      (mergedVisual.builder.createGeometry !== builder.createGeometry ||
        mergedVisual.builder.geometryUtils !== builder.geometryUtils)
    ) {
      mergedVisual.destroy();
      mergedByKind.delete(kind);
      mergedVisual = undefined;
    }
    const changed = !!mergedVisual && !!_structure && !Structure.areUnitIdsAndIndicesEqual(structure, _structure);
    // a new merged visual published (`mergedByKind`: picked and marked) once its first update merged, the per-group
    // visuals drawn and asked till then; one that throws destroyed
    const created = !mergedVisual;
    // the builder's props those of the representation's visuals (`P`), typed by its geometry's
    const mv =
      mergedVisual ??
      (UnitsMergedVisual(builder as UnitsVisualGeometryBuilder<any, any>, materialId) as UnitsMergedVisual<P>);
    if (created) mv.setState(_state, _state);
    try {
      await mv.createOrUpdate(getVisualContext(runtime), _theme, _props, structure, groups);
    } catch (e) {
      if (created) mv.destroy();
      throw e;
    }
    if (!mv.mergeable) {
      // the per-group visuals kept (updated, not made again): the merged one ends at a part found not to merge
      mv.destroy();
      mergedByKind.delete(kind);
      unmergeableAt.set(kind, { props: _props, theme: _theme, structure });
      return false;
    }
    mergedByKind.set(kind, mv);
    unmergeableAt.delete(kind);
    // remove highlight, as the per-group visuals do on a new structure
    if (changed) mv.mark(EveryLoci, MarkerAction.RemoveHighlight);
    return true;
  }

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
      if (!_structure && !structure) throw new Error('missing structure');
      const target = structure ?? _structure;
      const groups =
        structure && structure !== _structure ? StructureSymmetry.computeChildAwareTransformGroups(structure) : _groups;
      const byKind = getGroupsByKind(target, groups);

      // the kinds that merge: their merged visuals
      const merged = new Set<string>();
      for (const [kind, kindGroups] of byKind) {
        const builder = getMergedBuilder(kind, target, kindGroups);
        if (builder && (await createOrUpdateMerged(runtime, kind, builder, target, kindGroups))) merged.add(kind);
      }
      for (const [kind, mv] of [...mergedByKind]) {
        if (merged.has(kind)) continue;
        mv.destroy();
        mergedByKind.delete(kind);
      }

      // the groups of the other kinds: a visual each, of the previous one of the group where there is one
      const newStructure = !!structure && !!_structure && structure !== _structure;
      const equivalent =
        newStructure &&
        Structure.areUnitIdsAndIndicesEqual(structure!, _structure) &&
        structure!.child === _structure.child;
      const oldVisuals = visuals;
      visuals = createVisualsMap<P>();
      visualsByUnitId = undefined;
      let done = 0;
      const total = groups.length;
      for (const [kind, kindGroups] of byKind) {
        if (merged.has(kind)) continue;
        for (const group of kindGroups) {
          const entry = oldVisuals.get(group);
          // the same structure: by the props and theme alone (its group as it is); another: its group again
          const updated = createOrUpdateVisual(
            runtime,
            entry?.visual,
            { group, structure: target },
            !entry || newStructure,
          );
          const visual = isPromiseLike(updated) ? await updated : updated;
          visuals.set(group, { visual, group, kind });
          if (entry) {
            oldVisuals.delete(group);
            // Remove highlight, through `mark`: its status and average follow (they stood for the layer before)
            if (newStructure && !equivalent) visual.mark(EveryLoci, MarkerAction.RemoveHighlight);
          }
          if (runtime.shouldUpdate) {
            await runtime.update({ message: 'Creating or updating UnitsVisual', current: ++done, max: total });
          }
        }
      }
      oldVisuals.forEach(({ visual }) => visual.destroy());

      // a kind's per-group visuals' values merging again (by a theme or props that now do): its merged visual instead
      for (const [kind, kindGroups] of byKind) {
        if (merged.has(kind)) continue;
        const againBuilder = mergesAgain(kind, target, kindGroups);
        if (againBuilder && (await createOrUpdateMerged(runtime, kind, againBuilder, target, kindGroups))) {
          merged.add(kind);
          visuals.forEach((entry) => {
            if (entry.kind !== kind) return;
            entry.visual.destroy();
            visuals.delete(entry.group);
          });
          visualsByUnitId = undefined;
        }
      }

      // update list of renderObjects: the per-group visuals' (merged by type where they merge), the merged visuals'
      const entries: { visual: UnitsVisual<P>; renderObject: GraphicsRenderObject }[] = [];
      visuals.forEach(({ visual }) => {
        if (visual.renderObject) {
          entries.push({ visual, renderObject: visual.renderObject });
          geometryState.add(visual.renderObject.id, visual.geometryVersion);
        }
      });
      updateRenderObjectList(entries);
      for (const mv of mergedByKind.values()) {
        for (const ro of mv.renderObjects) {
          renderObjects.push(ro);
          geometryState.add(ro.id, mv.geometryVersion);
          if (_state.previewTransform) Visual.setPreviewTransform(ro, _state.previewTransform);
        }
      }
      geometryState.snapshot();
      _groups = groups;
      _structure = target;
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
      canCullOnGpu(graphics),
      graphics,
      _state.previewTransform,
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
    return pickMergedLoci(mergedByType, pickingId, (visual, _index, localInstanceId, objectId) =>
      visual.getLoci({ ...pickingId, objectId, instanceId: localInstanceId }),
    );
  }

  function getLoci(pickingId: PickingId) {
    for (const mv of mergedByKind.values()) {
      if (mv.renderObjects.some((o) => o.id === pickingId.objectId)) return mv.getLoci(pickingId);
    }
    const mergedLoci = getMergedLoci(pickingId);
    if (mergedLoci !== undefined) return mergedLoci;

    const visual = visualsByObjectId.get(pickingId.objectId);
    return visual ? visual.getLoci(pickingId) : EmptyLoci;
  }

  function eachLocation(cb: LocationCallback) {
    for (const mv of mergedByKind.values()) mv.eachLocation(cb);
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
    for (const mv of mergedByKind.values()) changed = mv.mark(loci, action) || changed;
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
      previewTransform,
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
    if (previewTransform !== undefined && !Mat4.areEqual(previewTransform, _state.previewTransform, 0)) {
      newState.previewTransform = previewTransform;
    }
    if (unitTransforms !== _state.unitTransforms || unitTransforms?.version !== _state.unitTransformsVersion) {
      newState.unitTransforms = unitTransforms;
      _state.unitTransformsVersion = unitTransforms ? unitTransforms?.version : -1;
    }
    if (syncManually !== undefined) newState.syncManually = syncManually;
    if (markerActions !== undefined) newState.markerActions = markerActions;

    visuals.forEach(({ visual, group }) => setVisualState(visual, group, newState));
    for (const mv of mergedByKind.values()) mv.setState(newState, _state);
    // on what is drawn: a merged object moves as one, its members' are not read (`mergeRenderObjectsByType`)
    const preview = newState.previewTransform;
    if (preview) for (const ro of renderObjects) Visual.setPreviewTransform(ro, preview);

    StructureRepresentationStateBuilder.update(_state, newState);
  }

  function setTheme(theme: Theme) {
    _theme = theme;
  }

  function destroy() {
    for (const mv of mergedByKind.values()) mv.destroy();
    mergedByKind.clear();
    unmergeableAt.clear();
    groupsByKind.clear();
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
      for (const mv of mergedByKind.values()) groupCount += mv.groupCount;
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
    // `mark` takes structure loci of an equivalent root only
    getStructureRootHashes: () => (_structure ? [_structure.root.hashCode] : []),
    setGraphics: (value: GraphicsContext) => {
      graphics = value;
    },
    getSignatures: (props: PD.Values<P>, ctx: SignatureContext) => {
      // render objects of the same type are merged by `renderMerged` (groups of one stay unmerged unless culled on the
      // GPU); 'auto' merges above `AutoRenderMergedThreshold` units, taken as unmerged: such objects compile while they
      // build
      const p = props as { renderMerged?: RenderMergedValue; gpuCulling?: GpuCullingValue };
      if (!(ctx.data instanceof Structure))
        return withMergedSignatures(getFactorySignatures(visualCtor, props, ctx), p);
      // with the structure: the visuals its groups select, merged as their counts decide
      const structure = ctx.data;
      const groups = StructureSymmetry.computeChildAwareTransformGroups(structure);
      const data = groups.map((group) => ({ structure, group }));
      const lodPixelSize = (props as { lodPixelSize?: number }).lodPixelSize ?? 0;
      const counts = {
        objects: groups.length,
        instances: structure.units.length,
        spheres: structure.elementCount,
        levels: getLodPixelLevelCount(lodPixelSize),
        graphics: ctx.graphics,
      };
      return withMergedSignatures(getFactorySignatures(visualCtor, props, ctx, data), p, counts);
    },
    destroy,
  };
}
