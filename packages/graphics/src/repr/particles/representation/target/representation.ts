/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  Representation,
  type RepresentationContext,
  type RepresentationParamsGetter,
  type RepresentationProvider,
} from '../../../representation.js';
import {
  type ParticleList,
  Particle,
  type ParticleTarget,
  getParticleTargetGroups,
} from '@molstar/model/model/particles/particle-list';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Theme } from '@molstar/graphics/theme/theme';
import { type GraphicsRenderObject, getNextMaterialId } from '@molstar/graphics/gl/render-object';
import { Subject } from 'rxjs';
import { type RuntimeContext, Task } from '@molstar/core/task';
import { type Loci as ModelLoci, EmptyLoci, isEveryLoci, EveryLoci } from '@molstar/model/model/loci';
import { type MarkerAction, MarkerActions } from '@molstar/core/util/marker-action';
import type { LocationCallback } from '../../../util.js';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { getLodPixelLevelCount } from '@molstar/graphics/geo/geometry/lod-levels';
import { OrderedSet } from '@molstar/core/data/int';
import { mergeRenderObjectsByType, pickMergedLoci, type MergedEntry } from '../../../merged.js';
import { Visual } from '../../../visual.js';
import { canCullOnGpu } from '@molstar/graphics/gl/renderable/merged';
import {
  ParticleTargetRepresentationParams,
  type ParticleTargetRepresentationProps,
  createTargetVisual,
  createParticleTargetPropsProvider,
  getTargetGeometryChoice,
  getTargetGeometryUtils,
  isPosesGeometry,
  type ParticleTargetProps,
  type TargetVisual,
  isTargetOccluderOff,
} from './common.js';
import type { ObjectSignature } from '@molstar/graphics/gl/precompile';
import {
  getGeometrySignatures,
  getSurfaceSignatureDefines,
  type SignatureContext,
  type SignatureDefines,
  uniqueSignatures,
  withMergedSignatures,
} from '../../../signatures.js';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { getSitesSignatureDefines } from './sites.js';
import { getStructureSignatureDefines, type SiteColorProps } from './site-color.js';
import { createTargetPagedObject, type TargetSlotsObject } from './paged.js';
import { ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';
import { ParticleTiles } from '@molstar/model/model/particles/particle-tiles';
import { ParticleDensityGrid } from '@molstar/model/model/particles/particle-tile-opaque';
import { type OccluderGrid, setOccluderGrid, usesOccluderGrid } from '@molstar/graphics/gl/occluder-grid';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';

export { ParticleTargetRepresentationParams };
export type { ParticleTargetRepresentationProps };

export interface ParticleTargetRepresentation
  extends Representation<ParticleList, ParticleTargetRepresentationParams> {}

export type ParticleTargetRepresentationProvider<Id extends string = string> = RepresentationProvider<
  ParticleList,
  ParticleTargetRepresentationParams,
  Representation.State,
  Id
>;

export function ParticleTargetRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<ParticleList, ParticleTargetRepresentationParams>,
): ParticleTargetRepresentation {
  let { graphics } = ctx;
  const updated = new Subject<number>();
  const geometryState = new Representation.GeometryState();
  const renderObjects: GraphicsRenderObject[] = [];
  const _state = Representation.createState();
  let _theme = Theme.createEmpty();

  // shared across all target visuals of this representation so mergeable ones can actually be merged
  // (canMergeRenderObjects requires equal materialId; a per-visual id would never match)
  const materialId = getNextMaterialId();

  let version = 0;
  let _particles: ParticleList | undefined;
  let _params: ParticleTargetRepresentationParams;
  let _props: ParticleTargetRepresentationProps;
  // Map from targetId → TargetVisual; re-used across updates.
  const targetVisuals = new Map<number, TargetVisual>();

  let mergedByType = new Map<string, MergedEntry<TargetVisual>>();
  /**
   * a slot list's object on GPU culling (`paged.ts`): its instances a page store's, culled and composed on the GPU,
   * updated by the list's changes
   */
  let slots: TargetSlotsObject | undefined;
  let slotsVersion = 0;
  /** what the visuals draw: the list, or a slot list flattened where not culled on the GPU (its keys its slots) */
  let drawn: ParticleList | undefined;

  /** What the last full update built from: a later one changing the clip alone sets it on the visuals (`isClipOnly`) */
  let built: { particles: ParticleList; theme: Theme; props: ParticleTargetRepresentationProps } | undefined;

  /**
   * Whether the update to `_props` changes the clip alone, its defines kept: a clip object moved (a clipping gizmo's
   * drag), nothing to rebuild and nothing merged otherwise; every target's full update and the merge check are skipped
   */
  function isClipOnly() {
    if (!built || built.particles !== _particles || built.theme !== _theme || targetVisuals.size === 0) return false;
    const before = built.props;
    if (before.clip === _props.clip) return false;
    for (const k in _props) {
      if (k !== 'clip' && (_props as any)[k] !== (before as any)[k]) return false;
    }
    return BaseGeometry.isSameClipDefines(before.clip, _props.clip);
  }

  /** updates one after another, see `Representation.createUpdateQueue` */
  const serialize = Representation.createUpdateQueue();
  function createOrUpdate(props: Partial<ParticleTargetRepresentationProps> = {}, input?: ParticleList) {
    if (input) {
      _particles = input;
      _params = getParams(ctx, _particles);
      if (!_props) _props = PD.getDefaultValues(_params);
    }
    // Quality-derived props are resolved per target (see `createParticleTargetPropsProvider`),
    // since each target has its own data and its own geometry params.
    Object.assign(_props, props);

    if (isClipOnly()) {
      const clip = _props.clip;
      built!.props = { ..._props };
      return serialize(
        Task.create('ParticleTargetRepresentation clip', async () => {
          // every member, also of a merged object (its first member's cells are the merged object's): one leaving its
          // group draws right
          for (const visual of targetVisuals.values()) {
            if (visual.renderObject) BaseGeometry.updateClip(visual.renderObject.values as any, clip);
          }
          updated.next(version++);
        }),
      );
    }

    const task = Task.create('ParticleTargetRepresentation', async (runtime) => {
      if (!_particles) return;

      let particles = _particles;
      // a tile-ordered list (a pack written with tiles) drawn whole as a slot list of all its tiles: a tile a range,
      // as a tile view draws them; not a target an index set (its particles a run a tile)
      if (!particles.pages && canCullOnGpu(graphics)) particles = ParticleTileView.ofWhole(particles) ?? particles;
      if (particles.pages && canCullOnGpu(graphics)) {
        // a slot list: the store's object alone
        for (const visual of targetVisuals.values()) visual.destroy();
        targetVisuals.clear();
        mergedByType.clear();
        slots ??= createTargetPagedObject(materialId, graphics);
        // the themes of the slot list where it is another (the list's index their arrays by its particles)
        const theme = particles === _particles ? _theme : getSlotsTheme(_theme, particles);
        await slots.update({ graphics, runtime }, particles, theme, _props);
        renderObjects.length = 0;
        const object = slots.renderObject;
        if (object) {
          slots.applyState((ro) => applyState(ro, _state));
          if (_state.previewTransform) Visual.setPreviewTransform(object, _state.previewTransform);
          renderObjects.push(object);
          geometryState.add(object.id, slotsVersion++);
        }
        geometryState.snapshot();
        drawn = particles;
        built = { particles, theme: _theme, props: { ..._props } };
        updated.next(version++);
        return;
      }
      if (slots) {
        slots.destroy();
        slots = undefined;
      }
      // a slot list not culled on the GPU: drawn flat, picked by its keys (its slots), not marked
      if (particles.pages) particles = ParticleSlots.flatten(particles);
      drawn = particles;
      const targets = particles.targetMapping ?? new Map<number, ParticleTarget>();

      // Grouping by target id is computed once per particle list and cached on it.
      const groups = getParticleTargetGroups(particles);
      const targetProps = createParticleTargetPropsProvider(_props);

      const instanced: { targetId: number; target: ParticleTarget; indices: OrderedSet<number> }[] = [];

      for (let g = 0; g < groups.targetIds.length; g++) {
        const targetId = groups.targetIds[g];
        const target = targets.get(targetId);
        // Particles without a target object are not shown here; use the spacefill representation for them.
        if (!target) continue;

        const indices = groups.sets[g];
        if (OrderedSet.size(indices) === 0) continue;

        instanced.push({ targetId, target, indices });
      }

      const presentTargets = new Set<number>();
      const entries: { visual: TargetVisual; renderObject: GraphicsRenderObject }[] = [];

      for (const { targetId, target, indices } of instanced) {
        presentTargets.add(targetId);

        let visual = targetVisuals.get(targetId);
        if (!visual) {
          visual = createTargetVisual(targetId, materialId, graphics);
          targetVisuals.set(targetId, visual);
        }

        await visual.createOrUpdate({ graphics, runtime }, _theme, targetProps(target), particles, indices, target);

        if (visual.renderObject) {
          // Render objects can be created here, so they need the current state applied.
          applyState(visual.renderObject, _state);
          visual.setTransform(_state.transform);
          entries.push({ visual, renderObject: visual.renderObject });
          geometryState.add(visual.renderObject.id, visual.geometryVersion);
        }
      }

      // Remove visuals for targets no longer present.
      for (const [tid, visual] of targetVisuals) {
        if (!presentTargets.has(tid)) {
          visual.destroy();
          targetVisuals.delete(tid);
        }
      }

      updateRenderObjectList(entries);
      await setOccluder(runtime);
      geometryState.snapshot();
      built = { particles, theme: _theme, props: { ..._props } };
      updated.next(version++);
    });
    return serialize(task);
  }

  function updateRenderObjectList(entries: { visual: TargetVisual; renderObject: GraphicsRenderObject }[]) {
    mergedByType = mergeRenderObjectsByType(
      entries,
      _props.renderMerged,
      mergedByType,
      renderObjects,
      _props.gpuCulling,
      canCullOnGpu(graphics),
      graphics,
      _state.previewTransform,
    );
  }

  /**
   * What the list drawn whole by its targets hides on its objects (`setOccluderGrid`): its density grid, a pack's or
   * computed once in `runtime`'s task (`occluder`; without one only a grid known); none moved (the grid in the list's
   * frame), with clip objects (the grid's rays march what they cut away too: it would hide what they reveal; a change
   * of their count a full update, `isSameClipDefines`) or where the backend does not use it
   */
  async function setOccluder(runtime?: RuntimeContext) {
    let grid: OccluderGrid | undefined;
    if (_particles && !isTargetOccluderOff(_props) && Mat4.isIdentity(_state.transform) && usesOccluderGrid(graphics)) {
      grid = ParticleDensityGrid.get(_particles);
      if (!grid && _props.occluder === 'compute' && runtime) {
        grid = await ParticleDensityGrid.computeTask(_particles).runInContext(runtime);
        if (grid) ParticleDensityGrid.set(_particles, grid);
      }
    }
    for (const ro of renderObjects) setOccluderGrid(ro.state, grid);
  }

  function applyState(renderObject: GraphicsRenderObject, state: Partial<Representation.State>) {
    if (state.visible !== undefined) renderObject.state.visible = state.visible;
    if (state.pickable !== undefined) renderObject.state.pickable = state.pickable;
  }

  function setState(state: Partial<Representation.State>) {
    Representation.updateState(_state, state);
    // apply per visual, not per (possibly merged) renderObject: a merged render object only
    // aliases its first member's state, so looping renderObjects would miss the other members
    for (const visual of targetVisuals.values()) {
      if (visual.renderObject) applyState(visual.renderObject, state);
      if (state.transform !== undefined) visual.setTransform(state.transform);
    }
    // a slot list's object: its members' state (its own their first's); no transform of its own yet
    slots?.applyState((ro) => applyState(ro, state));
    // on what is drawn: a merged object moves as one, its members' are not read (`mergeRenderObjectsByType`)
    const preview = state.previewTransform;
    if (preview !== undefined) for (const ro of renderObjects) Visual.setPreviewTransform(ro, preview);
    if (state.transform !== undefined && targetVisuals.size > 0) void setOccluder();
  }

  function setTheme(theme: Theme) {
    _theme = theme;
  }

  /** the particle of index `i` of what the visuals draw: a flattened slot list's its key (its slot) */
  function particleLoci(i: number) {
    const index = drawn && drawn !== _particles && drawn.keys ? drawn.keys[i] : i;
    return Particle.Loci(_particles!, OrderedSet.ofSingleton(index));
  }

  function getMergedLoci(pickingId: PickingId): ModelLoci | undefined {
    if (!_particles) return undefined;
    return pickMergedLoci(mergedByType, pickingId, (visual, index, localInstanceId) => {
      const indices = visual.particleIndices;
      return !indices || localInstanceId < 0 || localInstanceId >= OrderedSet.size(indices)
        ? EmptyLoci
        : particleLoci(OrderedSet.getAt(indices, localInstanceId));
    });
  }

  function getLoci(pickingId: PickingId): ModelLoci {
    if (!_particles) return EmptyLoci;
    if (slots) {
      const loci = slots.getLoci(pickingId) ?? EmptyLoci;
      // a list drawn as its tiles' slot list (`ParticleTileView.ofWhole`): its slots its indices
      if (Particle.isLoci(loci) && loci.particles !== _particles) return Particle.Loci(_particles, loci.indices);
      return loci;
    }

    const mergedLoci = getMergedLoci(pickingId);
    if (mergedLoci !== undefined) return mergedLoci;

    const { objectId, instanceId } = pickingId;
    for (const visual of targetVisuals.values()) {
      if (visual.renderObject && visual.renderObject.id === objectId) {
        // instanceId is local to this target visual; map to global particle index
        const indices = visual.particleIndices;
        if (!indices || instanceId < 0 || instanceId >= OrderedSet.size(indices)) return EmptyLoci;
        return particleLoci(OrderedSet.getAt(indices, instanceId));
      }
    }
    return EmptyLoci;
  }

  function getAllLoci(): ModelLoci[] {
    if (!_particles) return [];
    return [Particle.Loci(_particles, OrderedSet.ofRange(0, _particles.count - 1))];
  }

  function eachLocation(cb: LocationCallback) {
    if (!_particles) return;
    const particles = _particles;
    const loc = Particle.Location(particles, 0);
    for (let i = 0; i < particles.count; i++) {
      loc.index = i;
      cb(loc, false);
    }
  }

  function mark(loci: ModelLoci, action: MarkerAction): boolean {
    if (!_particles) return false;
    if (!MarkerActions.is(_state.markerActions, action)) return false;
    if (!isEveryLoci(loci) && (!Particle.isLoci(loci) || loci.particles !== _particles)) return false;
    if (Particle.isLoci(loci) && Particle.lociSize(loci) === _particles.count) {
      // Change to `EveryLoci` to allow for downstream optimizations
      loci = EveryLoci;
    }
    if (slots) return slots.mark(loci, action);
    // a flattened slot list: its loci by slot, its visuals' by index
    if (drawn !== _particles) return false;
    let changed = false;
    for (const visual of targetVisuals.values()) {
      changed = visual.mark(loci, action) || changed;
    }
    return changed;
  }

  function destroy() {
    for (const visual of targetVisuals.values()) {
      visual.destroy();
    }
    targetVisuals.clear();
    mergedByType.clear();
    renderObjects.length = 0;
    slots?.destroy();
    slots = undefined;
    drawn = undefined;
  }

  return {
    label: 'Particle Target',
    get updated() {
      return updated;
    },
    get renderObjects(): ReadonlyArray<GraphicsRenderObject> {
      return renderObjects;
    },
    get props(): Readonly<ParticleTargetRepresentationProps> {
      return _props;
    },
    get params(): Readonly<ParticleTargetRepresentationParams> {
      return _params;
    },
    get state(): Readonly<Representation.State> {
      return _state;
    },
    get theme(): Readonly<Theme> {
      return _theme;
    },
    get groupCount() {
      return renderObjects.reduce((sum, ro) => sum + ((ro.values as any).uGroupCount?.ref?.value || 0), 0);
    },
    get geometryVersion() {
      return geometryState.version;
    },
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
    getSignatures: (props: ParticleTargetRepresentationProps, ctx: SignatureContext) =>
      getParticleTargetSignatures(props, ctx),
    destroy,
  };
}

/**
 * The signatures of the targets of the particle list (`ctx.data`), as `createOrUpdate` draws them: per target with
 * particles, the geometry it is drawn with (`getTargetGeometryChoice`) with its props, colored per particle, sized as
 * the visual sizes it (a structure's impostors per element, volume dots uniformly); merged per geometry type as their
 * counts decide. None without data: every kind the params configure would compile many pipelines a list never draws,
 * and the representation prepares with its data when added and updated
 */
function getParticleTargetSignatures(
  props: ParticleTargetRepresentationProps,
  ctx: SignatureContext,
): ObjectSignature[] {
  const particles = ctx.data;
  // a slot list's object: compiled as it is first drawn (its pages' targets not grouped); a tile-ordered list's too
  // where drawn as its tiles' slot list (`ParticleTileView.ofWhole`, culled on the GPU)
  if (!Particle.is(particles) || !particles.targetMapping?.size || particles.pages) return [];
  if (ctx.graphics && canCullOnGpu(ctx.graphics) && ParticleTiles.getDirectory(particles)) return [];
  const groups = getParticleTargetGroups(particles);
  // the same targets, groups, props and context as last: a list moved as a whole (the release of a move), where each
  // target's signatures cost a visual's values (848 targets of a cellpack model: 80 ms)
  const key = JSON.stringify([props, ctx.color, ctx.size, ctx.layers]);
  const cached = signaturesCache.get(particles.targetMapping);
  if (cached && cached.groups === groups && cached.graphics === ctx.graphics && cached.key === key) {
    return cached.signatures.slice();
  }
  const signatures = computeParticleTargetSignatures(props, ctx, particles.targetMapping, groups);
  signaturesCache.set(particles.targetMapping, { groups, graphics: ctx.graphics, key, signatures });
  return signatures.slice();
}

/** `getParticleTargetSignatures` per target mapping, with what it was computed from */
const signaturesCache = new WeakMap<
  ReadonlyMap<number, ParticleTarget>,
  { groups: unknown; graphics: unknown; key: string; signatures: ObjectSignature[] }
>();

function computeParticleTargetSignatures(
  props: ParticleTargetRepresentationProps,
  ctx: SignatureContext,
  targetMapping: ReadonlyMap<number, ParticleTarget>,
  groups: ReturnType<typeof getParticleTargetGroups>,
): ObjectSignature[] {
  const targetProps = createParticleTargetPropsProvider(props);
  const uniformSize = ctx.size.types.length === 1 && ctx.size.types[0] === 'uniform';
  // per geometry type: the signatures, the objects (a target each) and their particles
  const byType = new Map<string, { signatures: ObjectSignature[]; objects: number; instances: number }>();
  // a target's signatures by what they are made from (its props, one object per kind and quality, its geometry choice,
  // sizes and defines): the yeast cell's 5,639 targets a few, not a visual's values each (300 ms); added once
  const memo = new Map<ParticleTargetProps, Map<string, ObjectSignature[]>>();
  const added = new Set<ObjectSignature[]>();
  for (let g = 0; g < groups.targetIds.length; g++) {
    const target = targetMapping.get(groups.targetIds[g]);
    const count = OrderedSet.size(groups.sets[g]);
    if (!target || count === 0) continue;
    const baseProps = targetProps(target);
    // marking is per particle (see `createTargetVisual`)
    const tProps = { ...baseProps, instanceGranularity: true };
    const choice = getTargetGeometryChoice(target, tProps, ctx.graphics);
    const size: SignatureContext['size'] =
      choice === 'volume-dot-spheres' || choice === 'volume-dot-mesh'
        ? { types: ['uniform'] }
        : choice === 'structure-spheres'
          ? { types: [uniformSize ? 'uniform' : 'group'] }
          : ctx.size;
    // surfaces keep varying groups where flat interpolation is supported, a shape's mesh as it is
    const visualDefines: SignatureDefines | undefined =
      choice === 'volume-isosurface' || choice === 'blob-surface'
        ? getSurfaceSignatureDefines(tProps, ctx.graphics)
        : target.kind === 'shape' && target.shape.geometry.kind === 'mesh'
          ? { dVaryingGroup: [(target.shape.geometry as Mesh).varyingGroup.ref.value] }
          : // a sites target's impostors from its quantized sites, their second color (by its data)
            target.kind === 'sites' && choice === 'structure-spheres'
            ? getSitesSignatureDefines(target.sites, tProps as unknown as SiteColorProps, ctx.graphics)
            : target.kind === 'structure' && choice === 'structure-spheres'
              ? getStructureSignatureDefines(tProps as unknown as SiteColorProps)
              : undefined;
    const utils = getTargetGeometryUtils(target, choice);
    const shapeKind = target.kind === 'shape' ? target.shape.geometry.kind : '';
    const memoKey = `${choice}|${shapeKind}|${size.types.join()}|${visualDefines ? JSON.stringify(visualDefines) : ''}`;
    let byKey = memo.get(baseProps);
    if (!byKey) memo.set(baseProps, (byKey = new Map()));
    let signatures = byKey.get(memoKey);
    if (!signatures) {
      // placed by poses where the geometry is (`createTargetParticleTransform`): their instances read from storage on
      // WebGPU (`dInstanceStorage`)
      signatures = getGeometrySignatures(utils, tProps, { ...ctx, size }, visualDefines).map((s) =>
        isPosesGeometry(s.type) ? { ...s, instancePoses: true } : s,
      );
      byKey.set(memoKey, signatures);
    }
    for (const s of signatures) {
      let entry = byType.get(s.type);
      if (!entry) byType.set(s.type, (entry = { signatures: [], objects: 0, instances: 0 }));
      if (!added.has(signatures)) entry.signatures.push(s);
    }
    added.add(signatures);
    const type = signatures[0]?.type;
    if (type) {
      const entry = byType.get(type)!;
      entry.objects += 1;
      entry.instances += count;
    }
  }
  const out: ObjectSignature[] = [];
  const levels = getLodPixelLevelCount(props.lodPixelSize);
  for (const { signatures, objects, instances } of byType.values()) {
    out.push(...withMergedSignatures(signatures, props, { objects, instances, levels, graphics: ctx.graphics }));
  }
  return uniqueSignatures(out);
}

function getParticleTargetParams(ctx: RepresentationContext, particles: ParticleList | undefined) {
  let params = ParticleTargetRepresentationParams;
  if (particles?.targetMapping?.size) {
    params = PD.clone(ParticleTargetRepresentationParams);
    let hasStructure = false;
    let hasVolume = false;
    let hasShape = false;
    for (const target of particles.targetMapping.values()) {
      if (target.kind === 'structure' || target.kind === 'sites') hasStructure = true;
      if (target.kind === 'volume') hasVolume = true;
      if (target.kind === 'shape') hasShape = true;
    }
    if (!hasStructure) params.structure.isHidden = true;
    if (!hasVolume) params.volume.isHidden = true;
    if (!hasShape) params.shape.isHidden = true;
  }
  return params;
}

export const ParticleTargetRepresentationProvider: ParticleTargetRepresentationProvider<'target'> = {
  name: 'target',
  label: 'Target',
  description: 'Displays each particle as an instanced reference structure or shape.',
  factory: ParticleTargetRepresentation,
  getParams: getParticleTargetParams,
  defaultValues: PD.getDefaultValues(ParticleTargetRepresentationParams),
  defaultColorTheme: { name: 'particle-entity' },
  defaultSizeTheme: { name: 'uniform', props: { value: 1.6 } },
  // a slot list of no tiles yet (the camera seeing none of them, or before the first loads) by its catalog's targets:
  // its params normalized against it (`State`), the type not applicable then dropped for good
  isApplicable: (data: any) => {
    const list = data as ParticleList;
    return !!list.targetMapping?.size || (!!list.pages && list.targetInfo.size > 0);
  },
};

/** per theme and slot list the theme made again for it (`getSlotsTheme`) */
const slotsThemes = new WeakMap<Theme, { list: ParticleList; theme: Theme }>();

/**
 * `theme` (a list's) made again for `list`, a slot list of its particles (`ParticleTileView.ofWhole`): its color
 * and size themes of the same props, read by the slot list's pages
 */
function getSlotsTheme(theme: Theme, list: ParticleList): Theme {
  const kept = slotsThemes.get(theme);
  if (kept?.list === list) return kept.theme;
  const ctx = { particles: list };
  const next: Theme = {
    ...theme,
    color: theme.color.factory(ctx, theme.color.props),
    size: theme.size.factory(ctx, theme.size.props),
  };
  slotsThemes.set(theme, { list, theme: next });
  return next;
}
