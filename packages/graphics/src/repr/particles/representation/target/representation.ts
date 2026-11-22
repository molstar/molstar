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
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Theme } from '@molstar/graphics/theme/theme';
import { type GraphicsRenderObject, getNextMaterialId, isMergedRenderObject } from '@molstar/graphics/gl/render-object';
import { Subject } from 'rxjs';
import { type RuntimeContext, Task } from '@molstar/core/task';
import { type Loci as ModelLoci, EmptyLoci, isEveryLoci, EveryLoci } from '@molstar/model/model/loci';
import { type MarkerAction, MarkerActions } from '@molstar/core/util/marker-action';
import { ValueCell } from '@molstar/core/util/value-cell';
import { createTemplateCopyLayers } from './copy-layers.js';
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
  mapTargetLayers,
  TargetLayerKeys,
  type TargetLayers,
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
import { createTargetPagedObject, isNestedPagedDrawable, type TargetSlotsObject } from './paged.js';
import {
  canExpandTemplates,
  getExpandedInstance,
  createExpandedMarkerState,
  createTemplateCopyMarkers,
  getLevelLoci,
  getTemplateCopies,
  getTemplateWideLocals,
  getTemplates,
  setExpansionLayers,
  setTargetExpansion,
  type TemplateCopies,
} from './expanded.js';
import { getCopyBlockBudget } from './copy-markers.js';
import { ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import { ParticleTileDrawings } from '@molstar/model/model/particles/particle-tile-drawings';
import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';
import { ParticleTiles } from '@molstar/model/model/particles/particle-tiles';
import { ParticleDensityGrid } from '@molstar/model/model/particles/particle-tile-opaque';
import { type OccluderGrid, setOccluderGrid, usesOccluderGrid } from '@molstar/graphics/gl/occluder-grid';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';

export { ParticleTargetRepresentationParams };

/** a template's coarse level's visual of a target keyed by its id plus the level times this (`updateExpanded`) */
const CoarseLevelKeys = 2 ** 32;
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
  /** the list drawn: `_list` flattened (`ParticleNesting.flatten`), `_list` itself without templates */
  let _particles: ParticleList | undefined;
  /** the list given: its loci what the representation gives and marks (a list with templates' by path) */
  let _list: ParticleList | undefined;
  let _params: ParticleTargetRepresentationParams;
  let _props: ParticleTargetRepresentationProps;
  // Map from targetId → TargetVisual; re-used across updates.
  const targetVisuals = new Map<number, TargetVisual>();

  let mergedByType = new Map<string, MergedEntry<TargetVisual>>();
  /**
   * a list with templates drawn by its copies' poses (`expanded.ts`): the visuals of its templates' leaf targets, by
   * target id (the flattened list's), each its template and copies
   */
  const expandedVisuals = new Map<
    number,
    {
      template: ParticleList;
      copies: TemplateCopies;
      level?: ParticleNesting.TemplateLevel;
      /** an instance's template particle by its instance id */
      toLocal: (instance: number) => number;
      /** its target (its leaf list's id, `getLeafList`) */
      target: number;
    }
  >();
  /** their copies' markers (a copy marked apart from the others), their objects' own marker state with them */
  const copyMarkers = createTemplateCopyMarkers(() => getCopyBlockBudget(graphics));
  const expandedMarkers = createExpandedMarkerState();
  /** their copies' own layers (a layer on a copy apart from the others) */
  const copyLayers = createTemplateCopyLayers(() => getCopyBlockBudget(graphics));
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
    // what is drawn: the targets' visuals, or the store's object (a slot list on the GPU)
    if (!built || built.particles !== _particles || built.theme !== _theme) return false;
    if (targetVisuals.size === 0 && !slots?.renderObject) return false;
    const before = built.props;
    return BaseGeometry.isClipChangeOnly(_props, before) && BaseGeometry.isSameClipDefines(before.clip, _props.clip);
  }

  /** updates one after another, see `Representation.createUpdateQueue` */
  const serialize = Representation.createUpdateQueue();
  function createOrUpdate(props: Partial<ParticleTargetRepresentationProps> = {}, input?: ParticleList) {
    if (input) {
      _list = input;
      // a list with templates by its leaves' targets (`getThemeList`: no copy expanded)
      _params = getParams(ctx, ParticleNesting.hasTemplates(input) ? ParticleNesting.getThemeList(input) : input);
      if (!_props) _props = PD.getDefaultValues(_params);
    }
    // Quality-derived props are resolved per target (see `createParticleTargetPropsProvider`),
    // since each target has its own data and its own geometry params.
    Object.assign(_props, props);
    // a list with templates drawn by the store's object as it is where it draws its every leaf target (culled and
    // composed on the GPU, a template's particles once), else flattened (`ParticleNesting.flatten`, cached on it): its
    // copies expanded on the CPU only where drawn so; again with the props (a geometry the store does not draw)
    if (_list) {
      // drawn by its copies' poses where the vertex stage composes them (`expanded.ts`, WebGL2)
      _particles =
        ParticleNesting.hasTemplates(_list) &&
        !(canCullOnGpu(graphics) && isNestedPagedDrawable(_list, _props, graphics)) &&
        !(!canCullOnGpu(graphics) && canExpandTemplates(graphics))
          ? ParticleNesting.flatten(_list)
          : _list;
    }
    // how it draws a tile view (its clip objects, its opacity): what the view's picks go by, a moved cut too (a
    // clipping gizmo's drag updates the representation alone, not the state)
    ParticleTileDrawings.set(self, _particles, _props);

    if (isClipOnly()) {
      const clip = _props.clip;
      built!.props = { ..._props };
      return serialize(
        Task.create('ParticleTargetRepresentation clip', async () => {
          // what is drawn alone: a visual not merged, a merged object's first member (its clip cells the object's);
          // not the thousands of other members (not read while merged): a full update gives them the clip
          // (`BaseGeometry.updateValues`) before they are merged again
          for (const ro of renderObjects) {
            if (ro === slots?.renderObject) continue;
            const values = isMergedRenderObject(ro) ? ro.members[0].values : ro.values;
            BaseGeometry.updateClip(values as any, clip);
          }
          // the store's object: its uniforms and its cull's clip, no page written again
          slots?.updateClip(clip);

          updated.next(version++);
        }),
      );
    }

    const task = Task.create('ParticleTargetRepresentation', async (runtime) => {
      if (!_particles) return;

      let particles = _particles;
      // a list with templates drawn as it is: by the store's object, or by its copies' poses (see `createOrUpdate`)
      const nested = ParticleNesting.hasTemplates(particles);
      const nestedPaged = nested && canCullOnGpu(graphics);
      // a tile-ordered list (a pack written with tiles) drawn whole as a slot list of all its tiles: a tile a range,
      // as a tile view draws them; not a target an index set (its particles a run a tile)
      if (!nested && !particles.pages && canCullOnGpu(graphics))
        particles = ParticleTileView.ofWhole(particles) ?? particles;
      if ((particles.pages || nestedPaged) && canCullOnGpu(graphics)) {
        // a slot list: the store's object alone
        for (const visual of targetVisuals.values()) visual.destroy();
        targetVisuals.clear();
        mergedByType.clear();
        if (!slots) {
          slots = createTargetPagedObject(materialId, graphics);
          // the state's layers: applied by its update, as it lays out its texels
          const layers = getStateLayers(_state);
          if (layers) slots.setLayers(layers);
        }
        // the themes of the slot list where it is another (the list's index their arrays by its particles)
        // a list with templates: by its leaves by level (`getThemeList`: a page a template's leaves, colored once for
        // all copies), not an index of the flattened list
        const theme = nestedPaged
          ? getSlotsTheme(_theme, ParticleNesting.getThemeList(particles))
          : particles === _particles
            ? _theme
            : getSlotsTheme(_theme, particles);
        await slots.update({ graphics, runtime }, particles, theme, _props);
        renderObjects.length = 0;
        const object = slots.renderObject;
        if (object) {
          slots.applyState((ro) => applyState(ro, _state));
          if (_state.previewTransform) Visual.setPreviewTransform(object, _state.previewTransform);
          renderObjects.push(object);
          geometryState.add(object.id, slotsVersion++);
        }
        // a list with templates drawn whole: what it hides (its templates' grids placed by their copies), as a flat
        // list's; a tile view's its store's object sets
        if (nestedPaged) await setOccluder(runtime);
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
      if (nested) {
        await updateExpanded(runtime, particles);
        return;
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

      // the layers on the visuals as made (before they are merged)
      applyLayers(getStateLayers(_state));
      updateRenderObjectList(entries);
      await setOccluder(runtime);
      geometryState.snapshot();
      built = { particles, theme: _theme, props: { ..._props } };
      updated.next(version++);
    });
    return serialize(task);
  }

  /**
   * A list with templates drawn by its copies' poses (`expanded.ts`): its own leaves' visuals as a flat list's (of its
   * top list, `getTopList`), a visual a leaf target of each template (of its leaf list, `getLeafList`: its particles in
   * the template's frame) drawn by the template's copies, not merged
   */
  async function updateExpanded(runtime: RuntimeContext, list: ParticleList) {
    drawn = list;
    const targetProps = createParticleTargetPropsProvider(_props);
    const present = new Set<number>();
    const entries: { visual: TargetVisual; renderObject: GraphicsRenderObject }[] = [];
    const expandedObjects: GraphicsRenderObject[] = [];
    const previous = new Map(expandedVisuals);
    expandedVisuals.clear();
    /** `leaves`' visuals; of a template's at its coarse level `level` (none: the template whole) by `copies` */
    const add = async (
      leaves: ParticleList,
      template?: ParticleList,
      copies?: TemplateCopies,
      level?: ParticleNesting.TemplateLevel,
    ) => {
      const groups = getParticleTargetGroups(leaves);
      const targets = leaves.targetMapping;
      for (let g = 0; g < groups.targetIds.length; g++) {
        const id = groups.targetIds[g];
        const target = targets?.get(id);
        const indices = groups.sets[g];
        if (!target || target.kind === 'particles' || OrderedSet.size(indices) === 0) continue;
        // a coarse level's visual of a target keyed past the target ids
        const targetId = id + (level ? level.level * CoarseLevelKeys : 0);
        present.add(targetId);
        let visual = targetVisuals.get(targetId);
        // a visual drawn flat before (a transform of its own) made again for the template's frame
        if (visual && template && !previous.has(targetId)) {
          visual.destroy();
          visual = undefined;
        }
        if (!visual) {
          visual = createTargetVisual(targetId, materialId, graphics);
          targetVisuals.set(targetId, visual);
        }
        await visual.createOrUpdate({ graphics, runtime }, _theme, targetProps(target), leaves, indices, target);
        const ro = visual.renderObject;
        if (!ro) continue;
        applyState(ro, _state);
        if (copies) {
          // a leaf list's particle its template particle, a coarse level's by its indices
          const toLocal = level
            ? (i: number) => level.indices[OrderedSet.getAt(indices, i)]
            : (i: number) => OrderedSet.getAt(indices, i);
          setTargetExpansion(ro, copies, level?.level ?? 0, copyMarkerOf(template!), toLocal, copyLayerOf(template!));
          expandedVisuals.set(targetId, { template: template!, copies, level, toLocal, target: id });
          expandedObjects.push(ro);
        } else {
          visual.setTransform(_state.transform);
          entries.push({ visual, renderObject: ro });
        }
        geometryState.add(ro.id, visual.geometryVersion);
      }
    };
    await add(ParticleNesting.getTopList(list));
    const templates = getTemplates(list).map((t) => [t, getTemplateCopies(list, t, _state.transform)] as const);
    // the copies' markers laid out (cleared), as the objects' own are made again
    copyMarkers.layout(list, templates);
    // and their layers (none: the state's applied once the objects are made, `applyLayers`)
    copyLayers.layout(list, templates);
    expandedMarkers.clear();
    for (const [template, copies] of templates) {
      await add(ParticleNesting.getLeafList(list, template), template, copies);
      // its far copies by its coarse levels (a copy one of them by the view, `getCopyLevels`)
      for (const level of ParticleNesting.getTemplateLevels(list, template)) {
        await add(level.list, template, copies, level);
      }
    }

    for (const [tid, visual] of targetVisuals) {
      if (!present.has(tid)) {
        visual.destroy();
        targetVisuals.delete(tid);
      }
    }
    applyLayers(getStateLayers(_state));
    updateRenderObjectList(entries);
    for (const ro of expandedObjects) {
      if (_state.previewTransform) Visual.setPreviewTransform(ro, _state.previewTransform);
      renderObjects.push(ro);
    }
    await setOccluder(runtime);
    geometryState.snapshot();
    built = { particles: list, theme: _theme, props: { ..._props } };
    updated.next(version++);
  }

  /** the copies' layers' image an object of `template` reads (by the markers' copies), whether a copy has one */
  function copyLayerOf(template: ParticleList) {
    return { image: copyLayers.image, has: copyLayers.has(template) };
  }

  /** the copies' markers' image an object of `template` reads, its first copy in it */
  function copyMarkerOf(template: ParticleList) {
    return { image: copyMarkers.image, base: copyMarkers.base(template) };
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
      // a list with templates': its templates' grids placed by their copies (given where theirs are: a pack's)
      grid =
        ParticleDensityGrid.get(_particles) ??
        (ParticleNesting.hasTemplates(_particles) ? ParticleNesting.getDensityGrid(_particles) : undefined);
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

  /** the layers `state` gives, none where it gives none (an empty layer one: what clears) */
  function getStateLayers(state: Partial<Representation.State>): TargetLayers | undefined {
    let layers: TargetLayers | undefined;
    for (const k of TargetLayerKeys) {
      if (state[k] === undefined) continue;
      layers ??= {};
      (layers as any)[k] = state[k];
    }
    return layers;
  }

  /**
   * `layers` on the visuals: a list with templates flattened by its flat loci; drawn by its copies its objects by
   * their template's particles every copy holds (`getTemplateWideLocals`), the rest by copy (`TemplateCopyLayers`:
   * an object's layer kept where it reaches none of its own, its copies' read by it); the store's object its own
   * (by store index and copy); a flattened slot list's visuals none
   */
  function applyLayers(layers: TargetLayers | undefined) {
    if (!layers) return;
    // the store's object: by store index and by copy (`TargetSlotsObject.setLayers`), kept for its updates
    if (slots) {
      slots.setLayers(layers);
      return;
    }
    if (!_particles || drawn !== _particles) return;
    const list = _particles;
    const flattened = _list !== _particles;
    if (expandedVisuals.size > 0) copyLayers.set(list, layers as Parameters<typeof copyLayers.set>[1]);
    for (const [targetId, visual] of targetVisuals) {
      const e = expandedVisuals.get(targetId);
      if (e) {
        visual.setLayers(
          mapTargetLayers(layers, (loci) => {
            if (isEveryLoci(loci)) return loci;
            const locals = getTemplateWideLocals(list, e.template, e.copies.tops, loci);
            return (locals && getLevelLoci(list, e.template, locals, e.level)) ?? EmptyLoci;
          }),
        );
        const ro = visual.renderObject;
        if (ro) {
          setExpansionLayers(ro, copyLayerOf(e.template));
          raiseCopyAverages(ro, e.template);
        }
      } else if (flattened) {
        visual.setLayers(
          mapTargetLayers(layers, (loci) =>
            Particle.isLoci(loci) && loci.particles === _list ? ParticleNesting.toFlatLoci(loci) : loci,
          ),
        );
      } else {
        visual.setLayers(layers);
      }
    }
  }

  /**
   * An expanded object's averages where its copies have transparency or emission of their own (its own none, or
   * less): what puts it into the passes that draw them (its own recomputed where its layer is set again)
   */
  function raiseCopyAverages(ro: GraphicsRenderObject, template: ParticleList) {
    const values = ro.values as unknown as {
      transparencyAverage: ValueCell<number>;
      emissiveAverage: ValueCell<number>;
      wiggleAverage: ValueCell<number>;
    };
    if (copyLayers.has(template, 'transparency') && values.transparencyAverage.ref.value < 0.5) {
      ValueCell.update(values.transparencyAverage, 0.5);
    }
    if (copyLayers.has(template, 'emissive') && values.emissiveAverage.ref.value < 0.5) {
      ValueCell.update(values.emissiveAverage, 0.5);
    }
    // what animates it (drawn anew each frame)
    if (copyLayers.has(template, 'wiggle') && values.wiggleAverage.ref.value < 0.5) {
      ValueCell.update(values.wiggleAverage, 0.5);
    }
  }

  function setState(state: Partial<Representation.State>) {
    Representation.updateState(_state, state);
    applyLayers(getStateLayers(state));
    // apply per visual, not per (possibly merged) renderObject: a merged render object only
    // aliases its first member's state, so looping renderObjects would miss the other members
    const moved = new Map<ParticleList, TemplateCopies>();
    for (const [targetId, visual] of targetVisuals) {
      if (visual.renderObject) applyState(visual.renderObject, state);
      if (state.transform === undefined) continue;
      // a template's leaf target: its copies placed by the transform (its instances in the template's frame)
      const e = expandedVisuals.get(targetId);
      if (e && _particles && visual.renderObject) {
        // a template's copies once for its visuals
        let copies = moved.get(e.template);
        if (!copies) moved.set(e.template, (copies = getTemplateCopies(_particles, e.template, state.transform)));
        e.copies = copies;
        setTargetExpansion(
          visual.renderObject,
          copies,
          e.level?.level ?? 0,
          copyMarkerOf(e.template),
          e.toLocal,
          copyLayerOf(e.template),
        );
      } else {
        visual.setTransform(state.transform);
      }
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

  /** a loci of the drawn list as the given list's (a flattened list's instances by path) */
  function toListLoci(loci: ModelLoci): ModelLoci {
    if (!_list || _list === _particles || !Particle.isLoci(loci) || loci.particles !== _particles) return loci;
    return ParticleNesting.fromFlatIndices(_list, loci.indices);
  }

  function getLoci(pickingId: PickingId): ModelLoci {
    return toListLoci(getDrawnLoci(pickingId));
  }

  function getDrawnLoci(pickingId: PickingId): ModelLoci {
    if (!_particles) return EmptyLoci;
    if (slots) {
      const loci = slots.getLoci(pickingId) ?? EmptyLoci;
      // a list drawn as its tiles' slot list (`ParticleTileView.ofWhole`): its slots its indices
      if (Particle.isLoci(loci) && loci.particles !== _particles) return Particle.Loci(_particles, loci.indices);
      return loci;
    }

    const mergedLoci = getMergedLoci(pickingId);
    if (mergedLoci !== undefined) return mergedLoci;

    // a template's leaf target: its draw slot the copy and the template particle (`getExpandedInstance`)
    for (const [targetId, e] of expandedVisuals) {
      const visual = targetVisuals.get(targetId);
      const ro = visual?.renderObject;
      if (!ro || ro.id !== pickingId.objectId || !visual!.particleIndices) continue;
      const instance = getExpandedInstance(
        ro,
        e.copies.tops,
        visual!.particleIndices,
        pickingId.instanceId,
        e.level?.indices,
      );
      return instance ? ParticleNesting.ofInstance(_particles, instance.top, instance.local) : EmptyLoci;
    }

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
    return [toListLoci(Particle.Loci(_particles, OrderedSet.ofRange(0, _particles.count - 1)))];
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
    // a loci of the given list with templates marked as its flattened list's
    if (_list !== _particles && Particle.isLoci(loci) && loci.particles === _list)
      loci = ParticleNesting.toFlatLoci(loci);
    if (!isEveryLoci(loci) && (!Particle.isLoci(loci) || loci.particles !== _particles)) return false;
    // all of it: a list with templates' by its instances (a loci's size counts them), not its top particles
    const all = ParticleNesting.hasTemplates(_particles)
      ? ParticleNesting.getExpandedCount(_particles)
      : _particles.count;
    if (Particle.isLoci(loci) && Particle.lociSize(loci) === all) {
      // Change to `EveryLoci` to allow for downstream optimizations
      loci = EveryLoci;
    }
    if (slots) return slots.mark(loci, action);
    // a flattened slot list: its loci by slot, its visuals' by index
    if (drawn !== _particles) return false;
    if (expandedVisuals.size > 0) return markExpanded(_particles, loci, action);
    let changed = false;
    for (const visual of targetVisuals.values()) changed = visual.mark(loci, action) || changed;
    return changed;
  }

  /**
   * A list with templates drawn by its copies: its copies' markers by copy (`TemplateCopyMarkers`), a template's
   * particles in every copy and every particle by its objects' own; its own leaves' visuals as a flat list's
   */
  function markExpanded(list: ParticleList, loci: ModelLoci, action: MarkerAction) {
    const { changed: copiesChanged, templateWide } = copyMarkers.mark(list, loci, action);
    const every = isEveryLoci(loci);
    let changed = copiesChanged;
    for (const [targetId, visual] of targetVisuals) {
      const e = expandedVisuals.get(targetId);
      const ro = visual.renderObject;
      if (!e || !ro) {
        changed = visual.mark(loci, action) || changed;
        continue;
      }
      // the image grown: another texture
      ValueCell.updateIfChanged((ro.values as any).tCopyMarkers, copyMarkers.image);
      expandedMarkers.restore(ro);
      const locals = templateWide.get(e.template);
      const own = every ? loci : locals && getLevelLoci(list, e.template, locals, e.level);
      if (own) changed = visual.mark(own, action) || changed;
      expandedMarkers.apply(ro, copyMarkers.isMarked(e.template, e.target));
    }
    return changed;
  }

  function destroy() {
    ParticleTileDrawings.set(self, undefined, undefined);
    for (const visual of targetVisuals.values()) {
      visual.destroy();
    }
    targetVisuals.clear();
    expandedVisuals.clear();
    mergedByType.clear();
    renderObjects.length = 0;
    slots?.destroy();
    slots = undefined;
    drawn = undefined;
  }

  const self: ParticleTargetRepresentation = {
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
  return self;
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

/**
 * the params by the kinds of targets present (a bit each: structure, volume, shape), made once each: asked at every
 * resolve of a list's params (a tile set's, each step), a clone of the whole tree each time otherwise
 */
const ParticleTargetParamsByKinds = new Map<number, typeof ParticleTargetRepresentationParams>();

function getParticleTargetParams(ctx: RepresentationContext, particles: ParticleList | undefined) {
  if (!particles?.targetMapping?.size) return ParticleTargetRepresentationParams;
  let kinds = 0;
  for (const target of particles.targetMapping.values()) {
    if (target.kind === 'structure' || target.kind === 'sites') kinds |= 1;
    else if (target.kind === 'volume') kinds |= 2;
    else if (target.kind === 'shape') kinds |= 4;
    if (kinds === 7) break;
  }
  let params = ParticleTargetParamsByKinds.get(kinds);
  if (!params) {
    params = PD.clone(ParticleTargetRepresentationParams);
    if (!(kinds & 1)) params.structure.isHidden = true;
    if (!(kinds & 2)) params.volume.isHidden = true;
    if (!(kinds & 4)) params.shape.isHidden = true;
    ParticleTargetParamsByKinds.set(kinds, params);
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
