/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { isPromiseLike } from '@molstar/core/util/type-helpers';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { HashMap } from '@molstar/core/util/map';
import { hash2 } from '@molstar/core/data/util';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ValueCell } from '@molstar/core/util/value-cell';
import { Bond, type Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import { EmptyLoci, isEveryLoci, type Loci } from '@molstar/model/model/loci';
import type { Geometry } from '@molstar/graphics/geo/geometry/geometry';
import {
  BaseGeometry,
  type GpuCullingValue,
  type RenderMergedValue,
  resolveInstanceGranularity,
} from '@molstar/graphics/geo/geometry/base';
import { PickingId } from '@molstar/graphics/geo/geometry/picking';
import {
  createMergedRecordsRenderObject,
  type GraphicsRenderObject,
  hasForeignResources,
  type MergedGraphicsRenderObject,
} from '@molstar/graphics/gl/render-object';
import {
  canCullOnGpu,
  canMergeValues,
  getMergeKey,
  getMergeRecordKeys,
  isMergeableType,
  type MergeableValues,
  type MergeRecord,
  type MergedType,
} from '@molstar/graphics/gl/renderable/merged';
import { RenderableState } from '@molstar/graphics/gl/renderable';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { Visual, type VisualContext } from '../visual.js';
import { type LocationCallback, VisualUpdateState } from '../util.js';
import { resolveMergedCulling } from '../merged.js';
import type { StructureParams } from './params.js';
import type { StructureRepresentationState } from './representation.js';
import { createUnitsTransform, includesUnitKind, type StructureGroup } from './visual/util/common.js';
import { getGroupLoci } from './visual/util/loci-index.js';
import {
  applyUnitsLoci,
  setUnitsUpdateState,
  type UnitsIterators,
  type UnitsVisualGeometryBuilder,
  updateUnitsValues,
} from './units-visual.js';

/**
 * A segment's cells beside its record (`MergeRecord`): its transform's, which the transform writers read back
 * (`createUnitsTransform`, `Visual.setTransform`: `aTransform` computed of them), not the merged object. `transform`
 * its own; the extra transform the shared empty one until unit transforms are set (`ExtraTransformKeys`). The matrix
 * (the representation's transform) and the poses (none for units) the shared values'
 */
const SegmentTransformKeys = ['transform', 'extraTransform', 'hasExtraTransform'] as const;
const ExtraTransformKeys = ['extraTransform', 'hasExtraTransform'] as const;

/**
 * The segments of one merge key (`getMergeKey`) and their merged object: groups of the same props and theme differ in
 * it where their counts resolve the instance granularity apart (`auto`: a large unit of many copies by instance, a
 * small one by group and instance), or their geometry's defines. Its shared values those of its first segment made
 */
type Part<G extends Geometry> = {
  readonly shared: MergeableValues;
  readonly recordKeys: readonly string[];
  /**
   * a record of `recordKeys` as a fast object to copy (`JSON.parse`'s): a record built key by key is a dictionary
   * (some 40 keys), its properties 1.6 KB a segment
   */
  readonly recordTemplate: { readonly [k: string]: unknown };
  /** the extra transform cells its segments share until unit transforms are set (`ownExtraTransforms`) */
  sharedExtra: { [k: string]: ValueCell<any> } | undefined;
  /** its segments, in the groups' order */
  segments: Segment<G>[];
  /** its object's state, another for each object made (a disposed object's state stays disposed) */
  renderableState: RenderableState;
  object: MergedGraphicsRenderObject | undefined;
  objectRecords: readonly MergeRecord[];
  /** the segments of the object's records, in their order: those with locations (as a units visual's render object) */
  drawn: Segment<G>[];
  objectKey: string;
  objectGraphics: GraphicsContext | undefined;
  /** the visual's `recordsVersion` and the `gpuCulling` prop the object was last updated with (`updateObject`) */
  objectRecordsVersion: number;
  objectCulling: GpuCullingValue | undefined;
};

/** A symmetry group of a merged units visual: its record and what its updates read */
type Segment<G extends Geometry> = {
  group: Unit.SymmetryGroup;
  geometry: G;
  readonly part: Part<G>;
  /**
   * its cells: the merged object's record of it, and its transform's (`SegmentTransformKeys`); its part's shared values
   * its prototype, so that the writers read and write it as a units visual's values (`UnitsMergedVisual`)
   */
  readonly record: MergeRecord;
  /** the record, as a render object's values for the writers that take one (`Visual.mark`, the layer setters) */
  readonly values: MergeableValues;
  /** of its location iterator at its last update */
  groupCount: number;
  instanceCount: number;
  previousMark?: Visual.PreviousMark;
  /**
   * its update state's info (`VisualUpdateState.info`: the versions its geometry was built of, the builder's
   * `initUpdateState`), its own: a group's, as a units visual's, not the visual's shared update state's
   */
  readonly info: { [k: string]: unknown };
};

/**
 * One visual over the symmetry groups of a units representation, merged without a visual, values or render object a
 * group: a group a segment, its record its own cells (`MergeRecord`), the rest the shared values of its part (the
 * groups of its merge key, a merged object each: one or a few). The writers of a units visual's values
 * (`setUnitsUpdateState`, `updateUnitsValues`, `Visual.mark`, the layer setters) write a segment's record, the shared
 * values its prototype: its cells change, the shared ones written alike by every segment of the part
 */
export interface UnitsMergedVisual<P extends StructureParams> {
  readonly builder: UnitsVisualGeometryBuilder<P, any>;
  /** its parts' merged objects */
  readonly renderObjects: readonly GraphicsRenderObject[];
  /** advanced where a segment's geometry was made again or a merged object is another */
  readonly geometryVersion: number;
  readonly groupCount: number;
  /**
   * whether its groups merge (`canMergeValues` of each part's shared values); false: per-group visuals instead, the
   * visual destroyed (an update that finds a part not to merge ends there, its groups not all updated)
   */
  readonly mergeable: boolean;
  createOrUpdate(
    ctx: VisualContext,
    theme: Theme,
    props: PD.Values<P>,
    structure: Structure,
    groups: ReadonlyArray<Unit.SymmetryGroup>,
  ): Promise<void>;
  getLoci(pickingId: PickingId): Loci;
  /** `loci` of the representation's structure (`UnitsRepresentation.mark` remaps it), or `EveryLoci` */
  mark(loci: Loci, action: MarkerAction): boolean;
  eachLocation(cb: LocationCallback): void;
  /** the representation's state: kept, and applied to every segment made later */
  setState(state: Partial<StructureRepresentationState>, current: StructureRepresentationState): void;
  destroy(): void;
}

/** What the loci of a segment are applied by: its group and its location iterator's counts */
type SegmentPlace = { readonly group: Unit.SymmetryGroup; readonly groupCount: number; readonly instanceCount: number };

function createSegmentMap<G extends Geometry>() {
  return new HashMap<Unit.SymmetryGroup, Segment<G>>(
    (group) => hash2(group.hashCode, group.transformHash),
    Unit.SymmetryGroup.areInvariantElementsEqual,
  );
}

/** What of a segment's own values an update makes again: its geometry, transforms, colors or sizes (not the props') */
function hasRecordUpdates(state: VisualUpdateState) {
  return state.createGeometry || state.updateTransform || state.updateMatrix || state.updateColor || state.updateSize;
}

/** The versions of the cells of `seg`'s record (its own keys, not the shared values beneath) */
function recordVersions(seg: { readonly record: object }) {
  const record = seg.record as { [k: string]: ValueCell<unknown> | undefined };
  return Object.keys(record).map((k) => record[k]?.ref?.version ?? -1);
}

function isSameVersions(a: number[], b: number[]) {
  return a.length === b.length && a.every((v, i) => v === b[i]);
}

export function UnitsMergedVisual<G extends Geometry, P extends StructureParams & Geometry.Params<G>>(
  builder: UnitsVisualGeometryBuilder<P, G>,
  materialId: number,
): UnitsMergedVisual<P> {
  const utils = builder.geometryUtils;
  const type = utils.createEmpty().kind as MergedType;
  if (!isMergeableType(type)) throw new Error(`'${type}' does not merge`);

  const updateState = VisualUpdateState.create();
  let segments: Segment<G>[] = [];
  let byGroup = createSegmentMap<G>();
  /** segments by the ids of their units, to mark only those a small loci has; reset with the segments */
  let byUnitId: Map<number, Segment<G>> | undefined;
  let parts: Part<G>[] = [];
  let renderObjects: GraphicsRenderObject[] = [];
  let geometryVersion = 0;
  /** its last update found a part not to merge and ended there (`mergeable`) */
  let unmergeable = false;
  /**
   * advanced where the records' cells are written but the observed ones (a mark's): the merged objects then read all
   * of their members at their next sync, else those of changed observed cells (`Merged.recordsUnchanged`)
   */
  const recordsVersion = ValueCell.create(0);
  function recordsWritten() {
    ValueCell.update(recordsVersion, recordsVersion.ref.value + 1);
  }

  /** the props, theme and structure of the segments' last update */
  let currentProps: PD.Values<P> | undefined;
  let currentTheme: Theme | undefined;
  let currentStructure: Structure | undefined;
  /** the groups of the segments' last update, in the segments' order */
  let currentGroups: ReadonlyArray<Unit.SymmetryGroup> | undefined;
  /** the structure and instance granularity the loci are applied with: the update's from its start */
  let lociStructure: Structure | undefined;
  let lociGranularity: PD.Values<P>['instanceGranularity'] = 'auto';

  /** the representation's state, applied to segments made later (`setState`) */
  const state: Partial<StructureRepresentationState> = {};

  /** `seg` as the render object the writers take (`Visual.mark`, the layer setters): its values its record */
  function asObject(seg: Segment<G>) {
    return seg as unknown as GraphicsRenderObject;
  }

  function structureGroupOf(seg: SegmentPlace): StructureGroup {
    return { group: seg.group, structure: lociStructure! };
  }

  function getVisualLoci(seg: SegmentPlace) {
    return (loci: Loci) =>
      StructureElement.Loci.is(loci) ? getGroupLoci(loci, seg.group.units, seg.group.unitIndexMap) : loci;
  }

  function lociApply(seg: SegmentPlace): Visual.LociApply {
    const structureGroup = structureGroupOf(seg);
    const visualLoci = getVisualLoci(seg);
    return (loci, apply, isMarking) =>
      applyUnitsLoci(
        builder.eachLocation,
        lociGranularity,
        structureGroup,
        seg.groupCount,
        seg.instanceCount,
        loci,
        visualLoci,
        apply,
        isMarking,
      );
  }

  function createGeometry(
    ctx: VisualContext,
    structureGroup: StructureGroup,
    props: PD.Values<P>,
    theme: Theme,
    geometry?: G,
  ) {
    const unit = structureGroup.group.units[0];
    return includesUnitKind(props.unitKinds, unit)
      ? builder.createGeometry(ctx, unit, structureGroup.structure, theme, props, geometry)
      : utils.createEmpty(geometry);
  }

  /** The part of `values`' merge key, a new one (`values` its shared values) where there is none */
  function getPart(values: MergeableValues, props: PD.Values<P>): Part<G> {
    const key = getMergeKey(values);
    for (const part of parts) if (getMergeKey(part.shared) === key) return part;
    const recordKeys = getMergeRecordKeys(type, values);
    for (const k of SegmentTransformKeys) if (values[k] !== undefined) recordKeys.push(k);
    const part: Part<G> = {
      shared: values,
      recordKeys,
      recordTemplate: JSON.parse(JSON.stringify(Object.fromEntries(recordKeys.map((k) => [k, null])))),
      sharedExtra: state.unitTransforms
        ? undefined
        : { extraTransform: ValueCell.create(new Float32Array(0)), hasExtraTransform: ValueCell.create(false) },
      segments: [],
      renderableState: utils.createRenderableState(props),
      object: undefined,
      objectRecords: [],
      drawn: [],
      objectKey: '',
      objectGraphics: undefined,
      objectRecordsVersion: -1,
      objectCulling: undefined,
    };
    parts.push(part);
    applyPartState(part);
    return part;
  }

  /** A segment of `structureGroup` with values of its own (`createValues`), picked to its record in its part */
  function createSegment(
    ctx: VisualContext,
    structureGroup: StructureGroup,
    geometry: G,
    props: PD.Values<P>,
    theme: Theme,
  ) {
    const locationIt = builder.createLocationIterator(structureGroup, props);
    const transform = createUnitsTransform(
      structureGroup,
      props.includeParent,
      geometry.boundingSphere,
      props.cellSize,
      props.batchSize,
    );
    const values = utils.createValues(geometry, transform, locationIt, theme, props) as unknown as MergeableValues;
    builder.processValues?.(values as any, geometry, props, theme, ctx.graphics);
    // the kept state on its values before its part is picked: their layers' defines as its part's (one made with them),
    // its merge key that of the segments made before it
    const place: SegmentPlace = {
      group: structureGroup.group,
      groupCount: locationIt.groupCount,
      instanceCount: locationIt.instanceCount,
    };
    applyStateTo({ values } as unknown as GraphicsRenderObject, place);
    const part = getPart(values, props);
    // the shared values its prototype: once a record, made as a fast object (a prototype set then, not keys added)
    const record: { [k: string]: unknown } = { ...part.recordTemplate };
    Object.setPrototypeOf(record, part.shared);
    for (const k of part.recordKeys) record[k] = values[k];
    if (part.sharedExtra) for (const k of ExtraTransformKeys) record[k] = part.sharedExtra[k];
    const info: { [k: string]: unknown } = {};
    if (builder.initUpdateState) {
      VisualUpdateState.reset(updateState);
      updateState.info = info;
      builder.initUpdateState(updateState, props, theme, structureGroup);
    }
    const seg: Segment<G> = {
      info,
      group: structureGroup.group,
      geometry,
      part,
      record: record as MergeRecord,
      values: record as unknown as MergeableValues,
      groupCount: locationIt.groupCount,
      instanceCount: locationIt.instanceCount,
    };
    return seg;
  }

  /** Updates `seg` to `structureGroup`, `props` and `theme` from what it was last; the geometry made again where needed */
  function updateSegment(
    ctx: VisualContext,
    seg: Segment<G>,
    structureGroup: StructureGroup,
    props: PD.Values<P>,
    theme: Theme,
  ): Promise<void> | void {
    const current: StructureGroup = { group: seg.group, structure: currentStructure! };
    VisualUpdateState.reset(updateState);
    updateState.info = seg.info;
    const values = seg.record as any;
    setUnitsUpdateState(
      updateState,
      builder,
      values,
      seg.geometry,
      props,
      currentProps!,
      theme,
      currentTheme!,
      structureGroup,
      current,
      ctx.graphics,
      false,
    );
    const apply = (newGeometry?: G) => {
      // the position iterator made where a writer reads it (vertex colors and sizes; not a recolor of a uniform color)
      let positionIt: LocationIterator | undefined;
      const its: UnitsIterators = {
        locationIt: builder.createLocationIterator(current, currentProps!),
        get positionIt() {
          return (positionIt ??= utils.createPositionIterator(seg.geometry, values));
        },
        set positionIt(it) {
          positionIt = it;
        },
      };
      updateUnitsValues(
        builder,
        updateState,
        values,
        its,
        newGeometry,
        seg.geometry,
        props,
        currentProps!,
        theme,
        structureGroup,
        current,
      );
      seg.group = structureGroup.group;
      if (newGeometry) {
        seg.geometry = newGeometry;
        geometryVersion += 1;
      }
      seg.groupCount = its.locationIt.groupCount;
      seg.instanceCount = its.locationIt.instanceCount;
      builder.processValues?.(values, seg.geometry, props, theme, ctx.graphics);
    };
    if (updateState.createGeometry) {
      const geometry = createGeometry(ctx, structureGroup, props, theme, seg.geometry);
      return Visual.withPreparedGeometry(geometry, utils, theme, ctx.graphics, apply);
    }
    return Visual.withPreparedGeometry(seg.geometry, utils, theme, ctx.graphics, () => apply());
  }

  function disposeObject(part: Part<G>) {
    if (part.object) {
      const { state } = part.object;
      state.disposed = true;
      // the next object a state of its own: made with the disposed one, it was never drawn (a backend switch)
      part.renderableState = RenderableState.create({
        disposed: false,
        visible: state.visible,
        pickable: state.pickable,
        colorOnly: state.colorOnly,
        pickOnly: state.pickOnly,
        opaque: state.opaque,
        writeDepth: state.writeDepth,
      });
    }
    part.object = undefined;
    part.objectRecords = [];
  }

  function reset() {
    for (const seg of segments) builder.dispose?.(seg.geometry);
    segments = [];
    byGroup = createSegmentMap<G>();
    byUnitId = undefined;
    for (const part of parts) disposeObject(part);
    parts = [];
    renderObjects = [];
    currentProps = undefined;
    currentTheme = undefined;
    currentStructure = undefined;
    currentGroups = undefined;
    lociStructure = undefined;
  }

  /**
   * Every segment's extra transform its own (from its part's shared empty ones), before unit transforms are set: not
   * cells of the merged objects' records, so replaced without a new object
   */
  function ownExtraTransforms() {
    for (const part of parts) {
      if (!part.sharedExtra) continue;
      for (const seg of part.segments) {
        const r = seg.record as { [k: string]: unknown };
        r.extraTransform = ValueCell.create(new Float32Array(0));
        r.hasExtraTransform = ValueCell.create(false);
      }
      part.sharedExtra = undefined;
    }
  }

  /**
   * The kept state on `o`, the values of a segment of `seg` made just now (as a units visual's `setVisualState` on a
   * new one)
   */
  function applyStateTo(o: GraphicsRenderObject, seg: SegmentPlace) {
    const { overpaint, transparency, emissive, substance, clipping, wiggle, transform, unitTransforms } = state;
    if (transform && !Mat4.isIdentity(transform)) Visual.setTransform(o, transform);
    if (unitTransforms) Visual.setTransform(o, undefined, unitTransforms.getSymmetryGroupTransforms(seg.group));
    if (overpaint?.layers.length) Visual.setOverpaint(o, overpaint, lociApply(seg), true);
    if (transparency?.layers.length) Visual.setTransparency(o, transparency, lociApply(seg), true);
    if (emissive?.layers.length) Visual.setEmissive(o, emissive, lociApply(seg), true);
    if (substance?.layers.length) Visual.setSubstance(o, substance, lociApply(seg), true);
    if (clipping?.layers.length) Visual.setClipping(o, clipping, lociApply(seg), true);
    if (wiggle?.layers.length) Visual.setWiggle(o, wiggle, lociApply(seg), true);
  }

  /** The kept state on `part`'s object and shared values: visibility, pickability, theme strength */
  function applyPartState(part: Part<G>) {
    const { renderableState } = part;
    if (state.visible !== undefined) renderableState.visible = state.visible;
    if (state.pickable !== undefined) renderableState.pickable = state.pickable;
    if (state.themeStrength) {
      Visual.setThemeStrength({ values: part.shared } as GraphicsRenderObject, state.themeStrength);
    }
  }

  function updateObject(part: Part<G>, props: PD.Values<P>, graphics: GraphicsContext | undefined) {
    // nothing of the records nor of how they are culled changed (a clip moved): not a pass over every record
    if (
      part.object &&
      graphics === part.objectGraphics &&
      part.objectRecordsVersion === recordsVersion.ref.value &&
      part.objectCulling === (props as { gpuCulling?: GpuCullingValue }).gpuCulling
    ) {
      return;
    }
    part.objectRecordsVersion = recordsVersion.ref.value;
    part.objectCulling = (props as { gpuCulling?: GpuCullingValue }).gpuCulling;
    part.drawn = part.segments.filter((s) => s.groupCount * s.instanceCount > 0);
    const records = part.drawn.map((s) => s.record);
    if (records.length === 0) {
      disposeObject(part);
      return;
    }
    const p = props as { renderMerged?: RenderMergedValue; gpuCulling?: GpuCullingValue };
    // merged here whatever `renderMerged` says (the representation chose this visual by it): how it is culled
    const { gpu, rasterOnly } = resolveMergedCulling(
      type,
      records,
      true,
      p.gpuCulling ?? false,
      canCullOnGpu(graphics),
      graphics,
    );
    const key = `${gpu}|${rasterOnly}`;
    const { object, objectRecords } = part;
    const same =
      !!object &&
      key === part.objectKey &&
      graphics === part.objectGraphics &&
      records.length === objectRecords.length &&
      records.every((r, i) => r === objectRecords[i]);
    if (same) return;
    disposeObject(part);
    part.object = createMergedRecordsRenderObject(type, part.shared as any, records, part.renderableState, materialId, {
      gpuCulling: gpu,
      rasterOnly,
      recordsVersion,
    });
    part.objectRecords = records;
    part.objectKey = key;
    part.objectGraphics = graphics;
    geometryVersion += 1;
  }

  /** the segments a loci reaches: by its units where it has fewer than there are segments, else all */
  function getMarkedSegments(loci: Loci): Iterable<Segment<G>> {
    if (!byUnitId) {
      byUnitId = new Map();
      for (const seg of segments) for (const u of seg.group.units) byUnitId.set(u.id, seg);
    }
    if (StructureElement.Loci.is(loci) && loci.elements.length < segments.length) {
      const marked = new Set<Segment<G>>();
      for (const e of loci.elements) {
        const seg = byUnitId.get(e.unit.id);
        if (seg) marked.add(seg);
      }
      return marked;
    }
    if (Bond.isLoci(loci) && loci.bonds.length < segments.length) {
      const marked = new Set<Segment<G>>();
      for (const b of loci.bonds) {
        const a = byUnitId.get(b.aUnit.id);
        if (a) marked.add(a);
        const c = byUnitId.get(b.bUnit.id);
        if (c) marked.add(c);
      }
      return marked;
    }
    return segments;
  }

  return {
    builder,
    get renderObjects() {
      return renderObjects;
    },
    get geometryVersion() {
      return geometryVersion;
    },
    get groupCount() {
      let count = 0;
      for (const seg of segments) count += seg.groupCount * seg.instanceCount;
      return count;
    },
    get mergeable() {
      return !unmergeable && parts.length > 0 && parts.every((part) => canMergeValues(type, [part.shared]));
    },
    async createOrUpdate(ctx, theme, props, structure, groups) {
      const { graphics } = ctx;
      unmergeable = false;
      // after a backend switch or restore: all made again, for the new context. A new instance granularity: the
      // segments' parts again (`auto` resolves by a group's counts)
      if (
        (graphics &&
          parts.some((p) => p.object && graphics !== p.objectGraphics && hasForeignResources(p.object, graphics))) ||
        (currentProps && props.instanceGranularity !== currentProps.instanceGranularity)
      ) {
        reset();
      }
      lociStructure = structure;
      lociGranularity = props.instanceGranularity;

      // the same groups (a theme or props update): the segments in their order, not looked up by group
      const sameGroups = groups === currentGroups && segments.length === groups.length;
      const previous = byGroup;
      const next = sameGroups ? byGroup : createSegmentMap<G>();
      const nextSegments: Segment<G>[] = sameGroups ? segments : [];
      // the props, theme and structure as each segment was made or last updated; a clip moved alone (a drag, each
      // frame): the parts' shared values' clip once, not each segment's values again (the clip cells shared, not a
      // record's: 94k segments 3 s a move)
      const sameData = theme === currentTheme && structure === currentStructure;
      const clipOnly = sameGroups && sameData && !!currentProps && BaseGeometry.isClipChangeOnly(props, currentProps);
      if (clipOnly) for (const part of parts) BaseGeometry.updateClip(part.shared as any, props.clip);
      const unchanged = clipOnly || (props === currentProps && sameData);
      // the props changed alone: a part's first segment updated, and where that wrote its shared values alone (no
      // cell of its record, nothing made again nor recolored) the part's other segments not (their values the same
      // shared ones: 94k segments seconds a change, e.g. a material's)
      const propsOnly = !unchanged && sameGroups && sameData && !!currentProps;
      /** by part, whether its first segment's update wrote its shared values alone (`propsOnly`) */
      const sharedOnly = new Map<Part<G>, boolean>();
      /**
       * the parts checked to merge in this update, once each after its first segment's update: one that does not ends
       * the update there (the visual left to be destroyed, `mergeable` false), not every group's geometry made first
       */
      const checked = new Set<Part<G>>();
      const merges = (seg: Segment<G>) => {
        if (checked.has(seg.part)) return true;
        checked.add(seg.part);
        return canMergeValues(type, [seg.part.shared]);
      };
      let written = !sameGroups;
      // what this update made and which of the previous segments it kept: the visual's own map and segments untouched
      // till it ends (one that does not, by a part not to merge, a throw or an abort, disposes what it made and leaves
      // the visual as it was)
      const made: Segment<G>[] = [];
      const kept = new Set<Segment<G>>();
      const partsBefore = new Set(parts);
      let completed = false;
      try {
        for (let i = 0, il = groups.length; i < il; ++i) {
          const group = groups[i];
          const structureGroup: StructureGroup = { group, structure };
          let seg = sameGroups ? segments[i] : previous.get(group);
          if (seg) {
            if (!sameGroups) kept.add(seg);
            if (propsOnly && seg.group === group) {
              const probed = sharedOnly.get(seg.part);
              if (probed) continue;
              if (probed === undefined) {
                const versions = recordVersions(seg);
                const updating = updateSegment(ctx, seg, structureGroup, props, theme);
                if (isPromiseLike(updating)) await updating;
                const clean = !hasRecordUpdates(updateState) && isSameVersions(versions, recordVersions(seg));
                sharedOnly.set(seg.part, clean);
                if (!clean) written = true;
                if (!merges(seg)) {
                  unmergeable = true;
                  return;
                }
                continue;
              }
            }
            if (!unchanged || seg.group !== group) {
              written = true;
              const updating = updateSegment(ctx, seg, structureGroup, props, theme);
              if (isPromiseLike(updating)) await updating;
            }
          } else {
            const geometry = createGeometry(ctx, structureGroup, props, theme);
            let segment: Segment<G> | undefined;
            const preparing = Visual.withPreparedGeometry(geometry, utils, theme, graphics, (g) => {
              segment = createSegment(ctx, structureGroup, g, props, theme);
            });
            if (isPromiseLike(preparing)) await preparing;
            seg = segment!;
            made.push(seg);
            geometryVersion += 1;
            // a new segment of a marked representation: its markers as made, clear (as a new units visual's)
          }
          if (!merges(seg)) {
            unmergeable = true;
            return;
          }
          if (!sameGroups) {
            next.set(group, seg);
            nextSegments.push(seg);
          }
          if (ctx.runtime.shouldUpdate) {
            await ctx.runtime.update({ message: 'Creating or updating merged units visual', current: i, max: il });
          }
        }
        completed = true;
      } finally {
        if (!completed) {
          for (const seg of made) builder.dispose?.(seg.geometry);
          // the parts made for them alone: no object yet (made with the render objects), dropped
          parts = parts.filter((part) => partsBefore.has(part));
        }
      }
      if (written) recordsWritten();
      if (!sameGroups) {
        previous.forEach((seg) => {
          if (!kept.has(seg)) builder.dispose?.(seg.geometry);
        });
        byGroup = next;
        segments = nextSegments;
        byUnitId = undefined;
        for (const part of parts) part.segments = [];
        for (const seg of segments) seg.part.segments.push(seg);
        for (const part of parts) if (part.segments.length === 0) disposeObject(part);
        parts = parts.filter((part) => part.segments.length > 0);
      }
      currentGroups = groups;

      for (const part of parts) {
        if (props !== currentProps) utils.updateRenderableState(part.renderableState, props);
        applyPartState(part);
      }
      currentProps = props;
      currentTheme = theme;
      currentStructure = structure;
      renderObjects = [];
      for (const part of parts) {
        updateObject(part, props, graphics);
        if (part.object) renderObjects.push(part.object);
      }
    },
    getLoci(pickingId) {
      const part = parts.find((p) => p.object?.id === pickingId.objectId);
      const object = part?.object;
      if (!part || !object) return EmptyLoci;
      const { instanceBases, instanceCounts, count } = object.merged.segments;
      const { instanceId } = pickingId;
      // the last segment starting at or before the instance
      let lo = 0;
      let hi = count - 1;
      while (lo < hi) {
        const mid = (lo + hi + 1) >> 1;
        if (instanceBases[mid] <= instanceId) lo = mid;
        else hi = mid - 1;
      }
      const base = instanceBases[lo];
      if (instanceId < base || instanceId >= base + instanceCounts[lo]) return EmptyLoci;
      const seg = part.drawn[lo];
      let id: PickingId = { ...pickingId, instanceId: instanceId - base };
      if (resolveInstanceGranularity(lociGranularity, seg.groupCount, seg.instanceCount)) {
        id = { ...id, groupId: PickingId.Null };
      }
      return builder.getLoci(id, structureGroupOf(seg), object.id);
    },
    mark(loci, action) {
      if (renderObjects.length === 0) return false;
      const every = isEveryLoci(loci);
      let changed = false;
      for (const seg of getMarkedSegments(loci)) {
        const visualLoci = getVisualLoci(seg)(loci);
        if (!visualLoci) continue;
        // a previous mark made for a segment's own loci (a hover's), not for each of all of them by a loci of all: one
        // kept is followed by every mark (a stale one would restore a status since changed)
        const previous = every
          ? seg.previousMark
          : (seg.previousMark ??= { loci: EmptyLoci, action: MarkerAction.None, status: -1 });
        if (Visual.mark(asObject(seg), visualLoci, action, lociApply(seg), previous)) changed = true;
      }
      return changed;
    },
    eachLocation(cb) {
      for (const seg of segments) {
        const it = builder.createLocationIterator(structureGroupOf(seg), currentProps!);
        it.reset();
        while (it.hasNext) {
          const { location, isSecondary } = it.move();
          cb(location, isSecondary);
        }
      }
    },
    setState(next, current) {
      const { visible, pickable, overpaint, transparency, emissive, substance, clipping, wiggle, themeStrength } = next;
      if (visible !== undefined) state.visible = visible;
      if (pickable !== undefined) state.pickable = pickable;
      if (themeStrength !== undefined) state.themeStrength = themeStrength;
      for (const part of parts) applyPartState(part);
      for (const [layer, set] of [
        [overpaint, Visual.setOverpaint],
        [transparency, Visual.setTransparency],
        [emissive, Visual.setEmissive],
        [substance, Visual.setSubstance],
        [clipping, Visual.setClipping],
        [wiggle, Visual.setWiggle],
      ] as const) {
        if (layer === undefined) continue;
        for (const seg of segments) (set as any)(asObject(seg), layer, lociApply(seg), true);
      }
      if (
        [overpaint, transparency, emissive, substance, clipping, wiggle, next.transform, next.unitTransforms].some(
          (v) => v !== undefined,
        )
      ) {
        recordsWritten();
      }
      if (overpaint !== undefined) state.overpaint = overpaint;
      if (transparency !== undefined) state.transparency = transparency;
      if (emissive !== undefined) state.emissive = emissive;
      if (substance !== undefined) state.substance = substance;
      if (clipping !== undefined) state.clipping = clipping;
      if (wiggle !== undefined) state.wiggle = wiggle;
      if (next.transform !== undefined) {
        state.transform = next.transform;
        for (const seg of segments) Visual.setTransform(asObject(seg), next.transform);
      }
      if ('unitTransforms' in next && next.unitTransforms !== undefined) {
        const unitTransforms = next.unitTransforms;
        if (unitTransforms) {
          ownExtraTransforms();
          for (const seg of segments) {
            Visual.setTransform(asObject(seg), undefined, unitTransforms.getSymmetryGroupTransforms(seg.group));
          }
        } else if (unitTransforms !== current.unitTransforms) {
          for (const seg of segments) Visual.setTransform(asObject(seg), undefined, null);
        }
        state.unitTransforms = unitTransforms;
      }
    },
    destroy() {
      reset();
    },
  };
}
