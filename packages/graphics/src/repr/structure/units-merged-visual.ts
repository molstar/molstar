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
import { isDebugMode } from '@molstar/core/util/debug';
import { Bond, type Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import { EmptyLoci, isEveryLoci, type Loci } from '@molstar/model/model/loci';
import type { Geometry } from '@molstar/graphics/geo/geometry/geometry';
import {
  BaseGeometry,
  getColorSmoothingProps,
  type GpuCullingValue,
  hasColorSmoothingProp,
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
  withMergeLodKey,
  createEmptyMergeMember,
  getMergeSegmentKeys,
  isMergeableType,
  type MergeableValues,
  type MergeRecord,
  type MergedType,
} from '@molstar/graphics/gl/renderable/merged';
import { RenderableState } from '@molstar/graphics/gl/renderable';
import { type DataArray, releaseGpuData } from '@molstar/graphics/gl/renderable/util';
import { prepareLayerSampling } from '@molstar/graphics/gl/compute/layer-sampling';
import { releaseGpuMeshVertices } from '@molstar/graphics/geo/geometry/mesh/gpu-vertices';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { Transparency } from '@molstar/graphics/theme/transparency';
import { Emissive } from '@molstar/graphics/theme/emissive';
import { Substance } from '@molstar/graphics/theme/substance';
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
 * A segment's extra transform cells (of `MergeSegmentKeys`): its part's shared empty ones until unit transforms are set
 * (`ownExtraTransforms`). The matrix (the representation's transform) and the poses (none for units) the shared values'
 */
const ExtraTransformKeys = ['extraTransform', 'hasExtraTransform'] as const;

/** A cell's value to compare writes of it by (`checkSharedWrites`): an array by its length and first entries */
function writeFingerprint(v: unknown): unknown {
  if (v === null || typeof v !== 'object' || !ArrayBuffer.isView(v) || v instanceof DataView) return v;
  const a = v as unknown as ArrayLike<number>;
  let s = `${a.length}:`;
  for (let i = 0, il = Math.min(a.length, 16); i < il; ++i) s += `${a[i]},`;
  return s;
}
/** the shared keys told written apart by segments (debug mode), once a key */
const warnedSharedWrites = new Set<string>();

/**
 * The layers smoothed on surfaces (`Visual.setOverpaint` and the others with a smoothing context), each its type's define:
 * a shared cell, in the merge key, so one type a part
 */
const SmoothedLayerKinds = ['overpaint', 'transparency', 'emissive', 'substance'] as const;
type SmoothedLayerKind = (typeof SmoothedLayerKinds)[number];
const SmoothedLayerDefines = {
  overpaint: 'dOverpaintType',
  transparency: 'dTransparencyType',
  emissive: 'dEmissiveType',
  substance: 'dSubstanceType',
} as const;
const SmoothedLayerTextures = {
  overpaint: 'tOverpaint',
  transparency: 'tTransparency',
  emissive: 'tEmissive',
  substance: 'tSubstance',
} as const;
type LayersOf = { readonly layers: ReadonlyArray<{ readonly loci: Loci }> };
type SmoothedLayers = { [K in SmoothedLayerKind]?: LayersOf };
type LayerSetter = (
  o: GraphicsRenderObject,
  layers: LayersOf,
  lociApply: Visual.LociApply,
  clear: boolean,
  smoothing?: Visual.SmoothingContext,
) => void;
/** whether two layers of a kind are alike (`Overpaint.areEqual`, ...) */
const SmoothedLayerEquals = {
  overpaint: Overpaint.areEqual,
  transparency: Transparency.areEqual,
  emissive: Emissive.areEqual,
  substance: Substance.areEqual,
} as unknown as { [K in SmoothedLayerKind]: (a: LayersOf, b: LayersOf) => boolean };

const SmoothedLayerSetters = {
  overpaint: Visual.setOverpaint,
  transparency: Visual.setTransparency,
  emissive: Visual.setEmissive,
  substance: Visual.setSubstance,
} as unknown as { [K in SmoothedLayerKind]: LayerSetter };

/** Whether `a` and `b` smooth alike (`smoothColors`: its name and params) */
function isSameSmoothing(a: PD.Values<any>, b: PD.Values<any>) {
  const x = a.smoothColors;
  const y = b.smoothColors;
  return x === y || (x?.name === y?.name && JSON.stringify(x?.params) === JSON.stringify(y?.params));
}

/** The layers of `state` that smoothing may reach (with layers), by kind */
function getSmoothedLayers(state: Partial<StructureRepresentationState>): SmoothedLayers {
  const layers: { [K in SmoothedLayerKind]?: LayersOf } = {};
  for (const kind of SmoothedLayerKinds) if (state[kind]?.layers.length) layers[kind] = state[kind];
  return layers;
}

/**
 * The segments of one merge key (`getMergeKey`) and their merged object: groups of the same props and theme differ in
 * it where their counts resolve the instance granularity apart (`auto`: a large unit of many copies by instance, a
 * small one by group and instance), or their geometry's defines or levels. Its shared values those of the first
 * segment made, its segments' own cells (`recordKeys`) those of its first segment now (`bindShared`): read through the
 * shared values (`granularity`, `getMergeKey`) they follow a live segment, not one gone
 */
type Part<G extends Geometry> = {
  readonly shared: MergeableValues;
  /** the keys its segments own (`getMergeSegmentKeys`): their records' and their own beside */
  readonly recordKeys: readonly string[];
  readonly recordKeySet: ReadonlySet<string>;
  /** the segment whose own cells the shared values hold (`bindShared`) */
  sharedFrom: Segment<G> | undefined;
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
  /** its merge key's part (`regroupByKey`) */
  part: Part<G>;
  /**
   * its cells: the merged object's record of it, and its own beside (`MergeSegmentKeys`); its part's shared values
   * its prototype, so that the writers read and write it as a units visual's values (`UnitsMergedVisual`)
   */
  readonly record: MergeRecord;
  /** the record, as a render object's values for the writers that take one (`Visual.mark`, the layer setters) */
  readonly values: MergeableValues;
  /** of its location iterator at its last update */
  groupCount: number;
  instanceCount: number;
  /**
   * the props, theme and structure it was last written with: its update from these (`updateSegment`), so that one
   * of an update that did not end (a throw, an abort) is brought from where it is
   */
  props: PD.Values<any>;
  theme: Theme;
  structure: Structure;
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
  /**
   * advanced where a segment's geometry was made again (the layers' nodes smooth again by it), not where a merged object
   * is another over the same segments (a layer regrouping them)
   */
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
  /**
   * the representation's state: kept, and applied to every segment made later; its layers may move segments to other
   * parts (smoothed layers, `writeLayers`), `renderObjects` then another list
   */
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
   * the updates running (`createOrUpdate`, suspended at its awaits): a `setState` meanwhile sets what its parts read
   * (`visible`, `pickable`, `themeStrength`) and keeps the rest (`pendingState`, `pendingCurrent`) for when they end;
   * written meanwhile, its layers would regroup the parts an update is making (`writeLayers`, `relistParts`) and miss the
   * segments it has yet to make
   */
  let updateDepth = 0;
  let pendingState: Partial<StructureRepresentationState> | undefined;
  let pendingCurrent: StructureRepresentationState | undefined;
  /**
   * advanced where the records' cells are written but the observed ones (a mark's): the merged objects then read all
   * of their members at their next sync, else those of changed observed cells (`Merged.recordsUnchanged`)
   */
  const recordsVersion = ValueCell.create(0);
  function recordsWritten() {
    ValueCell.update(recordsVersion, recordsVersion.ref.value + 1);
  }

  /** the props, theme, structure and context of the segments' last update */
  let currentProps: PD.Values<P> | undefined;
  let currentTheme: Theme | undefined;
  let currentStructure: Structure | undefined;
  let currentGraphics: GraphicsContext | undefined;
  /**
   * whether the last update ended (every segment written with its props, theme and structure): the fast paths of the
   * next one (not a segment's update state each) only then; one that did not leaves segments of either
   */
  let consistent = true;
  /** the next update a rebuild (a part's segments found of other merge keys, `keySplit`) */
  let rebuildNext = false;
  /** the segments an update wrote (`updateSegment`): the records then read again (`recordsWritten`) */
  let segmentWrites = 0;
  /** the segments' geometries made again by updates (`updateSegment`): their smoothed layers then written again */
  let geometriesRemade = 0;
  /** the groups of the segments' last update, in the segments' order */
  let currentGroups: ReadonlyArray<Unit.SymmetryGroup> | undefined;
  /** the structure and instance granularity the loci are applied with: the update's from its start */
  let lociStructure: Structure | undefined;
  let lociGranularity: PD.Values<P>['instanceGranularity'] = 'auto';

  /** the representation's state, applied to segments made later (`setState`) */
  const state: Partial<StructureRepresentationState> = {};

  /** `seg` dropped: its geometry, and what writes its layers on the GPU (`DataArray.gpu`) */
  function disposeSegment(seg: Segment<G>) {
    if (keptGeometries.delete(seg.geometry)) {
      // its geometry another segment's now (a rebuild's, `keepGeometry`): its vertices on the GPU with it
    } else {
      builder.dispose?.(seg.geometry);
      // its mesh's vertices on the GPU (the smoothed layers' sampling), kept by its position cell across layers
      releaseGpuMeshVertices(seg.values);
    }
    const values = seg.values as unknown as { [k: string]: ValueCell<DataArray<any>> | undefined };
    for (const kind of SmoothedLayerKinds) releaseGpuData(values[SmoothedLayerTextures[kind]]?.ref.value);
    // its colors smoothed on the GPU (`processValues` per vertex) as the layers
    releaseGpuData(values.tColor?.ref.value);
  }

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

  /**
   * the geometries a rebuild's segments took from the segments they replace (`keepGeometry`): disposed with neither,
   * the other one's then (the old segments' once the rebuild ends, the made ones' where it fails)
   */
  const keptGeometries = new Set<G>();

  /**
   * `old`'s geometry for a rebuild's segment of `structureGroup` where an update of `old` would not make it again (its
   * props, theme, structure and context as `setUnitsUpdateState` reads them): an instance granularity or a key split
   * remakes the parts, not the surfaces; another context the ones that depend on it
   */
  function keepGeometry(
    ctx: VisualContext,
    old: Segment<G> | undefined,
    structureGroup: StructureGroup,
    props: PD.Values<P>,
    theme: Theme,
  ): G | undefined {
    if (!old || keptGeometries.has(old.geometry)) return;
    VisualUpdateState.reset(updateState);
    updateState.info = old.info;
    setUnitsUpdateState(
      updateState,
      builder,
      old.record as any,
      old.geometry,
      props,
      old.props as PD.Values<P>,
      theme,
      old.theme,
      structureGroup,
      { group: old.group, structure: old.structure },
      ctx.graphics,
      false,
    );
    if (updateState.createGeometry) return;
    keptGeometries.add(old.geometry);
    return old.geometry;
  }

  /** The part of `values`' merge key, a new one (`values` its shared values) where there is none */
  function getPart(values: MergeableValues, props: PD.Values<P>): Part<G> {
    const key = getMergeKey(values);
    for (const part of parts) if (getMergeKey(part.shared) === key) return part;
    return createPart(values, getMergeSegmentKeys(type, values), props);
  }

  /** A part of the shared values `shared`, its segments' own keys `recordKeys` */
  function createPart(shared: MergeableValues, recordKeys: string[], props: PD.Values<P>): Part<G> {
    const part: Part<G> = {
      shared,
      recordKeys,
      recordKeySet: new Set(recordKeys),
      sharedFrom: undefined,
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
    // smoothing per vertex, merged (`UnitsVisualGeometryBuilder.processValues`): sampled on the GPU where it computes
    builder.processValues?.(values as any, geometry, props, theme, ctx.graphics, true);
    // the kept state on its values before its part is picked: their layers' defines as its part's (one made with them),
    // its merge key that of the segments made before it
    const place: SegmentPlace = {
      group: structureGroup.group,
      groupCount: locationIt.groupCount,
      instanceCount: locationIt.instanceCount,
    };
    applyStateTo({ values } as unknown as GraphicsRenderObject, place, {
      geometry,
      props,
      graphics: ctx.graphics,
      perVertexOnGpu: true,
    });
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
      props,
      theme,
      structure: structureGroup.structure,
    };
    // the first segment of a new part: its own cells the shared values' (its values those)
    part.sharedFrom ??= seg;
    return seg;
  }

  /**
   * `part`'s shared values holding its first segment's own cells, where that is another than they hold (the one they
   * held gone, or moved to another part): what is read of a segment through them (the granularity of a layer, the merge
   * key's levels) a live segment's
   */
  function bindShared(part: Part<G>) {
    const first = part.segments[0];
    if (!first || first === part.sharedFrom) return;
    const shared = part.shared as unknown as { [k: string]: unknown };
    const record = first.record as { [k: string]: unknown };
    for (const k of part.recordKeys) shared[k] = record[k];
    part.sharedFrom = first;
  }

  /**
   * Every segment in the part of its merge key, after an update that changed segments' levels (a size theme or a level
   * setting: each segment's own, its key's levels, `withMergeLodKey`): a segment whose key left its part's moves to a
   * part of its key, its record the same (its cells, markers and layers kept) over that part's shared values; parts of
   * one key become one. A key of no part a new part, its shared values a copy of the segment's part's (their cells of
   * their own, `createEmptyMergeMember`; its own cells its first segment's, `bindShared`). Whether a segment moved
   */
  function regroupByKey(props: PD.Values<P>) {
    const byKey = new Map<string, Part<G>>();
    const claimed = new Set<Part<G>>();
    let moved = false;
    for (const seg of segments) {
      const from = seg.part;
      const key = withMergeLodKey(getMergeKey(from.shared), seg.record.lodLevels);
      let to = byKey.get(key);
      if (!to) {
        if (!claimed.has(from)) {
          to = from;
        } else {
          const shared = createEmptyMergeMember(type, from.shared);
          to = createPart(shared, [...from.recordKeys], props);
        }
        byKey.set(key, to);
        claimed.add(to);
      }
      if (to !== from) {
        Object.setPrototypeOf(seg.record, to.shared);
        seg.part = to;
        moved = true;
      }
    }
    return moved;
  }

  /** The parts' segments again after segments moved: parts without segments gone, their shared values a live segment's */
  function relistParts() {
    byUnitId = undefined;
    for (const part of parts) part.segments = [];
    for (const seg of segments) seg.part.segments.push(seg);
    for (const part of parts) if (part.segments.length === 0) disposeObject(part);
    parts = parts.filter((part) => part.segments.length > 0);
    for (const part of parts) bindShared(part);
  }

  /** The parts' merged objects, made where their records changed, as the visual's render objects */
  function listObjects(props: PD.Values<P>, graphics: GraphicsContext | undefined) {
    renderObjects = [];
    for (const part of parts) {
      updateObject(part, props, graphics);
      if (part.object) renderObjects.push(part.object);
    }
  }

  /**
   * Whether `layers` smooth on `seg`: a mesh of a surface's resolution, `smoothColors` resolving on (`Visual.setOverpaint`
   * and the others: per vertex and instance), by group (not instance) and reached; else they are its plain type
   */
  /** whether a layer may smooth on a segment (`smoothsOn`: the props' smoothing), or one is smoothed (per vertex) */
  function maySmooth() {
    if (currentProps && hasColorSmoothingProp(currentProps as PD.Values<any>)) return true;
    for (const seg of segments) {
      const values = seg.values as unknown as { [k: string]: ValueCell<string> };
      for (const kind of SmoothedLayerKinds)
        if (values[SmoothedLayerDefines[kind]].ref.value === 'vertexInstance') return true;
    }
    return false;
  }

  function smoothsOn(seg: Segment<G>, layers: LayersOf) {
    if (layers.layers.length === 0 || seg.geometry.kind !== 'mesh') return false;
    const props = seg.props as PD.Values<any>;
    if (!hasColorSmoothingProp(props)) return false;
    if ((seg.values as MergeableValues).uInstanceGranularity.ref.value) return false;
    const resolution = (seg.geometry as { meta?: { resolution?: number } }).meta?.resolution;
    if (!getColorSmoothingProps(props.smoothColors, true, resolution)) return false;
    return Visual.reachesAny(layers.layers, lociApply(seg));
  }

  /**
   * per kind the layer last written to every segment, and the geometry it was written for (`geometryVersion`): the same
   * set again for that geometry not written again (`setState`). Why: a geometry made again smooths the layers again
   * here (`resmooth`), then the layers' nodes, recreated by the new `geometryVersion`, set them alike (a trajectory
   * frame, a surface's props: all units' layers twice, 33 ms on WebGPU, seconds on the CPU)
   */
  let layersWritten: { [K in SmoothedLayerKind]?: { readonly layer: LayersOf; readonly geometry: number } } = {};

  /**
   * Writes `layers` (the kinds given; the others kept) into every segment, smoothed where they smooth on it
   * (`smoothsOn`). The parts first: a part's layer types are its shared cells (in its merge key), so its segments of
   * other smoothed layers move to parts of their own (the most stay), each part's segments then writing the same types;
   * parts of one key after the writes joined (`regroupByKey`). Whether segments moved: the parts' objects then made again
   */
  function writeLayers(layers: SmoothedLayers): boolean {
    const kinds = SmoothedLayerKinds.filter((k) => layers[k] !== undefined);
    if (kinds.length === 0 || segments.length === 0) return false;
    for (const kind of kinds) layersWritten[kind] = { layer: layers[kind]!, geometry: geometryVersion };
    // smoothing per vertex: meshes alone (the surfaces merged on the CPU path), where a layer may smooth or a segment's
    // did (its types back); else not grouped and regrouped by it (a cartoon's every layer change: a pass over its segments)
    const flags = type === 'mesh' && maySmooth() ? new Int32Array(segments.length) : undefined;
    let split = false;
    if (flags) {
      const byPart = new Map<Part<G>, Map<number, Segment<G>[]>>();
      for (let i = 0, il = segments.length; i < il; ++i) {
        const seg = segments[i];
        let f = 0;
        for (let k = 0; k < SmoothedLayerKinds.length; ++k) {
          const kind = SmoothedLayerKinds[k];
          const layer = layers[kind];
          const smoothed = layer
            ? smoothsOn(seg, layer)
            : (seg.values as unknown as { [k: string]: ValueCell<string> })[SmoothedLayerDefines[kind]].ref.value ===
              'vertexInstance';
          if (smoothed) f |= 1 << k;
        }
        flags[i] = f;
        let groups = byPart.get(seg.part);
        if (!groups) byPart.set(seg.part, (groups = new Map()));
        const group = groups.get(f);
        if (group) group.push(seg);
        else groups.set(f, [seg]);
      }
      for (const [part, groups] of byPart) {
        if (groups.size < 2) continue;
        split = true;
        let keep = -1;
        let most = -1;
        for (const [f, segs] of groups) {
          if (segs.length > most) {
            most = segs.length;
            keep = f;
          }
        }
        for (const [f, segs] of groups) {
          if (f === keep) continue;
          const to = createPart(createEmptyMergeMember(type, part.shared), [...part.recordKeys], currentProps!);
          for (const seg of segs) {
            Object.setPrototypeOf(seg.record, to.shared);
            seg.part = to;
          }
        }
      }
    }
    for (const kind of kinds) {
      const set = SmoothedLayerSetters[kind];
      const layer = layers[kind]!;
      const bit = 1 << SmoothedLayerKinds.indexOf(kind);
      for (let i = 0, il = segments.length; i < il; ++i) {
        const seg = segments[i];
        const smoothing =
          flags && flags[i] & bit
            ? { geometry: seg.geometry, props: seg.props, graphics: currentGraphics, perVertexOnGpu: true }
            : undefined;
        set(asObject(seg), layer, lociApply(seg), true, smoothing);
      }
    }
    if (!flags) return false;
    const moved = regroupByKey(currentProps!);
    if (split || moved) relistParts();
    return split || moved;
  }

  /** the versions of `part`'s shared cells not its segments' own, before a segment's update (`checkSharedWrites`) */
  function sharedVersions(part: Part<G>) {
    const shared = part.shared as unknown as { [k: string]: ValueCell<unknown> | undefined };
    const versions = new Map<string, number>();
    for (const k of Object.keys(shared)) {
      if (!part.recordKeySet.has(k)) versions.set(k, shared[k]?.ref?.version ?? -1);
    }
    return versions;
  }

  /**
   * Debug mode: the shared cells `seg`'s update wrote (through its record's prototype) with another value than an
   * earlier segment's of this update: a cell that differs by segment, not its own (`MergeSegmentKeys`), every segment
   * reading the last one's
   */
  function checkSharedWrites(
    seg: Segment<G>,
    before: Map<string, number>,
    written: Map<Part<G>, Map<string, unknown>>,
  ) {
    const shared = seg.part.shared as unknown as { [k: string]: ValueCell<unknown> | undefined };
    let values = written.get(seg.part);
    if (!values) written.set(seg.part, (values = new Map()));
    for (const [k, version] of before) {
      const cell = shared[k];
      if (!cell?.ref || cell.ref.version === version) continue;
      const fingerprint = writeFingerprint(cell.ref.value);
      if (values.has(k) && values.get(k) !== fingerprint && !warnedSharedWrites.has(k)) {
        warnedSharedWrites.add(k);
        console.error(
          `merged units visual: '${k}' written by its segments with different values through their part's shared values (a segment's own cell: MergeSegmentKeys)`,
        );
      }
      values.set(k, fingerprint);
    }
  }

  /**
   * Updates `seg` to `structureGroup`, `props` and `theme` from what it was last written with (its own), by its update
   * state (`setUnitsUpdateState`, as a units visual's: its context's compatibility and its builder's versions too); the
   * geometry made again where needed; nothing where the state is none and what it was written with the same, or
   * another that writes the shared values alone (`sharedAlone`: its part's first segment's update wrote nothing of its
   * own, `createOrUpdate`)
   */
  function updateSegment(
    ctx: VisualContext,
    seg: Segment<G>,
    structureGroup: StructureGroup,
    props: PD.Values<P>,
    theme: Theme,
    sharedAlone = false,
  ): Promise<void> | void {
    const current: StructureGroup = { group: seg.group, structure: seg.structure };
    const currentProps = seg.props as PD.Values<P>;
    VisualUpdateState.reset(updateState);
    updateState.info = seg.info;
    const values = seg.record as any;
    setUnitsUpdateState(
      updateState,
      builder,
      values,
      seg.geometry,
      props,
      currentProps,
      theme,
      seg.theme,
      structureGroup,
      current,
      ctx.graphics,
      false,
    );
    if (
      (props === currentProps || sharedAlone) &&
      theme === seg.theme &&
      structureGroup.group === seg.group &&
      structureGroup.structure === seg.structure &&
      !hasRecordUpdates(updateState) &&
      !updateState.createNew &&
      !updateState.updateLocation
    ) {
      seg.props = props;
      return;
    }
    segmentWrites += 1;
    const apply = (newGeometry?: G) => {
      // the position iterator made where a writer reads it (vertex colors and sizes; not a recolor of a uniform color)
      let positionIt: LocationIterator | undefined;
      const its: UnitsIterators = {
        locationIt: builder.createLocationIterator(current, currentProps),
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
        currentProps,
        theme,
        structureGroup,
        current,
      );
      seg.group = structureGroup.group;
      seg.props = props;
      seg.theme = theme;
      seg.structure = structureGroup.structure;
      if (newGeometry) {
        seg.geometry = newGeometry;
        geometryVersion += 1;
        geometriesRemade += 1;
      }
      seg.groupCount = its.locationIt.groupCount;
      seg.instanceCount = its.locationIt.instanceCount;
      builder.processValues?.(values, seg.geometry, props, theme, ctx.graphics, true);
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
    for (const seg of segments) disposeSegment(seg);
    segments = [];
    byGroup = createSegmentMap<G>();
    byUnitId = undefined;
    layersWritten = {};
    for (const part of parts) disposeObject(part);
    parts = [];
    renderObjects = [];
    currentProps = undefined;
    currentTheme = undefined;
    currentStructure = undefined;
    currentGraphics = undefined;
    currentGroups = undefined;
    consistent = true;
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
  function applyStateTo(o: GraphicsRenderObject, seg: SegmentPlace, smoothing: Visual.SmoothingContext) {
    const { overpaint, transparency, emissive, substance, clipping, wiggle, transform, unitTransforms } = state;
    if (transform && !Mat4.isIdentity(transform)) Visual.setTransform(o, transform);
    if (unitTransforms) Visual.setTransform(o, undefined, unitTransforms.getSymmetryGroupTransforms(seg.group));
    // its values its own: smoothed where reached, its part then of those types (`getPart`, `writeLayers`)
    if (overpaint?.layers.length) Visual.setOverpaint(o, overpaint, lociApply(seg), true, smoothing);
    if (transparency?.layers.length) Visual.setTransparency(o, transparency, lociApply(seg), true, smoothing);
    if (emissive?.layers.length) Visual.setEmissive(o, emissive, lociApply(seg), true, smoothing);
    if (substance?.layers.length) Visual.setSubstance(o, substance, lociApply(seg), true, smoothing);
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
    const p = props as { gpuCulling?: GpuCullingValue };
    // merged here, also one record: how it is culled
    const { gpu, rasterOnly } = resolveMergedCulling(
      type,
      records,
      p.gpuCulling ?? false,
      canCullOnGpu(graphics),
      graphics,
      true,
    );
    const key = `${gpu}|${rasterOnly}`;
    const { object, objectRecords } = part;
    // another context alone keeps the object (the scene carries it over, as a plain one), unless it holds resources of
    // the previous one; how the device culls is in the key
    const foreign = !!object && !!graphics && graphics !== part.objectGraphics && hasForeignResources(object, graphics);
    const same =
      !!object &&
      !foreign &&
      key === part.objectKey &&
      records.length === objectRecords.length &&
      records.every((r, i) => r === objectRecords[i]);
    if (same) {
      part.objectGraphics = graphics;
      return;
    }
    disposeObject(part);
    part.object = createMergedRecordsRenderObject(type, part.shared as any, records, part.renderableState, materialId, {
      gpuCulling: gpu,
      rasterOnly,
      recordsVersion,
    });
    part.objectRecords = records;
    part.objectKey = key;
    part.objectGraphics = graphics;
    // not a new geometry (`geometryVersion`): the segments' own, made again where they are
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

  /** `createOrUpdate`'s update (`updateDepth` counts it, a key split's nested one too) */
  async function update(
    ctx: VisualContext,
    theme: Theme,
    props: PD.Values<P>,
    structure: Structure,
    groups: ReadonlyArray<Unit.SymmetryGroup>,
  ) {
    const { graphics } = ctx;
    unmergeable = false;
    segmentWrites = 0;
    const remadeBefore = geometriesRemade;
    // after a backend switch or restore: all made again, for the new context. A new instance granularity: the
    // segments' parts again (`auto` resolves by a group's counts). Made beside the current segments and parts, which
    // are disposed once the update ends (one that does not leaves them as they were)
    const rebuild =
      rebuildNext ||
      (!!graphics &&
        parts.some((p) => p.object && graphics !== p.objectGraphics && hasForeignResources(p.object, graphics))) ||
      (!!currentProps && props.instanceGranularity !== currentProps.instanceGranularity);
    // asked for by a key split (`keySplit`): asked again where this update ends early (the others' causes remain)
    const rebuildAsked = rebuildNext;
    rebuildNext = false;
    const oldParts = parts;
    const oldSegments = segments;
    if (rebuild) parts = [];
    const lociBefore = { structure: lociStructure, granularity: lociGranularity };
    lociStructure = structure;
    lociGranularity = props.instanceGranularity;

    // the same groups (a theme or props update): the segments in their order, not looked up by group
    const sameGroups = !rebuild && groups === currentGroups && segments.length === groups.length;
    const previous = rebuild ? createSegmentMap<G>() : byGroup;
    /** the segments a rebuild replaces, by group (their geometries kept where they stand, `keepGeometry`) */
    const oldByGroup = byGroup;
    const next = sameGroups ? byGroup : createSegmentMap<G>();
    const nextSegments: Segment<G>[] = sameGroups ? segments : [];
    // the fast paths, not a segment's update state each, where the last update ended (its segments written with
    // its props, theme and structure), in the same context, of a builder with no versions of its own a group
    // (`initUpdateState`: e.g. the secondary structure, the interactions, read by its update state)
    const ended = consistent && !rebuild && graphics === currentGraphics;
    const fast = ended && !builder.initUpdateState;
    // a clip moved alone (a drag, each frame): the parts' shared values' clip once, not each segment's values again
    // (the clip cells shared, not a record's: 94k segments 3 s a move)
    const sameData = theme === currentTheme && structure === currentStructure;
    const clipOnly =
      fast && sameGroups && sameData && !!currentProps && BaseGeometry.isClipChangeOnly(props, currentProps);
    if (clipOnly) for (const part of parts) BaseGeometry.updateClip(part.shared as any, props.clip);
    const unchanged = clipOnly || (fast && props === currentProps && sameData);
    // the props changed alone: a part's first segment updated, and where that wrote its shared values alone (no
    // cell of its record, nothing made again nor recolored) the part's other segments not (their values the same
    // shared ones: 94k segments seconds a change, e.g. a material's); of a builder with versions a group, the others'
    // update state each, written where it says
    const propsOnly = ended && !unchanged && sameGroups && sameData && !!currentProps;
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
    /** whether an update changed a segment's levels, its merge key's (`regroupByKey`) */
    let levelsChanged = false;
    const partsBefore = new Set(parts);
    /**
     * of a builder that writes its values by the segment (`processValues`: smoothing by a surface's resolution), each
     * part's merge key after its first segment's update: a later segment's another wrote its own key into the shared
     * values (`keySplit`), the visual then made again (each segment in the part of its key)
     */
    const firstKeys = builder.processValues ? new Map<Part<G>, string>() : undefined;
    let keySplit = false;
    const checkKey = (seg: Segment<G>) => {
      const key = getMergeKey(seg.part.shared);
      const first = firstKeys!.get(seg.part);
      if (first === undefined) firstKeys!.set(seg.part, key);
      else if (key !== first) keySplit = true;
    };
    /** debug mode: the shared cells written in this update, by part (`checkSharedWrites`) */
    const sharedWrites = isDebugMode ? new Map<Part<G>, Map<string, unknown>>() : undefined;
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
            if (probed && !builder.initUpdateState) {
              // written alike by its part's first segment: these props'
              seg.props = props;
              continue;
            }
            if (probed === undefined) {
              const versions = recordVersions(seg);
              const before = sharedWrites && sharedVersions(seg.part);
              const updating = updateSegment(ctx, seg, structureGroup, props, theme);
              if (isPromiseLike(updating)) await updating;
              if (before) checkSharedWrites(seg, before, sharedWrites!);
              if (firstKeys) checkKey(seg);
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
          if (unchanged && seg.group === group) {
            seg.props = props;
          } else {
            const lodVersion = seg.record.lodLevels?.ref.version ?? -1;
            const before = sharedWrites && sharedVersions(seg.part);
            const writes = segmentWrites;
            const sharedAlone = propsOnly && seg.group === group && sharedOnly.get(seg.part) === true;
            const updating = updateSegment(ctx, seg, structureGroup, props, theme, sharedAlone);
            if (isPromiseLike(updating)) await updating;
            if (segmentWrites !== writes) {
              written = true;
              if (firstKeys) checkKey(seg);
            }
            if (before) checkSharedWrites(seg, before, sharedWrites!);
            if ((seg.record.lodLevels?.ref.version ?? -1) !== lodVersion) levelsChanged = true;
          }
        } else {
          // a rebuild's: the replaced segment's geometry where it stands (a new one not made)
          const kept = rebuild ? keepGeometry(ctx, oldByGroup.get(group), structureGroup, props, theme) : undefined;
          const geometry = kept ?? createGeometry(ctx, structureGroup, props, theme);
          let segment: Segment<G> | undefined;
          const preparing = Visual.withPreparedGeometry(geometry, utils, theme, graphics, (g) => {
            segment = createSegment(ctx, structureGroup, g, props, theme);
          });
          if (isPromiseLike(preparing)) await preparing;
          seg = segment!;
          made.push(seg);
          if (!kept) geometryVersion += 1;
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
        if (rebuildAsked) rebuildNext = true;
        for (const seg of made) disposeSegment(seg);
        // the geometries kept for segments not made (the update ended before them) the old segments' alone
        keptGeometries.clear();
        // the parts made for them alone: no object yet (made with the render objects), dropped; a rebuild's the
        // current ones again
        parts = rebuild ? oldParts : parts.filter((part) => partsBefore.has(part));
        lociStructure = lociBefore.structure;
        lociGranularity = lociBefore.granularity;
        if (segmentWrites > 0) {
          // segments written before it ended: the next update by each one's own (no fast path), the merged objects
          // reading every record again (their layouts as written)
          consistent = false;
          recordsWritten();
        }
      }
    }
    if (written) recordsWritten();
    if (rebuild) {
      for (const seg of oldSegments) disposeSegment(seg);
      keptGeometries.clear();
      for (const part of oldParts) disposeObject(part);
    }
    if (!sameGroups) {
      previous.forEach((seg) => {
        if (!kept.has(seg)) disposeSegment(seg);
      });
      byGroup = next;
      segments = nextSegments;
    }
    if (levelsChanged) regroupByKey(props);
    if (!sameGroups || levelsChanged) relistParts();
    currentGroups = groups;

    for (const part of parts) {
      if (props !== currentProps) utils.updateRenderableState(part.renderableState, props);
      applyPartState(part);
    }
    // smoothed layers written again where what they were smoothed by changed: a geometry made again (its vertex data
    // the old one's), the smoothing props (on or off: the layers' types, and the segments' parts, with them)
    const smoothedLayers = type === 'mesh' ? getSmoothedLayers(state) : undefined;
    const resmooth =
      !!smoothedLayers &&
      Object.keys(smoothedLayers).length > 0 &&
      // not at the first update: its segments made with the layers (`createSegment`)
      !!currentProps &&
      (geometriesRemade !== remadeBefore || !isSameSmoothing(props as PD.Values<any>, currentProps as PD.Values<any>));
    currentProps = props;
    currentTheme = theme;
    currentStructure = structure;
    currentGraphics = graphics;
    consistent = true;
    // the per-vertex sampling's kernels ahead of a smoothed layer (`perVertexOnGpu`), not at its first
    if (type === 'mesh' && graphics && hasColorSmoothingProp(props) && props.smoothColors.name !== 'off') {
      prepareLayerSampling(graphics);
    }
    if (keySplit) {
      // a part's segments of other keys: made again, each in the part of its own (their geometries kept), the
      // layers and the objects by that rebuild alone (not smoothed and listed here first, then thrown away)
      rebuildNext = true;
      await visual.createOrUpdate(ctx, theme, props, structure, groups);
      return;
    }
    if (resmooth) {
      writeLayers(smoothedLayers!);
      recordsWritten();
    }
    listObjects(props, graphics);
  }

  const visual: UnitsMergedVisual<P> = {
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
      updateDepth += 1;
      try {
        await update(ctx, theme, props, structure, groups);
      } finally {
        updateDepth -= 1;
        // the state set while it updated, on what it made (or kept, where it ended early)
        if (updateDepth === 0 && pendingState) {
          const next = pendingState;
          pendingState = undefined;
          visual.setState(next, pendingCurrent!);
        }
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
      if (updateDepth > 0) {
        const { visible: _v, pickable: _p, themeStrength: _t, ...rest } = next;
        if (Object.keys(rest).length > 0) {
          pendingState = { ...pendingState, ...rest };
          pendingCurrent = current;
        }
        return;
      }
      // the layers smoothed on surfaces: segments by their layers' types (`writeLayers`), their parts' objects then
      // a layer as written for this geometry left as it is (`layersWritten`)
      const layers: SmoothedLayers = {};
      for (const kind of SmoothedLayerKinds) {
        const layer = next[kind];
        if (layer === undefined) continue;
        const written = layersWritten[kind];
        if (written && written.geometry === geometryVersion && SmoothedLayerEquals[kind](written.layer, layer))
          continue;
        layers[kind] = layer;
      }
      const regrouped = writeLayers(layers);
      for (const [layer, set] of [
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
      if (regrouped) listObjects(currentProps!, currentGraphics);
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
  return visual;
}
