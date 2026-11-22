/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { UserTiming } from '@molstar/core/task/util/user-timing';
import { Interval, OrderedSet, SortedArray } from '@molstar/core/data/int';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { Boundary } from '@molstar/core/math/geometry/boundary';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { applyMarkerAction, MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  type ParticleList,
  Particle,
  ParticleRotations,
  ParticleTarget,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { type ParticlePage, ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import type { Clip } from '@molstar/core/util/clip';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';
import { type ParticleDensityGrid, ParticleTileOpaque } from '@molstar/model/model/particles/particle-tile-opaque';
import { setOccluderGrid } from '@molstar/graphics/gl/occluder-grid';
import type { ParticleTileSet } from '@molstar/model/model/particles/particle-tiles';
import { arrayExtend } from '@molstar/core/util/array';

import { type Loci as ModelLoci, EmptyLoci, isEveryLoci } from '@molstar/model/model/loci';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { PagedCullChanges } from '@molstar/graphics/gl/compute/paged-cull';
import {
  createMergedRecordsRenderObject,
  type GraphicsRenderObject,
  type MergedGraphicsRenderObject,
} from '@molstar/graphics/gl/render-object';
import {
  createEmptyMergeMember,
  fillMergeRecord,
  getSourceRecordFillKeys,
  pickMergeRecord,
  type MergeableValues,
  type MergedInstanceSource,
  type MergedType,
  type MergeRecord,
} from '@molstar/graphics/gl/renderable/merged';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import { createDataArray, type DataArray } from '@molstar/graphics/gl/renderable/util';
import { TextureChanges } from '@molstar/graphics/gl/renderable/texture-changes';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { getMarkersAverage } from '@molstar/graphics/geo/geometry/marker-data';
import { createObjectColor } from '@molstar/graphics/geo/geometry/color-data';
import { Color } from '@molstar/core/util/color';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { fillParticleColors } from '@molstar/graphics/theme/color/particle-pages';
import { createCopyMarkers, getCopyBlockBudget } from './copy-markers.js';
import {
  type AnyLayers,
  type CopyLayerKey,
  type CopyLayerRoom,
  type CopyState,
  getCopyLayerSize,
  getCopyLayerStates,
  warnCopyLayerRoom,
  writeLayerItems,
} from './copy-layers.js';
import { getTemplateWideLocals } from './expanded.js';
import { CopyLayerBlock, CopyLayerKinds, CopyLayerValue } from '@molstar/graphics/gl/instance-expansion';
import { getTransparencyAverage } from '@molstar/graphics/geo/geometry/transparency-data';
import { getEmissiveAverage } from '@molstar/graphics/geo/geometry/emissive-data';
import { getWiggleAverage } from '@molstar/graphics/geo/geometry/wiggle-data';
import { SizeTheme } from '@molstar/graphics/theme/size';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Visual, type VisualContext } from '../../../visual.js';
import {
  ParticleTargetRepresentationParams,
  type ParticleTargetRepresentationProps,
  createParticleTargetPropsProvider,
  createTargetMember,
  getTargetGeometryChoice,
  type TargetMember,
  type TargetLayers,
  isTargetOccluderOff,
} from './common.js';
import {
  CopyPoseFloats,
  getStoreFormatOf,
  PageFlagsWord,
  TemplatePage,
  PageRecordWords,
  ParticlePageStore,
  StorePageBits,
  StorePageSize,
  type StoreBlock,
  type StorePlacement,
  writeCopyPose,
} from './page-store.js';

// A slot list's targets as one object of a page store (`ParticlePageStore`): a page a block of
// store pages, added and let go by the list's changes, its instances culled and composed on the GPU
// (`MergedInstanceSource.paged`). A target a member (a pool of records, `MergeRecord`), the last of no geometry
// (its sphere the object's: the pages' boxes and the members' reach); the colors and markers layers by store index
// (a page's particles at its store page's), a drawn instance's id its store index (its slot by its page's record).

/** The page store's object of a slot list's targets (`ParticleTargetRepresentation` on GPU culling) */
export interface TargetSlotsObject {
  readonly renderObject: MergedGraphicsRenderObject | undefined;
  readonly store: ParticlePageStore;
  /** `list` (a slot list) drawn: by its changes from the last list where it is that one's next version, else all */
  update(ctx: VisualContext, list: ParticleList, theme: Theme, props: ParticleTargetRepresentationProps): Promise<void>;
  /** the state of what is drawn and of its members (visible, pickable) */
  applyState(apply: (renderObject: GraphicsRenderObject) => void): void;
  /**
   * `clip` set on what is drawn, nothing else changed (a clip object moved: a clipping gizmo's drag; its defines the
   * same, `BaseGeometry.isSameClipDefines`): its uniforms, the cull's clip objects; no page written again
   */
  updateClip(clip: Clip.Props): void;
  /** the particle picked, undefined where `pickingId` is not this object's */
  getLoci(pickingId: PickingId): ModelLoci | undefined;
  mark(loci: ModelLoci, action: MarkerAction): boolean;
  /** the layers given (each its kind's again), kept for the updates to come (`TargetLayers`) */
  setLayers(layers: TargetLayers): void;
  destroy(): void;
}

/** The layers a store holds by store index: its texture, define and type cells, a kind's bytes an item */
const StoreLayers = [
  { key: 'overpaint', tex: 'tOverpaint', define: 'dOverpaint', type: 'overpaintType' },
  { key: 'transparency', tex: 'tTransparency', define: 'dTransparency', type: 'transparencyType' },
  { key: 'clipping', tex: 'tClipping', define: 'dClipping', type: 'clippingType' },
  { key: 'emissive', tex: 'tEmissive', define: 'dEmissive', type: 'emissiveType' },
  { key: 'substance', tex: 'tSubstance', define: 'dSubstance', type: 'substanceType' },
  { key: 'wiggle', tex: 'tWiggle', define: 'dWiggle', type: 'wiggleType' },
] as const;

/** The pool's members past the store's (the room it is made with), and its least */
const PoolHeadroom = 1.5;
const MinPoolMembers = 16;

/** a pool member's color before its target's (colors by member: drawn by none till then) */
const DefaultMemberColor = Color(0xcccccc);

/** The geometry a slot object's targets are drawn with */
const SlotsGeometry = 'structure-spheres';

/**
 * Whether a list with templates is drawn by a store's object (`addNested`): every leaf target of it (its own and its
 * templates', in the flattened list's ids) drawn as the store draws (`SlotsGeometry`); else flattened (a target the
 * store does not draw would be left out)
 */
export function isNestedPagedDrawable(
  list: ParticleList,
  props: ParticleTargetRepresentationProps,
  graphics: GraphicsContext | undefined,
) {
  if (!ParticleNesting.hasTemplates(list)) return false;
  const targetProps = createParticleTargetPropsProvider(props);
  for (const target of ParticleNesting.getIdSpace(list).targetMapping.values()) {
    if (getTargetGeometryChoice(target, targetProps(target), graphics) !== SlotsGeometry) return false;
  }
  return true;
}

/** a target's center where it has none */
const Origin = Vec3();

/** a member's bounds none of its own (the last member's the object's): skipped where the object's are made */
const NoBounds = Sphere3D.create(Vec3.create(NaN, NaN, NaN), NaN);

/** whether the targets `b` holds of `a`'s ids are `a`'s */
function isSameTargets(a: ReadonlyMap<number, ParticleTarget>, b: ReadonlyMap<number, ParticleTarget>) {
  if (a === b) return true;
  for (const [k, v] of b) {
    const before = a.get(k);
    if (before !== undefined && before !== v) return false;
  }
  return true;
}

function emptyChanges() {
  return {
    written: [] as number[],
    released: [] as number[],
    runs: [] as number[],
    extras: [] as number[],
    copies: [] as number[],
    members: new Set<number>(),
    grown: false,
  };
}

export function createTargetPagedObject(materialId: number, graphics: GraphicsContext | undefined): TargetSlotsObject {
  let store = new ParticlePageStore();

  let list: ParticleList | undefined;
  /**
   * the versions of a slot list since the store was laid out (`ParticleSlots.getChanges`): their slots the same
   * particles (a slot kept while it lives), a layer's loci by an earlier one's slots the store's (a view's step)
   */
  let slotVersions = new WeakSet<ParticleList>();
  let currentTheme: Theme | undefined;
  let currentProps: ParticleTargetRepresentationProps | undefined;
  let currentTargets: ReadonlyMap<number, ParticleTarget> | undefined;
  /** the pages held, by start (a block's key) */
  const pages = new Map<number, ParticlePage>();
  /**
   * a template's coarse levels' blocks (`addNested`), by key: the template's block (`of`) and each particle's index in
   * it (what colors and picks it as the template's particle)
   */
  const coarse = new Map<number, { readonly of: number; readonly indices: Int32Array }>();

  let byTarget = new Map<number, TargetMember>();
  let acceptedByTarget = new Map<number, boolean>();

  let object: MergedGraphicsRenderObject | undefined;
  /**
   * the object's members (a pool of records, filled with a target's values: the cells of `fillKeys` a record's own, the
   * others the empty record's, of what the source holds), the last of no geometry; no render object a member
   */
  let pool: MergeRecord[] = [];
  /** the values the members share (a copy of a target member's: defines, uniforms, props, clip), the object's */
  let shared: MergeableValues | undefined;
  let fillKeys: readonly string[] = [];
  let poolType: MergedType = 'spheres';
  let poolState: RenderableState | undefined;
  let poolTargets: number[] = [];
  /** per pool member its geometry's reach from its target's center (the culler's: a particle's sphere its scale times it) */
  let reaches = new Float32Array(0);
  /** per pool member its instance count as its cell holds it, and its bounds (none but the last's) */
  let poolCounts = new Int32Array(0);
  let poolBounds = new Float64Array(0);
  let changedMembers = new Set<number>();
  let countedMembers = new Set<number>();
  /** the store's changes since the culler's last sync (merged across updates) */
  let pending = emptyChanges();
  const objectSphere = Sphere3D();

  // the layers by store index (a texel a store page's particle: its `aInstance`), the merged object's own cells
  const colors = ValueCell.create(createDataArray(StorePageSize, 4, Uint8Array));
  const markers = ValueCell.create(createDataArray(StorePageSize, 1, Uint8Array));
  /**
   * the layers by store index (`StoreLayers`: a flat list's particles, a template's particles in every copy), the
   * merged object's own cells read where their define is on; a texel while off
   */
  const layerCells = new Map<string, ValueCell<DataArray<Uint8Array>>>(
    StoreLayers.map((l) => [l.tex, ValueCell.create(createDataArray(1, getCopyLayerSize(l.key), Uint8Array))]),
  );
  /** the layers given (`setLayers`): applied again where the store is updated */
  const layerState: TargetLayers = {};
  /**
   * the copies' own layers (a list with templates': `getCopyLayerStates`), by the copy markers' copies
   * (`CulledSchema.sCopyLayers`): `CopyLayerKinds` pairs of words a copy, then their bytes
   */
  let copyLayerBuffer = new Uint32Array(0);
  let copyLayerVersion = 0;
  const copyLayerStates = new Map<CopyLayerKey, Map<number, Map<number, CopyState>>>();
  /** warned of copies layered whole since the copies were laid out (`warnCopyLayerRoom`) */
  let copyLayersWarned = false;
  /** per template block's key its template (its own block's, not a coarse level's) */
  const templateLists = new Map<number, ParticleList>();
  const markerCells = {
    markerAverage: ValueCell.create(0),
    markerStatus: ValueCell.create(0),
    uMarker: ValueCell.create(0),
  };
  /**
   * the store's markers' own aggregates (what `Visual.mark` reads and sets), the object's (`markerCells`) by them and
   * the copies' markers (`setObjectMarkers`). Why: the object's status -1 while a copy was marked made every mark of the
   * store's scan its markers whole (a hover's clear, then the average: the yeast cell's 31.5 M template particles and
   * their levels a mark)
   */
  const storeCells = {
    markerAverage: ValueCell.create(0),
    markerStatus: ValueCell.create(0),
    uMarker: ValueCell.create(0),
  };
  const markValues = {
    tMarker: markers,
    ...storeCells,
    uGroupCount: ValueCell.create(1),
    instanceCount: ValueCell.create(0),
    uInstanceGranularity: ValueCell.create(true),
  };
  /**
   * the copies' markers (a list with templates': a template particle drawn once a copy from one store index, the
   * store's marker every copy's), laid out again where the templates or their copies are others (cleared); per top
   * particle placing a copy its copy (the cull's index)
   */
  const copyMarkers = createCopyMarkers(() => getCopyBlockBudget(graphics));
  let copyMarkersVersion = 0;
  /**
   * per top particle placing a copy its copies (the cull's indices: one a template block, its coarse levels' each) and
   * their templates' pages; a level's template particles by their indices in it (`coarse`'s inverted, made as marked)
   */
  let copyOfTop = new Map<number, { copy: number; key: number; pages: Int32Array }[]>();
  let levelLocals = new Map<number, Map<number, number>>();
  /**
   * per template (its block's key) its copies at every level and its blocks (its own and its coarse levels'): what a
   * template's marks (the store's markers of its particles: every copy's) are applied to, or made each copy's
   */
  let templateCopies = new Map<number, { copy: number; key: number; pages: Int32Array }[]>();
  let templateBlocks = new Map<number, { key: number; pages: Int32Array; tops: number }[]>();
  /** a block's template (a coarse level's its template's) */
  const templateOf = (key: number) => coarse.get(key)?.of ?? key;
  const copyMarkersSource = {
    get version() {
      return copyMarkersVersion;
    },
    get words() {
      return copyMarkers.words;
    },
    get bytes() {
      return copyMarkers.bytes;
    },
    get byteLength() {
      return copyMarkers.byteLength;
    },
    takeChanged: () => copyMarkers.takeChanged(),
  };
  function layoutCopyMarkers() {
    copyLayersWarned = false;
    const templates = store.templates;
    copyMarkers.layout(templates.map((t) => ({ pages: t.pages, copyCount: t.copies.tops.length })));
    copyOfTop = new Map();
    levelLocals = new Map();
    templateCopies = new Map();
    templateBlocks = new Map();
    let c = 0;
    for (const t of templates) {
      const of = templateOf(t.key);
      let blocks = templateBlocks.get(of);
      if (!blocks) templateBlocks.set(of, (blocks = []));
      blocks.push({ key: t.key, pages: t.pages, tops: t.copies.tops.length });
      let all = templateCopies.get(of);
      if (!all) templateCopies.set(of, (all = []));
      for (let k = 0; k < t.copies.tops.length; ++k) {
        const top = t.copies.tops[k];
        let copies = copyOfTop.get(top);
        if (!copies) copyOfTop.set(top, (copies = []));
        const entry = { copy: c++, key: t.key, pages: t.pages };
        copies.push(entry);
        all.push(entry);
      }
    }
    copyMarkersVersion += 1;
  }

  function createSource(): MergedInstanceSource {
    return {
      layers: new Map<string, ValueCell<DataArray<any>>>([['tColor', colors], ['tMarker', markers], ...layerCells]),
      markers: markerCells,
      takeMembers: () => {
        const taken = changedMembers;
        changedMembers = new Set();
        return taken;
      },
      counts: {
        instanceCounts: poolCounts,
        spheres: poolBounds,
        take: () => {
          const taken = countedMembers;
          countedMembers = new Set();
          return taken;
        },
      },
      paged: {
        source: store,
        takeChanges: (): PagedCullChanges => {
          const taken = pending;
          pending = emptyChanges();
          return taken;
        },
        get reaches() {
          return reaches;
        },
        copyMarkers: copyMarkersSource,
        copyLayers: {
          get version() {
            return copyLayerVersion;
          },
          get buffer() {
            return copyLayerBuffer;
          },
        },
      },
    };
  }

  /**
   * whether the colors are a member's (the theme's `targetColor` for every target: one color a target, read by the
   * member's segment, `object` colors), not a layer by store index filled a particle at a time; the other way the
   * object made again
   */
  let byMember = false;
  /**
   * `values`' colors one `color` (`createObjectColor`'s cells' values: the shared values' their color type, a record's
   * its texel alone)
   */
  function setMemberColor(values: MergeableValues | MergeRecord, color: Color) {
    const object = createObjectColor(color) as unknown as { [k: string]: ValueCell<any> };
    const cells = values as unknown as { [k: string]: ValueCell<any> | undefined };
    for (const k in object) {
      const cell = cells[k];
      if (cell) ValueCell.update(cell, object[k].ref.value);
    }
  }
  /** whether `theme` gives each of `targets` its color (`known`: those of them it gave one already, not asked again) */
  function isByMember(
    theme: Theme,
    targets: ReadonlyMap<number, ParticleTarget>,
    known?: ReadonlyMap<number, ParticleTarget>,
  ) {
    const targetColor = (theme.color as ColorTheme<any>).targetColor;
    if (!targetColor) return false;
    for (const t of targets.keys()) if (!known?.has(t) && targetColor(t) === undefined) return false;
    return true;
  }

  function reset() {
    store = new ParticlePageStore();
    pages.clear();
    coarse.clear();
    templateLists.clear();
    object = undefined;
    pool = [];
    shared = undefined;
    poolTargets = [];
    changedMembers = new Set();
    countedMembers = new Set();
    pending = emptyChanges();
    reaches = new Float32Array(0);
    poolCounts = new Int32Array(0);
    poolBounds = new Float64Array(0);
  }

  function placementOf(targets: ReadonlyMap<number, ParticleTarget>, props: ParticleTargetRepresentationProps) {
    const targetProps = createParticleTargetPropsProvider(props);
    const accepted = acceptedByTarget;
    return {
      centerOf: (t: number) => ParticleTarget.getCenter(targets.get(t)!) ?? Origin,
      accepts: (t: number) => {
        let a = accepted.get(t);
        if (a === undefined) {
          const target = targets.get(t);
          a = !!target && getTargetGeometryChoice(target, targetProps(target), graphics) === SlotsGeometry;
          if (target) accepted.set(t, a);
        }
        return a;
      },
    };
  }

  async function getMember(
    ctx: VisualContext,
    target: number,
    targets: ReadonlyMap<number, ParticleTarget>,
    theme: Theme,
    props: ParticleTargetRepresentationProps,
  ) {
    const kept = byTarget.get(target);
    if (kept) return kept;
    const t = targets.get(target)!;
    const member = await createTargetMember(
      ctx,
      t,
      theme,
      createParticleTargetPropsProvider(props)(t),
      graphics,
      materialId,
    );
    byTarget.set(target, member);
    return member;
  }

  /**
   * member `m`'s count, in the source's copy alone (`counts.instanceCounts`): the merged object reads a paged source's
   * counts from it, not its members' cells (thousands of them a step, a cell's box each)
   */
  function setMemberCount(m: number, count: number) {
    if (poolCounts[m] === count) return;
    poolCounts[m] = count;
    countedMembers.add(m);
  }

  /** the pool for the store's members and room, the last of no geometry; the object made again with it */
  async function createPool(
    ctx: VisualContext,
    targets: ReadonlyMap<number, ParticleTarget>,
    theme: Theme,
    props: ParticleTargetRepresentationProps,
    targetCount: number,
  ) {
    let like: TargetMember | undefined;
    for (let m = 0; m < store.pool.count && !like; ++m) {
      const t = store.pool.targets[m];
      if (t >= 0) like = await getMember(ctx, t, targets, theme, props);
    }
    object = undefined;
    pool = [];
    shared = undefined;
    poolTargets = [];
    changedMembers = new Set();
    countedMembers = new Set();
    if (!like) return;
    // room for every target of the list's catalog (a slot list's `targetInfo`: its tile set's, all of them): the pool
    // and the object made once, not again as tiles bring targets (the object's layers all uploaded again each time)
    const capacity = Math.max(MinPoolMembers, Math.ceil(store.pool.count * PoolHeadroom), targetCount + 1);
    const { type, state } = like.renderObject;
    poolType = type as MergedType;
    poolState = state;
    // the shared values: a target member's without geometry, in cells of their own (the object's), the clip as drawn
    // now (its member's from when it was made, a cut moved since)
    shared = createEmptyMergeMember(poolType, like.renderObject.values as unknown as MergeableValues);
    BaseGeometry.updateClip(shared as never, props.clip, true);
    // a color a member: the object's colors a table by segment, its members' all of one type from the start
    if (byMember) setMemberColor(shared, DefaultMemberColor);
    // an empty record: the shared values' cells (no geometry nor instances); a member's own cells those a fill writes
    const empty = pickMergeRecord(poolType, shared);
    fillKeys = getSourceRecordFillKeys(poolType, shared, new Set(['tColor', 'tMarker']));
    for (let m = 0; m <= capacity; ++m) {
      const record: { [k: string]: ValueCell<any> | undefined } = { ...empty };
      for (const k of fillKeys) record[k] = ValueCell.create(empty[k]!.ref.value);
      ValueCell.update(record.boundingSphere!, NoBounds);
      pool.push(record as MergeRecord);
      poolTargets.push(-1);
    }
    reaches = new Float32Array(pool.length);
    poolCounts = new Int32Array(pool.length);
    poolBounds = new Float64Array(pool.length * 4).fill(NaN);
    poolTargets.length = capacity;
    store.pool.capacity = capacity;
  }

  /**
   * the object's sphere (the last member's): where the list's slots' particles all are (a tile set's,
   * `ParticlePages.sphere`: the same from the first step, the camera reset to the whole), else the list's particles'
   * (`Particle.getBoundary`: a slot list's of its pages', each cached on its page's list), padded by the largest reach. Not the store pages' boxes' union box's sphere:
   * a round cell's tiles' box has corners √3 as far as its surface, a page's box a whole tile where its particles fill
   * a part
   */
  function setObjectBounds(list: ParticleList) {
    const sphere = list.pages?.sphere ?? Particle.getBoundary(list).sphere;
    if (!(sphere.radius >= 0) || ParticleSlots.getLiveCount(list) === 0) return;
    let reach = 0;
    for (let m = 0; m < reaches.length; ++m) reach = Math.max(reach, reaches[m]);
    const last = pool.length - 1;
    if (last < 0) return;
    Sphere3D.copy(objectSphere, sphere);
    // a coarse sample's instances drawn larger (`ParticleList.scale`); not by a template's coarse levels' (`coarse`:
    // drawn by far copies alone, up to 63 times larger: a sphere half again too large), a far copy's at its edge past
    // it by a little
    let scale = 0;
    for (const block of store.blockList) if (!coarse.has(block.key)) scale = Math.max(scale, block.maxScale);
    objectSphere.radius += reach * (scale || store.maxScale);
    const o = last * 4;
    const b = poolBounds;
    const c = objectSphere.center;
    if (b[o] === c[0] && b[o + 1] === c[1] && b[o + 2] === c[2] && b[o + 3] === objectSphere.radius) return;
    b[o] = c[0];
    b[o + 1] = c[1];
    b[o + 2] = c[2];
    b[o + 3] = objectSphere.radius;
    ValueCell.update(pool[last].boundingSphere, Sphere3D.clone(objectSphere));
    countedMembers.add(last);
  }

  /** the layers for the store's pages (made again where they grew: all written) */
  function ensureLayers() {
    const capacity = store.pageCapacity * StorePageSize;
    let grown = false;
    if (!byMember && colors.ref.value.count < capacity) {
      ValueCell.update(colors, createDataArray(capacity, 4, Uint8Array));
      grown = true;
    }
    if (markers.ref.value.count < capacity) {
      ValueCell.update(markers, createDataArray(capacity, 1, Uint8Array));
      grown = true;
    }
    return grown;
  }

  /** `block`'s store pages' texel ranges (its markers', colors by member) */
  function pushBlockRanges(block: StoreBlock, ranges: number[]) {
    for (let k = 0; k < block.pages.length; ++k) {
      const lo = k << StorePageBits;
      const hi = Math.min(block.count, lo + StorePageSize);
      const at = block.pages[k] << StorePageBits;
      ranges.push(at, at + hi - lo);
    }
  }

  let colorScratch = new Uint8Array(0);
  let templateScratch = new Uint8Array(0);
  /**
   * the template whose colors `templateScratch` holds in this pass of `fillColors` (cleared before each: a theme's
   * colors another's), filled once for its block and its coarse levels'. Why: filled a block each, the yeast cell's
   * 31.5 M particles four times (itself and its three levels: 0.9 s of its 6 s load)
   */
  let templateFilled: ParticleList | undefined;
  function getTemplateColors(template: ParticleList, theme: Theme) {
    if (templateFilled !== template) {
      if (templateScratch.length < template.count * 4) templateScratch = new Uint8Array(template.count * 4);
      fillParticleColors(theme.color as ColorTheme<any>, template, templateScratch, 0);
      templateFilled = template;
    }
    return templateScratch;
  }
  /**
   * `block`'s colors into its store pages' texels: the ranges written. A template's coarse level's its template's
   * particles' (`coarse`: the theme's of the template, not of a sample). `templateFilled` cleared before a pass
   */
  function fillColors(block: StoreBlock, theme: Theme, ranges: number[]) {
    const page = pages.get(block.key)!;
    const n = block.count;
    const level = coarse.get(block.key);
    let source: Uint8Array;
    if (level) {
      const from = new Uint32Array(getTemplateColors(pages.get(level.of)!.list, theme).buffer);
      if (colorScratch.length < n * 4) colorScratch = new Uint8Array(n * 4);
      // a color a word
      const to = new Uint32Array(colorScratch.buffer);
      const { indices } = level;
      for (let i = 0; i < n; ++i) to[i] = from[indices[i]];
      source = colorScratch;
    } else if (block.template) {
      source = getTemplateColors(page.list, theme);
    } else {
      if (colorScratch.length < n * 4) colorScratch = new Uint8Array(n * 4);
      fillParticleColors(theme.color as ColorTheme<any>, page.list, colorScratch, 0);
      source = colorScratch;
    }
    const array = colors.ref.value.array;
    for (let k = 0; k < block.pages.length; ++k) {
      const lo = k << StorePageBits;
      const hi = Math.min(n, lo + StorePageSize);
      const at = block.pages[k] << StorePageBits;
      array.set(source.subarray(lo * 4, hi * 4), at * 4);
      ranges.push(at, at + hi - lo);
    }
  }

  /** the tiles the object's occluder grid is of (`setOccluderGrid`), and the grid of them (none: none) */
  let occluderTiles: ArrayLike<number> | undefined;
  let occluder:
    | { readonly set: ParticleTileSet; readonly grid: ParticleDensityGrid; readonly tiles: ArrayLike<number> }
    | undefined;

  async function update(
    ctx: VisualContext,
    next: ParticleList,
    theme: Theme,
    props: ParticleTargetRepresentationProps,
  ) {
    // a list with templates: its targets by the flattened list's ids (its templates' moved past the top list's), what
    // its themes and members are keyed by
    const nested = ParticleNesting.hasTemplates(next);
    const targets = nested
      ? (ParticleNesting.getIdSpace(next).targetMapping as ReadonlyMap<number, ParticleTarget>)
      : (next.targetMapping ?? new Map<number, ParticleTarget>());
    const propsChanged = !currentProps || !PD.areEqual(ParticleTargetRepresentationParams, currentProps, props);
    const sizeChanged = !currentTheme || !SizeTheme.areEqual(theme.size, currentTheme.size);
    const colorChanged = !currentTheme || !ColorTheme.areEqual(theme.color, currentTheme.color);
    const newMembers = propsChanged || sizeChanged || !currentTargets || !isSameTargets(currentTargets, targets);
    if (newMembers) {
      byTarget = new Map();
      acceptedByTarget = new Map();
    }
    // the same colors and by member before: the targets new since alone asked (a step's view: thousands of targets)
    const known = !colorChanged && byMember && !newMembers ? currentTargets : undefined;
    // a slot list of no pages (an empty view) keeps the mode: its theme colors no target, made again for the next
    const memberColors = next.pages?.pages.length === 0 ? byMember : isByMember(theme, targets, known);
    const modeChanged = memberColors !== byMember;
    byMember = memberColors;
    // the same pages again (a theme, props or clipping changed): no page changed, not all of them written again
    const changes =
      !newMembers && !modeChanged && list
        ? list.pages && list.pages === next.pages
          ? { from: list.pages, added: [], removed: [], posed: [] }
          : ParticleSlots.getChanges(list, next)
        : undefined;
    // all made again (its pages, members and object): counted in timing mode, a cost of every page's
    if (!changes) UserTiming.span('store: reset', reset);

    const placement = placementOf(targets, props);
    const written: StoreBlock[] = [];
    const add = (page: ParticlePage) => {
      pages.set(page.start, page);
      // its box where known (a tile's list: its bounds from the tiling's pass), not a pass over its words
      const known = page.list._propertyData[Particle.BoundaryDescriptor.name] as Boundary | undefined;
      written.push(store.add(page.start, page.list, page.start, placement, known?.box));
    };
    UserTiming.span('store: pages', () => {
      if (changes) {
        for (const p of changes.removed) {
          store.remove(p.start);
          pages.delete(p.start);
        }
        for (const p of changes.posed) {
          store.remove(p.start);
          add(p);
        }
        for (const p of changes.added) add(p);
      } else if (nested) {
        addNested(next, placement, written);
      } else {
        ParticleSlots.forEachPage(next, add);
      }
    });
    const storeChanges = store.takeChanges();
    // appended one by one: a spread of a store's every page (a reset of a large view) is past the arguments' limit
    arrayExtend(pending.written, storeChanges.written);
    arrayExtend(pending.released, storeChanges.released);
    arrayExtend(pending.runs, storeChanges.runs);
    arrayExtend(pending.extras, storeChanges.extras);
    arrayExtend(pending.copies, storeChanges.copies);
    pending.grown = pending.grown || storeChanges.grown;

    // the pool: made where the store's members are more (the object made again); a member filled where the store
    // gave it another target (its geometry allocated again, `Merged.replaceGeometry`)
    const poolMade = store.pool.count > pool.length - 1;
    const targetCount = nested ? ParticleNesting.getIdSpace(next).targetInfo.size : (next.targetInfo?.size ?? 0);
    if (poolMade) await createPool(ctx, targets, theme, props, targetCount);
    for (let m = 0; m < store.pool.count; ++m) {
      const t = store.pool.targets[m];
      if (t < 0 || poolTargets[m] === t) continue;
      const member = await getMember(ctx, t, targets, theme, props);
      const values = pool[m];
      UserTiming.span('store: fill member', () =>
        fillMergeRecord(values, member.renderObject.values as unknown as MergeableValues, fillKeys),
      );
      if (byMember) setMemberColor(values, (theme.color as ColorTheme<any>).targetColor!(t)!);
      // a member's bounds none of its own: the object's the last member's
      ValueCell.update(values.boundingSphere, NoBounds);
      poolTargets[m] = t;
      UserTiming.span('store: replace geometry', () => object?.merged.replaceGeometry(m));
      changedMembers.add(m);
      storeChanges.members.add(m);
      // a pose places the geometry's point x at p + R (x - c): within |C - c| + r of its particle
      const ibs = values.invariantBoundingSphere.ref.value;
      reaches[m] = Vec3.distance(ibs.center, ParticleTarget.getCenter(targets.get(t)!) ?? Origin) + ibs.radius;
    }
    UserTiming.span('store: members', () => {
      for (const m of storeChanges.members) {
        pending.members.add(m);
        if (m >= pool.length - 1) continue;
        setMemberCount(m, m < store.pool.count ? store.pool.counts[m] : 0);
      }
      setObjectBounds(next);
    });

    // colors by member: the theme's changed, every member's again (its segment's texel)
    if (byMember && colorChanged && changes) {
      const targetColor = (theme.color as ColorTheme<any>).targetColor!;
      for (let m = 0; m < store.pool.count; ++m) {
        const t = poolTargets[m];
        if (t < 0) continue;
        setMemberColor(pool[m], targetColor(t)!);
        changedMembers.add(m);
      }
    }

    // the layers: all where they are new or the theme changed, else the pages written; colors by member none. The
    // markers all cleared where the store's texels are laid out again, else those of the pages written (a store page's
    // texels another tile's before, marked or not): a theme changed keeps what is marked
    const layersMade = ensureLayers();
    const all = layersMade || !changes || colorChanged || poolMade || storeChanges.grown;
    const clearMarkers = layersMade || !changes || poolMade || storeChanges.grown;
    const ranges: number[] = [];
    UserTiming.span('store: colors', () => {
      templateFilled = undefined;
      if (all) {
        if (!byMember) {
          for (const block of store.blockList) fillColors(block, theme, ranges);
          TextureChanges.update(colors, undefined);
        }
      } else if (written.length > 0) {
        for (const block of written) {
          if (byMember) pushBlockRanges(block, ranges);
          else fillColors(block, theme, ranges);
        }
        if (!byMember) TextureChanges.update(colors, ranges);
      }
    });
    // the copies' markers cleared with the store's, and where its templates or their copies changed
    if (clearMarkers || storeChanges.copies.length > 0 || written.some((b) => b.template)) layoutCopyMarkers();
    if (clearMarkers) {
      markers.ref.value.array.fill(0);
      setStoreMarkers(0);
      TextureChanges.update(markers, undefined);
    } else if (written.length > 0 && storeCells.markerStatus.ref.value !== 0) {
      const writtenRanges: number[] = [];
      for (const block of written) pushBlockRanges(block, writtenRanges);
      const array = markers.ref.value.array;
      for (let k = 0; k < writtenRanges.length; k += 2) array.fill(0, writtenRanges[k], writtenRanges[k + 1]);
      setStoreMarkers(getMarkersAverage(array, store.pageCapacity * StorePageSize));
      TextureChanges.update(markers, writtenRanges);
    }
    // the copies' markers laid out again (cleared) or the store's changed
    setObjectMarkers(false);

    if (!object && pool.length > 0 && shared && poolState) {
      object = createMergedRecordsRenderObject(poolType, shared as never, pool, poolState, materialId, {
        gpuCulling: true,
        source: createSource(),
      });
      occluderTiles = undefined;
    }
    // what its shown tiles hide (their opacities): the depth pyramid primed from it (`HiZPass.setOccluder`); a view's
    // step the last grid updated in place by the tiles it added and let go, not made again
    const view = isTargetOccluderOff(props) ? undefined : ParticleTileView.of(next);
    if (object && view?.tiles !== occluderTiles) {
      occluderTiles = view?.tiles;
      const last = occluder?.set === view?.set ? occluder : undefined;
      const grid = view
        ? UserTiming.span('store: occluder', () => ParticleTileOpaque.updateDensityGrid(view.set, view.tiles, last))
        : undefined;
      occluder = view && grid ? { set: view.set, grid, tiles: view.tiles } : undefined;
      setOccluderGrid(object.state, grid);
    }

    if (!changes || !list) slotVersions = new WeakSet();
    else if (next.pages) slotVersions.add(list);
    list = next;
    currentTheme = theme;
    currentProps = { ...props };
    currentTargets = targets;
    // the layers on the store's texels as now laid out: all where laid out again or its templates or their copies
    // changed (the copies' own too), else the pages written alone (the others' texels at the same store indices)
    if (StoreLayers.some((l) => (layerState[l.key]?.layers.length ?? 0) > 0)) {
      const relaid = all || storeChanges.copies.length > 0 || written.some((b) => b.template);
      if (relaid) applyLayers(layerState);
      else if (written.length > 0) {
        const writtenRanges: number[] = [];
        for (const block of written) pushBlockRanges(block, writtenRanges);
        applyLayers(layerState, writtenRanges);
      }
    }
  }

  /**
   * A list with templates in the store: its leaves a block (its key and slots 0, its particles' as a flat list's; its
   * copies' particles in no run), each template a block of its leaf list (`ParticleNesting.getLeafList`: its targets in
   * the flattened list's ids, a key below 0 a template, its pages drawn by its copies) and its copies (their top
   * particles' poses about the template's center: `ParticleTarget.getCenter`, the template's boundary's by default)
   */
  function addNested(next: ParticleList, placement: StorePlacement, written: StoreBlock[]) {
    // its leaves by the flattened list's tables (its themes'), its copies' particles of no target drawn
    const top: ParticlePage = { start: 0, list: ParticleNesting.getTopList(next) };
    // the store's format of all its lists (the top list's first: without rotations where its templates have them)
    const templates = [...ParticleNesting.getIdSpace(next).moves.keys()];
    store.format ??= getStoreFormatOf([top.list, ...templates]);
    pages.set(0, top);
    written.push(store.add(0, top.list, 0, placement));
    const p = Vec3();
    const q = Quat();
    let k = 0;
    for (const template of ParticleNesting.getIdSpace(next).moves.keys()) {
      const key = -++k;
      const leaves = ParticleNesting.getLeafList(next, template);
      pages.set(key, { start: 0, list: leaves });
      templateLists.set(key, template);
      written.push(store.add(key, leaves, 0, placement, undefined, true));
      const tops = Int32Array.from(ParticleNesting.getCopiesOf(next, template));
      const poses = new Float32Array(tops.length * CopyPoseFloats);
      let center: Vec3 | undefined;
      for (let c = 0; c < tops.length; ++c) {
        const t = tops[c];
        center = ParticleTarget.getCenter(ParticleNesting.getTemplate(next, ParticleValues.get(next.targets, t))!)!;
        GridPositions.get(next.coordinates, t, p);
        ParticleRotations.get(next.rotations, t, q);
        writeCopyPose(poses, c * CopyPoseFloats, p, q, center);
      }
      // its coarse levels each a template block of the same copies, a copy drawn at one by its distance: what picks
      // it its center and its particles' spacing (the cube root of its bounding sphere's volume a particle: a flat
      // template's too, its box of no volume)
      const levels = ParticleNesting.getTemplateLevels(next, template);
      // the template's own (its leaves its positions: cached on it, as its center's), not the leaf list's (a pass)
      const bounds = Particle.getBoundary(template).sphere;
      const r = bounds.radius;
      const spacing = Math.cbrt(((4 / 3) * Math.PI * r * r * r) / Math.max(1, leaves.count));
      // its sphere about its center (where a copy places it): a copy out of view culled whole
      const radius = r + (center ? Vec3.distance(bounds.center, center) : 0);
      const of = { tops, poses, levels: levels.length, center, spacing, radius };
      store.setCopies(key, { ...of, level: 0 });
      for (const l of levels) {
        const levelKey = -++k;
        pages.set(levelKey, { start: 0, list: l.list });
        coarse.set(levelKey, { of: key, indices: l.indices });
        written.push(store.add(levelKey, l.list, 0, placement, undefined, true));
        store.setCopies(levelKey, { ...of, level: l.level });
      }
    }
  }

  /** a store index's slot, -1 where its page holds none */
  function slotOf(index: number) {
    const p = index >> StorePageBits;
    const j = index & (StorePageSize - 1);
    if (p >= store.pageCount) return -1;
    const o = p * PageRecordWords;
    if (j >= store.records[o + 6]) return -1;
    return store.records[o + 7] + j;
  }

  /** a slot's store index, -1 where none (a list with templates: a top particle's, its leaves' block's) */
  function indexOf(slot: number) {
    if (!list!.pages) {
      const block = store.getBlock(0);
      if (!block || slot < 0 || slot >= block.count) return -1;
      return (block.pages[slot >> StorePageBits] << StorePageBits) + (slot & (StorePageSize - 1));
    }
    const page = ParticleSlots.findPage(list!, slot);
    if (page < 0) return -1;
    const block = store.getBlock(list!.pages!.pages[page].start);
    if (!block) return -1;
    const j = slot - block.start;
    if (j < 0 || j >= block.count) return -1;
    return (block.pages[j >> StorePageBits] << StorePageBits) + (j & (StorePageSize - 1));
  }

  /**
   * The path of template page `page`'s particle `j` in copy `copy` (the cull's: the templates' copies in the store's
   * order, `ParticlePageStore.templates`), undefined where none
   */
  function copyPathOf(page: number, j: number, copy: number) {
    if (copy < 0) return;
    let base = 0;
    for (const t of store.templates) {
      const n = t.copies.tops.length;
      if (copy >= base + n) {
        base += n;
        continue;
      }
      const k = t.pages.indexOf(page);
      if (k < 0) return;
      // a coarse level's particle its template's (`coarse`)
      const local = k * StorePageSize + j;
      const level = coarse.get(t.key);
      return { top: t.copies.tops[copy - base], local: level ? level.indices[local] : local };
    }
  }

  function getLoci(pickingId: PickingId): ModelLoci | undefined {
    if (!object || !list || pickingId.objectId !== object.id) return undefined;
    // a template's particle: its copy the pick's (gathered from the composed pose it drew: `paged-pick.ts`), its path
    // (the copy's top particle, its index in the template)
    const page = pickingId.instanceId >> StorePageBits;
    if (page < store.pageCount && store.records[page * PageRecordWords + PageFlagsWord] & TemplatePage) {
      const copy = pickingId.copyId ?? -1;
      const path = copyPathOf(page, pickingId.instanceId & (StorePageSize - 1), copy);
      return path ? ParticleNesting.ofInstance(list, path.top, path.local) : EmptyLoci;
    }
    const slot = slotOf(pickingId.instanceId);
    if (slot < 0) return EmptyLoci;
    return Particle.Loci(list, OrderedSet.ofSingleton(slot));
  }

  /** `locals` (template particles) as the store indices of a block of their template (a coarse level's: its sample's) */
  function storeIndicesOf(key: number, pages: Int32Array, locals: OrderedSet<number>) {
    const inLevel = getLevelLocals(key);
    const indices: number[] = [];
    OrderedSet.forEach(locals, (local) => {
      const j = inLevel ? inLevel.get(local) : local;
      if (j === undefined) return;
      const page = pages[j >> StorePageBits];
      if (page !== undefined) indices.push((page << StorePageBits) + (j & (StorePageSize - 1)));
    });
    indices.sort((a, b) => a - b);
    return indices;
  }

  /**
   * template `of`'s marks (the store's markers of its particles, every copy's) made each copy's and cleared: before
   * an action that takes away from a part of them (a copy's, one instance), not all of the template's
   */
  function spreadTemplateMarks(of: number) {
    const array = markers.ref.value.array;
    let marked = false;
    for (const { pages } of templateBlocks.get(of) ?? []) {
      for (let i = 0; i < pages.length && !marked; ++i) {
        const lo = pages[i] << StorePageBits;
        for (let k = lo; k < lo + StorePageSize && !marked; ++k) marked = array[k] !== 0;
      }
    }
    if (!marked) return false;
    for (const { copy } of templateCopies.get(of) ?? []) copyMarkers.include(copy, array);
    for (const { pages } of templateBlocks.get(of) ?? []) {
      for (let i = 0; i < pages.length; ++i) array.fill(0, pages[i] << StorePageBits, (pages[i] + 1) << StorePageBits);
    }
    return true;
  }

  /**
   * `loci`'s copies marked (a list with templates'): its top particles placing a copy that copy whole, its paths a
   * copy's particles by their template's store indices; the same particles of every copy of a template (a selection
   * by entity: `ParticleNesting.ofTemplateParticles`) the template's (the store's markers of its particles, read by
   * every copy: no block a copy); its other particles left to the store's markers. An action taking away from a part of
   * a template's marks makes them each copy's first. Whether the copies' markers changed, and the store's
   */
  function markCopies(loci: ModelLoci, action: MarkerAction) {
    const none = { copies: false, store: false };
    if (copyMarkers.copyCount === 0) return none;
    if (isEveryLoci(loci)) return { copies: copyMarkers.markAll(action), store: false };
    if (!Particle.isLoci(loci) || loci.particles !== list) return none;
    const adds = action === MarkerAction.Highlight || action === MarkerAction.Select;
    let copies = false;
    let storeChanged = false;
    const spread = (of: number) => {
      if (adds || !spreadTemplateMarks(of)) return;
      copies = true;
      storeChanged = true;
    };
    // a copy at every level it may be drawn at (its template's block and its coarse levels')
    OrderedSet.forEach(loci.indices, (top) => {
      const entries = copyOfTop.get(top);
      if (!entries) return;
      spread(templateOf(entries[0].key));
      for (const { copy } of entries) if (copyMarkers.markCopy(copy, action)) copies = true;
    });
    const partial = loci.partial;
    if (!partial) return { copies, store: storeChanged };
    // by template: the same particles of each of its copies (an entity's) its own
    const byTemplate = new Map<number, number[]>();
    for (let k = 0; k < partial.tops.length; ++k) {
      const entries = copyOfTop.get(partial.tops[k]);
      if (!entries) continue;
      const of = templateOf(entries[0].key);
      let ks = byTemplate.get(of);
      if (!ks) byTemplate.set(of, (ks = []));
      ks.push(k);
    }
    const array = markers.ref.value.array;
    for (const [of, ks] of byTemplate) {
      const blocks = templateBlocks.get(of) ?? [];
      const locals = partial.locals[ks[0]];
      const whole =
        ks.length === blocks[0]?.tops &&
        ks.every((k) => partial.locals[k] === locals || OrderedSet.areEqual(partial.locals[k], locals));
      if (whole && action !== MarkerAction.Toggle) {
        for (const { key, pages } of blocks) {
          const indices = storeIndicesOf(key, pages, locals);
          if (indices.length > 0) applyMarkerAction(array, SortedArray.ofSortedArray(indices), action);
        }
        storeChanged = true;
        // taken away: from a copy's own marks of them too
        if (!adds) {
          for (const k of ks) {
            for (const { copy, key, pages } of copyOfTop.get(partial.tops[k])!) {
              if (!copyMarkers.isMarked(copy)) continue;
              const indices = storeIndicesOf(key, pages, partial.locals[k]);
              if (copyMarkers.markIndices(copy, SortedArray.ofSortedArray(indices), action)) copies = true;
            }
          }
        }
        continue;
      }
      spread(of);
      for (const k of ks) {
        for (const { copy, key, pages } of copyOfTop.get(partial.tops[k])!) {
          // a level's particles those of its sample (by their indices in it), the others none
          const indices = storeIndicesOf(key, pages, partial.locals[k]);
          if (indices.length === 0) continue;
          if (copyMarkers.markIndices(copy, SortedArray.ofSortedArray(indices), action)) copies = true;
        }
      }
    }
    return { copies, store: storeChanged };
  }

  /** a coarse level's particles by their template particles' indices (none for a template's own block) */
  function getLevelLocals(key: number) {
    const level = coarse.get(key);
    if (!level) return undefined;
    let m = levelLocals.get(key);
    if (!m) {
      m = new Map();
      for (let i = 0; i < level.indices.length; ++i) m.set(level.indices[i], i);
      levelLocals.set(key, m);
    }
    return m;
  }

  function mark(loci: ModelLoci, action: MarkerAction) {
    if (!object || !list) return false;
    const capacity = store.pageCapacity * StorePageSize;
    ValueCell.updateIfChanged(markValues.instanceCount, capacity);
    const byCopies = markCopies(loci, action);
    const copiesChanged = byCopies.copies;
    // a template's marks in the store's markers: uploaded whole (a template's pages, its levels' too), its aggregates
    // counted (`Visual.mark` not told)
    if (byCopies.store) {
      TextureChanges.update(markers, undefined);
      setStoreMarkers(getMarkersAverage(markers.ref.value.array, capacity));
    }
    const storeChanged = markStore(loci, action, capacity) || byCopies.store;
    if (!copiesChanged && !storeChanged) return false;
    if (copiesChanged) copyMarkersVersion += 1;
    setObjectMarkers(copiesChanged);
    return true;
  }

  /**
   * `loci`'s store indices as runs (`[start, end)` to `f`): its slots' (a flat list's, a list with templates' own
   * leaves), and of a list with templates the particles of a template it holds in every copy
   * (`getTemplateWideLocals`), at each of its levels; every loci all of them
   */
  function forEachStoreRun(loci: ModelLoci, capacity: number, f: (start: number, end: number) => void) {
    if (isEveryLoci(loci)) {
      f(0, capacity);
      return;
    }
    if (!list || !Particle.isLoci(loci)) return;
    if (loci.particles !== list && !slotVersions.has(loci.particles)) return;
    const indices: number[] = [];
    OrderedSet.forEach(loci.indices, (slot) => {
      const i = indexOf(slot);
      if (i >= 0) indices.push(i);
    });
    for (const t of store.templates) {
      const template = templateLists.get(t.key);
      if (!template) continue;
      const locals = getTemplateWideLocals(list, template, t.copies.tops, loci);
      if (!locals) continue;
      for (const { key, pages } of templateBlocks.get(t.key) ?? []) {
        for (const i of storeIndicesOf(key, pages, locals)) indices.push(i);
      }
    }
    if (indices.length === 0) return;
    indices.sort((a, b) => a - b);
    let start = indices[0];
    let end = start + 1;
    for (let k = 1; k < indices.length; ++k) {
      const i = indices[k];
      if (i === end) end += 1;
      else if (i > end) {
        f(start, end);
        start = i;
        end = i + 1;
      }
    }
    f(start, end);
  }

  function setLayers(layers: TargetLayers) {
    for (const l of StoreLayers) {
      const value = layers[l.key];
      if (value !== undefined) (layerState as Record<string, unknown>)[l.key] = value;
    }
    applyLayers(layers);
  }

  /**
   * `layers` (each given its kind's again) on the store's texels and its copies': a kind with layers its define on
   * (by instance), its texels by store index, its copies' own (a list with templates'); the object's transparency and
   * emission averages by them (a copy's at least a half), its members' (what the merged object averages). `ranges`
   * (`[start, end)` store index pairs): those texels alone written again and uploaded, the copies' own kept (an
   * update's pages written, its templates and their copies the same)
   */
  function applyLayers(layers: TargetLayers, ranges?: readonly number[]) {
    if (!shared || !list) return;
    const capacity = store.pageCapacity * StorePageSize;
    const values = shared as unknown as Record<string, ValueCell<any>>;
    let copiesChanged = false;
    // the copies' blocks within the budget (`CopyLayerRoom`), the kinds kept (not given) holding theirs
    const budget = getCopyBlockBudget(graphics);
    const room: CopyLayerRoom = { bytes: budget, coarse: 0 };
    if (!ranges) {
      for (const [key, byTemplate] of copyLayerStates) {
        if (layers[key]) continue;
        const size = getCopyLayerSize(key);
        for (const [of, byCopy] of byTemplate) {
          const span = getCopyLayerSpan(of);
          for (const st of byCopy.values()) if (st.block) room.bytes -= span * size;
        }
      }
    }
    // sorted and merged (a gap's texels written again as well: all by the layers), undefined: most of them, all
    const sorted = ranges && TextureChanges.mergeRanges(ranges, capacity);
    for (const l of StoreLayers) {
      const value = layers[l.key] as AnyLayers | undefined;
      if (!value) continue;
      const on = value.layers.length > 0;
      const size = getCopyLayerSize(l.key);
      const cell = layerCells.get(l.tex)!;
      if (on && sorted && cell.ref.value.count >= capacity) {
        const array = cell.ref.value.array;
        for (let k = 0; k < sorted.length; k += 2) array.fill(0, sorted[k] * size, sorted[k + 1] * size);
        for (const layer of value.layers) {
          forEachStoreRun(layer.loci, capacity, (start, end) =>
            forEachRangeIn(sorted, start, end, (s, e) => writeLayerItems(l.key, array, s, e, layer)),
          );
        }
        TextureChanges.update(cell, sorted);
      } else if (on) {
        if (cell.ref.value.count < capacity) ValueCell.update(cell, createDataArray(capacity, size, Uint8Array));
        const array = cell.ref.value.array;
        array.fill(0);
        for (const layer of value.layers) {
          forEachStoreRun(layer.loci, capacity, (start, end) => writeLayerItems(l.key, array, start, end, layer));
        }
        TextureChanges.update(cell, undefined);
      } else if (cell.ref.value.count > 1) {
        ValueCell.update(cell, createDataArray(1, size, Uint8Array));
      }
      if (values[l.type]) ValueCell.updateIfChanged(values[l.type], 'instance');
      ValueCell.updateIfChanged(values[l.define], on);
      // the copies' own (a list with templates'), not by store index: as they were where texels alone are written
      if (ranges) continue;
      const byTemplate = new Map<number, Map<number, CopyState>>();
      if (on) {
        for (const t of store.templates) {
          const template = templateLists.get(t.key);
          if (!template) continue;
          const blockBytes = getCopyLayerSpan(t.key) * size;
          const states = getCopyLayerStates(list, template, t.copies.tops, l.key, value, room, blockBytes);
          if (states.size > 0) byTemplate.set(t.key, states);
        }
      }
      if (byTemplate.size > 0 || (copyLayerStates.get(l.key)?.size ?? 0) > 0) copiesChanged = true;
      copyLayerStates.set(l.key, byTemplate);
    }
    if (!copyLayersWarned) copyLayersWarned = warnCopyLayerRoom(room, budget);
    if (copiesChanged) writeCopyLayers();
    // the averages: the store's texels', at least a half where a copy has its own
    const average = (key: CopyLayerKey, get: (array: Uint8Array, count: number) => number) => {
      const cell = layerCells.get(StoreLayers.find((l) => l.key === key)!.tex)!;
      const own = cell.ref.value.count > 1 ? get(cell.ref.value.array, capacity) : 0;
      return (copyLayerStates.get(key)?.size ?? 0) > 0 ? Math.max(own, 0.5) : own;
    };
    if (layers.transparency || layers.emissive || layers.wiggle) {
      const transparency = average('transparency', getTransparencyAverage);
      const emissive = average('emissive', getEmissiveAverage);
      const wiggle = average('wiggle', getWiggleAverage);
      for (let m = 0; m < pool.length; ++m) {
        const record = pool[m] as unknown as Record<string, ValueCell<number> | undefined>;
        if (record.transparencyAverage) ValueCell.updateIfChanged(record.transparencyAverage, transparency);
        if (record.emissiveAverage) ValueCell.updateIfChanged(record.emissiveAverage, emissive);
        if (record.wiggleAverage) ValueCell.updateIfChanged(record.wiggleAverage, wiggle);
        changedMembers.add(m);
      }
    }
  }

  /** a template's copy's block of layers in store indices (`writeCopyLayers`): its blocks' pages' spans, its levels' */
  function getCopyLayerSpan(of: number) {
    let span = 0;
    for (const { pages: p } of templateBlocks.get(of) ?? []) {
      let lo = Infinity;
      let hi = -Infinity;
      for (let i = 0; i < p.length; ++i) {
        lo = Math.min(lo, p[i]);
        hi = Math.max(hi, p[i] + 1);
      }
      if (p.length > 0) span += (hi - lo) << StorePageBits;
    }
    return span;
  }

  /** the parts of `[start, end)` within sorted disjoint `ranges` (pairs) to `f` */
  function forEachRangeIn(ranges: readonly number[], start: number, end: number, f: (s: number, e: number) => void) {
    // the first range ending past start
    let lo = 0;
    let hi = ranges.length / 2;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (ranges[2 * mid + 1] <= start) lo = mid + 1;
      else hi = mid;
    }
    for (let k = 2 * lo; k < ranges.length && ranges[k] < end; k += 2) {
      f(Math.max(start, ranges[k]), Math.min(end, ranges[k + 1]));
    }
  }

  /**
   * The copies' own layers as the GPU reads them (`sCopyLayers`): per copy (the copy markers', a template's at each
   * level) `CopyLayerKinds` pairs of words, none, a value (`CopyLayerValue`, its first byte) or a block
   * (`CopyLayerBlock`, its first byte less its first store index times the size: by store index over its pages' span),
   * then their bytes
   */
  function writeCopyLayers() {
    let copyCount = 0;
    for (const t of store.templates) copyCount += t.copies.tops.length;
    const any = [...copyLayerStates.values()].some((m) => m.size > 0);
    if (!any) {
      copyLayerBuffer = new Uint32Array(0);
      copyLayerVersion += 1;
      return;
    }
    const wordCount = copyCount * CopyLayerKinds * 2;
    const words: number[] = new Array(wordCount).fill(0);
    const chunks: Uint8Array[] = [];
    let at = wordCount * 4;
    let c = 0;
    for (const t of store.templates) {
      const of = templateOf(t.key);
      const level = coarse.get(t.key);
      const count = pages.get(t.key)?.list.count ?? 0;
      // the block's span of store indices
      let lo = Infinity;
      let hi = -Infinity;
      for (let i = 0; i < t.pages.length; ++i) {
        lo = Math.min(lo, t.pages[i] << StorePageBits);
        hi = Math.max(hi, (t.pages[i] + 1) << StorePageBits);
      }
      for (let k = 0; k < t.copies.tops.length; ++k, ++c) {
        StoreLayers.forEach((l, kind) => {
          const state = copyLayerStates.get(l.key)?.get(of)?.get(k);
          if (!state) return;
          const size = getCopyLayerSize(l.key);
          const o = 2 * (c * CopyLayerKinds + kind);
          let bytes: Uint8Array;
          if (state.value) {
            bytes = state.value;
            words[o] = CopyLayerValue;
            words[o + 1] = at;
          } else {
            bytes = new Uint8Array((hi - lo) * size);
            const block = state.block!;
            for (let i = 0; i < count; ++i) {
              const local = level ? level.indices[i] : i;
              const index = (t.pages[i >> StorePageBits] << StorePageBits) + (i & (StorePageSize - 1));
              bytes.set(block.subarray(local * size, local * size + size), (index - lo) * size);
            }
            words[o] = CopyLayerBlock;
            words[o + 1] = (at - lo * size) | 0;
          }
          chunks.push(bytes);
          at += (bytes.length + 3) & ~3;
        });
      }
    }
    const buffer = new Uint8Array(at);
    const head = new Uint32Array(buffer.buffer, 0, wordCount);
    for (let i = 0; i < wordCount; ++i) head[i] = words[i] >>> 0;
    let offset = wordCount * 4;
    for (const bytes of chunks) {
      buffer.set(bytes, offset);
      offset += (bytes.length + 3) & ~3;
    }
    copyLayerBuffer = new Uint32Array(buffer.buffer);
    copyLayerVersion += 1;
  }

  /** the store's markers' aggregates by their average (`storeCells`; a status of 0 or -1) */
  function setStoreMarkers(average: number) {
    const status = average === 0 ? 0 : -1;
    ValueCell.updateIfChanged(storeCells.markerAverage, average);
    ValueCell.updateIfChanged(storeCells.markerStatus, status);
    ValueCell.updateIfChanged(storeCells.uMarker, status);
  }

  /**
   * The object's marker aggregates (`markerCells`) by the store's and the copies' markers: the store's where no copy
   * is marked, else read per instance (`uMarker` -1), its average not known but above 0; `moved`: the uniform's version
   * moved all the same (the object updated, the culler's sync uploading the copies')
   */
  function setObjectMarkers(moved: boolean) {
    const storeAverage = storeCells.markerAverage.ref.value;
    const storeStatus = storeCells.markerStatus.ref.value;
    const anyCopy = copyMarkers.markedCount > 0;
    const average = anyCopy ? (storeAverage > 0 ? storeAverage : 0.5) : storeAverage;
    const status = anyCopy ? -1 : storeStatus;
    // the store's markers stood for by its status (their changes' upload deferred, `Visual.mark`): read per instance
    // now, what was deferred uploaded
    if (status === -1 && storeStatus !== -1) TextureChanges.update(markers, []);
    ValueCell.updateIfChanged(markerCells.markerAverage, average);
    ValueCell.updateIfChanged(markerCells.markerStatus, status);
    if (moved) ValueCell.update(markerCells.uMarker, status);
    else ValueCell.updateIfChanged(markerCells.uMarker, status);
  }

  function markStore(loci: ModelLoci, action: MarkerAction, capacity: number) {
    const lociApply: Visual.LociApply = (loci, apply) => {
      if (isEveryLoci(loci)) return apply(Interval.ofBounds(0, capacity));
      if (!Particle.isLoci(loci) || Particle.isLociEmpty(loci)) return false;
      const indices: number[] = [];
      OrderedSet.forEach(loci.indices, (slot) => {
        const i = indexOf(slot);
        if (i >= 0) indices.push(i);
      });
      if (indices.length === 0) return false;
      indices.sort((a, b) => a - b);
      let changed = false;
      let start = indices[0];
      let end = start + 1;
      for (let k = 1; k < indices.length; ++k) {
        const i = indices[k];
        if (i === end) {
          end += 1;
          continue;
        }
        if (apply(Interval.ofBounds(start, end))) changed = true;
        start = i;
        end = i + 1;
      }
      if (apply(Interval.ofBounds(start, end))) changed = true;
      return changed;
    };
    return Visual.mark({ values: markValues } as unknown as GraphicsRenderObject, loci, action, lociApply);
  }

  return {
    get renderObject() {
      return object;
    },
    get store() {
      return store;
    },
    update,
    applyState(apply) {
      if (object) apply(object);
    },
    updateClip(clip) {
      // the object's alone (what draws and what the cull reads; its clip cells the shared values'): the members' records
      // have none
      if (object) BaseGeometry.updateClip(object.values as never, clip);
      // the next full update's props compared with this clip: not changed props, no reset
      if (currentProps) currentProps = { ...currentProps, clip };
    },
    getLoci,
    mark,
    setLayers,
    destroy() {
      reset();
      byTarget = new Map();
      acceptedByTarget = new Map();
      list = undefined;
      currentTheme = undefined;
      currentProps = undefined;
      currentTargets = undefined;
    },
  };
}
