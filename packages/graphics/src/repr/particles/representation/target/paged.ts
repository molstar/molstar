/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { Interval, OrderedSet } from '@molstar/core/data/int';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { Boundary } from '@molstar/core/math/geometry/boundary';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ParticleList, Particle, ParticleTarget } from '@molstar/model/model/particles/particle-list';
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
  createMergedRenderObject,
  createRenderObject,
  type GraphicsRenderObject,
  type MergedGraphicsRenderObject,
} from '@molstar/graphics/gl/render-object';
import {
  createEmptyMergeMember,
  fillMergeMember,
  type MergeableValues,
  type MergedInstanceSource,
  type MergedType,
} from '@molstar/graphics/gl/renderable/merged';
import { createDataArray, type DataArray } from '@molstar/graphics/gl/renderable/util';
import { TextureChanges } from '@molstar/graphics/gl/renderable/texture-changes';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { getMarkersAverage } from '@molstar/graphics/geo/geometry/marker-data';
import { createObjectColor } from '@molstar/graphics/geo/geometry/color-data';
import { Color } from '@molstar/core/util/color';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { fillParticleColors } from '@molstar/graphics/theme/color/particle-pages';
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
  isTargetOccluderOff,
} from './common.js';
import { PageRecordWords, ParticlePageStore, StorePageBits, StorePageSize, type StoreBlock } from './page-store.js';

// A slot list's targets as one object of a page store (`ParticlePageStore`): a page a block of
// store pages, added and let go by the list's changes, its instances culled and composed on the GPU
// (`MergedInstanceSource.paged`). A target a member (a pool, `TargetMemberPool`), the last of no geometry
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
  destroy(): void;
}

/** The pool's members past the store's (the room it is made with), and its least */
const PoolHeadroom = 1.5;
const MinPoolMembers = 16;

/** a pool member's color before its target's (colors by member: drawn by none till then) */
const DefaultMemberColor = Color(0xcccccc);

/** The geometry a slot object's targets are drawn with */
const SlotsGeometry = 'structure-spheres';

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
    members: new Set<number>(),
    grown: false,
  };
}

export function createTargetPagedObject(materialId: number, graphics: GraphicsContext | undefined): TargetSlotsObject {
  let store = new ParticlePageStore();

  let list: ParticleList | undefined;
  let currentTheme: Theme | undefined;
  let currentProps: ParticleTargetRepresentationProps | undefined;
  let currentTargets: ReadonlyMap<number, ParticleTarget> | undefined;
  /** the pages held, by start (a block's key) */
  const pages = new Map<number, ParticlePage>();

  let byTarget = new Map<number, TargetMember>();
  let acceptedByTarget = new Map<number, boolean>();

  let object: MergedGraphicsRenderObject | undefined;
  /** the object's members (a pool, filled with a target's values), the last of no geometry */
  let pool: GraphicsRenderObject[] = [];
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
  const markerCells = {
    markerAverage: ValueCell.create(0),
    markerStatus: ValueCell.create(0),
    uMarker: ValueCell.create(0),
  };
  const markValues = {
    tMarker: markers,
    ...markerCells,
    uGroupCount: ValueCell.create(1),
    instanceCount: ValueCell.create(0),
    uInstanceGranularity: ValueCell.create(true),
  };

  function createSource(): MergedInstanceSource {
    return {
      layers: new Map<string, ValueCell<DataArray<any>>>([
        ['tColor', colors],
        ['tMarker', markers],
      ]),
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
      },
    };
  }

  /**
   * whether the colors are a member's (the theme's `targetColor` for every target: one color a target, read by the
   * member's segment, `object` colors), not a layer by store index filled a particle at a time; the other way the
   * object made again
   */
  let byMember = false;
  /** member `values`' colors one `color` (`createObjectColor`'s cells' values: the member's own, not its target's) */
  function setMemberColor(values: MergeableValues, color: Color) {
    const object = createObjectColor(color) as unknown as { [k: string]: ValueCell<any> };
    const cells = values as unknown as { [k: string]: ValueCell<any> };
    for (const k in object) if (cells[k]) ValueCell.update(cells[k], object[k].ref.value);
  }
  /** whether `theme` gives each of `targets` its color */
  function isByMember(theme: Theme, targets: ReadonlyMap<number, ParticleTarget>) {
    const targetColor = (theme.color as ColorTheme<any>).targetColor;
    if (!targetColor) return false;
    for (const t of targets.keys()) if (targetColor(t) === undefined) return false;
    return true;
  }

  function reset() {
    store = new ParticlePageStore();
    pages.clear();
    object = undefined;
    pool = [];
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

  function setMemberCount(m: number, count: number) {
    if (poolCounts[m] === count) return;
    poolCounts[m] = count;
    ValueCell.update((pool[m].values as unknown as MergeableValues).instanceCount, count);
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
    poolTargets = [];
    changedMembers = new Set();
    countedMembers = new Set();
    if (!like) return;
    // room for every target of the list's catalog (a slot list's `targetInfo`: its tile set's, all of them): the pool
    // and the object made once, not again as tiles bring targets (the object's layers all uploaded again each time)
    const capacity = Math.max(MinPoolMembers, Math.ceil(store.pool.count * PoolHeadroom), targetCount + 1);
    const { type, state } = like.renderObject;
    const values = like.renderObject.values as unknown as MergeableValues;
    for (let m = 0; m <= capacity; ++m) {
      const member = createEmptyMergeMember(type as MergedType, values);
      ValueCell.update(member.boundingSphere, NoBounds);
      // a color a member: the object's colors a table by segment, its members' all of one type from the start
      if (byMember) setMemberColor(member, DefaultMemberColor);
      pool.push(createRenderObject(type, member as never, state, materialId));
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
    // a coarse sample's instances drawn larger (`ParticleList.scale`)
    objectSphere.radius += reach * store.maxScale;
    const o = last * 4;
    const b = poolBounds;
    const c = objectSphere.center;
    if (b[o] === c[0] && b[o + 1] === c[1] && b[o + 2] === c[2] && b[o + 3] === objectSphere.radius) return;
    b[o] = c[0];
    b[o + 1] = c[1];
    b[o + 2] = c[2];
    b[o + 3] = objectSphere.radius;
    ValueCell.update((pool[last].values as unknown as MergeableValues).boundingSphere, Sphere3D.clone(objectSphere));
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
  /** `block`'s colors into its store pages' texels: the ranges written */
  function fillColors(block: StoreBlock, theme: Theme, ranges: number[]) {
    const page = pages.get(block.key)!;
    const n = block.count;
    if (colorScratch.length < n * 4) colorScratch = new Uint8Array(n * 4);
    fillParticleColors(theme.color as ColorTheme<any>, page.list, colorScratch, 0);
    const array = colors.ref.value.array;
    for (let k = 0; k < block.pages.length; ++k) {
      const lo = k << StorePageBits;
      const hi = Math.min(n, lo + StorePageSize);
      const at = block.pages[k] << StorePageBits;
      array.set(colorScratch.subarray(lo * 4, hi * 4), at * 4);
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
    const targets = next.targetMapping ?? new Map<number, ParticleTarget>();
    const propsChanged = !currentProps || !PD.areEqual(ParticleTargetRepresentationParams, currentProps, props);
    const sizeChanged = !currentTheme || !SizeTheme.areEqual(theme.size, currentTheme.size);
    const colorChanged = !currentTheme || !ColorTheme.areEqual(theme.color, currentTheme.color);
    const newMembers = propsChanged || sizeChanged || !currentTargets || !isSameTargets(currentTargets, targets);
    if (newMembers) {
      byTarget = new Map();
      acceptedByTarget = new Map();
    }
    // a slot list of no pages (an empty view) keeps the mode: its theme colors no target, made again for the next
    const memberColors = next.pages?.pages.length === 0 ? byMember : isByMember(theme, targets);
    const modeChanged = memberColors !== byMember;
    byMember = memberColors;
    // the same pages again (a theme, props or clipping changed): no page changed, not all of them written again
    const changes =
      !newMembers && !modeChanged && list
        ? list.pages && list.pages === next.pages
          ? { from: list.pages, added: [], removed: [], posed: [] }
          : ParticleSlots.getChanges(list, next)
        : undefined;
    if (!changes) reset();

    const placement = placementOf(targets, props);
    const written: StoreBlock[] = [];
    const add = (page: ParticlePage) => {
      pages.set(page.start, page);
      // its box where known (a tile's list: its bounds from the tiling's pass), not a pass over its words
      const known = page.list._propertyData[Particle.BoundaryDescriptor.name] as Boundary | undefined;
      written.push(store.add(page.start, page.list, page.start, placement, known?.box));
    };
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
    } else {
      ParticleSlots.forEachPage(next, add);
    }
    const storeChanges = store.takeChanges();
    // appended one by one: a spread of a store's every page (a reset of a large view) is past the arguments' limit
    arrayExtend(pending.written, storeChanges.written);
    arrayExtend(pending.released, storeChanges.released);
    arrayExtend(pending.runs, storeChanges.runs);
    arrayExtend(pending.extras, storeChanges.extras);
    pending.grown = pending.grown || storeChanges.grown;

    // the pool: made where the store's members are more (the object made again); a member filled where the store
    // gave it another target (its geometry allocated again, `Merged.replaceGeometry`)
    const poolMade = store.pool.count > pool.length - 1;
    if (poolMade) await createPool(ctx, targets, theme, props, next.targetInfo?.size ?? 0);
    for (let m = 0; m < store.pool.count; ++m) {
      const t = store.pool.targets[m];
      if (t < 0 || poolTargets[m] === t) continue;
      const member = await getMember(ctx, t, targets, theme, props);
      const values = pool[m].values as unknown as MergeableValues;
      fillMergeMember(values, member.renderObject.values as never);
      if (byMember) setMemberColor(values, (theme.color as ColorTheme<any>).targetColor!(t)!);
      // a member's bounds none of its own: the object's the last member's
      ValueCell.update(values.boundingSphere, NoBounds);
      poolTargets[m] = t;
      object?.merged.replaceGeometry(m);
      changedMembers.add(m);
      storeChanges.members.add(m);
      // a pose places the geometry's point x at p + R (x - c): within |C - c| + r of its particle
      const ibs = values.invariantBoundingSphere.ref.value;
      reaches[m] = Vec3.distance(ibs.center, ParticleTarget.getCenter(targets.get(t)!) ?? Origin) + ibs.radius;
    }
    for (const m of storeChanges.members) {
      pending.members.add(m);
      if (m >= pool.length - 1) continue;
      setMemberCount(m, m < store.pool.count ? store.pool.counts[m] : 0);
    }
    setObjectBounds(next);

    // colors by member: the theme's changed, every member's again (its segment's texel)
    if (byMember && colorChanged && changes) {
      const targetColor = (theme.color as ColorTheme<any>).targetColor!;
      for (let m = 0; m < store.pool.count; ++m) {
        const t = poolTargets[m];
        if (t < 0) continue;
        setMemberColor(pool[m].values as unknown as MergeableValues, targetColor(t)!);
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
    if (clearMarkers) {
      markers.ref.value.array.fill(0);
      ValueCell.updateIfChanged(markerCells.markerAverage, 0);
      ValueCell.updateIfChanged(markerCells.markerStatus, 0);
      ValueCell.updateIfChanged(markerCells.uMarker, 0);
      TextureChanges.update(markers, undefined);
    } else if (written.length > 0 && markerCells.markerStatus.ref.value !== 0) {
      const writtenRanges: number[] = [];
      for (const block of written) pushBlockRanges(block, writtenRanges);
      const array = markers.ref.value.array;
      for (let k = 0; k < writtenRanges.length; k += 2) array.fill(0, writtenRanges[k], writtenRanges[k + 1]);
      const average = getMarkersAverage(array, store.pageCapacity * StorePageSize);
      ValueCell.updateIfChanged(markerCells.markerAverage, average);
      ValueCell.updateIfChanged(markerCells.markerStatus, average === 0 ? 0 : -1);
      ValueCell.updateIfChanged(markerCells.uMarker, average === 0 ? 0 : -1);
      TextureChanges.update(markers, writtenRanges);
    }

    if (!object && pool.length > 0) {
      object = createMergedRenderObject(pool, { gpuCulling: true, source: createSource() });
      occluderTiles = undefined;
    }
    // what its shown tiles hide (their opacities): the depth pyramid primed from it (`HiZPass.setOccluder`); a view's
    // step the last grid updated in place by the tiles it added and let go, not made again
    const view = isTargetOccluderOff(props) ? undefined : ParticleTileView.of(next);
    if (object && view?.tiles !== occluderTiles) {
      occluderTiles = view?.tiles;
      const last = occluder?.set === view?.set ? occluder : undefined;
      const grid = view ? ParticleTileOpaque.updateDensityGrid(view.set, view.tiles, last) : undefined;
      occluder = view && grid ? { set: view.set, grid, tiles: view.tiles } : undefined;
      setOccluderGrid(object.state, grid);
    }

    list = next;
    currentTheme = theme;
    currentProps = { ...props };
    currentTargets = targets;
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

  /** a slot's store index, -1 where none */
  function indexOf(slot: number) {
    const page = ParticleSlots.findPage(list!, slot);
    if (page < 0) return -1;
    const block = store.getBlock(list!.pages!.pages[page].start);
    if (!block) return -1;
    const j = slot - block.start;
    if (j < 0 || j >= block.count) return -1;
    return (block.pages[j >> StorePageBits] << StorePageBits) + (j & (StorePageSize - 1));
  }

  function getLoci(pickingId: PickingId): ModelLoci | undefined {
    if (!object || !list || pickingId.objectId !== object.id) return undefined;
    const slot = slotOf(pickingId.instanceId);
    if (slot < 0) return EmptyLoci;
    return Particle.Loci(list, OrderedSet.ofSingleton(slot));
  }

  function mark(loci: ModelLoci, action: MarkerAction) {
    if (!object || !list) return false;
    const capacity = store.pageCapacity * StorePageSize;
    ValueCell.updateIfChanged(markValues.instanceCount, capacity);
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
      for (const ro of pool) apply(ro);
    },
    updateClip(clip) {
      // the object's (the cull reads it), its members' (a merged object's cells a member's), the members kept to fill a
      // member again (`fillMergeMember` copies their values)
      if (object) BaseGeometry.updateClip(object.values as never, clip);
      for (const ro of pool) BaseGeometry.updateClip(ro.values as never, clip);
      for (const m of byTarget.values()) BaseGeometry.updateClip(m.renderObject.values as never, clip);
      // the next full update's props compared with this clip: not changed props, no reset
      if (currentProps) currentProps = { ...currentProps, clip };
    },
    getLoci,
    mark,
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
