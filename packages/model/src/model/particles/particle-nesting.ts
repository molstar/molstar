/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval, OrderedSet } from '@molstar/core/data/int';
import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { Sphere3D } from '@molstar/core/math/geometry';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { CustomProperties, CustomPropertyDescriptor } from '../custom-property.js';
import {
  Particle,
  type ParticleCompartmentInfo,
  type ParticleEntityInfo,
  type ParticleList,
  ParticleRotations,
  ParticleTarget,
  type ParticleTargetInfo,
  ParticleValues,
} from './particle-list.js';
import { DeferredPageBits } from './particle-page-cache.js';
import { type ParticlePage, ParticleSlots } from './particle-slots.js';
import { ParticleDensityGrid, ParticleTileOpaque } from './particle-tile-opaque.js';
import { ParticleTileCatalog, ParticleTiles } from './particle-tiles.js';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';

// Nested particle lists: a particle of a `particles` target (a template) places the template's whole list, its pose
// composed with each of the template's particles'. Stored once, expanded where drawn; an expanded instance's identity
// its top particle and its index in the template (`Path`)

/** A template's own particles' targets are leaves: a template places no template (two levels) */
export const MaxNestingDepth = 1;

/** An expanded instance: its particle of the top list, and its particle of that one's template (0 for a leaf) */
export interface ParticlePath {
  top: number;
  local: number;
}

/** How the top list's particles expand: `count` instances, top particle `i` the ones `[offsets[i], offsets[i + 1])` */
export interface ParticleExpansion {
  readonly count: number;
  /** none where the list has no template: an instance a particle */
  readonly offsets?: Float64Array;
}

const ExpansionDescriptor = CustomPropertyDescriptor({ name: 'particle-expansion' });
const FlattenDescriptor = CustomPropertyDescriptor({ name: 'particle-flatten' });

export namespace ParticleNesting {
  /** The template target `id` of `list` stands for, if it is one */
  export function getTemplate(
    list: ParticleList,
    id: number,
  ): Extract<ParticleTarget, { kind: 'particles' }> | undefined {
    const target = list.targetMapping?.get(id);
    return target?.kind === 'particles' ? target : undefined;
  }

  /** Whether some target of `list` is a template */
  export function hasTemplates(list: ParticleList): boolean {
    if (!list.targetMapping) return false;
    for (const target of list.targetMapping.values()) if (target.kind === 'particles') return true;
    return false;
  }

  /**
   * Throws where `list` nests deeper than `MaxNestingDepth` (a template with templates: a cycle among them), or a
   * template is a slot list (its pages not placed); what makes a nested list calls it
   */
  export function check(list: ParticleList) {
    if (!list.targetMapping) return;
    for (const [id, target] of list.targetMapping) {
      if (target.kind !== 'particles') continue;
      const template = target.particles;
      if (template === list) throw new Error(`particle target ${id} places its own list`);
      if (template.pages) throw new Error(`particle target ${id}: a slot list as a template`);
      if (hasTemplates(template)) {
        throw new Error(`particle target ${id}: a template with templates (nesting deeper than ${MaxNestingDepth})`);
      }
    }
  }

  /** How `list`'s particles expand (cached on it): an instance a leaf particle, a template's count a template's */
  export function getExpansion(list: ParticleList): ParticleExpansion {
    const cached = list._propertyData[ExpansionDescriptor.name] as ParticleExpansion | undefined;
    if (cached) return cached;
    let expansion: ParticleExpansion;
    if (!hasTemplates(list)) {
      expansion = { count: list.count };
    } else {
      const offsets = new Float64Array(list.count + 1);
      let o = 0;
      ParticleValues.forEachRun(list.targets, list.count, (start, end, id) => {
        const n = getTemplate(list, id)?.particles.count ?? 1;
        for (let i = start; i < end; ++i) {
          offsets[i] = o;
          o += n;
        }
      });
      offsets[list.count] = o;
      expansion = { count: o, offsets };
    }
    list.customProperties.add(ExpansionDescriptor);
    list._propertyData[ExpansionDescriptor.name] = expansion;
    return expansion;
  }

  /** The instance count of `list` expanded */
  export function getExpandedCount(list: ParticleList) {
    return getExpansion(list).count;
  }

  /** The index of the instance `path` in `list` expanded (`flatten`'s order) */
  export function getExpandedIndex(list: ParticleList, top: number, local: number) {
    const { offsets } = getExpansion(list);
    return offsets ? offsets[top] + local : top;
  }

  /** The path of instance `index` of `list` expanded into `out` (a binary search of the top particles) */
  export function getPath(list: ParticleList, index: number, out: ParticlePath): ParticlePath {
    const { offsets } = getExpansion(list);
    if (!offsets) {
      out.top = index;
      out.local = 0;
      return out;
    }
    let lo = 0;
    let hi = list.count - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (offsets[mid] <= index) lo = mid;
      else hi = mid - 1;
    }
    out.top = lo;
    out.local = index - offsets[lo];
    return out;
  }

  /** Ids of maps merged: a map's ids moved past those taken (-1, none, kept) */
  function createIdSpace<T>() {
    const merged = new Map<number, T>();
    let next = 0;
    return {
      merged,
      /** the ids of `map` as they are */
      keep(map: ReadonlyMap<number, T> | undefined) {
        if (!map) return;
        for (const [id, value] of map) {
          merged.set(id, value);
          if (id >= next) next = id + 1;
        }
      },
      /** the ids of `map` moved past those taken: the move (its first id) */
      move(map: ReadonlyMap<number, T> | undefined) {
        const base = next;
        if (!map) return base;
        for (const [id, value] of map) {
          merged.set(base + id, value);
          if (base + id >= next) next = base + id + 1;
        }
        return base;
      },
    };
  }

  const _p = Vec3();
  const _c = Vec3();
  const _q = Quat();
  const _qc = Quat();

  /** A template's ids in its list's flattened list: its target, entity and compartment ids moved by these */
  export interface IdMoves {
    readonly target: number;
    readonly entity: number;
    readonly compartment: number;
  }

  /** The id space of `list`'s flattened list (`getIdSpace`) */
  export interface IdSpace {
    /** per template its moves */
    readonly moves: ReadonlyMap<ParticleList, IdMoves>;
    readonly targetMapping: ReadonlyMap<number, ParticleTarget>;
    readonly targetInfo: ReadonlyMap<number, ParticleTargetInfo>;
    readonly entityInfo: ReadonlyMap<number, ParticleEntityInfo>;
    readonly compartmentInfo: ReadonlyMap<number, ParticleCompartmentInfo>;
  }

  const IdSpaceDescriptor = CustomPropertyDescriptor({ name: 'particle-nesting-ids' });

  /**
   * The ids of `list`'s flattened list (cached on it): the top list's leaves keep their target ids; a template's moved
   * past every id of the top list (its templates' too: no id of the nested list reused for another target), once a
   * template (shared by its targets), the templates in the order of their first copy; its infos at the same moved ids;
   * its entities and compartments alike. What `flatten` writes, and what draws a list's templates by their own lists
   * keys their targets by (the themes, of the flattened list, agree)
   */
  export function getIdSpace(list: ParticleList): IdSpace {
    const cached = list._propertyData[IdSpaceDescriptor.name] as IdSpace | undefined;
    if (cached) return cached;
    const targetMapping = new Map<number, ParticleTarget>();
    const targetInfo = new Map<number, ParticleTargetInfo>();
    let nextTarget = 0;
    const takeTarget = (id: number, target: ParticleTarget | undefined, info: ParticleTargetInfo | undefined) => {
      if (target) targetMapping.set(id, target);
      if (info) targetInfo.set(id, info);
      if (id >= nextTarget) nextTarget = id + 1;
    };
    for (const [id, target] of list.targetMapping ?? []) {
      if (target.kind !== 'particles') takeTarget(id, target, undefined);
      else if (id >= nextTarget) nextTarget = id + 1;
    }
    for (const [id, info] of list.targetInfo) {
      if (!getTemplate(list, id)) takeTarget(id, undefined, info);
      else if (id >= nextTarget) nextTarget = id + 1;
    }
    const entityIds = createIdSpace<ParticleEntityInfo>();
    const compartmentIds = createIdSpace<ParticleCompartmentInfo>();
    entityIds.keep(list.entityInfo);
    compartmentIds.keep(list.compartmentInfo);
    const moves = new Map<ParticleList, IdMoves>();
    ParticleValues.forEachRun(list.targets, list.count, (_start, _end, id) => {
      const template = getTemplate(list, id)?.particles;
      if (!template || moves.has(template)) return;
      const target = nextTarget;
      for (const [id, info] of template.targetInfo) takeTarget(target + id, undefined, info);
      if (template.targetMapping) for (const [id, t] of template.targetMapping) takeTarget(target + id, t, undefined);
      // a template without infos or targets still moves its ids past the top list's
      nextTarget = Math.max(nextTarget, target + 1);
      moves.set(template, {
        target,
        entity: entityIds.move(template.entityInfo),
        compartment: compartmentIds.move(template.compartmentInfo),
      });
    });
    const space: IdSpace = {
      moves,
      targetMapping,
      targetInfo,
      entityInfo: entityIds.merged,
      compartmentInfo: compartmentIds.merged,
    };
    list.customProperties.add(IdSpaceDescriptor);
    list._propertyData[IdSpaceDescriptor.name] = space;
    return space;
  }

  const DensityDescriptor = CustomPropertyDescriptor({ name: 'particle-nesting-density' });
  const _dc = Vec3();
  const _dp = Vec3();
  const _dq = Quat();

  /**
   * The optical densities of `list` (with templates) over the centers' box of its instances, `longest` blocks along its
   * longest side (cached on it; undefined where none): each template's grid (`ParticleDensityGrid`: a pack's, else
   * computed once on the template) placed by each copy: a list's block the copy's template grid covers whole the average
   * of the template blocks its box spans in the template's frame (a table of the grid's sums: eight reads a block), a
   * block on a copy's edge none (the gaps between copies open, what lies behind seen through them; the blocks coarser
   * than the gaps would close them); its own leaves' molecules as a flat list's. No copy expanded: the list's blocks a
   * copy covers once a copy (at 4x4x4 isg copies about 15³ each, not the template's 230 thousand blocks). What primes
   * the GPU's occlusion of a list drawn whole (`setOccluderGrid`), as a flat list's grid does (without it what lies
   * behind dense content is culled only by the last frame's depth: rotating, the second phase draws it)
   */
  export function getDensityGrid(list: ParticleList, longest = 64): ParticleDensityGrid | undefined {
    if (list._propertyData[DensityDescriptor.name] !== undefined) {
      return list._propertyData[DensityDescriptor.name] || undefined;
    }
    const templates = new Map<ParticleList, ParticleDensityGrid | undefined>();
    for (const t of getIdSpace(list).moves.keys()) {
      templates.set(t, ParticleDensityGrid.get(t) ?? ParticleDensityGrid.compute(t));
    }
    // the box: the leaves' centers and the copies' grids' corners placed
    const box = Box3D.setEmpty(Box3D());
    const forEachCopy = (f: (grid: ParticleDensityGrid, i: number, center: Vec3) => void) => {
      ParticleValues.forEachRun(list.targets, list.count, (start, end, id) => {
        const template = getTemplate(list, id);
        const grid = template && templates.get(template.particles);
        if (!grid) return;
        const center = ParticleTarget.getCenter(template)!;
        for (let i = start; i < end; ++i) f(grid, i, center);
      });
    };
    const place = (i: number, center: Vec3, x: Vec3, out: Vec3) => {
      GridPositions.get(list.coordinates, i, _dp);
      ParticleRotations.get(list.rotations, i, _dq);
      Vec3.transformQuat(out, Vec3.sub(out, x, center), _dq);
      return Vec3.add(out, out, _dp);
    };
    ParticleValues.forEachRun(list.targets, list.count, (start, end, id) => {
      if (getTemplate(list, id)) return;
      for (let i = start; i < end; ++i) Box3D.add(box, GridPositions.get(list.coordinates, i, _dp));
    });
    forEachCopy((grid, i, center) => {
      for (let k = 0; k < 8; ++k) {
        Vec3.set(
          _dc,
          k & 1 ? grid.max[0] : grid.min[0],
          k & 2 ? grid.max[1] : grid.min[1],
          k & 4 ? grid.max[2] : grid.min[2],
        );
        Box3D.add(box, place(i, center, _dc, _dc));
      }
    });
    if (Box3D.isEmpty(box)) {
      list._propertyData[DensityDescriptor.name] = null;
      return undefined;
    }
    const dimensions = ParticleTileOpaque.getGridDimensions(box, longest);
    const [bx, by, bz] = dimensions;
    const o = box.min;
    const s = [(box.max[0] - o[0]) / bx || 1, (box.max[1] - o[1]) / by || 1, (box.max[2] - o[2]) / bz || 1];
    // the leaves' molecules (their targets' volumes, the catalog's; the copies' particles of none)
    const leaves = ParticleTileOpaque.computeDensity(box, dimensions, (sink) => {
      const targets = ParticleTileCatalog.fromList(list).targets;
      ParticleTileOpaque.forEachMolecule(list, (id) => targets.get(id)?.volume, sink);
    });
    const density = leaves ?? new Float32Array(bx * by * bz);
    let any = !!leaves;
    // per template its grid's sums (a block's the sum of those before it on every axis): a box of blocks' sum by eight
    // reads, a list's block's average density whatever the template blocks it spans
    const sums = new Map<ParticleDensityGrid, Float64Array>();
    const getSums = (grid: ParticleDensityGrid) => {
      let t = sums.get(grid);
      if (t) return t;
      const [tx, ty, tz] = grid.dimensions;
      const sx = tx + 1;
      const sxy = sx * (ty + 1);
      t = new Float64Array(sxy * (tz + 1));
      const d = grid.density;
      for (let z = 0, b = 0; z < tz; ++z) {
        for (let y = 0; y < ty; ++y) {
          for (let x = 0; x < tx; ++x, ++b) {
            const o1 = x + 1 + sx * (y + 1) + sxy * (z + 1);
            t[o1] = d[b] + t[o1 - 1] + t[o1 - sx] + t[o1 - sxy] - t[o1 - 1 - sx] - t[o1 - 1 - sxy] - t[o1 - sx - sxy];
            t[o1] += t[o1 - 1 - sx - sxy];
          }
        }
      }
      sums.set(grid, t);
      return t;
    };
    const inverse = Quat();
    const base = Vec3();
    const ex = Vec3();
    const ey = Vec3();
    const ez = Vec3();
    const lo = [0, 0, 0];
    const hi = [0, 0, 0];
    const at = [0, 0, 0];
    forEachCopy((grid, i, center) => {
      const [tx, ty, tz] = grid.dimensions;
      const sx = tx + 1;
      const sxy = sx * (ty + 1);
      const table = getSums(grid);
      const ts = [(grid.max[0] - grid.min[0]) / tx, (grid.max[1] - grid.min[1]) / ty, (grid.max[2] - grid.min[2]) / tz];
      const eps = 1e-3 * Math.max(ts[0], ts[1], ts[2]);
      GridPositions.get(list.coordinates, i, _dp);
      Quat.invert(inverse, ParticleRotations.get(list.rotations, i, _dq));
      // the list's block (x, y, z)'s corner in the template's frame: base + x ex + y ey + z ez (the copy's pose undone),
      // its eight corners' box that plus the unit block's corners' extremes (lo, hi)
      Vec3.add(base, Vec3.transformQuat(base, Vec3.sub(base, o as Vec3, _dp), inverse), center);
      Vec3.transformQuat(ex, Vec3.set(ex, s[0], 0, 0), inverse);
      Vec3.transformQuat(ey, Vec3.set(ey, 0, s[1], 0), inverse);
      Vec3.transformQuat(ez, Vec3.set(ez, 0, 0, s[2]), inverse);
      for (let a = 0; a < 3; ++a) {
        lo[a] = Math.min(0, ex[a]) + Math.min(0, ey[a]) + Math.min(0, ez[a]);
        hi[a] = Math.max(0, ex[a]) + Math.max(0, ey[a]) + Math.max(0, ez[a]);
      }
      // the list's blocks the copy's placed grid may cover: its corners' placed box
      const range = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
      for (let k = 0; k < 8; ++k) {
        Vec3.set(
          _dc,
          k & 1 ? grid.max[0] : grid.min[0],
          k & 2 ? grid.max[1] : grid.min[1],
          k & 4 ? grid.max[2] : grid.min[2],
        );
        place(i, center, _dc, _dc);
        for (let a = 0; a < 3; ++a) {
          const c = (_dc[a] - o[a]) / s[a];
          if (c < range[a]) range[a] = c;
          if (c > range[a + 3]) range[a + 3] = c;
        }
      }
      const x0 = Math.max(0, Math.floor(range[0]));
      const y0 = Math.max(0, Math.floor(range[1]));
      const z0 = Math.max(0, Math.floor(range[2]));
      const x1 = Math.min(bx, Math.ceil(range[3]));
      const y1 = Math.min(by, Math.ceil(range[4]));
      const z1 = Math.min(bz, Math.ceil(range[5]));
      for (let z = z0; z < z1; ++z) {
        for (let y = y0; y < y1; ++y) {
          for (let x = x0; x < x1; ++x) {
            // the block's box in the template's frame; a block not within the template's grid whole none (a block
            // straddling a copy's edge holds a gap between copies, what lies behind seen through it: the blocks
            // coarser than the gaps would close them)
            let inside = true;
            for (let a = 0; a < 3; ++a) {
              at[a] = base[a] + x * ex[a] + y * ey[a] + z * ez[a];
              if (at[a] + lo[a] < grid.min[a] - eps || at[a] + hi[a] > grid.max[a] + eps) inside = false;
            }
            if (!inside) continue;
            // its average: the sum over the template blocks the box spans, by their count
            const i0 = Math.max(0, Math.floor((at[0] + lo[0] - grid.min[0]) / ts[0]));
            const j0 = Math.max(0, Math.floor((at[1] + lo[1] - grid.min[1]) / ts[1]));
            const k0 = Math.max(0, Math.floor((at[2] + lo[2] - grid.min[2]) / ts[2]));
            const i1 = Math.min(tx, Math.max(i0 + 1, Math.ceil((at[0] + hi[0] - grid.min[0]) / ts[0])));
            const j1 = Math.min(ty, Math.max(j0 + 1, Math.ceil((at[1] + hi[1] - grid.min[1]) / ts[1])));
            const k1 = Math.min(tz, Math.max(k0 + 1, Math.ceil((at[2] + hi[2] - grid.min[2]) / ts[2])));
            const sum =
              table[i1 + sx * j1 + sxy * k1] -
              table[i0 + sx * j1 + sxy * k1] -
              table[i1 + sx * j0 + sxy * k1] -
              table[i1 + sx * j1 + sxy * k0] +
              table[i0 + sx * j0 + sxy * k1] +
              table[i0 + sx * j1 + sxy * k0] +
              table[i1 + sx * j0 + sxy * k0] -
              table[i0 + sx * j0 + sxy * k0];
            const v = sum / ((i1 - i0) * (j1 - j0) * (k1 - k0));
            if (!(v > 0)) continue;
            density[x + bx * (y + by * z)] += v;
            any = true;
          }
        }
      }
    });
    const grid: ParticleDensityGrid | undefined = any
      ? { min: Vec3.clone(box.min), max: Vec3.clone(box.max), dimensions, density }
      : undefined;
    list.customProperties.add(DensityDescriptor);
    list._propertyData[DensityDescriptor.name] = grid ?? null;
    return grid;
  }

  const LeafListsDescriptor = CustomPropertyDescriptor({ name: 'particle-nesting-leaf-lists' });
  const TopListDescriptor = CustomPropertyDescriptor({ name: 'particle-nesting-top-list' });
  const ThemeListDescriptor = CustomPropertyDescriptor({ name: 'particle-nesting-theme-list' });
  const TemplateLevelsDescriptor = CustomPropertyDescriptor({ name: 'particle-nesting-template-levels' });

  /**
   * A slot list of `list`'s leaves by level (cached on it): its own particles (`getTopList`) from slot 0, each template's
   * (`getLeafList`, the templates in their ids' order) after; the flattened list's tables. What themes a list with
   * templates without expanding it: a theme of a slot list colors a page by its own list (a template particle once for
   * all its copies, by its leaf's data), not by an index of the flattened list
   */
  export function getThemeList(list: ParticleList): ParticleList {
    const cached = list._propertyData[ThemeListDescriptor.name] as ParticleList | undefined;
    if (cached) return cached;
    const ids = getIdSpace(list);
    const pages: ParticlePage[] = [{ start: 0, list: getTopList(list) }];
    let count = list.count;
    for (const template of ids.moves.keys()) {
      pages.push({ start: count, list: getLeafList(list, template) });
      count += template.count;
    }
    const themed = ParticleSlots.create(
      {
        label: list.label,
        entryId: list.entryId,
        targetInfo: ids.targetInfo,
        targetMapping: ids.targetMapping,
        entityInfo: ids.entityInfo.size > 0 ? ids.entityInfo : undefined,
        compartmentInfo: ids.compartmentInfo.size > 0 ? ids.compartmentInfo : undefined,
        sourceData: list.sourceData,
      },
      count,
      pages,
    );
    list.customProperties.add(ThemeListDescriptor);
    list._propertyData[ThemeListDescriptor.name] = themed;
    return themed;
  }

  /**
   * `list`'s own particles with the flattened list's tables (cached on it): its ids as they are (the flattened list
   * keeps the top list's), its infos and targets the flattened list's (a template's target none: its particles
   * placing copies, no target of their own drawn); what colors and sizes its leaves by the flattened list's themes
   */
  export function getTopList(list: ParticleList): ParticleList {
    const cached = list._propertyData[TopListDescriptor.name] as ParticleList | undefined;
    if (cached) return cached;
    const ids = getIdSpace(list);
    const top: ParticleList = {
      ...list,
      targetInfo: ids.targetInfo as Map<number, ParticleTargetInfo>,
      targetMapping: ids.targetMapping as Map<number, ParticleTarget>,
      entityInfo: ids.entityInfo.size > 0 ? (ids.entityInfo as Map<number, ParticleEntityInfo>) : undefined,
      compartmentInfo:
        ids.compartmentInfo.size > 0 ? (ids.compartmentInfo as Map<number, ParticleCompartmentInfo>) : undefined,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    list.customProperties.add(TopListDescriptor);
    list._propertyData[TopListDescriptor.name] = top;
    return top;
  }

  /**
   * `template`'s particles (a template of `list`) in `list`'s flattened ids (cached on `list`): its targets, entities and
   * compartments moved (`getIdSpace`; one without an entity or compartment of its own none, not a copy's top
   * particle's: that is by copy), the flattened list's infos and targets; its positions, rotations and radii its own
   * (views, in its own frame). What colors and sizes a template particle once for all its copies by the flattened
   * list's themes, and what a store of its particles keys their targets by; no copy expanded
   */
  export function getLeafList(list: ParticleList, template: ParticleList): ParticleList {
    let byTemplate = list._propertyData[LeafListsDescriptor.name] as Map<ParticleList, ParticleList> | undefined;
    if (!byTemplate) {
      list.customProperties.add(LeafListsDescriptor);
      list._propertyData[LeafListsDescriptor.name] = byTemplate = new Map();
    }
    const cached = byTemplate.get(template);
    if (cached) return cached;
    const ids = getIdSpace(list);
    const m = ids.moves.get(template);
    if (!m) throw new Error('particle nesting: not a template of the list');
    const n = template.count;
    const moved = (v: ParticleValues<Int32Array> | undefined, base: number) =>
      v && ParticleValues.map(v, n, (x) => (x >= 0 ? base + x : -1));
    const leaves: ParticleList = {
      ...template,
      label: template.label,
      targets: moved(template.targets, m.target)!,
      targetInfo: ids.targetInfo as Map<number, ParticleTargetInfo>,
      targetMapping: ids.targetMapping as Map<number, ParticleTarget>,
      entities: moved(template.entities, m.entity),
      entityInfo: ids.entityInfo.size > 0 ? (ids.entityInfo as Map<number, ParticleEntityInfo>) : undefined,
      compartments: moved(template.compartments, m.compartment),
      compartmentInfo:
        ids.compartmentInfo.size > 0 ? (ids.compartmentInfo as Map<number, ParticleCompartmentInfo>) : undefined,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    // the template's positions and radii: its boundary and its order (a pack's, set where read)
    if (Particle.hasBoundary(template)) Particle.setBoundary(leaves, Particle.getBoundary(template));
    if (Particle.isSampleOrdered(template)) Particle.setSampleOrdered(leaves);
    // its particles not read (a pack opened by range: read by page when drawn, by its tiles' boxes before)
    const deferred = Particle.getDeferred(template);
    if (deferred) Particle.setDeferred(leaves, deferred);
    const directory = ParticleTiles.getDirectory(template);
    if (directory) ParticleTiles.setDirectory(leaves, directory);
    byTemplate.set(template, leaves);
    return leaves;
  }

  /**
   * a template's coarse level kept where it keeps this many particles a target on average (`getTemplateLevels`). Why:
   * a level of a particle or two a target draws each of them as large as the template (the isg cell's fourth: 199
   * particles of 100 targets, up to 63 times larger, a bound padded by it; its third 960, up to 22.6)
   */
  const MinLevelKeptPerTarget = 8;
  /** a sample's targets (its runs) */
  const runsOf = (s: Particle.SampleLevel) => {
    const scales = s.list.scales;
    return scales && ParticleValues.isRuns(scales) ? scales.offsets.length - 1 : 1;
  };

  /** a template's coarse levels at most (`getTemplateLevels`) */
  export const MaxTemplateLevels = 4;

  /** A coarse level of a template (`getTemplateLevels`) */
  export interface TemplateLevel {
    /** from 1: of each target every `8^level`-th particle */
    readonly level: number;
    /** its particles, as `getLeafList`'s (its ids the flattened list's), drawn larger (`ParticleList.scales`) */
    readonly list: ParticleList;
    /** a particle's index in the template */
    readonly indices: Int32Array;
  }

  /**
   * `template`'s coarse levels (a template of `list`, cached on `list`): its leaf list (`getLeafList`) sampled by target
   * (`Particle.sampleByTarget`), level by level while a level keeps at most half the one before and eight particles a
   * target on average (`MinLevelKeptPerTarget`; at most
   * `MaxTemplateLevels`; none for a template of a few particles a target). What its far copies draw: a copy a level by
   * its distance, its sample drawn larger so it fills about the space the template does
   */
  export function getTemplateLevels(list: ParticleList, template: ParticleList): readonly TemplateLevel[] {
    let byTemplate = list._propertyData[TemplateLevelsDescriptor.name] as
      | Map<ParticleList, readonly TemplateLevel[]>
      | undefined;
    if (!byTemplate) {
      list.customProperties.add(TemplateLevelsDescriptor);
      list._propertyData[TemplateLevelsDescriptor.name] = byTemplate = new Map();
    }
    const cached = byTemplate.get(template);
    if (cached) return cached;
    const leaves = getLeafList(list, template);
    const levels: TemplateLevel[] = [];
    let before = leaves.count;
    // stored with it (a pack's): read, not sampled (the yeast cell's 461 ms at load, its gathers scattered
    // tile-ordered); kept by the same rules
    const stored = Particle.getStoredLevels(template);
    if (stored) {
      const m = getIdSpace(list).moves.get(template)!;
      const listScale = leaves.scale !== undefined && leaves.scale > 0 ? leaves.scale : 1;
      for (const s of stored) {
        const runCount = s.targets.length;
        const targets = new Int32Array(runCount);
        const scales = new Float32Array(runCount);
        for (let d = 0; d < runCount; ++d) {
          const t = s.targets[d];
          targets[d] = t >= 0 ? m.target + t : -1;
          scales[d] = s.scales[d] * listScale;
        }
        if (s.keys.length * 2 > before || s.keys.length < MinLevelKeptPerTarget * runCount) break;
        const sample = Particle.levelOfRuns(leaves, s.keys, s.offsets, scales, {
          targets,
          coordinates: s.coordinates,
          rotations: s.rotations,
        });
        levels.push({ level: s.level, list: sample.list, indices: sample.indices });
        before = s.keys.length;
      }
      byTemplate.set(template, levels);
      return levels;
    }
    const order = Particle.getSampleOrder(leaves);
    let finer: Particle.SampleLevel | undefined;
    for (let level = 1; level <= MaxTemplateLevels; ++level) {
      // a level from the one before (every eighth of a target's kept), the first from all
      const sample = Particle.sampleByTarget(leaves, level, order, finer);
      if (!sample || sample.list.count * 2 > before) break;
      if (sample.list.count < MinLevelKeptPerTarget * runsOf(sample)) break;
      levels.push({ level, list: sample.list, indices: sample.indices });
      before = sample.list.count;
      finer = sample;
    }
    byTemplate.set(template, levels);
    return levels;
  }

  /**
   * `list` expanded (cached on it): an instance a particle, in `getExpandedIndex`'s order, a template's particle placed
   * by its top particle's pose (`q (p - c) + t`, its rotation `q qc`; `c` the template's center); a list without
   * templates as it is. The templates' targets, entities and compartments taken past the top list's ids; an instance
   * with no entity or compartment of its own takes its top particle's. Radii where the top list and every template
   * have them; fibers, attributes, keys and cells not kept for a list with templates. The reference every consumer
   * expanding on its own is tested against
   */
  export function flatten(list: ParticleList): ParticleList {
    if (!hasTemplates(list)) return list;
    const cached = list._propertyData[FlattenDescriptor.name] as ParticleList | undefined;
    if (cached) return cached;
    check(list);
    // every copy's particles: its templates' read (not by range: `Particle.DeferredParticles`)
    for (const t of getIdSpace(list).moves.keys()) Particle.assertRead(t, 'flattening a list with templates');
    const { count: n } = getExpansion(list);
    const count = list.count;

    const ids = getIdSpace(list);
    const { targetMapping, targetInfo } = ids;
    const getMoves = (template: ParticleList) => ids.moves.get(template)!;

    const templates: (Extract<ParticleTarget, { kind: 'particles' }> | undefined)[] = [];
    const topTemplate = (i: number) => getTemplate(list, ParticleValues.get(list.targets, i));
    let withRotations = !!list.rotations;
    let withRadii = !!list.radii;
    let withEntities = !!list.entities;
    let withCompartments = !!list.compartments;
    for (let i = 0; i < count; ++i) templates[i] = topTemplate(i);
    for (const target of list.targetMapping!.values()) {
      if (target.kind !== 'particles') continue;
      const t = target.particles;
      withRotations ||= !!t.rotations;
      withRadii &&= !!t.radii;
      withEntities ||= !!t.entities;
      withCompartments ||= !!t.compartments;
    }

    const coordinates = new Float32Array(n * 3);
    const rotations = withRotations ? new Float32Array(n * 4) : undefined;
    const targets = new Int32Array(n);
    const radii = withRadii ? new Float32Array(n) : undefined;
    const entities = withEntities ? new Int32Array(n) : undefined;
    const compartments = withCompartments ? new Int32Array(n) : undefined;
    const value = (v: ParticleValues<Int32Array> | undefined, i: number) => (v ? ParticleValues.get(v, i) : -1);

    let o = 0;
    for (let i = 0; i < count; ++i) {
      const template = templates[i];
      GridPositions.get(list.coordinates, i, _p);
      ParticleRotations.get(list.rotations, i, _q);
      const entity = value(list.entities, i);
      const compartment = value(list.compartments, i);
      if (!template) {
        coordinates.set(_p, o * 3);
        if (rotations) Quat.toArray(_q, rotations, o * 4);
        targets[o] = ParticleValues.get(list.targets, i);
        if (radii) radii[o] = ParticleValues.get(list.radii!, i);
        if (entities) entities[o] = entity;
        if (compartments) compartments[o] = compartment;
        o += 1;
        continue;
      }
      const t = template.particles;
      const m = getMoves(t);
      const center = ParticleTarget.getCenter(template)!;
      const tx = _p[0],
        ty = _p[1],
        tz = _p[2];
      for (let j = 0, jl = t.count; j < jl; ++j, ++o) {
        GridPositions.get(t.coordinates, j, _p);
        Vec3.transformQuat(_p, Vec3.sub(_c, _p, center), _q);
        coordinates[o * 3] = _p[0] + tx;
        coordinates[o * 3 + 1] = _p[1] + ty;
        coordinates[o * 3 + 2] = _p[2] + tz;
        if (rotations)
          Quat.toArray(Quat.multiply(_qc, _q, ParticleRotations.get(t.rotations, j, _qc)), rotations, o * 4);
        targets[o] = m.target + ParticleValues.get(t.targets, j);
        if (radii) radii[o] = ParticleValues.get(t.radii!, j);
        if (entities) {
          const e = value(t.entities, j);
          entities[o] = e >= 0 ? m.entity + e : entity;
        }
        if (compartments) {
          const c = value(t.compartments, j);
          compartments[o] = c >= 0 ? m.compartment + c : compartment;
        }
      }
    }

    const path: ParticlePath = { top: 0, local: 0 };
    const flat: ParticleList = {
      entryId: list.entryId,
      label: list.label,
      count: n,
      targets: ParticleValues.compact(targets),
      targetInfo,
      targetMapping,
      entities: entities ? ParticleValues.compact(entities) : undefined,
      entityInfo: withEntities ? ids.entityInfo : undefined,
      compartments: compartments ? ParticleValues.compact(compartments) : undefined,
      compartmentInfo: withCompartments ? ids.compartmentInfo : undefined,
      coordinates,
      rotations,
      radii: radii ? ParticleValues.compact(radii) : undefined,
      getParticleLabel: (index) => {
        getPath(list, index, path);
        const template = templates[path.top];
        const top = list.getParticleLabel(path.top);
        return template ? `${top} › ${template.particles.getParticleLabel(path.local)}` : top;
      },
      sourceData: list.sourceData,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    list.customProperties.add(FlattenDescriptor);
    list._propertyData[FlattenDescriptor.name] = flat;
    return flat;
  }

  // Loci by path sets: `Particle.Loci`'s `indices` the top particles held whole (a copy its every instance), its
  // `partial` the copies held in part, each its template particles. Canonical (`createLoci`): a copy's set neither
  // empty nor all, no top in both; so two loci of the same instances are equal by their parts

  /** The instances top particle `top` of `list` expands to: its template's count, 1 a leaf */
  export function getCopyCount(list: ParticleList, top: number) {
    const { offsets } = getExpansion(list);
    return offsets ? offsets[top + 1] - offsets[top] : 1;
  }

  /**
   * The loci of `list` holding the top particles `whole` whole and, of each of `tops` (ascending), the template
   * particles `locals` (a copy's set past its template's count clipped; a set of all moved to `whole`, an empty one
   * dropped, a top already whole ignored)
   */
  export function createLoci(
    list: ParticleList,
    whole: OrderedSet<number>,
    tops: ArrayLike<number> = [],
    locals: readonly OrderedSet<number>[] = [],
  ): Particle.Loci {
    const wholeTops: number[] = [];
    const partTops: number[] = [];
    const partLocals: OrderedSet<number>[] = [];
    for (let k = 0; k < tops.length; ++k) {
      const t = tops[k];
      if (OrderedSet.has(whole, t)) continue;
      const n = getCopyCount(list, t);
      let l = locals[k];
      if (OrderedSet.size(l) > 0 && OrderedSet.max(l) >= n) l = OrderedSet.intersect(l, Interval.ofLength(n));
      const size = OrderedSet.size(l);
      if (size === 0) continue;
      if (size === n) wholeTops.push(t);
      else {
        partTops.push(t);
        partLocals.push(l);
      }
    }
    const indices = wholeTops.length > 0 ? OrderedSet.union(whole, OrderedSet.ofSortedArray(wholeTops)) : whole;
    return Particle.Loci(list, indices, { tops: Int32Array.from(partTops), locals: partLocals });
  }

  /** A copy (or a leaf) whole: top particle `top`'s every instance */
  export function ofCopy(list: ParticleList, top: number) {
    return Particle.Loci(list, OrderedSet.ofSingleton(top));
  }

  /** One instance: template particle `local` in copy `top` (a leaf's `local` 0) */
  export function ofInstance(list: ParticleList, top: number, local: number) {
    return createLoci(list, OrderedSet.Empty, [top], [OrderedSet.ofSingleton(local)]);
  }

  /** The top particles placing `template` (ascending), by their targets' runs: no pass over a template's particles */
  export function getCopiesOf(list: ParticleList, template: ParticleList): number[] {
    const tops: number[] = [];
    ParticleValues.forEachRun(list.targets, list.count, (start, end, id) => {
      if (getTemplate(list, id)?.particles !== template) return;
      for (let i = start; i < end; ++i) tops.push(i);
    });
    return tops;
  }

  /** Template particles `locals` of `template` in every copy of it: `[*, locals]` (a template's set, not a copy's) */
  export function ofTemplateParticles(list: ParticleList, template: ParticleList, locals: OrderedSet<number>) {
    const tops = getCopiesOf(list, template);
    return createLoci(
      list,
      OrderedSet.Empty,
      tops,
      tops.map(() => locals),
    );
  }

  /** The top particles `test` holds, whole (a query of the top level: by its targets, entities, compartments) */
  export function selectTop(list: ParticleList, test: (top: number) => boolean) {
    const tops: number[] = [];
    for (let i = 0, il = list.count; i < il; ++i) if (test(i)) tops.push(i);
    return Particle.Loci(list, OrderedSet.ofSortedArray(tops));
  }

  /**
   * The particles of `template` `test` holds, in every copy of it (a query of the template's level: its particles
   * tested once, not a copy each)
   */
  export function selectInTemplate(list: ParticleList, template: ParticleList, test: (local: number) => boolean) {
    const locals: number[] = [];
    for (let j = 0, jl = template.count; j < jl; ++j) if (test(j)) locals.push(j);
    return ofTemplateParticles(list, template, OrderedSet.ofSortedArray(locals));
  }

  /** `loci`'s set of copy `top`: all where whole, its partial set, else none */
  function getLocals(loci: Particle.Loci, top: number): OrderedSet<number> {
    if (OrderedSet.has(loci.indices, top)) return Interval.ofLength(getCopyCount(loci.particles, top));
    const p = loci.partial;
    if (!p) return OrderedSet.Empty;
    const k = binarySearch(p.tops, top);
    return k >= 0 ? p.locals[k] : OrderedSet.Empty;
  }

  function binarySearch(xs: Int32Array, x: number) {
    let lo = 0;
    let hi = xs.length - 1;
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      const v = xs[mid];
      if (v === x) return mid;
      if (v < x) lo = mid + 1;
      else hi = mid - 1;
    }
    return -1;
  }

  const NoTops = new Int32Array(0);

  /** The tops of `a`'s and `b`'s partial copies, ascending, distinct */
  function partialTops(a: Particle.Loci, b: Particle.Loci): number[] {
    const x = a.partial?.tops ?? NoTops;
    const y = b.partial?.tops ?? NoTops;
    const out: number[] = [];
    let i = 0;
    let j = 0;
    while (i < x.length || j < y.length) {
      const v = j >= y.length || (i < x.length && x[i] <= y[j]) ? x[i] : y[j];
      out.push(v);
      if (i < x.length && x[i] === v) ++i;
      if (j < y.length && y[j] === v) ++j;
    }
    return out;
  }

  function checkSameList(a: Particle.Loci, b: Particle.Loci) {
    if (a.particles !== b.particles) throw new Error('particle loci of different lists');
  }

  /** The instances `a` or `b` holds (of the same list) */
  export function union(a: Particle.Loci, b: Particle.Loci): Particle.Loci {
    checkSameList(a, b);
    if (!a.partial && !b.partial) return Particle.Loci(a.particles, OrderedSet.union(a.indices, b.indices));
    const tops = partialTops(a, b);
    return createLoci(
      a.particles,
      OrderedSet.union(a.indices, b.indices),
      tops,
      tops.map((t) => OrderedSet.union(getLocals(a, t), getLocals(b, t))),
    );
  }

  /** The instances `a` and `b` both hold */
  export function intersect(a: Particle.Loci, b: Particle.Loci): Particle.Loci {
    checkSameList(a, b);
    if (!a.partial && !b.partial) return Particle.Loci(a.particles, OrderedSet.intersect(a.indices, b.indices));
    const tops = partialTops(a, b);
    return createLoci(
      a.particles,
      OrderedSet.intersect(a.indices, b.indices),
      tops,
      tops.map((t) => OrderedSet.intersect(getLocals(a, t), getLocals(b, t))),
    );
  }

  /** The instances `a` holds and `b` does not */
  export function subtract(a: Particle.Loci, b: Particle.Loci): Particle.Loci {
    checkSameList(a, b);
    if (!a.partial && !b.partial) return Particle.Loci(a.particles, OrderedSet.subtract(a.indices, b.indices));
    const bTops = b.partial ? OrderedSet.ofSortedArray(b.partial.tops) : OrderedSet.Empty;
    const tops = partialTops(a, b);
    return createLoci(
      a.particles,
      OrderedSet.subtract(a.indices, OrderedSet.union(b.indices, bTops)),
      tops,
      tops.map((t) => OrderedSet.subtract(getLocals(a, t), getLocals(b, t))),
    );
  }

  /** Whether `a` and `b` hold the same partial copies (canonical: by their parts) */
  export function arePartialEqual(a: Particle.PartialCopies | undefined, b: Particle.PartialCopies | undefined) {
    if (a === b) return true;
    const n = a?.tops.length ?? 0;
    if (n !== (b?.tops.length ?? 0)) return false;
    for (let k = 0; k < n; ++k) {
      if (a!.tops[k] !== b!.tops[k] || !OrderedSet.areEqual(a!.locals[k], b!.locals[k])) return false;
    }
    return true;
  }

  /** The instances `loci` holds: a copy held whole its template's count */
  export function getLociSize(loci: Particle.Loci) {
    const { offsets } = getExpansion(loci.particles);
    const { indices, partial } = loci;
    let size = 0;
    if (!offsets) size = OrderedSet.size(indices);
    else if (OrderedSet.isInterval(indices)) {
      if (Interval.size(indices) > 0) size = offsets[Interval.end(indices)] - offsets[Interval.start(indices)];
    } else {
      for (let k = 0, kl = indices.length; k < kl; ++k) size += offsets[indices[k] + 1] - offsets[indices[k]];
    }
    if (partial) for (const l of partial.locals) size += OrderedSet.size(l);
    return size;
  }

  /** Calls `f` for each instance of `loci` by its path, in the expansion's order (`getExpandedIndex`) */
  export function forEachInstance(loci: Particle.Loci, f: (top: number, local: number) => void) {
    const { particles: list, indices, partial } = loci;
    const tops = partial?.tops ?? NoTops;
    const n = OrderedSet.size(indices);
    let i = 0;
    let p = 0;
    while (i < n || p < tops.length) {
      const w = i < n ? OrderedSet.getAt(indices, i) : Infinity;
      if (p >= tops.length || w < tops[p]) {
        for (let j = 0, jl = getCopyCount(list, w); j < jl; ++j) f(w, j);
        ++i;
      } else {
        const t = tops[p];
        OrderedSet.forEach(partial!.locals[p], (j) => f(t, j));
        ++p;
      }
    }
  }

  /**
   * `loci` (of a list with templates) as a loci of its flattened list (`flatten`): an instance its expanded index.
   * Where a list's representation draws the flattened list: what it marks
   */
  export function toFlatLoci(loci: Particle.Loci): Particle.Loci {
    const list = loci.particles;
    const { offsets } = getExpansion(list);
    if (!offsets) return loci;
    const flat = flatten(list);
    const { indices, partial } = loci;
    if (!partial && OrderedSet.isInterval(indices)) {
      if (Interval.size(indices) === 0) return Particle.Loci(flat, OrderedSet.Empty);
      return Particle.Loci(flat, Interval.ofBounds(offsets[Interval.start(indices)], offsets[Interval.end(indices)]));
    }
    const out = new Int32Array(getLociSize(loci));
    let o = 0;
    forEachInstance(loci, (top, local) => {
      out[o++] = offsets[top] + local;
    });
    return Particle.Loci(flat, OrderedSet.ofSortedArray(out));
  }

  const _path: ParticlePath = { top: 0, local: 0 };

  /**
   * The loci of `list` (with templates) holding the instances `indices` of its flattened list: a copy whose every
   * instance is held whole, else its held template particles. By the indices' runs a copy (a search a copy, not an
   * instance): what a flattened list's pick or query gives, as a path set
   */
  export function fromFlatIndices(list: ParticleList, indices: OrderedSet<number>): Particle.Loci {
    const { offsets } = getExpansion(list);
    if (!offsets) return Particle.Loci(list, indices);
    const whole: number[] = [];
    const tops: number[] = [];
    const locals: OrderedSet<number>[] = [];
    const n = OrderedSet.size(indices);
    let k = 0;
    while (k < n) {
      const v = OrderedSet.getAt(indices, k);
      if (v >= offsets[list.count]) break;
      const { top } = getPath(list, v, _path);
      const start = offsets[top];
      const end = offsets[top + 1];
      // the first index of the next copy's run
      const kEnd = Math.min(n, OrderedSet.findPredecessorIndex(indices, end));
      if (kEnd - k === end - start) whole.push(top);
      else {
        const l = new Int32Array(kEnd - k);
        for (let m = k; m < kEnd; ++m) l[m - k] = OrderedSet.getAt(indices, m) - start;
        tops.push(top);
        locals.push(OrderedSet.ofSortedArray(l));
      }
      k = kEnd;
    }
    return createLoci(list, OrderedSet.ofSortedArray(whole), tops, locals);
  }

  /** `loci` of another list (`Particle.remapLoci`): the tops past `list`'s count dropped, a copy's set clipped */
  export function remapLoci(loci: Particle.Loci, list: ParticleList): Particle.Loci {
    const count = list.count;
    const indices = OrderedSet.intersect(loci.indices, Interval.ofLength(count));
    const p = loci.partial;
    if (!p) return createLoci(list, indices);
    const tops: number[] = [];
    const locals: OrderedSet<number>[] = [];
    for (let k = 0; k < p.tops.length && p.tops[k] < count; ++k) {
      tops.push(p.tops[k]);
      locals.push(p.locals[k]);
    }
    return createLoci(list, indices, tops, locals);
  }

  const _sp = Vec3();
  const _sc = Vec3();
  const _sq = Quat();

  /** template particles a copy holds from which its part is bounded by their sphere once a set (`forEachLociSphere`) */
  const PartSphereMinCount = 16;

  /**
   * Calls `f` with a sphere of each part of `loci`: a leaf its particle's, a copy held whole its template's boundary
   * placed by the copy's pose, a copy held in part its instances' (their poses composed), or where it holds many
   * (`PartSphereMinCount`) the sphere of its set in its template (computed once a template and set, in `spheres`: a
   * template particle set in every copy, one pass over the set, not a copy each) placed; `center` reused
   */
  export function forEachLociSphere(
    loci: Particle.Loci,
    f: (center: Vec3, radius: number) => void,
    spheres: Map<ParticleList, Map<OrderedSet<number>, Sphere3D>> = new Map(),
  ) {
    const { particles: list, indices, partial } = loci;
    OrderedSet.forEach(indices, (top) => {
      const template = getTemplate(list, ParticleValues.get(list.targets, top));
      GridPositions.get(list.coordinates, top, _sp);
      if (!template) {
        f(_sp, list.radii ? ParticleValues.get(list.radii, top) : 0);
        return;
      }
      const sphere = Particle.getBoundary(template.particles).sphere;
      ParticleRotations.get(list.rotations, top, _sq);
      Vec3.transformQuat(_sc, Vec3.sub(_sc, sphere.center, ParticleTarget.getCenter(template)!), _sq);
      f(Vec3.add(_sc, _sc, _sp), sphere.radius);
    });
    if (!partial) return;
    for (let k = 0; k < partial.tops.length; ++k) {
      const top = partial.tops[k];
      const template = getTemplate(list, ParticleValues.get(list.targets, top))!;
      const t = template.particles;
      const center = ParticleTarget.getCenter(template)!;
      GridPositions.get(list.coordinates, top, _sp);
      ParticleRotations.get(list.rotations, top, _sq);
      const locals = partial.locals[k];
      // a template not read (by range: `Particle.DeferredParticles`): its set's sphere by its tiles' boxes alone
      if (OrderedSet.size(locals) >= PartSphereMinCount || Particle.getDeferred(t)) {
        let byTemplate = spheres.get(t);
        if (!byTemplate) spheres.set(t, (byTemplate = new Map()));
        let sphere = byTemplate.get(locals);
        if (!sphere)
          byTemplate.set(locals, (sphere = Particle.getBoundingSphere(Particle.Loci(t, locals), Sphere3D())));
        Vec3.transformQuat(_sc, Vec3.sub(_sc, sphere.center, center), _sq);
        f(Vec3.add(_sc, _sc, _sp), sphere.radius);
        continue;
      }
      OrderedSet.forEach(locals, (j) => {
        GridPositions.get(t.coordinates, j, _sc);
        Vec3.transformQuat(_sc, Vec3.sub(_sc, _sc, center), _sq);
        f(Vec3.add(_sc, _sc, _sp), t.radii ? ParticleValues.get(t.radii, j) : 0);
      });
    }
  }

  const _boundaryHelper = new BoundaryHelper('98');
  const _deferredHelper = new BoundaryHelper('98');
  const _tileCenter = Vec3();
  const _keptPosition = Vec3();

  /**
   * a page not kept holding at most this many of a loci's particles: each looked for in its list's stored levels (a
   * coarse copy's click), the rest by their tiles; more: by their tiles alone (a tile's box then about theirs)
   */
  const LevelLookupMax = 64;

  /** A list's stored levels' particles by page (`Particle.getStoredLevels`): what a page not kept has in memory */
  interface StoredLevelIndex {
    /** per page its entries' first (`DeferredPageBits`), the entries' particle, level and index in it */
    readonly offsets: Int32Array;
    readonly keys: Int32Array;
    readonly levels: Uint8Array;
    readonly at: Int32Array;
  }
  const storedLevelIndices = new WeakMap<ParticleList, StoredLevelIndex | null>();

  /** `list`'s stored levels' particles by page (made once a list, a pass over their keys), null where it has none */
  function getStoredLevelIndex(list: ParticleList): StoredLevelIndex | null {
    let index = storedLevelIndices.get(list);
    if (index !== undefined) return index;
    const stored = Particle.getStoredLevels(list);
    if (!stored?.length) {
      storedLevelIndices.set(list, null);
      return null;
    }
    const pages = (list.count >> DeferredPageBits) + 1;
    const offsets = new Int32Array(pages + 1);
    let total = 0;
    for (const level of stored) {
      for (let i = 0; i < level.keys.length; ++i) offsets[(level.keys[i] >> DeferredPageBits) + 1] += 1;
      total += level.keys.length;
    }
    for (let p = 0; p < pages; ++p) offsets[p + 1] += offsets[p];
    const keys = new Int32Array(total);
    const levels = new Uint8Array(total);
    const at = new Int32Array(total);
    const fill = offsets.slice(0, pages);
    for (let l = 0; l < stored.length; ++l) {
      const k = stored[l].keys;
      for (let i = 0; i < k.length; ++i) {
        const e = fill[k[i] >> DeferredPageBits]++;
        keys[e] = k[i];
        levels[e] = l;
        at[e] = i;
      }
    }
    index = { offsets, keys, levels, at };
    storedLevelIndices.set(list, index);
    return index;
  }

  /** A tile directory's tiles by their first particle (made once a directory) */
  const tileOrders = new WeakMap<object, { readonly tiles: Int32Array; readonly firsts: Int32Array }>();
  function getTileOrder(directory: { readonly ranges: Int32Array }) {
    let order = tileOrders.get(directory);
    if (order) return order;
    const { ranges } = directory;
    const all = Array.from({ length: ranges.length / 3 }, (_, t) => t).filter((t) => ranges[t * 3 + 2] > 0);
    all.sort((a, b) => ranges[a * 3 + 1] - ranges[b * 3 + 1]);
    const tiles = Int32Array.from(all);
    order = { tiles, firsts: Int32Array.from(all, (t) => ranges[t * 3 + 1]) };
    tileOrders.set(directory, order);
    return order;
  }

  /**
   * The sphere bounding particles `indices` of `list` whose particles are not read (`Particle.DeferredParticles`, a
   * pack opened by range), at once (no read begun): a particle by its position and radius where it is in memory (its
   * page kept, `DeferredParticles.getPage`: its pack's cache; or in its stored levels, `LevelLookupMax`: a coarse
   * copy's click), else by the box of its tile (its tile directory: conservative) grown by the largest radius of those
   * in it; its boundary's where it has no directory. Tighter once a page is read
   */
  export function getDeferredBoundingSphere(list: ParticleList, indices: OrderedSet<number>, out?: Sphere3D): Sphere3D {
    const helper = _deferredHelper;
    helper.reset();
    forEachDeferredSphere(list, indices, (c, r) => helper.includePositionRadius(c, r));
    helper.finishedIncludeStep();
    forEachDeferredSphere(list, indices, (c, r) => helper.radiusPositionRadius(c, r));
    return helper.getSphere(out);
  }

  /** `getDeferredBoundingSphere`'s spheres: a particle in memory its own, the others' their tiles' boxes */
  function forEachDeferredSphere(
    list: ParticleList,
    indices: OrderedSet<number>,
    f: (center: Vec3, radius: number) => void,
  ) {
    const size = OrderedSet.size(indices);
    if (size === 0) return;
    const deferred = Particle.getDeferred(list)!;
    const { radii } = list;
    const radiusOf = (i: number) => (radii ? ParticleValues.get(radii, i) : 0);
    const grid = list.coordinates as GridPositions;
    const stored = Particle.getStoredLevels(list);
    const levelIndex = stored ? getStoredLevelIndex(list) : null;
    // the loci's particles not in memory: ranges of their places in `indices` (`[a, b)` pairs, ascending)
    const missing: number[] = [];
    const miss = (a: number, b: number) => {
      const n = missing.length;
      if (n > 0 && missing[n - 1] === a) missing[n - 1] = b;
      else missing.push(a, b);
    };
    for (let a = 0; a < size; ) {
      const page = OrderedSet.getAt(indices, a) >> DeferredPageBits;
      const lo = page << DeferredPageBits;
      const b = OrderedSet.findPredecessorIndex(indices, lo + (1 << DeferredPageBits));
      const kept = deferred.getPage(page);
      if (kept) {
        const view: GridPositions = { ...grid, data: kept.positions };
        for (let k = a; k < b; ++k) {
          const i = OrderedSet.getAt(indices, k);
          GridPositions.get(view, i - lo, _keptPosition);
          f(_keptPosition, radiusOf(i));
        }
      } else if (levelIndex && b - a <= LevelLookupMax) {
        const { offsets, keys, levels, at } = levelIndex;
        const e0 = offsets[page];
        const e1 = offsets[page + 1];
        for (let k = a; k < b; ++k) {
          const i = OrderedSet.getAt(indices, k);
          let e = e0;
          while (e < e1 && keys[e] !== i) ++e;
          if (e === e1) {
            miss(k, k + 1);
            continue;
          }
          GridPositions.get(stored![levels[e]].coordinates, at[e], _keptPosition);
          f(_keptPosition, radiusOf(i));
        }
      } else {
        miss(a, b);
      }
      a = b;
    }
    if (missing.length === 0) return;
    const directory = ParticleTiles.getDirectory(list);
    if (!directory) {
      const { sphere } = Particle.getBoundary(list);
      f(sphere.center, sphere.radius);
      return;
    }
    const { ranges, boxes } = directory;
    const { tiles, firsts } = getTileOrder(directory);
    // per missing range, per tile it spans: the tile's box's center and half diagonal grown by their largest radius
    for (let m = 0; m < missing.length; m += 2) {
      const end = missing[m + 1];
      for (let a = missing[m]; a < end; ) {
        const i = OrderedSet.getAt(indices, a);
        // the tile holding particle `i`: the last whose first is at most `i`
        let lo = 0;
        let hi = firsts.length;
        while (hi - lo > 1) {
          const mid = (lo + hi) >> 1;
          if (firsts[mid] <= i) lo = mid;
          else hi = mid;
        }
        const t = tiles[lo];
        const tileEnd = ranges[t * 3 + 1] + ranges[t * 3 + 2];
        if (firsts.length === 0 || firsts[lo] > i || i >= tileEnd) {
          // in no tile: the list's boundary
          const { sphere } = Particle.getBoundary(list);
          f(sphere.center, sphere.radius);
          return;
        }
        const b = Math.min(end, OrderedSet.findPredecessorIndex(indices, tileEnd));
        let reach = 0;
        if (radii) for (let k = a; k < b; ++k) reach = Math.max(reach, radiusOf(OrderedSet.getAt(indices, k)));
        const o = t * 6;
        const dx = boxes[o + 3] - boxes[o];
        const dy = boxes[o + 4] - boxes[o + 1];
        const dz = boxes[o + 5] - boxes[o + 2];
        Vec3.set(
          _tileCenter,
          (boxes[o] + boxes[o + 3]) / 2,
          (boxes[o + 1] + boxes[o + 4]) / 2,
          (boxes[o + 2] + boxes[o + 5]) / 2,
        );
        f(_tileCenter, Math.sqrt(dx * dx + dy * dy + dz * dz) / 2 + reach);
        a = b;
      }
    }
  }

  /** The sphere bounding `loci`'s instances (`forEachLociSphere`): all of them the list's cached boundary */
  export function getLociBoundingSphere(loci: Particle.Loci, out?: Sphere3D): Sphere3D {
    if (!loci.partial && OrderedSet.size(loci.indices) === loci.particles.count) {
      const sphere = Particle.getBoundary(loci.particles).sphere;
      return out ? Sphere3D.copy(out, sphere) : Sphere3D.clone(sphere);
    }
    const helper = _boundaryHelper;
    const spheres = new Map<ParticleList, Map<OrderedSet<number>, Sphere3D>>();
    helper.reset();
    forEachLociSphere(loci, (c, r) => helper.includePositionRadius(c, r), spheres);
    helper.finishedIncludeStep();
    forEachLociSphere(loci, (c, r) => helper.radiusPositionRadius(c, r), spheres);
    return helper.getSphere(out);
  }

  /**
   * The label of `loci`: an instance by its path (`Particle 12 › Particle 3,401`), a copy held whole by its top
   * particle and its count, else the instances' count
   */
  export function getLociLabel(loci: Particle.Loci): string {
    const { particles: list, indices, partial } = loci;
    const size = getLociSize(loci);
    if (size === 0) return 'None';
    if (!partial && OrderedSet.size(indices) === 1) {
      const label = list.getParticleLabel(OrderedSet.start(indices));
      return size === 1 ? label : `${label} · ${size.toLocaleString('en')} particles`;
    }
    if (size === 1 && partial) {
      const top = partial.tops[0];
      const template = getTemplate(list, ParticleValues.get(list.targets, top))!;
      const local = OrderedSet.start(partial.locals[0]);
      return `${list.getParticleLabel(top)} › ${template.particles.getParticleLabel(local)}`;
    }
    return `${size.toLocaleString('en')} particles`;
  }
}
