/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import type { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Particle, type ParticleList, ParticleValues } from '@molstar/model/model/particles/particle-list';
import { ParticleTiles } from '@molstar/model/model/particles/particle-tiles';
import { createRangeAllocator, type RangeAllocator } from '@molstar/graphics/gl/range-allocator';
import {
  PagedCullPageBits,
  PagedFlagsWord,
  PagedPerParticleMembers,
  PagedRecordWords,
  PagedTemplatePage,
} from '@molstar/graphics/gl/compute/paged-cull';
import { TargetMemberPool } from './member-pool.js';

// A slot list's particles in pages of a fixed size, for a GPU that culls them and composes their poses (no instance
// data on the CPU): a slot list's page (a tile) its store pages, each its first slot, its values as stored (positions,
// rotations), a record (its grid, its runs, its count, its slots, its box, its extras), its runs of targets a block of
// a runs pool; a page's member and scale a particle where its runs are many or a target scales by radii given a
// particle (its extras, a block of an extras pool). A tile let go its pages free-listed: fixed size, no compaction.
// What changed (`takeChanges`) what a GPU uploads.

export const StorePageBits = PagedCullPageBits;
/** particles a store page */
export const StorePageSize = 1 << StorePageBits;
/**
 * words of a page's record: its grid's origin and step (f32; a step 0 for float positions), its runs' first and count
 * (u32; `PerParticleMembers`: its members its extras'), its count and slot base (u32), its box's min (f32), its extras'
 * first (u32, in particles; `NoExtras` for none), its box's max (f32), its flags (u32: `TemplatePage`)
 */
export const PageRecordWords = PagedRecordWords;
/** a record's flags word */
export const PageFlagsWord = PagedFlagsWord;
/**
 * a page of a template's particles (a nested list's, `ParticleNesting`): in its template's frame, drawn by its copies
 * alone (not an active page of its own)
 */
export const TemplatePage = PagedTemplatePage;
/** floats of a copy's pose (`StoreCopies.poses`): `t = p - R c` (`c` the template's center) and its scale 1, `q` */
export const CopyPoseFloats = 8;
/** a page's members a particle (the runs' count word): its extras' */
export const PerParticleMembers = PagedPerParticleMembers;
/** a page of no extras (its record's extras' first) */
export const NoExtras = 0xffffffff;
/** words of a particle's extras: its member (u32), its scale (f32) */
export const ExtraWords = 2;
/** words of a run: its first in its page, its count, its member (u32), its scale (f32; `ScaleByParticle`: its particles') */
export const RunWords = 4;
/** a run's scale read a particle (the page's scales: a target scaling by radii given a particle) */
export const ScaleByParticle = -1;
/** a member a particle of no member (left out by the placement) */
export const NoMember = 0xffffffff;
/** a page whose runs are more than its particles divided by this: its members a particle instead */
const RunsShare = 8;

/**
 * The bytes of a deferred list's page values kept in memory once given up (`releaseValues`: a culler wrote them), least
 * recently given up let go first: read again by range past it
 */
export const PageValueCacheBytes = 256 * 1024 * 1024;

/**
 * an add's scratch, kept between adds (grown, never shrunk): its runs (first, count, member, scale), its members'
 * instances by member (all 0 between adds), the members it counted
 */
let RunScratch = new Float64Array(256);
let MemberCountScratch = new Float64Array(64);
let TouchedScratch = new Int32Array(64);

/** The values of a store: one format for all its pages (the first list's) */
export interface PageStoreFormat {
  readonly positions: 'grid' | 'float32';
  readonly rotations: 'none' | 'float32' | 'snorm16';
}

/** A page's values as uploaded: views of its list's arrays (their format the store's) or converted copies */
export interface PageValues {
  /** two words a particle (a grid's) or three (float32) */
  readonly positions: Uint32Array | Float32Array;
  readonly rotations?: Float32Array | Int16Array;
}

/**
 * A template's copies (a nested list's top particles placing it): their top particles and poses (`CopyPoseFloats`
 * each: a copy's instance of template particle `x` at `t + R_c (t_x + s_x R_x c_x)`, its rotation `q_c q_x`)
 */
export interface StoreCopies {
  readonly tops: Int32Array;
  readonly poses: Float32Array;
  /**
   * a template's coarse level of its own (`ParticleNesting.getTemplateLevels`: from 1; 0 or none the template), its
   * copies the template's: a copy drawn at one level by its distance (`PagedCullView.coarse`); `levels` the
   * template's coarse levels, `center` its center (in its frame) and `spacing` its particles' (what picks a level)
   */
  readonly level?: number;
  readonly levels?: number;
  readonly center?: ArrayLike<number>;
  readonly spacing?: number;
  /** its template's particles' centers' distance from `center` at most: a copy out of view culled whole */
  readonly radius?: number;
}

/** A slot list's page in the store: its store pages, and its members' instances (what a let go takes back) */
export interface StoreBlock {
  readonly key: number;
  /** a template's particles (`TemplatePage`): drawn by its copies (`ParticlePageStore.getCopies`) */
  readonly template?: boolean;
  /** its first slot */
  readonly start: number;
  /** its particles (its list's) */
  readonly count: number;
  readonly pages: Int32Array;
  /** member, instances: two a member */
  readonly members: Int32Array;
  /** its runs' largest scale (`ParticleList.scale`, `scales`: a coarse sample's particles drawn larger) */
  readonly maxScale: number;
}

/** What changed since the last `takeChanges`: what a GPU uploads */
export interface PageStoreChanges {
  /** the pages written (their values and records) */
  readonly written: number[];
  /** the pages let go (their records cleared) */
  readonly released: number[];
  /** the templates' copies set or let go (by their blocks' keys) */
  readonly copies: number[];
  /** the runs pool's ranges written, `[start, end)` in runs (two numbers a range) */
  readonly runs: number[];
  /** the extras pool's ranges written, `[start, end)` in particles (two numbers a range) */
  readonly extras: number[];
  /** the members whose instances came or went */
  readonly members: Set<number>;
  /** the members freed (their last instance gone) */
  readonly freed: number[];
  /** the pages' capacity, the runs pool's, the members pool's grown (their buffers larger) */
  grown: boolean;
}

export interface PageStoreOptions {
  /** the pages there is room for at first */
  readonly pageCapacity?: number;
  /** a deferred list's page values kept once given up (`PageValueCacheBytes` by default) */
  readonly valueCacheBytes?: number;
  /** a deferred list's read failing: told once a list (by default a console warning), its page read again later */
  readonly onReadError?: (list: ParticleList, error: unknown) => void;
  /** the clock a failed read's retry waits by (ms, `performance.now` by default) */
  readonly now?: () => number;
}

/** a failed page read's first wait before it is read again (ms), doubled each failure up to `PageReadRetryMaxMs` */
export const PageReadRetryMs = 250;
export const PageReadRetryMaxMs = 30_000;

function valuesBytes(v: PageValues) {
  return v.positions.byteLength + (v.rotations?.byteLength ?? 0);
}

function emptyChanges(): PageStoreChanges {
  return { written: [], released: [], runs: [], extras: [], copies: [], members: new Set(), freed: [], grown: false };
}

/** A tile's placement of its targets: a target's center (where a particle places it) and whether it scales by radius */
export interface StorePlacement {
  centerOf(target: number): Vec3;
  scaleByRadius?(target: number): boolean;
  /** whether the store holds `target`'s particles (all without): the others left out */
  accepts?(target: number): boolean;
}

export class ParticlePageStore {
  /** the store's values' format: its first list's */
  format: PageStoreFormat | undefined;
  /** the targets' members */
  readonly pool = new TargetMemberPool();
  /** the largest scale of a list added (`ParticleList.scale`, 1 at least): what pads the object's bounds */
  maxScale = 1;
  /** the page records (`PageRecordWords` words a page), as words and as floats */
  records: Uint32Array;
  recordFloats: Float32Array;
  /** the runs pool (`RunWords` words a run), as words and as floats */
  runs: Uint32Array;
  runFloats: Float32Array;
  /** the extras pool (`ExtraWords` words a particle: member, scale), as words and as floats */
  extras: Uint32Array;
  extraFloats: Float32Array;
  /** per page its values as uploaded (none: a free page, or a deferred list's page not read: `requestValues`) */
  values: (PageValues | undefined)[] = [];
  /** per page of a deferred list (`Particle.DeferredParticles`): its list and range, its read in flight */
  private readonly deferredPages = new Map<
    number,
    {
      readonly list: ParticleList;
      readonly lo: number;
      readonly hi: number;
      reading?: Promise<void>;
      /** its reads failed in a row, and when it is read again (`now`): not every frame the cull wants it */
      failures?: number;
      retryAt?: number;
    }
  >();
  private readonly onReadError: (list: ParticleList, error: unknown) => void;
  private readonly now: () => number;
  /** the lists whose read failure was told (`onReadError`) */
  private readonly toldFailures = new WeakSet<ParticleList>();
  /** a deferred list's page values given up (`releaseValues`), by when (the oldest first), their bytes */
  private readonly valueCache = new Map<number, PageValues>();
  private readonly valueCacheLimit: number;
  /** the bytes of the values given up and kept (`releaseValues`) */
  cachedValueBytes = 0;
  /** the bytes of a deferred list's particles read (`requestValues`: its positions and rotations as read) */
  readValueBytes = 0;
  /** the pages there are (used and free): `[0, pageCount)` */
  pageCount = 0;

  private readonly runAlloc: RangeAllocator;
  private readonly extraAlloc: RangeAllocator;
  /** per page its runs' first in the pool (where its members are runs), its extras' count (0: none) */
  private readonly pageRuns: number[] = [];
  private readonly pageExtras: number[] = [];
  private readonly freePages: number[] = [];
  private readonly blocks = new Map<number, StoreBlock>();
  private readonly copies = new Map<number, StoreCopies>();
  /** per template block its members' instances counted this many times (its copies'; 1 as added, before any) */
  private readonly copyFactors = new Map<number, number>();
  private changes = emptyChanges();

  constructor(options: PageStoreOptions = {}) {
    const pages = Math.max(4, options.pageCapacity ?? 16);
    this.valueCacheLimit = options.valueCacheBytes ?? PageValueCacheBytes;
    this.onReadError =
      options.onReadError ??
      ((list, error) => console.warn(`particle list '${list.label}': a page's read failed, read again later`, error));
    this.now = options.now ?? (() => performance.now());
    this.records = new Uint32Array(pages * PageRecordWords);
    this.recordFloats = new Float32Array(this.records.buffer);
    this.runs = new Uint32Array(pages * 64 * RunWords);
    this.runFloats = new Float32Array(this.runs.buffer);
    this.runAlloc = createRangeAllocator(pages * 64, 1);
    this.extras = new Uint32Array(1024 * ExtraWords);
    this.extraFloats = new Float32Array(this.extras.buffer);
    this.extraAlloc = createRangeAllocator(1024, 1);
  }

  get extraCapacity() {
    return this.extras.length / ExtraWords;
  }

  get pageCapacity() {
    return this.records.length / PageRecordWords;
  }

  get runCapacity() {
    return this.runs.length / RunWords;
  }

  getBlock(key: number) {
    return this.blocks.get(key);
  }

  get blockList(): Iterable<StoreBlock> {
    return this.blocks.values();
  }

  /** What changed since the last look, cleared */
  takeChanges(): PageStoreChanges {
    const c = this.changes;
    this.changes = emptyChanges();
    return c;
  }

  /**
   * The flat `list`'s particles as a block of `key` (a slot list's page: its slots from `start`), its targets' members
   * made where they have none (their centers `placement`'s); the particles of targets the placement does not accept in
   * no run (their pages' slots, drawn by none). `template`: a template's particles (a nested list's leaf list,
   * `ParticleNesting.getLeafList`: its targets in the flattened list's ids), its pages `TemplatePage`, drawn by its
   * copies (`setCopies`)
   */
  add(
    key: number,
    list: ParticleList,
    start: number,
    placement: StorePlacement,
    box?: Box3D,
    template = false,
  ): StoreBlock {
    if (this.blocks.has(key)) throw new Error(`particle page store: a block of ${key} already`);
    const format = (this.format ??= getFormat(list));
    const n = list.count;
    const pageCount = Math.ceil(n / StorePageSize);
    const accepts = placement.accepts;
    const radii = list.radii;
    const radiiByParticle = radii && !ParticleValues.isRuns(radii) ? (radii as Float32Array) : undefined;
    // a coarse sample's particles drawn larger (`ParticleList.scale`, or a run's own: `scales`): its runs' scales by it
    const listScale = list.scale !== undefined && list.scale > 0 ? list.scale : 1;
    const scales = list.scales;
    if (!scales && listScale > this.maxScale) this.maxScale = listScale;
    let blockScale = scales ? 0 : listScale;

    // the runs by page: first (in the page), count, member, scale; split at the pages' ends. In particle order, so a
    // page's runs consecutive (`pageRunStarts`) in one scratch; a member's instances counted in an array by member, not
    // a map: a run's cost (the yeast drawn whole by its tiles: 2.29 M runs)
    const pageRunStarts = new Int32Array(pageCount + 1);
    const pool = this.pool;
    let runScratch = RunScratch;
    let runN = 0;
    let memberCounts = MemberCountScratch;
    let touched = TouchedScratch;
    let touchedN = 0;
    let byParticle = false;
    let lastPage = 0;
    // a target's member and its scaling by its radius once a target (`ofTarget`: twice the member, one more where
    // scaled; -1 not accepted), not once a run: a tile-ordered template's runs a target a tile (the yeast's 2.3 M of
    // 5,639 targets: about 200 ms of its add in the callbacks a run)
    const ofTarget = new Map<number, number>();
    const visit = (s: number, e: number, target: number, runScale: number) => {
      let o = ofTarget.get(target);
      if (o === undefined) {
        const scaledTarget = !!radii && !!placement.scaleByRadius?.(target);
        o = accepts && !accepts(target) ? -1 : pool.acquire(target, placement.centerOf) * 2 + (scaledTarget ? 1 : 0);
        ofTarget.set(target, o);
      }
      if (o < 0) return;
      const m = o >> 1;
      if (m >= memberCounts.length) {
        const grown = new Float64Array(Math.max(m + 1, memberCounts.length * 2));
        grown.set(memberCounts);
        memberCounts = MemberCountScratch = grown;
      }
      if (memberCounts[m] === 0) {
        if (touchedN === touched.length) {
          const grown = new Int32Array(Math.max(64, touchedN * 2));
          grown.set(touched);
          touched = TouchedScratch = grown;
        }
        touched[touchedN++] = m;
      }
      memberCounts[m] += e - s;
      const scaled = (o & 1) === 1;
      if (runScale > this.maxScale) this.maxScale = runScale;
      if (runScale > blockScale) blockScale = runScale;
      let scale = runScale;
      if (scaled && radiiByParticle) {
        scale = ScaleByParticle;
        byParticle = true;
      } else if (scaled) {
        scale = ParticleValues.get(radii!, s) * runScale;
      }
      for (let a = s; a < e; ) {
        const k = a >> StorePageBits;
        const b = Math.min(e, (k + 1) << StorePageBits);
        while (lastPage < k) pageRunStarts[++lastPage] = runN;
        if ((runN + 1) * 4 > runScratch.length) {
          const grown = new Float64Array(Math.max(256, runScratch.length * 2));
          grown.set(runScratch);
          runScratch = RunScratch = grown;
        }
        const q = runN * 4;
        runScratch[q] = a - (k << StorePageBits);
        runScratch[q + 1] = b - a;
        runScratch[q + 2] = m;
        runScratch[q + 3] = scale;
        runN += 1;
        a = b;
      }
    };
    const targets = list.targets;
    if (ParticleValues.isRuns(targets)) {
      // the runs' arrays walked (no callback a run); a run's scale by the same runs where its scales are (a level's:
      // its targets' and its scales' runs one), else by its first particle
      const { offsets, values } = targets;
      const runScales =
        scales && ParticleValues.isRuns(scales) && scales.offsets === offsets ? scales.values : undefined;
      for (let r = 0, rl = values.length; r < rl; ++r) {
        const s = offsets[r];
        const e = offsets[r + 1];
        if (e <= s) continue;
        visit(s, e, values[r], runScales ? runScales[r] : scales ? ParticleValues.get(scales, s) : listScale);
      }
    } else {
      ParticleValues.forEachRun(targets, n, (s, e, target) =>
        visit(s, e, target, scales ? ParticleValues.get(scales, s) : listScale),
      );
    }
    while (lastPage < pageCount) pageRunStarts[++lastPage] = runN;

    // a deferred list's pages (a pack opened by range): their boxes its tiles' (the directory's, the union of those a
    // page holds: conservative), their values read when asked (`requestValues`)
    const deferred = Particle.getDeferred(list);
    const directory = deferred ? ParticleTiles.getDirectory(list) : undefined;
    if (deferred && !box && !directory) throw new Error('particle page store: a deferred list without its tiles');
    let tile = 0;
    const pages = new Int32Array(pageCount);
    for (let k = 0; k < pageCount; ++k) {
      const p = this.allocatePage();
      pages[k] = p;
      const lo = k << StorePageBits;
      const hi = Math.min(n, lo + StorePageSize);
      const count = hi - lo;
      const run0 = pageRunStarts[k] * 4;
      const runEnd = pageRunStarts[k + 1] * 4;
      const runCount = (runEnd - run0) / 4;
      const o = p * PageRecordWords;
      const r = this.records;
      const f = this.recordFloats;
      const values = deferred ? undefined : getPageValues(list, format, lo, hi);
      writeGrid(list, f, o);
      r[o + 6] = count;
      r[o + 7] = start + lo;
      r[o + PageFlagsWord] = template ? TemplatePage : 0;
      // the list's box where given (a tile's: its tiling's), each of its pages' (larger than a page's own: what predicts
      // its members' draws, never fewer); else its particles' (a pass over the page's words)
      if (box) writeGivenBox(box, f, o);
      else if (values) writeBox(values.positions, format, f, o, count);
      else tile = writeTilesBox(directory!, tile, lo, hi, f, o);
      this.values[p] = values;
      if (deferred) this.deferredPages.set(p, { list, lo, hi });
      const perParticle = runCount > count / RunsShare;
      // a member and a scale a particle: its members a particle (a flat list's targets alternating), or a run scaled
      // by its particles' radii
      let pageScaled = false;
      for (let i = run0 + 3; i < runEnd && byParticle; i += 4) if (runScratch[i] === ScaleByParticle) pageScaled = true;
      if (perParticle || pageScaled) {
        const e = this.allocateExtras(count);
        const w = this.extras;
        const wf = this.extraFloats;
        for (let j = 0; j < count; ++j) {
          w[(e + j) * ExtraWords] = NoMember;
          wf[(e + j) * ExtraWords + 1] = 1;
        }
        for (let i = run0; i < runEnd; i += 4) {
          const first = runScratch[i];
          const end = first + runScratch[i + 1];
          const m = runScratch[i + 2];
          const scale = runScratch[i + 3];
          for (let j = first; j < end; ++j) {
            w[(e + j) * ExtraWords] = m;
            wf[(e + j) * ExtraWords + 1] =
              scale === ScaleByParticle
                ? radiiByParticle![lo + j] * (scales ? ParticleValues.get(scales, lo + j) : listScale)
                : scale;
          }
        }
        r[o + 11] = e;
        this.pageExtras[p] = count;
        this.changes.extras.push(e, e + count);
      } else {
        r[o + 11] = NoExtras;
        this.pageExtras[p] = 0;
      }
      if (perParticle) {
        r[o + 4] = 0;
        r[o + 5] = PerParticleMembers;
        this.pageRuns[p] = -1;
      } else {
        const first = this.allocateRuns(runCount);
        const w = this.runs;
        const wf = this.runFloats;
        for (let i = 0; i < runCount; ++i) {
          const q = (first + i) * RunWords;
          const u = run0 + i * 4;
          w[q] = runScratch[u];
          w[q + 1] = runScratch[u + 1];
          w[q + 2] = runScratch[u + 2];
          wf[q + 3] = runScratch[u + 3];
        }
        r[o + 4] = first;
        r[o + 5] = runCount;
        this.pageRuns[p] = first;
        if (runCount > 0) this.changes.runs.push(first, first + runCount);
      }
      this.changes.written.push(p);
    }

    const members = new Int32Array(touchedN * 2);
    for (let i = 0; i < touchedN; ++i) {
      const m = touched[i];
      const c = memberCounts[m];
      memberCounts[m] = 0;
      members[i * 2] = m;
      members[i * 2 + 1] = c;
      pool.counts[m] += c;
      this.changes.members.add(m);
    }
    const block: StoreBlock = template
      ? { key, template: true, start, count: n, pages, members, maxScale: blockScale }
      : { key, start, count: n, pages, members, maxScale: blockScale };
    this.blocks.set(key, block);
    return block;
  }

  /** The block of `key` let go: its pages and runs free, its members' counts less, a member of none freed */
  /**
   * Page `p`'s values read where they are not (a deferred list's page: `Particle.DeferredParticles`), once while in
   * flight: the read (its values set when it arrives, where the page is still that list's; set at once where given up
   * and kept: `releaseValues`); undefined where the page has its values or is not deferred, or its read failed and is
   * waited for (`PageReadRetryMs`, doubled a failure: told once a list, `onReadError`)
   */
  requestValues(p: number): Promise<void> | undefined {
    const d = this.deferredPages.get(p);
    if (!d || this.values[p]) return undefined;
    const kept = this.valueCache.get(p);
    if (kept) {
      this.uncache(p, kept);
      this.values[p] = kept;
      return Promise.resolve();
    }
    if (d.reading) return d.reading;
    if (d.retryAt !== undefined && this.now() < d.retryAt) return undefined;
    const format = this.format!;
    const deferred = Particle.getDeferred(d.list)!;
    const grid = d.list.coordinates;
    d.reading = deferred.read(d.lo, d.hi - d.lo).then(
      (read) => {
        d.reading = undefined;
        d.failures = d.retryAt = undefined;
        if (this.deferredPages.get(p) !== d) return;
        this.readValueBytes += read.positions.byteLength + (read.rotations?.byteLength ?? 0);
        // the read's words on the list's grid, as `getPageValues` takes a list's
        const page = {
          coordinates: GridPositions.is(grid) ? { ...grid, data: read.positions } : read.positions,
          rotations: read.rotations,
        } as ParticleList;
        this.values[p] = getPageValues(page, format, 0, d.hi - d.lo);
      },
      (error) => {
        // read again when asked, after a wait
        d.reading = undefined;
        d.failures = (d.failures ?? 0) + 1;
        d.retryAt = this.now() + Math.min(PageReadRetryMaxMs, PageReadRetryMs * 2 ** (d.failures - 1));
        if (this.deferredPages.get(p) !== d || this.toldFailures.has(d.list)) return;
        this.toldFailures.add(d.list);
        this.onReadError(d.list, error);
      },
    );
    return d.reading;
  }

  /**
   * Page `p`'s values given up where they are a deferred list's (a culler wrote them: its copy on the GPU), kept within
   * the cache's bytes (the least recently given up let go first): read again when asked (`requestValues`)
   */
  releaseValues(p: number) {
    const values = this.values[p];
    if (!values || !this.deferredPages.has(p)) return;
    this.values[p] = undefined;
    this.valueCache.set(p, values);
    this.cachedValueBytes += valuesBytes(values);
    for (const [q, v] of this.valueCache) {
      if (this.cachedValueBytes <= this.valueCacheLimit) break;
      this.uncache(q, v);
    }
  }

  private uncache(p: number, values: PageValues) {
    this.valueCache.delete(p);
    this.cachedValueBytes -= valuesBytes(values);
  }

  remove(key: number) {
    const block = this.blocks.get(key);
    if (!block) return;
    this.blocks.delete(key);
    if (this.copies.delete(key)) this.changes.copies.push(key);
    const factor = this.copyFactors.get(key) ?? 1;
    this.copyFactors.delete(key);
    for (const p of block.pages) {
      const o = p * PageRecordWords;
      const runFirst = this.pageRuns[p];
      if (runFirst >= 0 && this.records[o + 5] > 0) this.runAlloc.free(runFirst, this.records[o + 5]);
      if (this.pageExtras[p] > 0) this.extraAlloc.free(this.records[o + 11], this.pageExtras[p]);
      this.records.fill(0, o, o + PageRecordWords);
      this.values[p] = undefined;
      this.deferredPages.delete(p);
      const kept = this.valueCache.get(p);
      if (kept) this.uncache(p, kept);
      this.freePages.push(p);
      this.changes.released.push(p);
    }
    const { members } = block;
    for (let i = 0; i < members.length; i += 2) {
      const m = members[i];
      this.pool.counts[m] -= members[i + 1] * factor;
      this.changes.members.add(m);
      if (this.pool.counts[m] === 0) {
        this.pool.release(m);
        this.changes.freed.push(m);
      }
    }
  }

  /**
   * The copies of template block `key` (`StoreCopies`): what draws it, replacing any it had; its members' instances
   * its particles' times the copies (a member of none of them freed)
   */
  setCopies(key: number, copies: StoreCopies) {
    const block = this.blocks.get(key);
    if (!block?.template) throw new Error(`particle page store: no template block of ${key}`);
    if (copies.poses.length !== copies.tops.length * CopyPoseFloats)
      throw new Error('particle page store: a copy pose a top particle');
    // a template of no copies is none: not placed (`addNested` makes a template's block where a copy places it)
    if (copies.tops.length === 0) throw new Error('particle page store: a template of no copies');
    this.copies.set(key, copies);
    this.changes.copies.push(key);
    const before = this.copyFactors.get(key) ?? 1;
    const after = copies.tops.length;
    this.copyFactors.set(key, after);
    if (before === after) return;
    const { members } = block;
    for (let i = 0; i < members.length; i += 2) {
      const m = members[i];
      this.pool.counts[m] += members[i + 1] * (after - before);
      this.changes.members.add(m);
      if (this.pool.counts[m] === 0) {
        this.pool.release(m);
        this.changes.freed.push(m);
      }
    }
  }

  /** the times template block `key`'s members' instances are counted (its copies; 1 before any) */
  getCopyFactor(key: number) {
    return this.copyFactors.get(key) ?? 1;
  }

  getCopies(key: number): StoreCopies | undefined {
    return this.copies.get(key);
  }

  /** the templates' blocks' keys with copies */
  get copyKeys(): Iterable<number> {
    return this.copies.keys();
  }

  /** the templates with copies: their pages and copies (`PagedCullSource.templates`), by their blocks' keys' order */
  get templates(): { readonly key: number; readonly pages: Int32Array; readonly copies: StoreCopies }[] {
    const out: { key: number; pages: Int32Array; copies: StoreCopies }[] = [];
    for (const [key, copies] of this.copies) {
      const block = this.blocks.get(key);
      if (block) out.push({ key, pages: block.pages, copies });
    }
    return out;
  }

  /**
   * The poses of template block `key`'s instances in each of its copies (`compose`'s composed with each copy's: `t_c +
   * R_c t_x`, the scale `s_x`, `q_c q_x`), copy by copy, a copy's by its template's order; their copies' top particles,
   * slots (the template particle's) and members: the CPU's expansion (tests, a CPU reference of the GPU's)
   */
  composeCopies(key: number) {
    const block = this.blocks.get(key)!;
    const copies = this.copies.get(key);
    const base = this.compose(block);
    const n = base.slots.length;
    const k = copies ? copies.tops.length : 0;
    const poses = new Float32Array(n * k * 8);
    const tops = new Int32Array(n * k);
    const slots = new Int32Array(n * k);
    const members = new Int32Array(n * k);
    for (let c = 0; c < k; ++c) {
      const cp = copies!.poses;
      const o = c * CopyPoseFloats;
      const tx = cp[o];
      const ty = cp[o + 1];
      const tz = cp[o + 2];
      const qx = cp[o + 4];
      const qy = cp[o + 5];
      const qz = cp[o + 6];
      const qw = cp[o + 7];
      for (let i = 0; i < n; ++i) {
        const b = i * 8;
        const out = (c * n + i) * 8;
        // R_c t_x: t + 2 w (q x t) + 2 q x (q x t)
        const vx = base.poses[b];
        const vy = base.poses[b + 1];
        const vz = base.poses[b + 2];
        const ax = 2 * (qy * vz - qz * vy);
        const ay = 2 * (qz * vx - qx * vz);
        const az = 2 * (qx * vy - qy * vx);
        poses[out] = tx + vx + qw * ax + (qy * az - qz * ay);
        poses[out + 1] = ty + vy + qw * ay + (qz * ax - qx * az);
        poses[out + 2] = tz + vz + qw * az + (qx * ay - qy * ax);
        poses[out + 3] = base.poses[b + 3];
        // q_c q_x
        const bx = base.poses[b + 4];
        const by = base.poses[b + 5];
        const bz = base.poses[b + 6];
        const bw = base.poses[b + 7];
        poses[out + 4] = qw * bx + qx * bw + qy * bz - qz * by;
        poses[out + 5] = qw * by - qx * bz + qy * bw + qz * bx;
        poses[out + 6] = qw * bz + qx * by - qy * bx + qz * bw;
        poses[out + 7] = qw * bw - qx * bx - qy * by - qz * bz;
        tops[c * n + i] = copies!.tops[c];
        slots[c * n + i] = base.slots[i];
        members[c * n + i] = base.members[i];
      }
    }
    return { tops, slots, members, poses };
  }

  /**
   * The poses of `block`'s instances (`t = p - s R c` and `s`, the unit quaternion: eight floats each) with their
   * slots and members, by page and particle (those of no member left out): the CPU's composition (a CPU cull, tests)
   */
  compose(block: StoreBlock) {
    const format = this.format!;
    let total = 0;
    for (let i = 1; i < block.members.length; i += 2) total += block.members[i];
    const slots = new Int32Array(total);
    const members = new Int32Array(total);
    const poses = new Float32Array(total * 8);
    const p = Vec3();
    const c = Vec3();
    let o = 0;
    for (const page of block.pages) {
      const r = page * PageRecordWords;
      const f = this.recordFloats;
      const values = this.values[page]!;
      const count = this.records[r + 6];
      const slotBase = this.records[r + 7];
      const extras = this.records[r + 11];
      const each = (j: number, m: number, scale: number) => {
        if (m === NoMember) return;
        readPosition(values.positions, format, f, r, j, p);
        let s = scale === ScaleByParticle ? this.extraFloats[(extras + j) * ExtraWords + 1] : scale;
        if (!(s > 0)) s = 1;
        const q = readRotation(values.rotations, format, j);
        const centers = this.pool.centers;
        Vec3.set(c, centers[m * 3], centers[m * 3 + 1], centers[m * 3 + 2]);
        writePose(poses, o * 8, p, s, q, c);
        slots[o] = slotBase + j;
        members[o] = m;
        o += 1;
      };
      if (this.records[r + 5] === PerParticleMembers) {
        for (let j = 0; j < count; ++j) {
          each(j, this.extras[(extras + j) * ExtraWords], this.extraFloats[(extras + j) * ExtraWords + 1]);
        }
      } else {
        const first = this.records[r + 4];
        for (let k = 0; k < this.records[r + 5]; ++k) {
          const q = (first + k) * RunWords;
          const s = this.runs[q];
          const e = s + this.runs[q + 1];
          for (let j = s; j < e; ++j) each(j, this.runs[q + 2], this.runFloats[q + 3]);
        }
      }
    }
    return { slots, members, poses };
  }

  private allocatePage() {
    const free = this.freePages.pop();
    if (free !== undefined) return free;
    const p = this.pageCount++;
    if (p >= this.pageCapacity) {
      const records = new Uint32Array(this.pageCapacity * 2 * PageRecordWords);
      records.set(this.records);
      this.records = records;
      this.recordFloats = new Float32Array(records.buffer);
      this.changes.grown = true;
    }
    return p;
  }

  private allocateExtras(count: number) {
    let first = this.extraAlloc.allocate(count);
    if (first < 0) {
      const capacity = Math.max(this.extraCapacity * 2, this.extraCapacity + count);
      this.extraAlloc.grow(capacity);
      const extras = new Uint32Array(capacity * ExtraWords);
      extras.set(this.extras);
      this.extras = extras;
      this.extraFloats = new Float32Array(extras.buffer);
      this.changes.grown = true;
      first = this.extraAlloc.allocate(count);
    }
    return first;
  }

  private allocateRuns(count: number) {
    if (count === 0) return 0;
    let first = this.runAlloc.allocate(count);
    if (first < 0) {
      const capacity = Math.max(this.runCapacity * 2, this.runCapacity + count);
      this.runAlloc.grow(capacity);
      const runs = new Uint32Array(capacity * RunWords);
      runs.set(this.runs);
      this.runs = runs;
      this.runFloats = new Float32Array(runs.buffer);
      this.changes.grown = true;
      first = this.runAlloc.allocate(count);
    }
    return first;
  }
}

function getFormat(list: ParticleList): PageStoreFormat {
  return {
    positions: GridPositions.is(list.coordinates) ? 'grid' : 'float32',
    rotations: !list.rotations ? 'none' : list.rotations instanceof Int16Array ? 'snorm16' : 'float32',
  };
}

/**
 * The format of a store of all of `lists` (a list with templates': its top list's and its templates', set before its
 * first is added): the positions the largest list's (the others converted), the rotations the first list's with
 * rotations (a list without them its identities)
 */
export function getStoreFormatOf(lists: readonly ParticleList[]): PageStoreFormat {
  let largest = lists[0];
  for (const l of lists) if (l.count > largest.count) largest = l;
  const rotated = lists.find((l) => l.rotations);
  return {
    positions: getFormat(largest).positions,
    rotations: rotated ? getFormat(rotated).rotations : 'none',
  };
}

const Snorm = 32767;

/** particles `[lo, hi)`'s values in `format`: views of the list's where its format is the store's, else converted */
function getPageValues(list: ParticleList, format: PageStoreFormat, lo: number, hi: number): PageValues {
  const { coordinates, rotations } = list;
  let positions: Uint32Array | Float32Array;
  if (format.positions === 'grid') {
    // a list of float positions on a grid of its own (its record's origin and step)
    const grid = GridPositions.is(coordinates) ? coordinates : getOwnGrid(list);
    positions = grid.data.subarray(lo * 2, hi * 2);
  } else if (GridPositions.is(coordinates)) {
    positions = new Float32Array((hi - lo) * 3);
    const v = Vec3();
    for (let i = lo; i < hi; ++i) {
      GridPositions.get(coordinates, i, v);
      positions.set(v, (i - lo) * 3);
    }
  } else {
    positions = coordinates.subarray(lo * 3, hi * 3);
  }
  let rot: Float32Array | Int16Array | undefined;
  if (format.rotations === 'snorm16') {
    rot = rotations
      ? rotations instanceof Int16Array
        ? rotations.subarray(lo * 4, hi * 4)
        : Int16Array.from(rotations.subarray(lo * 4, hi * 4), (x) => Math.round(Math.max(-1, Math.min(1, x)) * Snorm))
      : identityRotations(hi - lo, true);
  } else if (format.rotations === 'float32') {
    rot = rotations
      ? rotations instanceof Float32Array
        ? rotations.subarray(lo * 4, hi * 4)
        : Float32Array.from(rotations.subarray(lo * 4, hi * 4), (x) => x / Snorm)
      : identityRotations(hi - lo, false);
  } else if (rotations) {
    throw new Error('particle page store: a list with rotations in a store without');
  }
  return { positions, rotations: rot };
}

/** a float list's positions on a grid of its own (made once a list) */
const ownGrids = new WeakMap<ParticleList, GridPositions>();
function getOwnGrid(list: ParticleList) {
  let grid = ownGrids.get(list);
  if (!grid) ownGrids.set(list, (grid = GridPositions.fromFloat(list.coordinates as Float32Array)));
  return grid;
}

function identityRotations(n: number, snorm: boolean) {
  const out = snorm ? new Int16Array(n * 4) : new Float32Array(n * 4);
  for (let i = 0; i < n; ++i) out[i * 4 + 3] = snorm ? Snorm : 1;
  return out;
}

/** the list's grid's origin and step into record `o` (a step 0 for float positions) */
function writeGrid(list: ParticleList, f: Float32Array, o: number) {
  const coordinates = list.coordinates;
  const grid = GridPositions.is(coordinates) ? coordinates : ownGrids.get(list);
  f[o] = grid ? grid.origin[0] : 0;
  f[o + 1] = grid ? grid.origin[1] : 0;
  f[o + 2] = grid ? grid.origin[2] : 0;
  f[o + 3] = grid ? grid.step : 0;
}

/** the extent of a page's positions into its record's box (a grid's steps' extent, then placed) */
/** `box` as a record's box (at `o`) */
/**
 * The box of a deferred list's page `[lo, hi)` (its particles not read): its tiles' centers' boxes' union (the
 * directory's, those whose particles it holds: conservative); `tile` the first of them to look at, the next page's
 * returned
 */
function writeTilesBox(
  directory: NonNullable<ReturnType<typeof ParticleTiles.getDirectory>>,
  tile: number,
  lo: number,
  hi: number,
  f: Float32Array,
  o: number,
) {
  const { ranges, boxes } = directory;
  const tiles = ranges.length / 3;
  let t = tile;
  while (t < tiles && ranges[t * 3 + 1] + ranges[t * 3 + 2] <= lo) t += 1;
  const next = t;
  const box = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
  for (; t < tiles && ranges[t * 3 + 1] < hi; ++t) {
    for (let a = 0; a < 3; ++a) {
      box[a] = Math.min(box[a], boxes[t * 6 + a]);
      box[3 + a] = Math.max(box[3 + a], boxes[t * 6 + 3 + a]);
    }
  }
  for (let a = 0; a < 3; ++a) {
    f[o + 8 + a] = box[a];
    f[o + 12 + a] = box[3 + a];
  }
  return next;
}

function writeGivenBox(box: Box3D, f: Float32Array, o: number) {
  f[o + 8] = box.min[0];
  f[o + 9] = box.min[1];
  f[o + 10] = box.min[2];
  f[o + 12] = box.max[0];
  f[o + 13] = box.max[1];
  f[o + 14] = box.max[2];
}

function writeBox(
  positions: Uint32Array | Float32Array,
  format: PageStoreFormat,
  f: Float32Array,
  o: number,
  n: number,
) {
  const min = [Infinity, Infinity, Infinity];
  const max = [-Infinity, -Infinity, -Infinity];
  if (format.positions === 'grid') {
    const w = positions as Uint32Array;
    let x0 = Infinity,
      y0 = Infinity,
      z0 = Infinity,
      x1 = -Infinity,
      y1 = -Infinity,
      z1 = -Infinity;
    for (let j = 0; j < n; ++j) {
      const w0 = w[j * 2];
      const w1 = w[j * 2 + 1];
      const x = w0 & 0x1fffff;
      const y = (w0 >>> 21) | ((w1 & 0x3ff) << 11);
      const z = (w1 >>> 10) & 0x1fffff;
      if (x < x0) x0 = x;
      if (y < y0) y0 = y;
      if (z < z0) z0 = z;
      if (x > x1) x1 = x;
      if (y > y1) y1 = y;
      if (z > z1) z1 = z;
    }
    const step = f[o + 3];
    min[0] = f[o] + step * x0;
    min[1] = f[o + 1] + step * y0;
    min[2] = f[o + 2] + step * z0;
    max[0] = f[o] + step * x1;
    max[1] = f[o + 1] + step * y1;
    max[2] = f[o + 2] + step * z1;
  } else {
    const v = positions as Float32Array;
    for (let j = 0; j < n; ++j) {
      for (let a = 0; a < 3; ++a) {
        const x = v[j * 3 + a];
        if (x < min[a]) min[a] = x;
        if (x > max[a]) max[a] = x;
      }
    }
  }
  f[o + 8] = min[0];
  f[o + 9] = min[1];
  f[o + 10] = min[2];
  f[o + 12] = max[0];
  f[o + 13] = max[1];
  f[o + 14] = max[2];
}

function readPosition(
  positions: Uint32Array | Float32Array,
  format: PageStoreFormat,
  f: Float32Array,
  o: number,
  j: number,
  out: Vec3,
) {
  if (format.positions === 'grid') {
    const w = positions as Uint32Array;
    const w0 = w[j * 2];
    const w1 = w[j * 2 + 1];
    const step = f[o + 3];
    out[0] = f[o] + step * (w0 & 0x1fffff);
    out[1] = f[o + 1] + step * ((w0 >>> 21) | ((w1 & 0x3ff) << 11));
    out[2] = f[o + 2] + step * ((w1 >>> 10) & 0x1fffff);
  } else {
    out[0] = positions[j * 3];
    out[1] = positions[j * 3 + 1];
    out[2] = positions[j * 3 + 2];
  }
  return out;
}

const _q = [0, 0, 0, 1];
/** particle `j`'s unit quaternion (the identity without rotations, or of length 0) */
function readRotation(rotations: Float32Array | Int16Array | undefined, format: PageStoreFormat, j: number) {
  if (!rotations || format.rotations === 'none') {
    _q[0] = _q[1] = _q[2] = 0;
    _q[3] = 1;
    return _q;
  }
  let x = rotations[j * 4],
    y = rotations[j * 4 + 1],
    z = rotations[j * 4 + 2],
    w = rotations[j * 4 + 3];
  const n = Math.sqrt(x * x + y * y + z * z + w * w);
  if (n > 0) {
    x /= n;
    y /= n;
    z /= n;
    w /= n;
  } else {
    x = y = z = 0;
    w = 1;
  }
  _q[0] = x;
  _q[1] = y;
  _q[2] = z;
  _q[3] = w;
  return _q;
}

/**
 * A copy's pose (`CopyPoseFloats`) at `o` of `out`: a top particle at `p` turned by `q` (unit) placing a template about
 * its center `c` (`t = p - R c`, scale 1)
 */
export function writeCopyPose(out: Float32Array, o: number, p: Vec3, q: ArrayLike<number>, c: Vec3) {
  writePose(out, o, p, 1, [q[0], q[1], q[2], q[3]], c);
}

/** `t = p - s R c` and `s`, then the quaternion, at `o` of `out` (`R c`: `c + 2 w (q × c) + 2 q × (q × c)`) */
function writePose(out: Float32Array, o: number, p: Vec3, s: number, q: number[], c: Vec3) {
  const [qx, qy, qz, qw] = q;
  const tx = 2 * (qy * c[2] - qz * c[1]);
  const ty = 2 * (qz * c[0] - qx * c[2]);
  const tz = 2 * (qx * c[1] - qy * c[0]);
  const rx = c[0] + qw * tx + (qy * tz - qz * ty);
  const ry = c[1] + qw * ty + (qz * tx - qx * tz);
  const rz = c[2] + qw * tz + (qx * ty - qy * tx);
  out[o] = p[0] - s * rx;
  out[o + 1] = p[1] - s * ry;
  out[o + 2] = p[2] - s * rz;
  out[o + 3] = s;
  out[o + 4] = qx;
  out[o + 5] = qy;
  out[o + 6] = qz;
  out[o + 7] = qw;
}
