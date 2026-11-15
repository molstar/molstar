/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { Column } from '@molstar/core/data/db';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { GridLookup3D } from '@molstar/core/math/geometry/lookup3d/grid';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { CentroidHelper } from '@molstar/core/math/geometry/centroid-helper';
import { Grid } from '../volume/grid.js';
import type { ModelFormat } from '@molstar/model/formats/format';
import { CustomProperties, CustomPropertyDescriptor } from '../custom-property.js';
import type { Boundary } from '@molstar/core/math/geometry/boundary';
import { fillIdentityTransform } from '@molstar/core/math/linear-algebra/3d/transform-array';
import { GridPositions, type Positions } from '@molstar/core/math/geometry/grid-positions';
import type { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { type ParticlePages, type ParticlePagesChange, type ParticleSlot, ParticleSlots } from './particle-slots.js';

/**
 * A single named per-particle scalar attribute (e.g. `'cc'`, `'class'`, `'score'`).
 *
 * `column` provides the per-particle values (indexed by particle position, length = `count`).
 * Formats should reuse an existing `Column` (e.g. via `Column.view`) or provide a lambda-based
 * accessor instead of copying data whenever possible.
 */
export interface ParticleAttribute {
  readonly column: Column<any>;
  readonly label: string;
  readonly min: number;
  readonly max: number;
}

export interface ParticleCompartmentInfo {
  /** Compartment name/path (e.g. `"root.mge.surface.proteins"`). */
  readonly name: string;
  /**
   * The entities of the compartment's particles (`entityInfo` ids, in its order), where the whole is known though not
   * all its particles are at hand (a tile set's catalog, `ParticleTileCatalog.fromList`): what colors a view of some
   * tiles as the whole (`particle-hierarchy`)
   */
  readonly entities?: readonly number[];
  /** whether some of the compartment's particles have no entity (with `entities`) */
  readonly withoutEntity?: boolean;
}

export interface ParticleEntityInfo {
  /** Entity name (e.g. `"MG_191_192_NAP"`). */
  readonly name: string;
  /** Optional functional annotation of the entity (e.g. `"Lipid"`). */
  readonly function?: string;
}

/** Lightweight metadata for a single target id, always available (unlike `targetMapping`). */
/**
 * The box of a range of particles' positions and their largest distance from its center (Å): what a target's placement
 * needs of its particles' positions without a pass over them (`Particle.getRangeBounds`)
 */
export interface ParticleRangeBounds {
  readonly min: Vec3;
  readonly max: Vec3;
  readonly distance: number;
}

/** Values a range's bounds take in `Particle.setTargetBounds`: min xyz, max xyz, the distance */
export const ParticleRangeBoundsValues = 7;

export interface ParticleTargetInfo {
  /** Optional display name for this target (e.g. chain id, tomogram id, agent type name). */
  readonly name?: string;
  /** Index into `entityInfo`; states/enforces that every particle of this target shares this entity. */
  readonly entity?: number;
  /**
   * Index into `compartmentInfo`, -1 for none, where every particle of this target shares it (a pack's: stored a
   * target); none where not known so: a theme by compartment colors by particle
   */
  readonly compartment?: number;
}

/**
 * A value a particle (`ParticleList.targets`, `entities`, `compartments`, `radii`): an array, or runs, a value a run of
 * consecutive particles. A list sorted by target (a particle pack's) holds a target's values once, as a run, where an
 * array would repeat them a particle. Read by `ParticleValues`: `get` a particle,
 * `forEachRun` all of them (an array's runs of equal values)
 */
export type ParticleValues<A extends Int32Array | Float32Array = Int32Array> = A | ParticleRuns<A>;

/** Runs of particles and a value each: run `r` the particles `[offsets[r], offsets[r + 1])` */
export interface ParticleRuns<A extends Int32Array | Float32Array = Int32Array> {
  readonly kind: 'runs';
  /** a run's first particle, the particle count last (run count + 1 values) */
  readonly offsets: Int32Array;
  readonly values: A;
}

export namespace ParticleValues {
  export function isRuns<A extends Int32Array | Float32Array>(v: ParticleValues<A>): v is ParticleRuns<A> {
    return (v as ParticleRuns<A>).kind === 'runs';
  }

  export function ofRuns<A extends Int32Array | Float32Array>(offsets: Int32Array, values: A): ParticleRuns<A> {
    return { kind: 'runs', offsets, values };
  }

  /** the run found last by `get`: a particle after another of the same run found without a search */
  let lastRuns: ParticleRuns<Int32Array | Float32Array> | undefined;
  let lastRun = 0;

  /** The run of `runs` holding particle `i` (a binary search, the run found last tried first) */
  export function runOf(runs: ParticleRuns<Int32Array | Float32Array>, i: number) {
    const { offsets } = runs;
    if (runs === lastRuns && offsets[lastRun] <= i && i < offsets[lastRun + 1]) return lastRun;
    let lo = 0;
    let hi = offsets.length - 2;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (offsets[mid] <= i) lo = mid;
      else hi = mid - 1;
    }
    lastRuns = runs;
    lastRun = lo;
    return lo;
  }

  /** Particle `i`'s value */
  export function get(v: ParticleValues<Int32Array | Float32Array>, i: number): number {
    if (!isRuns(v)) return v[i];
    return v.values[runOf(v, i)];
  }

  /** `f` of each run of `v` over `count` particles: a run's own, an array's runs of equal values */
  export function forEachRun(
    v: ParticleValues<Int32Array | Float32Array>,
    count: number,
    f: (start: number, end: number, value: number) => void,
  ) {
    if (isRuns(v)) {
      const { offsets, values } = v;
      for (let r = 0, rl = values.length; r < rl; ++r) {
        if (offsets[r + 1] > offsets[r]) f(offsets[r], offsets[r + 1], values[r]);
      }
      return;
    }
    let start = 0;
    for (let i = 1; i <= count; ++i) {
      if (i === count || v[i] !== v[start]) {
        if (count > 0) f(start, i, v[start]);
        start = i;
      }
    }
  }

  /** `f` of each run where both `a` and `b` are the same over `count` particles */
  export function forEachRun2(
    a: ParticleValues<Int32Array | Float32Array>,
    b: ParticleValues<Int32Array | Float32Array>,
    count: number,
    f: (start: number, end: number, va: number, vb: number) => void,
  ) {
    if (!isRuns(a) && !isRuns(b)) {
      let start = 0;
      for (let i = 1; i <= count; ++i) {
        if (i === count || a[i] !== a[start] || b[i] !== b[start]) {
          if (count > 0) f(start, i, a[start], b[start]);
          start = i;
        }
      }
      return;
    }
    // the runs of each with their values, merged: a run ends where either's does (no search a run: a tile-ordered
    // list's runs are millions)
    const runs = (v: ParticleValues<Int32Array | Float32Array>) => {
      const ends: number[] = [];
      const values: number[] = [];
      forEachRun(v, count, (_, end, value) => {
        ends.push(end);
        values.push(value);
      });
      return { ends, values };
    };
    const { ends: ea, values: va } = runs(a);
    const { ends: eb, values: vb } = runs(b);
    let start = 0;
    let ia = 0;
    let ib = 0;
    while (start < count) {
      const end = Math.min(ea[ia], eb[ib]);
      f(start, end, va[ia], vb[ib]);
      if (ea[ia] === end) ia++;
      if (eb[ib] === end) ib++;
      start = end;
    }
  }

  /**
   * `f` of each run of `runs` holding particles of `indices` and how many of them, in order: from run to run by a
   * search (not the runs between, a tile-ordered list's target a run a tile it is in); stops where `f` returns false
   */
  export function forEachRunOf(
    runs: ParticleRuns<Int32Array | Float32Array>,
    indices: OrderedSet<number>,
    f: (run: number, count: number) => boolean | void,
  ) {
    const n = OrderedSet.size(indices);
    for (let k = 0; k < n; ) {
      const r = runOf(runs, OrderedSet.getAt(indices, k));
      const next = OrderedSet.findPredecessorIndex(indices, runs.offsets[r + 1]);
      if (f(r, next - k) === false) return;
      k = next;
    }
  }

  /** `array` as runs where they are fewer than half its values, else as it is */
  export function compact<A extends Int32Array | Float32Array>(array: A): ParticleValues<A> {
    const n = array.length;
    let runs = n > 0 ? 1 : 0;
    for (let i = 1; i < n; ++i) if (array[i] !== array[i - 1]) runs += 1;
    if (runs * 2 > n) return array;
    const offsets = new Int32Array(runs + 1);
    const values = new (array.constructor as { new (n: number): A })(runs);
    for (let i = 0, r = -1; i < n; ++i) {
      if (i === 0 || array[i] !== array[i - 1]) {
        r += 1;
        offsets[r] = i;
        values[r] = array[i];
      }
    }
    offsets[runs] = n;
    return ofRuns(offsets, values);
  }

  /** `v` as an array, a value a particle (`count`): where one is needed whole */
  export function toArray<A extends Int32Array | Float32Array>(
    v: ParticleValues<A>,
    count: number,
    ctor: { new (n: number): A },
  ): A {
    if (!isRuns(v)) return v;
    const out = new ctor(count);
    forEachRun(v, count, (start, end, value) => out.fill(value, start, end));
    return out;
  }

  /**
   * Particles `[first, first + count)` of `v` from 0: an array's a view (`subarray`), runs' the runs over the range
   * (their values a view, their offsets made: a run a value, not a particle)
   */
  export function slice<A extends Int32Array | Float32Array>(
    v: ParticleValues<A>,
    first: number,
    count: number,
  ): ParticleValues<A> {
    if (!isRuns(v)) return v.subarray(first, first + count) as A;
    if (count === 0) return ofRuns(Int32Array.of(0), v.values.subarray(0, 0) as A);
    const r0 = runOf(v, first);
    const r1 = runOf(v, first + count - 1);
    const offsets = new Int32Array(r1 - r0 + 2);
    for (let r = r0; r <= r1; ++r) offsets[r - r0] = Math.max(0, v.offsets[r] - first);
    offsets[r1 - r0 + 1] = count;
    return ofRuns(offsets, v.values.subarray(r0, r1 + 1) as A);
  }

  /**
   * `f` of each value of `v`, in its form: runs mapped run by run (a value a run made, not a particle), an array value
   * by value
   */
  export function map(
    v: ParticleValues<Int32Array | Float32Array>,
    count: number,
    f: (value: number) => number,
  ): ParticleValues<Int32Array> {
    if (isRuns(v)) {
      // a loop, not `Int32Array.from(values, f)` (its mapping by the iterator protocol: 60 ns a run, 2.29 M in the yeast
      // tile-ordered); a run's value as the last one's reused (a tile's neighboring targets often of one entity)
      const values = v.values;
      const out = new Int32Array(values.length);
      let last = Number.NaN;
      let y = 0;
      for (let i = 0; i < values.length; ++i) {
        const x = values[i];
        if (x !== last) {
          y = f(x);
          last = x;
        }
        out[i] = y;
      }
      return ofRuns(v.offsets, out);
    }
    const out = new Int32Array(count);
    const cache = new Map<number, number>();
    for (let i = 0; i < count; ++i) {
      const x = v[i];
      let y = cache.get(x);
      if (y === undefined) cache.set(x, (y = f(x)));
      out[i] = y;
    }
    return out;
  }
}

export interface ParticleList {
  readonly entryId?: string;
  readonly label?: string;

  readonly count: number;

  /**
   * Optional periodic/bounding cell for the particle system (e.g. the simulation box
   * size of a Simularium trajectory).
   */
  readonly cell?: Cell;

  /** Unique keys for each particle for mapping to source data; none: a particle's key is its index. */
  readonly keys?: Int32Array;

  /**
   * Per-particle target index (`ParticleValues`: an array, or runs of particles). Each value identifies which
   * target structure/volume/shape this particle belongs to. Use 0 for
   * single-target data. The distinct values correspond to the
   * keys of `targetMapping` when present.
   */
  readonly targets: ParticleValues<Int32Array>;

  /**
   * Metadata for each unique target ID in `targets`. Always present (unlike `targetMapping`),
   * so target ids can be enumerated cheaply without scanning `targets`.
   */
  readonly targetInfo: ReadonlyMap<number, ParticleTargetInfo>;

  /**
   * Optional mapping from each unique target ID in `targets` to the reference object
   * instanced at every particle with that target ID. The objects are owned by the
   * particle list itself and are not part of the plugin state tree.
   */
  readonly targetMapping?: ReadonlyMap<number, ParticleTarget>;

  /**
   * Optional per-particle compartment index (`ParticleValues`). Each value identifies
   * which compartment this particle belongs to. A value of -1 means "no compartment".
   * The distinct non-negative values correspond to keys of `compartmentInfo`.
   */
  readonly compartments?: ParticleValues<Int32Array>;

  /**
   * Optional mapping from each unique compartment index in `compartments` to the
   * compartment info.
   */
  readonly compartmentInfo?: ReadonlyMap<number, ParticleCompartmentInfo>;

  /**
   * Optional per-particle entity index (`ParticleValues`). Each value identifies
   * which entity (molecule type) this particle belongs to. A value of -1 means "no entity".
   * The distinct non-negative values correspond to keys of `entityInfo`.
   */
  readonly entities?: ParticleValues<Int32Array>;

  /**
   * Optional mapping from each unique entity index in `entities` to the entity info.
   */
  readonly entityInfo?: ReadonlyMap<number, ParticleEntityInfo>;

  /** Particle positions in angstrom, packed as `[x0, y0, z0, x1, y1, z1, ...]`. */
  /**
   * Particle positions in angstrom: float32 `[x0, y0, z0, x1, ...]`, or on a grid (`GridPositions`: 21 bits an axis,
   * a pack's). Read either by `GridPositions.get`; loops over all branch on the kind once (`instance-transforms.ts`)
   */
  readonly coordinates: Positions;
  /**
   * Optional per-particle orientations as quaternions, packed as `[x0, y0, z0, w0, ...]` (`ParticleRotations`): unit
   * float32, or snorm16 (a pack's)
   */
  readonly rotations?: ParticleRotations;
  /** Optional per-particle bounding sphere radii in angstrom (`ParticleValues`). */
  readonly radii?: ParticleValues<Float32Array>;

  /**
   * Optional polyline (fiber) connectivity over particles, stored in compressed-sparse-row
   * form. Fiber `f` (for `0 <= f < count`) is the ordered polyline through the particles
   * `indices[offsets[f]]` .. `indices[offsets[f + 1]) - 1]`. Used by formats such as
   * Simularium where an agent expands into a chain of particles.
   */
  readonly fibers?: {
    readonly count: number;
    readonly offsets: Int32Array;
    readonly indices: Int32Array;
  };

  /**
   * Named per-particle scalar attributes, indexed by particle position (length = count).
   * Keys are short identifiers (e.g. 'cc', 'class', 'score').
   */
  readonly attributes?: ReadonlyMap<string, ParticleAttribute>;

  /** A particle in words, its index from 1 first: "Particle 165 · tomo 1 · class 2 · cc 0.4213" */
  readonly getParticleLabel: (index: number) => string;

  readonly sourceData: ModelFormat;

  /**
   * The list this one is a copy of moved as a whole (`Particle.transform`) and the move from it: the first of a line
   * of moves, the moves composed, so two lists of a line give the move between them (`Particle.getMove`) and what was
   * built for one can be moved for the other
   */
  readonly movedFrom?: { readonly particles: ParticleList; readonly transform: Mat4 };

  /**
   * The list this one is a copy of with its particles' poses changed alone (`Particle.withPoses`: their positions and
   * rotations; the same count, targets, entities, compartments, radii, fibers and attributes): the first of a line of
   * such copies, so two lists of a line are the same particles posed otherwise (`Particle.isSamePoses`: a trajectory's
   * frames, a drag), what was built of their identity holds for both (colors by entity, the instances' order)
   */
  readonly posesFrom?: ParticleList;

  /**
   * A slot list's pages (`ParticleSlots`): a particle's index its slot, `[0, count)` the slots, a page a range of them
   * with a flat list of its own (a tile's), the slots no page holds gaps. Its own per-particle fields are empty: what
   * reads its particles reads its pages (`ParticleSlots.forEachPage`); none for a flat list (every slot a particle)
   */
  readonly pages?: ParticlePages;

  /**
   * The slot list this one is the next version of, and the pages added and removed between them
   * (`ParticleSlots.getChanges`): what its consumers update by; one step only (a version holds the one before it alone)
   */
  readonly pagesFrom?: ParticlePagesChange;

  /**
   * its particles drawn scaled by this, 1 without: a coarse level's node (`ParticleTileLevels`), a sample of its tiles'
   * particles drawn larger so they fill about the space all fill
   */
  readonly scale?: number;

  customProperties: CustomProperties;
  _propertyData: { [name: string]: any };
}

/**
 * A reference object instanced at each particle of a given target id. Each particle has
 * exactly one target id (see `ParticleList.targets`); the distinct target ids map to these
 * targets via `ParticleList.targetMapping`.
 */
export type ParticleTarget =
  | {
      readonly kind: 'structure';
      readonly structure: import('../structure/structure.js').Structure;
      /**
       * The point of the structure placed at a particle's coordinates: the centroid of its atoms (the default,
       * `ParticleTarget.getCenter`), given when the format has it already (the centroid it placed its particles by)
       */
      readonly center?: Vec3;
    }
  | {
      /** typed spheres (atoms, coarse beads) without a structure: what a target's geometry needs, read from its file */
      readonly kind: 'sites';
      readonly sites: ParticleSites;
      /** as a structure target's: the point placed at a particle's coordinates, the sites' centroid by default */
      readonly center?: Vec3;
    }
  | { readonly kind: 'shape'; readonly shape: import('../shape/shape.js').Shape }
  | { readonly kind: 'volume'; readonly volume: import('../volume/volume.js').Volume };

/**
 * The sites of a target (`ParticleTarget` kind `sites`): positions and radii, a site's index its group in the target's
 * geometry. A few bytes a site, not a model and a structure of them (millions of atoms in a mesoscale model's targets)
 */
export interface ParticleSites {
  readonly count: number;
  /** x, y, z a site (decoded on each read where the sites are `quantized`: read a site by `ParticleSites.getPosition`) */
  readonly positions: Float32Array;
  /** the physical radius a site (an atom's van der Waals radius, a bead's own; decoded where `quantized`) */
  readonly radii: Float32Array;
  /**
   * the atomic number a site (`AtomNumber`, 0 for none: a bead), where the sites are atoms of known elements (decoded
   * where `quantized`, from its `typeElements`); read a site by `ParticleSites.getElement`
   */
  readonly elements?: Uint8Array;
  readonly label?: string;
  /** the sites' chains, their sites consecutive (as a structure's units): runs, a few bytes a chain, not an index a site */
  readonly chains?: ParticleSiteChains;
  /** the sites as a particle pack holds them: positions quantized, radii by type (8 bytes a site, not 16) */
  readonly quantized?: ParticleSitesQuantized;
  /**
   * the sites every `stride`-th of their target's (`ParticleSites.coarsen`), 1 without: drawn as that stride's level
   * of detail is near the camera too (fewer, larger spheres), a tenth of the bytes at 10
   */
  readonly stride?: number;
}

/**
 * Quantized sites: a site's position `origin + q * scale` (`q` int16 x, y, z), its radius `typeRadii[type]` and its
 * element `typeElements[type]` where known (the type a uint16, as int16 bits): `data` x, y, z, type a site
 */
export interface ParticleSitesQuantized {
  readonly data: Int16Array;
  readonly origin: Vec3;
  readonly scale: number;
  readonly typeRadii: Float32Array;
  /** the atomic number a type (`ParticleSites.elements`), where the sites' elements are known */
  readonly typeElements?: Uint8Array;
  /**
   * the extremes of the quantized positions of each chain (`ParticleSites.chains`; without, of all sites): min x, y, z,
   * max x, y, z, six a chain (`getQuantizedExtremes`); a pack's, where it has them: a chain's box with no pass
   */
  readonly extremes?: Int16Array;
}

/**
 * Particles' rotations, quaternions `[x, y, z, w]` a particle: unit as float32, or snorm16 (`Int16Array`: a unit
 * component times 32767, rounded; a pack's, 8 B a particle where float32 is 16 B, within about 1e-4 rad). Readers take
 * either as it is and rotate by `q` normalized (`R(q)` with `2 / |q|²` in place of 2), so no decoded copy is made
 */
export type ParticleRotations = Float32Array | Int16Array;

export namespace ParticleRotations {
  /** A snorm16 component's integer for 1 */
  export const Snorm16One = 32767;

  /** Particle `i`'s rotation of `rotations` (either kind) into `out`, normalized; the identity where it has none */
  export function get(rotations: ParticleRotations | undefined, i: number, out: Quat): Quat {
    if (!rotations) {
      Quat.setIdentity(out);
      return out;
    }
    const o = i * 4;
    const x = rotations[o],
      y = rotations[o + 1],
      z = rotations[o + 2],
      w = rotations[o + 3];
    const n = Math.sqrt(x * x + y * y + z * z + w * w);
    if (!(n > 0)) {
      Quat.setIdentity(out);
      return out;
    }
    return Quat.set(out, x / n, y / n, z / n, w / n);
  }

  /**
   * `rotations` as snorm16 (`Int16Array` as it is): each normalized, its sign so that `w >= 0` (q and -q the same
   * rotation), its components times 32767, rounded
   */
  export function toSnorm16(rotations: ParticleRotations): Int16Array {
    if (rotations instanceof Int16Array) return rotations;
    const out = new Int16Array(rotations.length);
    for (let o = 0; o < rotations.length; o += 4) {
      const x = rotations[o],
        y = rotations[o + 1],
        z = rotations[o + 2],
        w = rotations[o + 3];
      const n = Math.sqrt(x * x + y * y + z * z + w * w);
      if (!(n > 0)) {
        out[o + 3] = Snorm16One;
        continue;
      }
      const k = (w < 0 ? -Snorm16One : Snorm16One) / n;
      out[o] = Math.round(x * k);
      out[o + 1] = Math.round(y * k);
      out[o + 2] = Math.round(z * k);
      out[o + 3] = Math.round(w * k);
    }
    return out;
  }
}

export namespace ParticleSites {
  /** Site `i`'s position into `out`: read in place, quantized or not (no decoded array made) */
  export function getPosition(sites: ParticleSites, i: number, out: Vec3): Vec3 {
    const q = sites.quantized;
    if (!q) return Vec3.fromArray(out, sites.positions, i * 3);
    const { data, origin, scale } = q;
    const o = i * 4;
    out[0] = origin[0] + data[o] * scale;
    out[1] = origin[1] + data[o + 1] * scale;
    out[2] = origin[2] + data[o + 2] * scale;
    return out;
  }

  /** Site `i`'s radius: read in place, by its type where quantized */
  export function getRadius(sites: ParticleSites, i: number): number {
    const q = sites.quantized;
    if (!q) return sites.radii[i];
    return q.typeRadii[q.data[i * 4 + 3] & 0xffff];
  }

  /**
   * The extremes of the quantized positions `[start, end)` of `data` (x, y, z, type a site) into `out` at `offset`:
   * min x, y, z, max x, y, z (`ParticleSitesQuantized.extremes`); none (`start === end`): zeros
   */
  export function getQuantizedExtremes(
    data: Int16Array,
    start: number,
    end: number,
    out: Int16Array,
    offset = 0,
  ): Int16Array {
    if (start >= end) {
      out.fill(0, offset, offset + 6);
      return out;
    }
    let x0 = 32767;
    let y0 = 32767;
    let z0 = 32767;
    let x1 = -32768;
    let y1 = -32768;
    let z1 = -32768;
    for (let o = start * 4, ol = end * 4; o < ol; o += 4) {
      const x = data[o];
      const y = data[o + 1];
      const z = data[o + 2];
      if (x < x0) x0 = x;
      if (x > x1) x1 = x;
      if (y < y0) y0 = y;
      if (y > y1) y1 = y;
      if (z < z0) z0 = z;
      if (z > z1) z1 = z;
    }
    out[offset] = x0;
    out[offset + 1] = y0;
    out[offset + 2] = z0;
    out[offset + 3] = x1;
    out[offset + 4] = y1;
    out[offset + 5] = z1;
    return out;
  }

  /** Whether the sites' elements are known (`elements`, a quantized type's `typeElements`) */
  export function hasElements(sites: ParticleSites): boolean {
    return sites.quantized ? !!sites.quantized.typeElements : !!sites.elements;
  }

  /** Site `i`'s atomic number (0 for none), by its type where quantized; 0 where the elements are not known */
  export function getElement(sites: ParticleSites, i: number): number {
    const q = sites.quantized;
    if (!q) return sites.elements ? sites.elements[i] : 0;
    return q.typeElements ? q.typeElements[q.data[i * 4 + 3] & 0xffff] : 0;
  }

  /**
   * Every `stride`-th of `sites` (their `stride` times `stride`): a copy, the sites' order as the levels of detail
   * take it (a level every stride-th site), so the coarse sites are a level's. A chain's first site the first of its
   * kept ones (a chain of fewer than `stride` sites may keep none); its extremes as they were (a box around the kept
   * ones). `sites` itself at a stride of 1
   */
  export function coarsen(sites: ParticleSites, stride: number): ParticleSites {
    stride = Math.floor(stride);
    if (!(stride > 1)) return sites;
    const count = Math.ceil(sites.count / stride);
    const chains = sites.chains && {
      offsets: Int32Array.from(sites.chains.offsets, (o) => Math.ceil(o / stride)),
      labels: sites.chains.labels,
    };
    const base = (sites.stride ?? 1) * stride;
    const q = sites.quantized;
    if (q) {
      const data = new Int16Array(count * 4);
      for (let i = 0, j = 0; i < count; ++i, j += stride) data.set(q.data.subarray(j * 4, j * 4 + 4), i * 4);
      return ofQuantized(count, { ...q, data }, sites.label, chains, base);
    }
    const positions = new Float32Array(count * 3);
    const radii = new Float32Array(count);
    const elements = sites.elements && new Uint8Array(count);
    for (let i = 0, j = 0; i < count; ++i, j += stride) {
      positions.set(sites.positions.subarray(j * 3, j * 3 + 3), i * 3);
      radii[i] = sites.radii[j];
      if (elements) elements[i] = sites.elements![j];
    }
    return { count, positions, radii, elements, label: sites.label, chains, stride: base };
  }

  /**
   * Sites of `quantized`: their `positions` and `radii` decoded on each read (arrays made where a reader takes them
   * whole); the readers of a site at a time read the quantized data (`getPosition`, `getRadius`)
   */
  export function ofQuantized(
    count: number,
    quantized: ParticleSitesQuantized,
    label?: string,
    chains?: ParticleSiteChains,
    stride?: number,
  ): ParticleSites {
    const sites: ParticleSites = {
      count,
      quantized,
      label,
      chains,
      stride,
      get positions() {
        const out = new Float32Array(count * 3);
        const p = Vec3();
        for (let i = 0; i < count; ++i) Vec3.toArray(getPosition(sites, i, p), out, i * 3);
        return out;
      },
      get radii() {
        const out = new Float32Array(count);
        for (let i = 0; i < count; ++i) out[i] = getRadius(sites, i);
        return out;
      },
      get elements() {
        if (!quantized.typeElements) return undefined;
        const out = new Uint8Array(count);
        for (let i = 0; i < count; ++i) out[i] = getElement(sites, i);
        return out;
      },
    };
    return sites;
  }
}

/** Chain `c` of `ParticleSites`: its sites `[offsets[c], offsets[c + 1])`, its label (`label_asym_id`) */
export interface ParticleSiteChains {
  /** a chain's first site, the site count last (chain count + 1 values) */
  readonly offsets: Int32Array;
  readonly labels: readonly string[];
}

export namespace ParticleTarget {
  export function data(target: ParticleTarget) {
    switch (target.kind) {
      case 'structure':
        return target.structure;
      case 'sites':
        return target.sites;
      case 'shape':
        return target.shape;
      case 'volume':
        return target.volume;
    }
  }

  const centroids = new WeakMap<object, Vec3>();
  const position = Vec3();

  /**
   * The point of the target placed at a particle's coordinates: a structure's centroid (`center` when given, cached
   * per structure otherwise); a volume's box center (voxel N/2, the subtomogram formats' convention), not its
   * isosurface's, which moves with the iso value; none for shapes, placed by their geometry's bounding sphere center
   */
  export function getCenter(target: ParticleTarget): Vec3 | undefined {
    if (target.kind === 'shape') return undefined;
    if (target.kind === 'volume') return Grid.getBoundingSphere(target.volume.grid).center;
    if (target.center) return target.center;
    if (target.kind === 'sites') {
      const { sites } = target;
      let center = centroids.get(sites);
      if (!center) {
        const helper = new CentroidHelper();
        for (let i = 0; i < sites.count; ++i) helper.includeStep(ParticleSites.getPosition(sites, i, position));
        helper.finishedIncludeStep();
        center = Vec3.clone(helper.center);
        centroids.set(sites, center);
      }
      return center;
    }
    const { structure } = target;
    let center = centroids.get(structure);
    if (!center) {
      const helper = new CentroidHelper();
      for (const unit of structure.units) {
        const { elements, conformation } = unit;
        for (let i = 0, il = elements.length; i < il; ++i) {
          helper.includeStep(conformation.position(elements[i], position));
        }
      }
      helper.finishedIncludeStep();
      center = Vec3.clone(helper.center);
      centroids.set(structure, center);
    }
    return center;
  }
}

const ParticleSpacingDescriptor = CustomPropertyDescriptor({ name: 'particle-spacing-radius' });

/** How many particles the spacing is measured over (evenly spread): the median settles well before */
const SpacingSamples = 2000;

/**
 * A radius for particles that carry none (picks, oriented points): half the median distance from a particle to its
 * nearest neighbor, so neighbors about touch, as cryo-ET picks of one kind do. Measured on up to 2,000 particles,
 * cached on the list; 1 for fewer than two particles or all at one point.
 */
/** A position read (`GridPositions.get`) by the functions here */
const _position = Vec3();

export function getParticleSpacingRadius(data: ParticleList): number {
  const cached = data._propertyData[ParticleSpacingDescriptor.name] as number | undefined;
  if (cached !== undefined) return cached;
  const { count, coordinates } = data;
  let radius = 1;
  if (count > 1) {
    const x = new Float32Array(count);
    const y = new Float32Array(count);
    const z = new Float32Array(count);
    for (let i = 0; i < count; i++) {
      GridPositions.get(coordinates, i, _position);
      x[i] = _position[0];
      y[i] = _position[1];
      z[i] = _position[2];
    }
    const lookup = GridLookup3D({ x, y, z, indices: OrderedSet.ofBounds(0, count) }, Particle.getBoundary(data));
    const step = Math.max(1, Math.floor(count / SpacingSamples));
    const nearest: number[] = [];
    for (let i = 0; i < count; i += step) {
      const r = lookup.nearest(x[i], y[i], z[i], 2);
      let best = Infinity;
      for (let k = 0; k < r.count; k++) {
        if (r.indices[k] !== i && r.squaredDistances[k] > 0) best = Math.min(best, r.squaredDistances[k]);
      }
      if (best < Infinity) nearest.push(Math.sqrt(best));
    }
    if (nearest.length) {
      nearest.sort((a, b) => a - b);
      radius = nearest[nearest.length >> 1] / 2;
    }
  }
  data.customProperties.add(ParticleSpacingDescriptor);
  data._propertyData[ParticleSpacingDescriptor.name] = radius;
  return radius;
}

/**
 * Particle `i`'s transform (column-major mat4, its rotation then its position) into `out` at `offset`, all 16 values:
 * as `getParticleTransforms` holds it, without the list's copy of all (a target's own transforms, millions of
 * particles)
 */
export function writeParticleTransform(data: ParticleList, i: number, out: Float32Array, offset: number) {
  const { rotations, coordinates } = data;
  const o = offset;
  if (rotations) {
    // quat → mat4 (see Mat4.fromQuat), by 2 / |q|² (`ParticleRotations`: either kind); the identity gives the identity
    const q = i * 4;
    const x = rotations[q],
      y = rotations[q + 1],
      z = rotations[q + 2],
      w = rotations[q + 3];
    const nq = x * x + y * y + z * z + w * w;
    const f = nq > 0 ? 2 / nq : 0;
    const x2 = x * f,
      y2 = y * f,
      z2 = z * f;
    const xx = x * x2,
      yx = y * x2,
      yy = y * y2;
    const zx = z * x2,
      zy = z * y2,
      zz = z * z2;
    const wx = w * x2,
      wy = w * y2,
      wz = w * z2;
    out[o] = 1 - yy - zz;
    out[o + 1] = yx + wz;
    out[o + 2] = zx - wy;
    out[o + 4] = yx - wz;
    out[o + 5] = 1 - xx - zz;
    out[o + 6] = zy + wx;
    out[o + 8] = zx + wy;
    out[o + 9] = zy - wx;
    out[o + 10] = 1 - xx - yy;
  } else {
    out[o] = 1;
    out[o + 1] = 0;
    out[o + 2] = 0;
    out[o + 4] = 0;
    out[o + 5] = 1;
    out[o + 6] = 0;
    out[o + 8] = 0;
    out[o + 9] = 0;
    out[o + 10] = 1;
  }
  out[o + 3] = 0;
  out[o + 7] = 0;
  out[o + 11] = 0;
  GridPositions.get(coordinates, i, _position);
  out[o + 12] = _position[0];
  out[o + 13] = _position[1];
  out[o + 14] = _position[2];
  out[o + 15] = 1;
}

const ParticleTransformsDescriptor = CustomPropertyDescriptor({ name: 'particle-transforms' });

/**
 * Per-particle transforms as a flat array, 16 consecutive floats (column-major mat4) per
 * particle. Computed once and cached on the `ParticleList`.
 */
export function getParticleTransforms(data: ParticleList): Float32Array {
  let cached = data._propertyData[ParticleTransformsDescriptor.name] as Float32Array | undefined;
  if (!cached) {
    const particleCount = data.count;
    const transformArray = new Float32Array(particleCount * 16);
    const { rotations, coordinates } = data;

    let hasRotations = false;
    if (rotations) {
      for (let i = 0; i < particleCount; ++i) {
        const o = i * 4;
        // the identity: x, y, z zero (w of either sign and kind)
        if (rotations[o + 0] !== 0 || rotations[o + 1] !== 0 || rotations[o + 2] !== 0) {
          hasRotations = true;
          break;
        }
      }
    }

    if (rotations && hasRotations) {
      for (let i = 0; i < particleCount; ++i) writeParticleTransform(data, i, transformArray, i * 16);
    } else {
      fillIdentityTransform(transformArray, particleCount);
      for (let i = 0; i < particleCount; ++i) {
        GridPositions.get(coordinates, i, _position);
        transformArray[i * 16 + 12] = _position[0];
        transformArray[i * 16 + 13] = _position[1];
        transformArray[i * 16 + 14] = _position[2];
      }
    }

    data.customProperties.add(ParticleTransformsDescriptor);
    data._propertyData[ParticleTransformsDescriptor.name] = transformArray;
    cached = transformArray;
  }
  return cached;
}

const ParticleTargetGroupsDescriptor = CustomPropertyDescriptor({ name: 'particle-target-groups' });

/** The particles of a `ParticleList` grouped by their target id. */
export interface ParticleTargetGroups {
  /** The distinct target ids, ascending. */
  readonly targetIds: Int32Array;
  /** The particles of group `g`, as indices into the unfiltered `ParticleList`. */
  readonly sets: ReadonlyArray<OrderedSet<number>>;
  readonly groupOfTarget: ReadonlyMap<number, number>;
}

/**
 * Group the particles by target id. Computed once and cached on the `ParticleList`.
 *
 * When `targets` is already non-decreasing - the usual case, since formats emit particles
 * grouped by target - each group is a contiguous range and is returned as an `Interval`,
 * which needs no per-group storage and makes `OrderedSet.getAt`/`indexOf` arithmetic rather
 * than a binary search. Otherwise the groups are built with a counting sort.
 */
export function getParticleTargetGroups(data: ParticleList): ParticleTargetGroups {
  if (!data._propertyData[ParticleTargetGroupsDescriptor.name]) {
    data.customProperties.add(ParticleTargetGroupsDescriptor);
    data._propertyData[ParticleTargetGroupsDescriptor.name] = computeParticleTargetGroups(data);
  }
  return data._propertyData[ParticleTargetGroupsDescriptor.name];
}

function computeParticleTargetGroups(data: ParticleList): ParticleTargetGroups {
  const { count } = data;

  const groupOfTarget = new Map<number, number>();

  if (count === 0) {
    return { targetIds: new Int32Array(0), sets: [], groupOfTarget };
  }
  if (ParticleValues.isRuns(data.targets)) return getRunTargetGroups(data.targets, groupOfTarget);
  const targets = data.targets;

  let sorted = true;
  let minId = targets[0];
  let maxId = targets[0];
  let runCount = 1;
  for (let i = 1; i < count; ++i) {
    const t = targets[i];
    const p = targets[i - 1];
    if (t < p) sorted = false;
    if (t !== p) runCount += 1;
    if (t < minId) minId = t;
    if (t > maxId) maxId = t;
  }

  if (sorted) {
    // runs of equal ids are the distinct ids and are already ascending
    const targetIds = new Int32Array(runCount);
    const sets: OrderedSet<number>[] = new Array(runCount);
    let g = 0;
    let start = 0;
    for (let i = 1; i <= count; ++i) {
      if (i === count || targets[i] !== targets[start]) {
        targetIds[g] = targets[start];
        sets[g] = OrderedSet.ofBounds(start, i);
        groupOfTarget.set(targets[start], g);
        g += 1;
        start = i;
      }
    }
    return { targetIds, sets, groupOfTarget };
  }

  // Counting sort keyed on the target id, offset so that negative ids are handled too.
  let targetIds: Int32Array;
  const span = maxId - minId + 1;
  const useHistogram = span <= 4 * count + 1024;
  if (useHistogram) {
    const histogram = new Int32Array(span);
    for (let i = 0; i < count; ++i) histogram[targets[i] - minId] += 1;
    let distinctCount = 0;
    for (let s = 0; s < span; ++s) {
      if (histogram[s] !== 0) distinctCount += 1;
    }
    targetIds = new Int32Array(distinctCount);
    let g = 0;
    for (let s = 0; s < span; ++s) {
      if (histogram[s] !== 0) {
        targetIds[g] = s + minId;
        groupOfTarget.set(s + minId, g);
        g += 1;
      }
    }
  } else {
    const distinct = new Set<number>();
    for (let i = 0; i < count; ++i) distinct.add(targets[i]);
    const ids = Array.from(distinct).sort((a, b) => a - b);
    targetIds = new Int32Array(ids.length);
    for (let g = 0; g < ids.length; ++g) {
      targetIds[g] = ids[g];
      groupOfTarget.set(ids[g], g);
    }
  }

  const groupCount = targetIds.length;
  const sets: OrderedSet<number>[] = new Array(groupCount);
  const offsets = new Int32Array(groupCount + 1);
  for (let i = 0; i < count; ++i) offsets[groupOfTarget.get(targets[i])! + 1] += 1;
  for (let g = 0; g < groupCount; ++g) offsets[g + 1] += offsets[g];

  const indices = new Int32Array(count);
  const cursor = Int32Array.from(offsets.subarray(0, groupCount));
  for (let i = 0; i < count; ++i) {
    const g = groupOfTarget.get(targets[i])!;
    indices[cursor[g]] = i;
    cursor[g] += 1;
  }
  for (let g = 0; g < groupCount; ++g) {
    sets[g] = OrderedSet.ofSortedArray(indices.subarray(offsets[g], offsets[g + 1]));
  }

  return { targetIds, sets, groupOfTarget };
}

/**
 * The groups of targets as runs (`ParticleValues`): a target of one run (a pack's) an `Interval`, no particle visited;
 * one of several their particles, in order
 */
function getRunTargetGroups(runs: ParticleRuns<Int32Array>, groupOfTarget: Map<number, number>): ParticleTargetGroups {
  const { offsets, values } = runs;
  const byTarget = new Map<number, number[]>();
  for (let r = 0; r < values.length; ++r) {
    if (offsets[r + 1] === offsets[r]) continue;
    let list = byTarget.get(values[r]);
    if (!list) byTarget.set(values[r], (list = []));
    list.push(r);
  }
  const targetIds = Int32Array.from(byTarget.keys()).sort();
  const sets: OrderedSet<number>[] = new Array(targetIds.length);
  for (let g = 0; g < targetIds.length; ++g) {
    const id = targetIds[g];
    groupOfTarget.set(id, g);
    const list = byTarget.get(id)!;
    if (list.length === 1) {
      sets[g] = OrderedSet.ofBounds(offsets[list[0]], offsets[list[0] + 1]);
      continue;
    }
    let n = 0;
    for (const r of list) n += offsets[r + 1] - offsets[r];
    const indices = new Int32Array(n);
    let k = 0;
    for (const r of list) for (let i = offsets[r]; i < offsets[r + 1]; ++i) indices[k++] = i;
    sets[g] = OrderedSet.ofSortedArray(indices);
  }
  return { targetIds, sets, groupOfTarget };
}

const FiberParticleMaskDescriptor = CustomPropertyDescriptor({ name: 'particle-fiber-mask' });

/**
 * Per-particle mask (length = `count`) marking fiber particles, i.e. particles referenced
 * by `fibers.indices`. Returns `undefined` if the particle list has no `fibers` data.
 * Computed once and cached on the `ParticleList`.
 */
export function getFiberParticleMask(data: ParticleList): Uint8Array | undefined {
  const { fibers } = data;
  if (!fibers) return undefined;
  if (!data._propertyData[FiberParticleMaskDescriptor.name]) {
    const mask = new Uint8Array(data.count);
    const { offsets, indices } = fibers;
    for (let f = 0; f < fibers.count; ++f) {
      const start = offsets[f];
      const end = offsets[f + 1];
      for (let i = start + 0; i < end; ++i) {
        mask[indices[i]] = 1;
      }
    }
    data.customProperties.add(FiberParticleMaskDescriptor);
    data._propertyData[FiberParticleMaskDescriptor.name] = mask;
  }
  return data._propertyData[FiberParticleMaskDescriptor.name];
}

export namespace Particle {
  export function is(x: any): x is ParticleList {
    return (
      typeof x?.count === 'number' &&
      (x?.keys === undefined || x.keys instanceof Int32Array) &&
      (x?.targets instanceof Int32Array || x?.targets?.kind === 'runs') &&
      (x?.coordinates instanceof Float32Array || GridPositions.is(x?.coordinates)) &&
      !!x?.sourceData &&
      !!x?.customProperties &&
      !!x?._propertyData
    );
  }

  /** A single particle within a `ParticleList`. */
  export interface Location {
    readonly kind: 'particle-location';
    particles: ParticleList;
    /** Particle index in the list. */
    index: number;
  }

  export function Location(particles?: ParticleList, index = 0): Location {
    return { kind: 'particle-location', particles: particles!, index };
  }
  export function isLocation(x: any): x is Location {
    return !!x && x.kind === 'particle-location';
  }
  const _slot: ParticleSlot = { list: undefined as unknown as ParticleList, index: 0 };
  /** Write the particle's position into `out` (a slot list's: its page's; a gap's the origin). */
  export function position(out: Vec3, location: Location): Vec3 {
    const { particles, index } = location;
    if (!particles.pages) return GridPositions.get(particles.coordinates, index, out);
    const s = ParticleSlots.locate(particles, index, _slot);
    return s ? GridPositions.get(s.list.coordinates, s.index, out) : Vec3.set(out, 0, 0, 0);
  }

  /** A loci over one or more particles in a `ParticleList`. */
  export interface Loci {
    readonly kind: 'particle-loci';
    readonly particles: ParticleList;
    readonly indices: OrderedSet<number>;
  }
  export function Loci(particles: ParticleList, indices: OrderedSet<number>): Loci {
    return { kind: 'particle-loci', particles, indices };
  }
  export function isLoci(x: any): x is Loci {
    return !!x && x.kind === 'particle-loci';
  }
  export function areLociEqual(a: Loci, b: Loci) {
    return a.particles === b.particles && OrderedSet.areEqual(a.indices, b.indices);
  }
  export function isLociEmpty(loci: Loci) {
    return OrderedSet.isEmpty(loci.indices);
  }
  export function lociSize(loci: Loci) {
    return OrderedSet.size(loci.indices);
  }
  /**
   * Remap a loci to a new `ParticleList`; indices outside the new range are dropped (a slot list's gaps too: a
   * particle kept by its slot across versions, gone with its page)
   */
  export function remapLoci(loci: Loci, particles: ParticleList): Loci {
    if (loci.particles === particles) return loci;
    const { count } = particles;
    if (count === 0) return Loci(particles, OrderedSet.Empty);
    const filtered: number[] = [];
    const paged = !!particles.pages;
    OrderedSet.forEach(loci.indices, (v) => {
      if (v < count && (!paged || ParticleSlots.findPage(particles, v) >= 0)) filtered.push(v);
    });
    return Loci(particles, OrderedSet.ofSortedArray(filtered));
  }

  const boundaryHelperCoarse = new BoundaryHelper('14');
  const boundaryHelperFine = new BoundaryHelper('98');
  function getBoundaryHelper(count: number) {
    return count > 10_000 ? boundaryHelperCoarse : boundaryHelperFine;
  }
  const _tmpPos = Vec3();

  export function getBoundingSphere(loci: Loci, boundingSphere?: Sphere3D): Sphere3D {
    const { particles, indices } = loci;
    const { coordinates, radii } = particles;
    if (isLociEmpty(loci)) {
      return boundingSphere ? Sphere3D.setZero(boundingSphere) : Sphere3D();
    } else if (lociSize(loci) === ParticleSlots.getLiveCount(particles)) {
      // all of them: the cached boundary (a slot list's of its pages', each cached on its list)
      const sphere = getBoundary(particles).sphere;
      return boundingSphere ? Sphere3D.copy(boundingSphere, sphere) : sphere;
    } else if (particles.pages) {
      // a slot list's: a page the loci holds whole by its list's cached boundary, the others' particles read from their
      // page's list (a gap none); not a search a particle: a view's loci are millions of particles
      const { pages } = particles.pages;
      const helper = getBoundaryHelper(OrderedSet.size(indices));
      helper.reset();
      for (let pass = 0; pass < 2; ++pass) {
        for (const page of pages) {
          const { start, list } = page;
          const lo = OrderedSet.findPredecessorIndex(indices, start);
          const hi = OrderedSet.findPredecessorIndex(indices, start + list.count);
          if (hi <= lo) continue;
          if (hi - lo === list.count) {
            const sphere = getBoundary(list).sphere;
            if (pass === 0) helper.includeSphere(sphere);
            else helper.radiusSphere(sphere);
            continue;
          }
          for (let k = lo; k < hi; ++k) {
            const i = OrderedSet.getAt(indices, k) - start;
            GridPositions.get(list.coordinates, i, _tmpPos);
            const r = list.radii ? ParticleValues.get(list.radii, i) : 0;
            if (pass === 0) helper.includePositionRadius(_tmpPos, r);
            else helper.radiusPositionRadius(_tmpPos, r);
          }
        }
        if (pass === 0) helper.finishedIncludeStep();
      }
      return helper.getSphere(boundingSphere);
    } else if (lociSize(loci) === particles.count) {
      const sphere = getBoundary(particles).sphere;
      return boundingSphere ? Sphere3D.copy(boundingSphere, sphere) : sphere;
    }
    const boundaryHelper = getBoundaryHelper(OrderedSet.size(indices));
    boundaryHelper.reset();
    if (radii) {
      OrderedSet.forEach(indices, (v) => {
        GridPositions.get(coordinates, v, _tmpPos);
        boundaryHelper.includePositionRadius(_tmpPos, ParticleValues.get(radii, v));
      });
      boundaryHelper.finishedIncludeStep();
      OrderedSet.forEach(indices, (v) => {
        GridPositions.get(coordinates, v, _tmpPos);
        boundaryHelper.radiusPositionRadius(_tmpPos, ParticleValues.get(radii, v));
      });
    } else {
      OrderedSet.forEach(indices, (v) => {
        GridPositions.get(coordinates, v, _tmpPos);
        boundaryHelper.includePosition(_tmpPos);
      });
      boundaryHelper.finishedIncludeStep();
      OrderedSet.forEach(indices, (v) => {
        GridPositions.get(coordinates, v, _tmpPos);
        boundaryHelper.radiusPosition(_tmpPos);
      });
    }
    return boundaryHelper.getSphere(boundingSphere);
  }

  export function getLabel(loci: Loci): string {
    const size = OrderedSet.size(loci.indices);
    if (size === 0) return 'None';
    if (size === 1) {
      const index = OrderedSet.start(loci.indices);
      return loci.particles.getParticleLabel(index);
    }
    return `${size.toLocaleString('en')} particles`;
  }

  /**
   * A copy of `particles` moved by `transform`: the coordinates transformed, the orientations turned by its rotation
   * (given where it turns and the list has none); the radii, the cell and the targets as they were. What is cached on
   * the list by position (the boundary, the transforms) is computed again, what is not kept (`keptByMove`). The copy
   * records the move (`movedFrom`). For rigid moves: the coordinates take any affine `transform`, but an orientation is
   * a unit quaternion, turned by the rotation part alone (`Mat4.extractRotation`), so a mirroring transform (negative
   * determinant) or a shear leaves the orientations turned, not mirrored or sheared, and a scale leaves the radii.
   */
  export function transform(particles: ParticleList, transform: Mat4): ParticleList {
    const { count, coordinates: positions } = particles;
    // float32 out (a grid's decoded: moves are rare)
    const src = GridPositions.toFloat(positions);
    const coordinates = new Float32Array(src.length);
    // `Vec3.transformMat4` and `Quat.multiply` inlined (millions of particles); the transform affine, w 1
    const [m0, m1, m2, , m4, m5, m6, , m8, m9, m10, , m12, m13, m14, m15] = transform;
    const w = m15 || 1;
    for (let o = 0, ol = count * 3; o < ol; o += 3) {
      const x = src[o],
        y = src[o + 1],
        z = src[o + 2];
      coordinates[o] = (m0 * x + m4 * y + m8 * z + m12) / w;
      coordinates[o + 1] = (m1 * x + m5 * y + m9 * z + m13) / w;
      coordinates[o + 2] = (m2 * x + m6 * y + m10 * z + m14) / w;
    }
    const turn = Quat.normalize(Quat(), Quat.fromMat4(Quat(), Mat4.extractRotation(Mat4(), transform)));
    let rotations = particles.rotations;
    // q and −q are the same turn
    if (rotations || 1 - Math.abs(turn[3]) > 1e-9) {
      const turned = new Float32Array(count * 4);
      const [ax, ay, az, aw] = turn;
      for (let o = 0, ol = count * 4; o < ol; o += 4) {
        // the identity where the list has none
        let bx = 0,
          by = 0,
          bz = 0,
          bw = 1;
        if (rotations) {
          // normalized: either kind (`ParticleRotations`), the turned float32
          bx = rotations[o];
          by = rotations[o + 1];
          bz = rotations[o + 2];
          bw = rotations[o + 3];
          const nb = Math.sqrt(bx * bx + by * by + bz * bz + bw * bw);
          if (nb > 0) {
            bx /= nb;
            by /= nb;
            bz /= nb;
            bw /= nb;
          } else {
            bw = 1;
          }
        }
        turned[o] = ax * bw + aw * bx + ay * bz - az * by;
        turned[o + 1] = ay * bw + aw * by + az * bx - ax * bz;
        turned[o + 2] = az * bw + aw * bz + ax * by - ay * bx;
        turned[o + 3] = aw * bw - ax * bx - ay * by - az * bz;
      }
      rotations = turned;
    }
    const from = particles.movedFrom;
    const customProperties = new CustomProperties();
    const propertyData = Object.create(null);
    const rigid = Math.abs(Mat4.getMaxScaleOnAxis(transform) - 1) < 1e-6;
    for (const descriptor of keptByMove) {
      const value = particles._propertyData[descriptor.name];
      if (value === undefined || (descriptor === ParticleSpacingDescriptor && !rigid)) continue;
      customProperties.add(descriptor);
      propertyData[descriptor.name] = value;
    }
    return {
      ...particles,
      coordinates,
      rotations,
      movedFrom: {
        particles: from?.particles ?? particles,
        transform: from ? Mat4.mul(Mat4(), transform, from.transform) : Mat4.copy(Mat4(), transform),
      },
      customProperties,
      _propertyData: propertyData,
    };
  }

  /**
   * A copy of `particles` with other poses: `coordinates` (and `rotations` where given) in place of theirs, all else
   * as it is; the copy records the line it is of (`posesFrom`). What is cached on the list by position (the boundary,
   * the transforms, the targets' bounds) is computed again, what is not kept (`keptByMove`, the spacing too: the
   * distances change). Not a move (`movedFrom` dropped): the instances' transforms built again, not moved
   */
  export function withPoses(
    particles: ParticleList,
    coordinates: ParticleList['coordinates'],
    rotations?: ParticleList['rotations'],
  ): ParticleList {
    if (GridPositions.count(coordinates) !== particles.count) throw new Error('withPoses: another count of particles');
    const customProperties = new CustomProperties();
    const propertyData = Object.create(null);
    for (const descriptor of keptByMove) {
      const value = particles._propertyData[descriptor.name];
      if (value === undefined || descriptor === ParticleSpacingDescriptor) continue;
      customProperties.add(descriptor);
      propertyData[descriptor.name] = value;
    }
    return {
      ...particles,
      coordinates,
      rotations: rotations ?? particles.rotations,
      posesFrom: particles.posesFrom ?? particles,
      movedFrom: undefined,
      customProperties,
      _propertyData: propertyData,
    };
  }

  /** Whether `a` and `b` are the same particles, posed alike or not (one line of `withPoses`, either its first) */
  export function isSamePoses(a: ParticleList, b: ParticleList) {
    return (a.posesFrom ?? a) === (b.posesFrom ?? b);
  }

  /**
   * What `transform` keeps of the list's cached data: what does not depend on positions (the target groups, the fiber
   * mask: the arrays they are of are shared), the spacing of a move without scale (distances kept); the release of a
   * move handle moves a list of millions, the groups alone cost 60 ms again
   */
  const keptByMove = [ParticleTargetGroupsDescriptor, FiberParticleMaskDescriptor, ParticleSpacingDescriptor];

  /**
   * The move from `from` to `to` where both are of one line of moves (`movedFrom`, either the first of it): `to` is
   * `from` moved by it, its particles the same; undefined otherwise
   */
  export function getMove(from: ParticleList, to: ParticleList): Mat4 | undefined {
    if ((from.movedFrom?.particles ?? from) !== (to.movedFrom?.particles ?? to)) return undefined;
    const t = to.movedFrom?.transform ?? Mat4.identity();
    const f = from.movedFrom?.transform;
    return f ? Mat4.mul(Mat4(), t, Mat4.invert(Mat4(), f)) : Mat4.copy(Mat4(), t);
  }

  /**
   * The particles `indices` of `particles` as a list of their own, in that order: their positions (on the same grid as
   * they are, the words copied: exact), rotations, radii, targets, entities, compartments (runs where fewer than half
   * the particles) and attributes; the targets, their info, entities, compartments and the cell the same objects (one
   * for all subsets: what is keyed by a target is one); a particle's label its own (by its index in `particles`); no
   * fibers (a fiber's particles may be in other subsets); `keys` where given (`ParticleList.keys`). A tile's particles
   * (`particle-tiles.ts`)
   */
  export function subset(
    particles: ParticleList,
    indices: Int32Array,
    options: { readonly label?: string; readonly keys?: Int32Array } = {},
  ): ParticleList {
    const n = indices.length;
    let attributes: Map<string, ParticleAttribute> | undefined;
    if (particles.attributes) {
      attributes = new Map();
      for (const [name, a] of particles.attributes) {
        const values = new Float32Array(n);
        for (let i = 0; i < n; ++i) values[i] = a.column.value(indices[i]);
        attributes.set(name, { ...a, column: Column.ofFloatArray(values) });
      }
    }
    return {
      entryId: particles.entryId,
      label: options.label ?? particles.label,
      count: n,
      cell: particles.cell,
      keys: options.keys,
      targets: gatherValues(particles.targets, indices, Int32Array),
      targetInfo: particles.targetInfo,
      targetMapping: particles.targetMapping,
      compartments: particles.compartments && gatherValues(particles.compartments, indices, Int32Array),
      compartmentInfo: particles.compartmentInfo,
      entities: particles.entities && gatherValues(particles.entities, indices, Int32Array),
      entityInfo: particles.entityInfo,
      coordinates: gatherPositions(particles.coordinates, indices),
      rotations: particles.rotations && gatherRotations(particles.rotations, indices),
      radii: particles.radii && gatherValues(particles.radii, indices, Float32Array),
      attributes,
      getParticleLabel: (index) => particles.getParticleLabel(indices[index]),
      sourceData: particles.sourceData,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
  }

  /**
   * Particles `[first, first + count)` of `particles` as a list of views of its arrays (no copy, no pass over them:
   * the positions' words and the rotations as stored, the values' runs over the range, `ParticleValues.slice`), its
   * targets the same objects; a tile of a list ordered by tile (`ParticleTiles`)
   */
  export function range(
    particles: ParticleList,
    first: number,
    count: number,
    options: { readonly label?: string; readonly keys?: Int32Array } = {},
  ): ParticleList {
    const end = first + count;
    let attributes: Map<string, ParticleAttribute> | undefined;
    if (particles.attributes) {
      attributes = new Map();
      for (const [name, a] of particles.attributes) {
        attributes.set(name, { ...a, column: Column.window(a.column, first, end) });
      }
    }
    const { coordinates } = particles;
    return {
      entryId: particles.entryId,
      label: options.label ?? particles.label,
      count,
      cell: particles.cell,
      keys: options.keys ?? particles.keys?.subarray(first, end),
      targets: ParticleValues.slice(particles.targets, first, count),
      targetInfo: particles.targetInfo,
      targetMapping: particles.targetMapping,
      compartments: particles.compartments && ParticleValues.slice(particles.compartments, first, count),
      compartmentInfo: particles.compartmentInfo,
      entities: particles.entities && ParticleValues.slice(particles.entities, first, count),
      entityInfo: particles.entityInfo,
      coordinates: GridPositions.is(coordinates)
        ? { ...coordinates, data: coordinates.data.subarray(first * 2, end * 2) }
        : coordinates.subarray(first * 3, end * 3),
      rotations: particles.rotations?.subarray(first * 4, end * 4) as ParticleRotations | undefined,
      radii: particles.radii && ParticleValues.slice(particles.radii, first, count),
      attributes,
      getParticleLabel: (index) => particles.getParticleLabel(first + index),
      sourceData: particles.sourceData,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
  }

  /** `values` into `out` (a typed array of their count) */
  function copyInto<A extends Int32Array | Float32Array>(out: A, values: readonly number[]): A {
    out.set(values);
    return out;
  }

  /** `coordinates` of `indices`, in their order: a grid's words copied (exact), else float32's */
  function gatherPositions(coordinates: Positions, indices: Int32Array): Positions {
    const n = indices.length;
    if (GridPositions.is(coordinates)) {
      const data = new Uint32Array(n * 2);
      for (let i = 0; i < n; ++i) {
        data[i * 2] = coordinates.data[indices[i] * 2];
        data[i * 2 + 1] = coordinates.data[indices[i] * 2 + 1];
      }
      return { kind: 'grid', data, origin: coordinates.origin, step: coordinates.step };
    }
    const out = new Float32Array(n * 3);
    for (let i = 0; i < n; ++i) out.set(coordinates.subarray(indices[i] * 3, indices[i] * 3 + 3), i * 3);
    return out;
  }

  /** `rotations` of `indices`, in their order, as they are stored (float32 or snorm16) */
  function gatherRotations(rotations: ParticleRotations, indices: Int32Array): ParticleRotations {
    const n = indices.length;
    const out = new (rotations.constructor as { new (n: number): ParticleRotations })(n * 4);
    // element by element: a `subarray` a particle is an allocation each (31M, a list's whole order)
    for (let i = 0; i < n; ++i) {
      const o = indices[i] * 4;
      const t = i * 4;
      out[t] = rotations[o];
      out[t + 1] = rotations[o + 1];
      out[t + 2] = rotations[o + 2];
      out[t + 3] = rotations[o + 3];
    }
    return out;
  }

  /** `values` of `indices`, in their order: as runs where they are fewer than half the values, else an array */
  function gatherValues<A extends Int32Array | Float32Array>(
    values: ParticleValues<A>,
    indices: Int32Array,
    ctor: { new (n: number): A },
  ): ParticleValues<A> {
    const n = indices.length;
    // runs gathered as runs (a value a run, not an array of all: a list's whole order, 31M particles a value kind)
    if (ParticleValues.isRuns(values)) {
      const offsets: number[] = [];
      const runValues: number[] = [];
      for (let i = 0; i < n; ++i) {
        const v = ParticleValues.get(values, indices[i]);
        if (i === 0 || v !== runValues[runValues.length - 1]) {
          offsets.push(i);
          runValues.push(v);
        }
      }
      if (runValues.length * 2 <= n) {
        offsets.push(n);
        return ParticleValues.ofRuns(Int32Array.from(offsets), copyInto(new ctor(runValues.length), runValues));
      }
    }
    const array = new ctor(n);
    for (let i = 0; i < n; ++i) array[i] = ParticleValues.get(values, indices[i]);
    return ParticleValues.compact(array);
  }

  /**
   * Return a copy of `particles` with the given per-target reference objects attached; of the same line of moves
   * (`movedFrom`: the coordinates the same), so a list moved under it gives its move
   */
  export function withTargets(
    particles: ParticleList,
    targetMapping: ReadonlyMap<number, ParticleTarget>,
  ): ParticleList {
    return { ...particles, targetMapping, movedFrom: particles.movedFrom ?? { particles, transform: Mat4.identity() } };
  }

  /**
   * Per run of a list's targets (`ParticleRuns`: a pack's, sorted by target) its particles' bounds
   * (`ParticleRangeBoundsValues` values a run): a pack's where it has them (`setTargetBounds`), else computed for all
   * runs at the first use (`getRangeBounds`); by position, so not kept by a move (`keptByMove`)
   */
  const TargetBoundsDescriptor = CustomPropertyDescriptor({ name: 'particle-target-bounds' });
  export function setTargetBounds(particles: ParticleList, values: Float64Array): void {
    particles.customProperties.add(TargetBoundsDescriptor);
    particles._propertyData[TargetBoundsDescriptor.name] = values;
  }

  /**
   * The bounds of positions `[start, end)` of `coordinates` into `out` from `offset` (`ParticleRangeBoundsValues`): the
   * box, then the largest distance from its center (a second pass); zeros for none
   */
  export function computeRangeBounds(
    coordinates: Positions,
    start: number,
    end: number,
    out: Float64Array,
    offset = 0,
  ): void {
    if (end <= start) {
      out.fill(0, offset, offset + ParticleRangeBoundsValues);
      return;
    }
    if (GridPositions.is(coordinates)) {
      computeGridRangeBounds(coordinates, start, end, out, offset);
      return;
    }
    const p = coordinates;
    let minX = Infinity,
      minY = Infinity,
      minZ = Infinity,
      maxX = -Infinity,
      maxY = -Infinity,
      maxZ = -Infinity;
    for (let o = start * 3, ol = end * 3; o < ol; o += 3) {
      const x = p[o],
        y = p[o + 1],
        z = p[o + 2];
      if (x < minX) minX = x;
      if (y < minY) minY = y;
      if (z < minZ) minZ = z;
      if (x > maxX) maxX = x;
      if (y > maxY) maxY = y;
      if (z > maxZ) maxZ = z;
    }
    const cx = (minX + maxX) / 2,
      cy = (minY + maxY) / 2,
      cz = (minZ + maxZ) / 2;
    let d2 = 0;
    for (let o = start * 3, ol = end * 3; o < ol; o += 3) {
      const dx = p[o] - cx,
        dy = p[o + 1] - cy,
        dz = p[o + 2] - cz;
      const d = dx * dx + dy * dy + dz * dz;
      if (d > d2) d2 = d;
    }
    setRangeBounds(out, offset, minX, minY, minZ, maxX, maxY, maxZ, Math.sqrt(d2));
  }

  /**
   * `computeRangeBounds` on a grid: in its steps (integers, the center doubled: exact), decoded once at the end; the
   * grid's words read inline, not an element at a time (`GridPositions.get`: 31.5 M particles of a pack)
   */
  function computeGridRangeBounds(g: GridPositions, start: number, end: number, out: Float64Array, offset: number) {
    const { data, origin, step } = g;
    let minX = Infinity,
      minY = Infinity,
      minZ = Infinity,
      maxX = -Infinity,
      maxY = -Infinity,
      maxZ = -Infinity;
    for (let o = start * 2, ol = end * 2; o < ol; o += 2) {
      const a = data[o],
        b = data[o + 1];
      const x = a & 0x1fffff,
        y = (a >>> 21) | ((b & 0x3ff) << 11),
        z = (b >>> 10) & 0x1fffff;
      if (x < minX) minX = x;
      if (y < minY) minY = y;
      if (z < minZ) minZ = z;
      if (x > maxX) maxX = x;
      if (y > maxY) maxY = y;
      if (z > maxZ) maxZ = z;
    }
    // twice the center: integers, exact
    const cx = minX + maxX,
      cy = minY + maxY,
      cz = minZ + maxZ;
    let d2 = 0;
    for (let o = start * 2, ol = end * 2; o < ol; o += 2) {
      const a = data[o],
        b = data[o + 1];
      const dx = 2 * (a & 0x1fffff) - cx,
        dy = 2 * ((a >>> 21) | ((b & 0x3ff) << 11)) - cy,
        dz = 2 * ((b >>> 10) & 0x1fffff) - cz;
      const d = dx * dx + dy * dy + dz * dz;
      if (d > d2) d2 = d;
    }
    const [ox, oy, oz] = origin;
    setRangeBounds(
      out,
      offset,
      ox + step * minX,
      oy + step * minY,
      oz + step * minZ,
      ox + step * maxX,
      oy + step * maxY,
      oz + step * maxZ,
      (step * Math.sqrt(d2)) / 2,
    );
  }

  function setRangeBounds(
    out: Float64Array,
    o: number,
    minX: number,
    minY: number,
    minZ: number,
    maxX: number,
    maxY: number,
    maxZ: number,
    distance: number,
  ) {
    out[o] = minX;
    out[o + 1] = minY;
    out[o + 2] = minZ;
    out[o + 3] = maxX;
    out[o + 4] = maxY;
    out[o + 5] = maxZ;
    out[o + 6] = distance;
  }

  /**
   * The bounds of particles `[start, end)` where they are one run of the list's targets (a target's particles of a
   * pack, `ParticleRuns`): the pack's, else all runs' computed once (`TargetBoundsDescriptor`); undefined otherwise
   */
  export function getRangeBounds(particles: ParticleList, start: number, end: number): ParticleRangeBounds | undefined {
    const { targets } = particles;
    if (!ParticleValues.isRuns(targets) || end <= start) return undefined;
    const r = ParticleValues.runOf(targets, start);
    const { offsets } = targets;
    if (offsets[r] !== start || offsets[r + 1] !== end) return undefined;
    const values = getRunBounds(particles, targets);
    const o = r * ParticleRangeBoundsValues;
    return {
      min: Vec3.create(values[o], values[o + 1], values[o + 2]),
      max: Vec3.create(values[o + 3], values[o + 4], values[o + 5]),
      distance: values[o + 6],
    };
  }

  /** the bounds of each run of `targets` (`particles`' runs): the pack's, else all computed once */
  function getRunBounds(particles: ParticleList, targets: ParticleRuns) {
    let values = particles._propertyData[TargetBoundsDescriptor.name] as Float64Array | undefined;
    if (!values) {
      const { offsets } = targets;
      const runs = offsets.length - 1;
      values = new Float64Array(runs * ParticleRangeBoundsValues);
      for (let k = 0; k < runs; ++k) {
        computeRangeBounds(particles.coordinates, offsets[k], offsets[k + 1], values, k * ParticleRangeBoundsValues);
      }
      setTargetBounds(particles, values);
    }
    return values;
  }

  /**
   * The bounds of particles `indices` where they are whole runs of the list's targets (`getRangeBounds`' of one; a
   * tile-ordered pack's target a run a tile it is in: its runs' united, no pass over its particles); undefined otherwise
   */
  export function getIndicesBounds(
    particles: ParticleList,
    indices: OrderedSet<number>,
  ): ParticleRangeBounds | undefined {
    const { targets } = particles;
    const n = OrderedSet.size(indices);
    if (!ParticleValues.isRuns(targets) || n === 0) return undefined;
    const first = OrderedSet.getAt(indices, 0);
    const last = OrderedSet.getAt(indices, n - 1);
    if (last - first + 1 === n) return getRangeBounds(particles, first, last + 1);
    const { offsets } = targets;
    const values = getRunBounds(particles, targets);
    const min = Vec3.create(Infinity, Infinity, Infinity);
    const max = Vec3.create(-Infinity, -Infinity, -Infinity);
    const runs: number[] = [];
    let whole = true;
    ParticleValues.forEachRunOf(targets, indices, (r, count) => {
      if (count !== offsets[r + 1] - offsets[r]) {
        whole = false;
        return false;
      }
      runs.push(r);
      const o = r * ParticleRangeBoundsValues;
      for (let a = 0; a < 3; ++a) {
        if (values[o + a] < min[a]) min[a] = values[o + a];
        if (values[o + 3 + a] > max[a]) max[a] = values[o + 3 + a];
      }
    });
    if (!whole) return undefined;
    // the largest distance from the united box's center: a run's from its box's center, that center's from it
    const cx = (min[0] + max[0]) / 2;
    const cy = (min[1] + max[1]) / 2;
    const cz = (min[2] + max[2]) / 2;
    let distance = 0;
    for (const r of runs) {
      const o = r * ParticleRangeBoundsValues;
      const dx = (values[o] + values[o + 3]) / 2 - cx;
      const dy = (values[o + 1] + values[o + 4]) / 2 - cy;
      const dz = (values[o + 2] + values[o + 5]) / 2 - cz;
      distance = Math.max(distance, Math.sqrt(dx * dx + dy * dy + dz * dz) + values[o + 6]);
    }
    return { min, max, distance };
  }

  export const BoundaryDescriptor: CustomPropertyDescriptor<Boundary> = CustomPropertyDescriptor({
    name: 'particle-boundary',
  });
  export function setBoundary(particles: ParticleList, boundary: Boundary): void {
    particles.customProperties.add(BoundaryDescriptor);
    particles._propertyData[BoundaryDescriptor.name] = boundary;
  }
  export function getBoundary(particles: ParticleList): Boundary {
    if (!particles._propertyData[BoundaryDescriptor.name] && particles.pages) {
      // a slot list's: of its pages' (each cached on its list: a version's by its pages, not its particles), all of
      // them computed first: a page's own computing resets the helpers this one shares
      const boundaries = particles.pages.pages.map((page) => getBoundary(page.list));
      const helper = getBoundaryHelper(boundaries.length * 1000);
      helper.reset();
      const box = Box3D.setEmpty(Box3D());
      for (const b of boundaries) {
        helper.includeSphere(b.sphere);
        Box3D.add(box, b.box.min);
        Box3D.add(box, b.box.max);
      }
      helper.finishedIncludeStep();
      for (const b of boundaries) helper.radiusSphere(b.sphere);

      setBoundary(particles, { box, sphere: helper.getSphere() });
    }
    if (!particles._propertyData[BoundaryDescriptor.name]) {
      // Compute boundary from particle positions and radii, and store it in the particle list for later retrieval.
      // loop over positions and radii to compute the boundary
      const { count, coordinates, radii } = particles;
      const boundaryHelper = getBoundaryHelper(count);
      boundaryHelper.reset();
      if (radii) {
        // by the radii's runs: a value a run, not a particle (`ParticleValues`)
        ParticleValues.forEachRun(radii, count, (start, end, radius) => {
          for (let i = start; i < end; i++) {
            GridPositions.get(coordinates, i, _tmpPos);
            boundaryHelper.includePositionRadius(_tmpPos, radius);
          }
        });
        boundaryHelper.finishedIncludeStep();
        ParticleValues.forEachRun(radii, count, (start, end, radius) => {
          for (let i = start; i < end; i++) {
            GridPositions.get(coordinates, i, _tmpPos);
            boundaryHelper.radiusPositionRadius(_tmpPos, radius);
          }
        });
      } else {
        for (let i = 0; i < count; i++) {
          GridPositions.get(coordinates, i, _tmpPos);
          boundaryHelper.includePosition(_tmpPos);
        }
        boundaryHelper.finishedIncludeStep();
        for (let i = 0; i < count; i++) {
          GridPositions.get(coordinates, i, _tmpPos);
          boundaryHelper.radiusPosition(_tmpPos);
        }
      }
      setBoundary(particles, { box: boundaryHelper.getBox(), sphere: boundaryHelper.getSphere() });
    }
    return particles._propertyData[BoundaryDescriptor.name];
  }
}
