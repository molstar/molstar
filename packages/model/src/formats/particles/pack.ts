/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '@molstar/core/data/db';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import type { ModelFormat } from '../format.js';
import { CustomProperties } from '../../model/custom-property.js';
import {
  Particle,
  ParticleRangeBoundsValues,
  ParticleRotations,
  ParticleSites,
  ParticleTarget,
  ParticleValues,
} from '../../model/particles/particle-list.js';
import {
  ParticleTileCatalog,
  type ParticleTargetMeasures,
  ParticleTileDirectoryBox,
  ParticleTileDirectoryRange,
  ParticleTileGrid,
  ParticleTiles,
  type ParticleTileDirectory,
} from '../../model/particles/particle-tiles.js';
import { mortonOrder3d } from '@molstar/core/data/util/hash-functions';
import { radixSortIndices } from '@molstar/core/data/util/sort';
import { ChunkedArray } from '@molstar/core/data/util/chunked-array';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { ParticleDensityGrid, ParticleTileOpaque } from '../../model/particles/particle-tile-opaque.js';
import { ParticleNesting } from '../../model/particles/particle-nesting.js';
import type {
  ParticleAttribute,
  ParticleCompartmentInfo,
  ParticleEntityInfo,
  ParticleList,
  ParticleTargetInfo,
} from '../../model/particles/particle-list.js';

/**
 * A particle pack (`formats.md` §3): a `ParticleList` and its targets' sites, nothing else, as a header, a JSON
 * document of the small tables and raw arrays read as views of the file's bytes:
 *
 * - offset 0, 64 bytes: `MSPK`, the version (u32), the byte order mark (u32 `0x01020304` as written), the JSON's
 *   offset and length (u64 each), zeros
 * - the blocks, each at a multiple of 256 bytes from the file's start (WebGPU's copy alignment, a typed array's many
 *   times), little-endian, no encoding
 * - the JSON (UTF-8, `ParticlePackDocument`) after them: the list's label, entry id and cell, its entities,
 *   compartments and attributes, a row a target, the blocks
 *
 * The particles are sorted by target, a target's particles a range (`particleOffset`, `particleCount`), or, tile-ordered
 * (version 6, `tiles`), by tile of a grid over their box, the tiles in Morton order, a tile's by target: a target's
 * particles runs (`target.runs`), the tiles a directory (`tile.ranges`, `tile.boxes`: `ParticleTileDirectory`), a tile
 * set made from it without a pass over the particles, a tile a range of every block; by target, a density grid over
 * their box where written (`density`: what hides what is behind as the whole is drawn, the tiles' opacities a
 * tile-ordered pack's). A particle's id, entity and
 * compartment are not stored per particle, nor a key (a particle's key is its index in the pack). A target is `sites`
 * (positions and radii, chains as runs), a template (version 7, `template`: a list of `templates`, a document of its own
 * whose blocks are of the same file, placed whole at each of the target's particles) or has none.
 */
export const ParticlePackMagic = 'MSPK';
/** The version written: 7 where the list has templates (a reader before it refuses what it would draw wrong), else 6 */
export const ParticlePackVersion = 7;
const ParticlePackVersionWithoutTemplates = 6;
/**
 * The versions read: 2's rotations float32, 2's and 3's positions float32 (a block's type tells); before 5 no
 * `siteElements` (the sites' elements not known); before 6 sorted by target; before 7 no templates
 */
const ReadVersions = [2, 3, 4, 5, 6, 7];
const HeaderBytes = 64;
const radiusBits = new Float32Array(1);
const radiusBitsU32 = new Uint32Array(radiusBits.buffer);
const BlockAlignment = 256;
const ByteOrderMark = 0x01020304;

/**
 * `snorm16`: int16 values of `v * 32767` (`ParticleRotations`), read as `Int16Array`; `grid21`: positions on the
 * document's `positionGrid`, two u32s an element (`GridPositions`), read as `Uint32Array`
 */
export type ParticlePackBlockType = 'float32' | 'int32' | 'int16' | 'snorm16' | 'grid21';

/** An array of the pack: `count` items of `components` values of `type`, at `byteOffset` from the file's start */
export interface ParticlePackBlock {
  readonly name: string;
  readonly type: ParticlePackBlockType;
  readonly components: number;
  readonly count: number;
  readonly byteOffset: number;
  readonly byteLength: number;
}

/**
 * A target's sites: those `[offset, offset + count)` of the `site` block (x, y, z, type as int16 a site: a position
 * `origin + q * scale`, a radius `siteRadii[type]`, an element `siteElements[type]`, the type as uint16), their chains
 * as runs
 */
export interface ParticlePackSites {
  readonly offset: number;
  readonly count: number;
  readonly origin: readonly [number, number, number];
  /** Å a unit: the target's extent over 32767 (0.01 Å for ±327 Å) */
  readonly scale: number;
  readonly label?: string;
  /** the point a particle's position places (`ParticleTarget.getCenter`) */
  readonly center: readonly [number, number, number];
  /** a chain a run of consecutive sites (`ParticleSiteChains`) */
  readonly chains?: { readonly labels: readonly string[]; readonly counts: readonly number[] };
  /**
   * the quantized positions' extremes of each chain (without chains, of all sites): min x, y, z, max x, y, z
   * (`ParticleSitesQuantized.extremes`); a pack without them (older, the same version) has its boxes computed where read
   */
  readonly extremes?: readonly number[];
}

export interface ParticlePackTarget {
  /** the target's id (`ParticleList.targets`) */
  readonly id: number;
  readonly name?: string;
  /** an index of `entities`, -1 for none */
  readonly entity: number;
  /** an index of `compartments`, -1 for none */
  readonly compartment: number;
  /** its particles' first, -1 tile-ordered (its particles runs of `target.runs`) */
  readonly particleOffset: number;
  readonly particleCount: number;
  /** its particles' radius (a target's particles' are the same; the largest where a list's differ in rounding) */
  readonly radius?: number;
  readonly sites?: ParticlePackSites;
  /**
   * how far its sites reach from its center and its molecule's volume (version 6: `ParticleTileCatalog.getTargetRadius`,
   * `getTargetVolume` of its sites as written): a tile set's catalog without a pass over the sites; computed where none
   */
  readonly reach?: number;
  readonly volume?: number;
  /**
   * its particles' positions' box and largest distance from its center (`Particle.getRangeBounds`), rounded outward
   * to `BoundsSteps`: a target placed without a pass over its particles; packs without them computed at read
   */
  readonly bounds?: {
    readonly min: readonly [number, number, number];
    readonly max: readonly [number, number, number];
    readonly distance: number;
  };
  /**
   * a template (version 7): its particles each place the list of `templates[index]` whole, about `center` (the point a
   * particle's position places: `ParticleTarget.getCenter`, the template's boundary's where written)
   */
  readonly template?: { readonly index: number; readonly center: readonly [number, number, number] };
}

/** Steps an Å the targets' bounds are rounded outward to (`ParticlePackTarget.bounds`): 0.001 Å */
const BoundsSteps = 1000;

export interface ParticlePackDocument {
  readonly label?: string;
  readonly entryId?: string;
  readonly count: number;
  readonly cell?: { readonly size: readonly number[]; readonly anglesInRadians: readonly number[] };
  readonly entities: readonly { readonly index: number; readonly name: string; readonly function?: string }[];
  readonly compartments: readonly { readonly index: number; readonly name: string }[];
  readonly targets: readonly ParticlePackTarget[];
  /** numeric attributes, a `attribute.<name>` block each */
  readonly attributes: readonly {
    readonly name: string;
    readonly label: string;
    readonly min: number;
    readonly max: number;
  }[];
  readonly fiberCount?: number;
  /** the sites' radii by type (`ParticlePackSites`) */
  readonly siteRadii: readonly number[];
  /**
   * the sites' atomic numbers by type (version 5), where every target's sites have elements: a type a distinct radius
   * and element
   */
  readonly siteElements?: readonly number[];
  /** the grid of the `position` block where `grid21` (version 4) */
  readonly positionGrid?: { readonly origin: readonly number[]; readonly step: number };
  /**
   * tile-ordered (version 6): the tiles' grid, `dimensions` cells an axis (a power of two) over the particles' box
   * `min`-`max`; the blocks `target.runs` (first particle, target id), `tile.ranges` (cell id, first particle,
   * count) and `tile.boxes` (the particles' centers' box)
   */
  readonly tiles?: {
    readonly dimensions: number;
    readonly min: readonly number[];
    readonly max: readonly number[];
    /**
     * the tiles' opacities' blocks an axis, where written (the block `tile.opacity`: a tile's `blocks`³ optical
     * densities, `ParticleTileOpaque.computeOpacity` of its particles as written, its targets' volumes)
     */
    readonly opacityBlocks?: number;
    /**
     * the coarse levels written (`WriteParticlePackOptions.levels`): level `L` (1 to `levels`) a node a cell of the
     * grid of `dimensions >> L` cells an axis, of each target every `8^L`-th of its particles in its `2^L`³ tiles (at
     * least one; in Morton order of a fine grid over them: spread evenly), copies of theirs in the blocks `node.*`
     * (`node.position`, `node.rotation`, `node.attribute.<name>`; `node.keys` a node particle's index in the pack,
     * `node.runs` its targets' runs, `node.scales` a run's scale: the cube root of its target's particles over those
     * kept, so each target's volume kept; `node.ranges` a node's level, cell, first and count, `node.boxes` its
     * centers' box); a reader of none ignores them, one of no `node.scales` draws a node `2^L` larger
     */
    readonly levels?: number;
  };
  /**
   * by target, where written (version 6): the particles' optical densities on a grid of `dimensions` blocks over their
   * box `min`-`max` (the block `density`, x fastest: `ParticleDensityGrid`, of the particles as written, the targets'
   * volumes); a tile-ordered pack's are its tiles' opacities
   */
  readonly density?: {
    readonly dimensions: readonly number[];
    readonly min: readonly number[];
    readonly max: readonly number[];
  };
  /**
   * the lists the targets' `template`s place (version 7): a document each, its blocks of the same file, by target, no
   * templates, tiles or density of its own
   */
  readonly templates?: readonly ParticlePackDocument[];
  /**
   * a template's (version 7, optional: a reader without it computes them): its particles by target (tile-ordered:
   * by tile, then by target, `WriteParticlePackOptions.templateTiles`), each target's in the sample order of its levels
   * (`Particle.getSampleOrder`, `setSampleOrdered`), and its boundary of the positions and radii as written, rounded
   * outward (`Particle.setBoundary`). Why: the yeast cell's levels sorted its 31.5 M
   * particles at every load (1.9 s), its sphere a pass over them (0.5 s)
   */
  readonly sampleOrdered?: boolean;
  /**
   * a template's coarse levels written (version 7, optional: a reader without them samples them,
   * `ParticleNesting.getTemplateLevels`): level `L` (1 to `templateLevels`) of each target every `8^L`-th of its
   * particles as written (at least one), a run a target (as first seen), in the blocks `level.*` (`level.position`,
   * `level.rotation`, `level.keys` a level particle's index in the template, `level.runs` its targets' runs,
   * `level.scales` a run's scale: the square root of its target's particles over those kept, the area a far copy covers
   * kept, `level.ranges` a level's level, first and count). Why: sampled at load, the yeast cell's 461 ms (its gathers
   * scattered tile-ordered). Not by tile (a tile set's `node.*`: a run a target a node, at least one a target a node):
   * the yeast's levels 11-42% larger and 724 k runs, not 5,639 (its frames rotating 1.70 to 1.97 ms), the isg cell's 2
   * levels by its 4 tiles an axis, not 3 (its far frames 2.40 to 3.93 ms)
   */
  readonly templateLevels?: number;
  readonly boundary?: {
    readonly min: readonly [number, number, number];
    readonly max: readonly [number, number, number];
    readonly center: readonly [number, number, number];
    readonly radius: number;
  };
  readonly blocks: readonly ParticlePackBlock[];
}

export interface WriteParticlePackOptions {
  /**
   * tile-ordered, the coarse levels (`ParticlePackDocument.tiles.levels`): a tile set's far tiles drawn as samples of
   * their particles, read as a few ranges; about a seventh more particles written (copies). False by default
   */
  readonly levels?: boolean;
  /**
   * tile-ordered (`ParticlePackDocument.tiles`): the particles by tile of a grid of `tiles` cells an axis (a power of
   * two, at most 1024) over their box, in Morton order of the tiles, a tile's by target; none: by target
   */
  readonly tiles?: number;
  /** tile-ordered, the tiles' opacities' blocks an axis (`ParticleTileOpaque`): 4 by default, 0 none written */
  readonly opacityBlocks?: number;
  /**
   * by target, the density grid's blocks along the box's longest side (`ParticlePackDocument.density`): 64 by default,
   * 0 none written; tile-ordered with the tiles' opacities none (they are its)
   */
  readonly densityBlocks?: number;
  /**
   * a list with templates, its templates' tiles an axis (tile-ordered as `tiles`, a tile's by target, each target's in
   * its sample order): by default about a store page's particles a tile (`templateTilesOf`), 0 or 1 by target alone
   */
  readonly templateTiles?: number;
  /** a template's coarse levels written (`ParticlePackDocument.templateLevels`; set for a list's templates) */
  readonly templateLevels?: number;
}

/** a template's particles a tile by default (a store page's: its pages then about a tile's, by space) */
const TemplateTileParticles = 8192;

/**
 * A template's tiles an axis by default: a power of two whose cube is about its particles' store pages (the yeast
 * cell's 31.5 M: 16, the isg cell's 460 k: 4). Why by tile: a page by target holds a target spread over the template
 * (its box the template's), every page reached by a near view; by tile its particles near each other
 */
export function templateTilesOf(count: number) {
  const tiles = Math.cbrt(count / TemplateTileParticles);
  if (tiles < 1.5) return 1;
  return Math.min(1024, 2 ** Math.round(Math.log2(tiles)));
}

/** A pack read: its document and the bytes its blocks are views of */
export interface ParticlePack {
  readonly document: ParticlePackDocument;
  readonly buffer: ArrayBuffer;
  /** the file's first byte in `buffer` */
  readonly byteOffset: number;
  /** its blocks read apart, not views of `buffer` (a pack opened by range: `ParticlePackBlocks`) */
  readonly blocks?: ParticlePackBlocks;
}

export type ParticlePackFormat = ModelFormat<ParticlePack>;

export namespace ParticlePackFormat {
  export function is(x?: ModelFormat): x is ParticlePackFormat {
    return x?.kind === 'particle-pack';
  }

  export function create(pack: ParticlePack): ParticlePackFormat {
    return { kind: 'particle-pack', name: 'particle-pack', data: pack };
  }
}

/** Whether `data` starts as a pack does (its magic) */
export function isParticlePack(data: Uint8Array) {
  return data.length >= HeaderBytes && String.fromCharCode(data[0], data[1], data[2], data[3]) === ParticlePackMagic;
}

/** a pack's header (its first `ParticlePackHeaderBytes`): its version and where its document is */
export interface ParticlePackHeader {
  readonly version: number;
  readonly jsonOffset: number;
  readonly jsonLength: number;
}

/** bytes of a pack's header */
export const ParticlePackHeaderBytes = HeaderBytes;

/** The header of a pack starting `data` (at least its header's bytes), checked: its magic, version and byte order */
export function readParticlePackHeader(data: Uint8Array): ParticlePackHeader {
  if (!isParticlePack(data)) throw new Error('not a particle pack (no MSPK magic)');
  if (new Uint8Array(new Uint16Array([1]).buffer)[0] !== 1)
    throw new Error('particle packs read on little-endian hosts only');
  const view = new DataView(data.buffer, data.byteOffset, HeaderBytes);
  const version = view.getUint32(4, true);
  if (!ReadVersions.includes(version))
    throw new Error(`particle pack version ${version}, ${ReadVersions.join(' and ')} read`);
  if (view.getUint32(8, true) !== ByteOrderMark) throw new Error('particle pack of another byte order');
  return { version, jsonOffset: Number(view.getBigUint64(12, true)), jsonLength: Number(view.getBigUint64(20, true)) };
}

/** bytes a value of a block's type */
const BlockTypeBytes: Record<ParticlePackBlockType, number> = { float32: 4, int32: 4, int16: 2, snorm16: 2, grid21: 4 };

/** blocks of a value a particle (`count` the pack's): their count checked against it */
function isPerParticleBlock(name: string) {
  return name === 'position' || name === 'rotation' || name.startsWith('attribute.');
}

/**
 * The document of a pack (its JSON's bytes), its blocks checked within `byteLength` (the file's), aligned, and of the
 * size their items say (`count` items of `components` values of `type`: what a view of them reads), those of a value a
 * particle of the pack's count. Why: a block whose sizes disagree would be read into the next, or fail unnamed
 */
export function parseParticlePackDocument(json: Uint8Array, byteLength: number): ParticlePackDocument {
  const document = JSON.parse(new TextDecoder().decode(json)) as ParticlePackDocument;
  checkParticlePackDocument(document, byteLength, false);
  return document;
}

/** `document`'s blocks and its templates' checked (`parseParticlePackDocument`), a template's none of its own */
function checkParticlePackDocument(document: ParticlePackDocument, byteLength: number, isTemplate: boolean) {
  if (document.templates?.length) {
    if (isTemplate) throw new Error('particle pack template with templates of its own');
    for (const t of document.templates) checkParticlePackDocument(t, byteLength, true);
    for (const t of document.targets) {
      if (t.template && !document.templates[t.template.index])
        throw new Error(`particle pack target ${t.id} of template ${t.template.index}, not a template`);
    }
  } else if (document.targets.some((t) => t.template)) {
    throw new Error('particle pack target of a template, but the pack has none');
  }
  for (const b of document.blocks) {
    if (b.byteOffset % BlockAlignment !== 0) throw new Error(`particle pack block '${b.name}' not aligned`);
    if (b.byteOffset + b.byteLength > byteLength) throw new Error(`particle pack truncated (block '${b.name}')`);
    const bytes = BlockTypeBytes[b.type];
    if (bytes === undefined) throw new Error(`particle pack block '${b.name}' of an unknown type '${b.type}'`);
    if (b.byteLength !== b.count * b.components * bytes) {
      throw new Error(
        `particle pack block '${b.name}' of ${b.byteLength} bytes, not its ${b.count} × ${b.components} ${b.type}`,
      );
    }
    if (isPerParticleBlock(b.name) && b.count !== document.count) {
      throw new Error(
        `particle pack block '${b.name}' of ${b.count} items, not the pack's ${document.count} particles`,
      );
    }
  }
}

/** bytes of a pack at an offset not keeping its blocks aligned copied at most (`readParticlePack`) */
export const ParticlePackMaxCopyBytes = 16 * 1024 * 1024;

/**
 * The pack in `data`: its header and document read, its blocks left as they are (views made where read, no copy).
 * `data` at an offset not a multiple of 4 (its blocks' typed arrays not aligned) copied where small (at most
 * `maxCopyBytes`, `ParticlePackMaxCopyBytes` by default), else refused. Why: a pack of GBs copied whole would take its
 * size again in memory. Bytes read from a file, a download or an inflate have a buffer of their own (offset 0), and
 * Node's pooled buffers are 8-byte aligned: a view at another offset is a caller's (an entry stored in a zip, a pack
 * within a larger buffer; the plugin's zips are read as files of their own)
 */
export function readParticlePack(data: Uint8Array, options: { readonly maxCopyBytes?: number } = {}): ParticlePack {
  const header = readParticlePackHeader(data);
  const aligned = data.byteOffset % 4 === 0;
  const maxCopy = options.maxCopyBytes ?? ParticlePackMaxCopyBytes;
  if (!aligned && data.byteLength > maxCopy) {
    throw new Error(
      `particle pack of ${data.byteLength} bytes at offset ${data.byteOffset} (not a multiple of 4): not copied ` +
        `(past ${maxCopy} bytes); give its bytes in a buffer of their own`,
    );
  }
  const bytes = aligned ? data : data.slice();
  const { jsonOffset, jsonLength } = header;
  if (jsonOffset + jsonLength > bytes.byteLength) throw new Error('particle pack truncated (its document)');
  const document = parseParticlePackDocument(bytes.subarray(jsonOffset, jsonOffset + jsonLength), bytes.byteLength);
  return { document, buffer: bytes.buffer as ArrayBuffer, byteOffset: bytes.byteOffset };
}

function getBlock(pack: ParticlePack, name: string) {
  return pack.document.blocks.find((b) => b.name === name);
}

/**
 * A pack's blocks read apart (a pack opened by range: `openParticleListFromPackStream`), by their byte offset, and
 * the reads of those not read (a template's particles: its positions, rotations and attributes, read by range when
 * asked, `Particle.DeferredParticles`)
 */
export interface ParticlePackBlocks {
  readonly arrays: ReadonlyMap<number, Float32Array | Int32Array | Int16Array | Uint32Array>;
  readonly read: (
    block: ParticlePackBlock,
    first: number,
    count: number,
    signal?: AbortSignal,
  ) => Promise<Float32Array | Int32Array | Int16Array | Uint32Array>;
}

/** The view of block `name` of `pack` (its bytes, no copy), none where the pack has no such block */
export function getParticlePackArray(
  pack: ParticlePack,
  name: string,
): Float32Array | Int32Array | Int16Array | Uint32Array | undefined {
  const b = getBlock(pack, name);
  if (!b) return;
  // read apart: the blocks read (a template's particles' none, `Particle.DeferredParticles`)
  if (pack.blocks) return pack.blocks.arrays.get(b.byteOffset);
  const length = b.count * b.components;
  const offset = pack.byteOffset + b.byteOffset;
  if (b.type === 'float32') return new Float32Array(pack.buffer, offset, length);
  if (b.type === 'int16' || b.type === 'snorm16') return new Int16Array(pack.buffer, offset, length);
  if (b.type === 'grid21') return new Uint32Array(pack.buffer, offset, length);
  return new Int32Array(pack.buffer, offset, length);
}

function getFloatArray(pack: ParticlePack, name: string) {
  const a = getParticlePackArray(pack, name);
  if (a && !(a instanceof Float32Array)) throw new Error(`particle pack block '${name}' not float32`);
  return a;
}

function getIntArray(pack: ParticlePack, name: string) {
  const a = getParticlePackArray(pack, name);
  if (a && !(a instanceof Int32Array)) throw new Error(`particle pack block '${name}' not int32`);
  return a;
}

/**
 * A warning where template `index` of a pack an older writer wrote: of more than a tile's particles and not
 * tile-ordered (its copies culled particle by particle: the isg cell's 16x16x16 grid near, frames 39.9 ms, 29.2
 * tile-ordered), or its coarse levels not stored (sampled at every load: the yeast cell's 461 ms). Drawn either way;
 * the pack written again (`writeParticlePack`, `scripts/particle-grid-pack.mjs`) is faster
 */
function warnOlderTemplate(template: ParticlePackDocument, index: number) {
  const slower: string[] = [];
  if (!template.tiles && templateTilesOf(template.count) > 1)
    slower.push('not tile-ordered (its copies culled slower)');
  if (!template.templateLevels) slower.push('its coarse levels not stored (sampled at every load)');
  if (slower.length === 0) return;
  console.warn(
    `particle pack: template ${index} (${template.label}) written by an older writer: ${slower.join(', ')}; ` +
      'write the pack again for faster frames and loads',
  );
}

/**
 * A template pack's coarse levels (`ParticlePackDocument.templateLevels`, the `level.*` blocks) as a list's stored
 * levels (`Particle.StoredLevel`), finest first, views of the pack's bytes; undefined where a block is missing. A
 * template's levels read, not sampled where loaded (`ParticleNesting.getTemplateLevels`)
 */
function getStoredLevels(pack: ParticlePack, grid: GridPositions): Particle.StoredLevel[] | undefined {
  const position = getParticlePackArray(pack, 'level.position');
  const rotation = getParticlePackArray(pack, 'level.rotation');
  const keys = getIntArray(pack, 'level.keys');
  const runs = getIntArray(pack, 'level.runs');
  const scales = getFloatArray(pack, 'level.scales');
  const ranges = getIntArray(pack, 'level.ranges');
  if (!(position instanceof Uint32Array) || !keys || !runs || !scales || !ranges) return undefined;
  if (rotation && !(rotation instanceof Int16Array))
    throw new Error("particle pack block 'level.rotation' not snorm16");
  const runCount = runs.length / 2;
  const levels: Particle.StoredLevel[] = [];
  let r0 = 0;
  for (let k = 0; k < ranges.length; k += 3) {
    const level = ranges[k];
    const lo = ranges[k + 1];
    const hi = lo + ranges[k + 2];
    // its runs: those from its first to its end (the runs ascending, a level's consecutive)
    if (r0 >= runCount || runs[r0 * 2] !== lo)
      throw new Error(`particle pack level ${level}'s runs not from its first`);
    let r1 = r0;
    while (r1 < runCount && runs[r1 * 2] < hi) r1 += 1;
    const offsets = new Int32Array(r1 - r0 + 1);
    const targets = new Int32Array(r1 - r0);
    for (let r = r0; r < r1; ++r) {
      offsets[r - r0] = runs[r * 2] - lo;
      targets[r - r0] = runs[r * 2 + 1];
    }
    offsets[r1 - r0] = hi - lo;
    levels.push({
      level,
      keys: keys.subarray(lo, hi),
      offsets,
      targets,
      scales: scales.subarray(r0, r1),
      coordinates: { kind: 'grid', data: position.subarray(lo * 2, hi * 2), origin: grid.origin, step: grid.step },
      rotations: rotation?.subarray(lo * 4, hi * 4),
    });
    r0 = r1;
  }
  return levels;
}

/**
 * The list of `pack`: its positions, rotations, sites, fibers and attributes views of the pack's bytes; its targets,
 * entities, compartments and radii runs (`ParticleValues`), a value a target, not a particle; no keys (a particle's key
 * its index)
 */
export function createParticleListFromPack(pack: ParticlePack): ParticleList {
  const doc = pack.document;
  const n = doc.count;
  // on a grid (version 4), float32 (versions 2 and 3)
  const positions = getParticlePackArray(pack, 'position');
  const grid = doc.positionGrid;
  let coordinates: GridPositions | Float32Array | undefined;
  if (positions instanceof Uint32Array && grid) {
    const origin = Vec3.create(grid.origin[0], grid.origin[1], grid.origin[2]);
    coordinates = { kind: 'grid', data: positions, origin, step: grid.step };
  } else if (positions instanceof Float32Array) {
    coordinates = positions;
  }
  // a template's particles not read (a pack opened by range): an empty grid of its origin and step, read by range
  const positionBlock = getBlock(pack, 'position');
  const deferred = !!pack.blocks && !coordinates && !!positionBlock && !!grid && !!doc.templateLevels;
  if (deferred) {
    coordinates = {
      kind: 'grid',
      data: new Uint32Array(0),
      origin: Vec3.create(grid!.origin[0], grid!.origin[1], grid!.origin[2]),
      step: grid!.step,
    };
  }
  if (!coordinates || (!deferred && GridPositions.count(coordinates) !== n))
    throw new Error('particle pack without positions');
  // snorm16 (version 3), float32 (version 2); a deferred list's empty, of its type
  const rotationBlock = getBlock(pack, 'rotation');
  const rotations = deferred
    ? rotationBlock && (rotationBlock.type === 'float32' ? new Float32Array(0) : new Int16Array(0))
    : getParticlePackArray(pack, 'rotation');
  if (rotations instanceof Int32Array || rotations instanceof Uint32Array)
    throw new Error("particle pack block 'rotation' not snorm16 or float32");

  // a run a target, in the pack's order (its particles' ranges); tile-ordered a run a target in a tile (its runs block)
  const tiled = doc.tiles ? getIntArray(pack, 'target.runs') : undefined;
  if (doc.tiles && !tiled) throw new Error("particle pack tile-ordered but without a 'target.runs' block");
  const sorted = tiled ? [] : [...doc.targets].sort((a, b) => a.particleOffset - b.particleOffset);
  const runCount = tiled ? tiled.length / 2 : sorted.length;
  const offsets = new Int32Array(runCount + 1);
  const targetIds = new Int32Array(runCount);
  const entityValues = new Int32Array(runCount);
  const compartmentValues = new Int32Array(runCount);
  const radiusValues = new Float32Array(runCount);
  const byId = tiled ? new Map(doc.targets.map((t) => [t.id, t])) : undefined;
  for (let r = 0; r < runCount; ++r) {
    let t: ParticlePackTarget | undefined;
    if (tiled) {
      offsets[r] = tiled[r * 2];
      t = byId!.get(tiled[r * 2 + 1]);
      if (!t) throw new Error(`particle pack run of target ${tiled[r * 2 + 1]}, not a target`);
      if (r > 0 && offsets[r] <= offsets[r - 1]) throw new Error('particle pack runs not ascending');
    } else {
      t = sorted[r];
      if (t.particleOffset !== offsets[r]) throw new Error('particle pack targets not consecutive ranges');
      offsets[r + 1] = t.particleOffset + t.particleCount;
    }
    targetIds[r] = t.id;
    entityValues[r] = t.entity;
    compartmentValues[r] = t.compartment;
    radiusValues[r] = t.radius ?? 0;
  }
  if (tiled) {
    if (runCount > 0 && offsets[0] !== 0) throw new Error('particle pack runs not from its first particle');
    offsets[runCount] = n;
  }
  if (offsets[runCount] !== n) throw new Error('particle pack targets not covering its particles');
  // the targets' bounds where the pack has them all (`Particle.setTargetBounds`), a run's values in its order (a
  // tile-ordered pack's runs' none: computed where asked)
  let bounds: Float64Array | undefined;
  if (!tiled && sorted.every((t) => t.bounds)) {
    bounds = new Float64Array(runCount * ParticleRangeBoundsValues);
    for (let r = 0; r < runCount; ++r) {
      const { min, max, distance } = sorted[r].bounds!;
      bounds.set([min[0], min[1], min[2], max[0], max[1], max[2], distance], r * ParticleRangeBoundsValues);
    }
  }
  const targets = ParticleValues.ofRuns(offsets, targetIds);
  const entities = doc.entities.length > 0 ? ParticleValues.ofRuns(offsets, entityValues) : undefined;
  const compartments = doc.compartments.length > 0 ? ParticleValues.ofRuns(offsets, compartmentValues) : undefined;
  const radii = doc.targets.some((t) => t.radius !== undefined)
    ? ParticleValues.ofRuns(offsets, radiusValues)
    : undefined;

  const siteData = getParticlePackArray(pack, 'site');
  if (siteData && !(siteData instanceof Int16Array)) throw new Error("particle pack block 'site' not int16");
  const { targetInfo, targetName, measures, entityInfo, compartmentInfo } = getParticlePackInfo(doc);
  const siteTypes = getParticlePackSiteTypes(doc);
  const targetMapping = new Map<number, ParticleTarget>();
  for (const t of doc.targets) {
    const s = t.sites;
    if (!s) continue;
    if (!siteData) throw new Error('particle pack target with sites but no site block');
    targetMapping.set(
      t.id,
      createParticlePackSitesTarget(s, siteData.subarray(s.offset * 4, (s.offset + s.count) * 4), siteTypes),
    );
  }
  // the templates: a list each, read once for the targets placing it (views of the same bytes)
  const templates: ParticleList[] = [];
  for (const t of doc.targets) {
    if (!t.template) continue;
    const { index, center } = t.template;
    if (!templates[index]) {
      warnOlderTemplate(doc.templates![index], index);
      templates[index] = createParticleListFromPack({ ...pack, document: doc.templates![index] });
    }
    targetMapping.set(t.id, {
      kind: 'particles',
      particles: templates[index],
      center: Vec3.create(center[0], center[1], center[2]),
    });
  }

  const fiberOffsets = getIntArray(pack, 'fiber.offsets');
  const fiberIndices = getIntArray(pack, 'fiber.indices');
  const fibers =
    doc.fiberCount !== undefined && fiberOffsets && fiberIndices
      ? { count: doc.fiberCount, offsets: fiberOffsets, indices: fiberIndices }
      : undefined;

  let attributes: Map<string, ParticleAttribute> | undefined;
  for (const a of doc.attributes) {
    const values = getFloatArray(pack, `attribute.${a.name}`);
    if (!values) continue;
    attributes ??= new Map();
    attributes.set(a.name, { column: Column.ofFloatArray(values), label: a.label, min: a.min, max: a.max });
  }

  const cell = doc.cell ? Cell.create(Vec3.ofArray(doc.cell.size), Vec3.ofArray(doc.cell.anglesInRadians)) : undefined;

  const list: ParticleList = {
    label: doc.label,
    entryId: doc.entryId,
    count: n,
    cell,
    targets,
    targetInfo,
    targetMapping: targetMapping.size > 0 ? targetMapping : undefined,
    entities,
    entityInfo: entityInfo.size > 0 ? entityInfo : undefined,
    compartments,
    compartmentInfo: compartmentInfo.size > 0 ? compartmentInfo : undefined,
    coordinates,
    rotations,
    radii,
    fibers,
    attributes,
    getParticleLabel: (index: number) => {
      const parts = [`Particle ${index + 1}`];
      const entity = entities ? entityInfo.get(ParticleValues.get(entities, index))?.name : undefined;
      if (entity) parts.push(entity);
      const name = targetName.get(ParticleValues.get(targets, index));
      if (name !== undefined) parts.push(`target ${name}`);
      return parts.join(' · ');
    },
    sourceData: ParticlePackFormat.create(pack),
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
  if (templates.length > 0) ParticleNesting.check(list);
  if (bounds) Particle.setTargetBounds(list, bounds);
  if (doc.sampleOrdered) Particle.setSampleOrdered(list);
  if (doc.boundary) {
    const { min, max, center, radius } = doc.boundary;
    Particle.setBoundary(list, {
      box: Box3D.create(Vec3.create(min[0], min[1], min[2]), Vec3.create(max[0], max[1], max[2])),
      sphere: Sphere3D.create(Vec3.create(center[0], center[1], center[2]), radius),
    });
  }
  if (measures.size > 0) ParticleTileCatalog.setTargetMeasures(list, measures);
  if (doc.tiles) ParticleTiles.setDirectory(list, getTileDirectory(pack, doc.tiles));
  if (deferred) {
    const blocks = pack.blocks!;
    Particle.setDeferred(list, {
      read: async (first, count, signal) => {
        const [positions, rotations] = await Promise.all([
          blocks.read(positionBlock!, first, count, signal),
          rotationBlock ? blocks.read(rotationBlock, first, count, signal) : Promise.resolve(undefined),
        ]);
        if (!(positions instanceof Uint32Array)) throw new Error("particle pack block 'position' not on a grid");
        if (rotations && !(rotations instanceof Int16Array || rotations instanceof Float32Array))
          throw new Error("particle pack block 'rotation' not snorm16 or float32");
        return { positions, rotations };
      },
    });
  }
  if (doc.templateLevels && GridPositions.is(coordinates)) {
    const levels = getStoredLevels(pack, coordinates);
    if (levels) Particle.setStoredLevels(list, levels);
  }
  const density = doc.density ? getFloatArray(pack, 'density') : undefined;
  if (doc.density && density) {
    const [x, y, z] = doc.density.dimensions;
    if (density.length !== x * y * z) throw new Error("particle pack 'density' not of its grid's blocks");
    ParticleDensityGrid.set(list, {
      min: Vec3.ofArray(doc.density.min),
      max: Vec3.ofArray(doc.density.max),
      dimensions: [x, y, z],
      density,
    });
  }
  return list;
}

/**
 * The info of `doc`'s targets (a target's name, entity and compartment, stored once), its entities and compartments,
 * its targets' names and measures (reach, volume where written)
 */
export function getParticlePackInfo(doc: ParticlePackDocument) {
  const targetInfo = new Map<number, ParticleTargetInfo>();
  const targetName = new Map<number, string>();
  const measures = new Map<number, ParticleTargetMeasures>();
  for (const t of doc.targets) {
    targetInfo.set(t.id, {
      name: t.name,
      entity: t.entity >= 0 ? t.entity : undefined,
      // a target's compartment stored once: its particles' (-1 none), where the pack has compartments
      compartment: doc.compartments.length > 0 ? t.compartment : undefined,
    });
    if (t.name !== undefined) targetName.set(t.id, t.name);
    if (t.sites && (t.reach !== undefined || t.volume !== undefined)) {
      measures.set(t.id, { radius: t.reach, volume: t.volume });
    }
  }
  const entityInfo = new Map<number, ParticleEntityInfo>();
  for (const e of doc.entities) entityInfo.set(e.index, { name: e.name, function: e.function });
  const compartmentInfo = new Map<number, ParticleCompartmentInfo>();
  for (const c of doc.compartments) compartmentInfo.set(c.index, { name: c.name });
  return { targetInfo, targetName, measures, entityInfo, compartmentInfo };
}

/** the sites' radii and elements (where every target's sites have them) by type, of `doc` */
export function getParticlePackSiteTypes(doc: ParticlePackDocument) {
  const typeRadii = Float32Array.from(doc.siteRadii ?? []);
  const typeElements =
    doc.siteElements && doc.siteElements.length === typeRadii.length ? Uint8Array.from(doc.siteElements) : undefined;
  return { typeRadii, typeElements };
}

/** A target of sites `s` (a pack target's), `data` its sites' words (four int16 a site: theirs alone) */
export function createParticlePackSitesTarget(
  s: ParticlePackSites,
  data: Int16Array,
  types: ReturnType<typeof getParticlePackSiteTypes>,
): ParticleTarget {
  let chains: ParticleSites['chains'];
  if (s.chains) {
    const offsets = new Int32Array(s.chains.counts.length + 1);
    for (let c = 0; c < s.chains.counts.length; ++c) offsets[c + 1] = offsets[c] + s.chains.counts[c];
    chains = { offsets, labels: s.chains.labels };
  }
  // the quantized sites as they are (views of the bytes): decoded where read
  const sites = ParticleSites.ofQuantized(
    s.count,
    {
      data,
      origin: Vec3.create(s.origin[0], s.origin[1], s.origin[2]),
      scale: s.scale,
      typeRadii: types.typeRadii,
      typeElements: types.typeElements,
      // six a chain (one run without chains)
      extremes:
        s.extremes && s.extremes.length === (s.chains ? s.chains.counts.length : 1) * 6
          ? Int16Array.from(s.extremes)
          : undefined,
    },
    s.label,
    chains,
  );
  return { kind: 'sites', sites, center: Vec3.create(s.center[0], s.center[1], s.center[2]) };
}

/** The tile directory of tile-ordered `pack` (`tiles` its document's), views of its blocks */
function getTileDirectory(
  pack: ParticlePack,
  tiles: NonNullable<ParticlePackDocument['tiles']>,
): ParticleTileDirectory {
  const ranges = getIntArray(pack, 'tile.ranges');
  const boxes = getFloatArray(pack, 'tile.boxes');
  if (!ranges || !boxes) throw new Error("particle pack tile-ordered but without 'tile.ranges' and 'tile.boxes'");
  if (ranges.length / ParticleTileDirectoryRange !== boxes.length / ParticleTileDirectoryBox)
    throw new Error("particle pack 'tile.ranges' and 'tile.boxes' of other tile counts");
  const blocks = tiles.opacityBlocks;
  const density = blocks ? getFloatArray(pack, 'tile.opacity') : undefined;
  if (blocks && density?.length !== (ranges.length / ParticleTileDirectoryRange) * blocks ** 3)
    throw new Error("particle pack 'tile.opacity' not a tile's blocks a tile");
  return {
    min: Vec3.create(tiles.min[0], tiles.min[1], tiles.min[2]),
    max: Vec3.create(tiles.max[0], tiles.max[1], tiles.max[2]),
    dimensions: tiles.dimensions,
    ranges,
    boxes,
    opacity: blocks && density ? { blocks, density } : undefined,
  };
}

/** `coordinates` on a grid in `order` (an old index a new position): a grid's words as they are, else float32's box */
function gatherPositions(coordinates: ParticleList['coordinates'], order: Int32Array): GridPositions {
  const n = order.length;
  if (coordinates instanceof Float32Array) {
    const sorted = new Float32Array(n * 3);
    for (let p = 0; p < n; ++p) sorted.set(coordinates.subarray(order[p] * 3, order[p] * 3 + 3), p * 3);
    return GridPositions.fromFloat(sorted);
  }
  const data = new Uint32Array(n * 2);
  for (let p = 0; p < n; ++p) {
    data[p * 2] = coordinates.data[order[p] * 2];
    data[p * 2 + 1] = coordinates.data[order[p] * 2 + 1];
  }
  return { kind: 'grid', data, origin: coordinates.origin, step: coordinates.step };
}

/** A block to write: its name, type, components and values */
interface BlockData {
  readonly name: string;
  readonly type: ParticlePackBlockType;
  readonly components: number;
  readonly values: Float32Array | Int32Array | Int16Array | Uint32Array;
}

/**
 * `list` as a pack: its particles sorted by target (each target's in their order; the targets in the order of their
 * first particle), its targets' sites (those of kind `sites`; other kinds not packed). A target's entity and
 * compartment must be the same for its particles (they are stored once); its radius is too, by definition (the
 * largest of its particles' kept, where a list's differ in rounding). Keys are not kept (a particle's key is its index in
 * the pack). Tile-ordered with `options.tiles` (`WriteParticlePackOptions`). A list with templates (`ParticleNesting`):
 * its templates once each (a document of `templates`, by target), its targets placing them by index; not tile-ordered
 * (refused: a tile set of copies is not one of their particles), no density grid of its own (a copy's volume its
 * template's particles')
 */
export function writeParticlePack(list: ParticleList, options: WriteParticlePackOptions = {}): Uint8Array<ArrayBuffer> {
  const nested = ParticleNesting.hasTemplates(list);
  if (nested) {
    ParticleNesting.check(list);
    for (const t of ParticleNesting.getIdSpace(list).moves.keys()) Particle.assertRead(t, 'writing a pack');
    if (options.tiles !== undefined) throw new Error('particle pack: a list with templates not tile-ordered');
  }
  const templates: ParticleList[] = [];
  const top = buildPack(list, nested ? { ...options, densityBlocks: 0 } : options, templates);
  const parts = templates.map((t) => {
    const tiles = options.templateTiles ?? templateTilesOf(t.count);
    return buildPack(
      t,
      // its coarse levels written (`level.*`): read, not sampled where loaded (`Particle.StoredLevel`)
      {
        densityBlocks: options.densityBlocks,
        templateLevels: ParticleNesting.MaxTemplateLevels,
        ...(tiles > 1 && { tiles, opacityBlocks: 0 }),
      },
      undefined,
    );
  });

  // the layout: the header, the blocks aligned (the top list's, then each template's), the document after them (their
  // offsets known before it is written)
  let end = HeaderBytes;
  const place = (blocks: readonly BlockData[]) =>
    blocks.map((b): ParticlePackBlock => {
      const byteOffset = align(end);
      end = byteOffset + b.values.byteLength;
      return {
        name: b.name,
        type: b.type,
        components: b.components,
        count: b.values.length / b.components,
        byteOffset,
        byteLength: b.values.byteLength,
      };
    });
  const placed = place(top.blocks);
  const templateDocuments = parts.map((part) => ({ ...part.document, blocks: place(part.blocks) }));
  const document: ParticlePackDocument = {
    ...top.document,
    templates: templateDocuments.length > 0 ? templateDocuments : undefined,
    blocks: placed,
  };
  const json = new TextEncoder().encode(JSON.stringify(document));

  const out = new Uint8Array(end + json.byteLength);
  const view = new DataView(out.buffer);
  for (let k = 0; k < 4; ++k) out[k] = ParticlePackMagic.charCodeAt(k);
  view.setUint32(4, nested ? ParticlePackVersion : ParticlePackVersionWithoutTemplates, true);
  view.setUint32(8, ByteOrderMark, true);
  view.setBigUint64(12, BigInt(end), true);
  view.setBigUint64(20, BigInt(json.byteLength), true);
  const write = (blocks: readonly BlockData[], at: readonly ParticlePackBlock[]) => {
    for (let k = 0; k < blocks.length; ++k) {
      const v = blocks[k].values;
      out.set(new Uint8Array(v.buffer, v.byteOffset, v.byteLength), at[k].byteOffset);
    }
  };
  write(top.blocks, placed);
  for (let k = 0; k < parts.length; ++k) write(parts[k].blocks, templateDocuments[k].blocks);
  out.set(json, end);
  return out;
}

/**
 * The blocks and document (but its blocks' places) of `list` (`writeParticlePack`); its targets' templates gathered
 * into `templates` (an index a list, once each), none where `templates` is undefined (a template's own list)
 */
function buildPack(
  list: ParticleList,
  options: WriteParticlePackOptions,
  templates: ParticleList[] | undefined,
): { blocks: BlockData[]; document: Omit<ParticlePackDocument, 'blocks'> } {
  const n = list.count;
  // runs (`ParticleValues`) expanded: the writer reads a value a particle
  const listTargets = ParticleValues.toArray(list.targets, n, Int32Array);
  const listEntities = list.entities && ParticleValues.toArray(list.entities, n, Int32Array);
  const listCompartments = list.compartments && ParticleValues.toArray(list.compartments, n, Int32Array);
  const listRadii = list.radii && ParticleValues.toArray(list.radii, n, Float32Array);
  // the targets in the order of their first particle, their particles in order
  const targetOrder: number[] = [];
  const targetCounts = new Map<number, number>();
  for (let i = 0; i < n; ++i) {
    const t = listTargets[i];
    const c = targetCounts.get(t);
    if (c === undefined) {
      targetOrder.push(t);
      targetCounts.set(t, 1);
    } else {
      targetCounts.set(t, c + 1);
    }
  }
  for (const id of list.targetInfo.keys()) {
    if (!targetCounts.has(id)) {
      targetOrder.push(id);
      targetCounts.set(id, 0);
    }
  }
  const targetOffset = new Map<number, number>();
  let offset = 0;
  for (const t of targetOrder) {
    targetOffset.set(t, offset);
    offset += targetCounts.get(t)!;
  }
  // old index of each new position, and new position of each old index
  let order = new Int32Array(n);
  const position = new Int32Array(n);
  const next = new Map(targetOffset);
  // a template's (`templates` none) each target's in its levels' sample order (`ParticlePackDocument.sampleOrdered`)
  const sampleOrdered = !templates && n > 0;
  const visit = sampleOrdered ? Particle.getSampleOrder(list) : undefined;
  for (let s = 0; s < n; ++s) {
    const i = visit ? visit[s] : s;
    const p = next.get(listTargets[i])!;
    next.set(listTargets[i], p + 1);
    order[p] = i;
  }

  // the positions on a grid, in the pack's order: the list's grid as it is, else one of its box
  let positions = gatherPositions(list.coordinates, order);
  // tile-ordered: by tile (a stable sort of the order by target), the tiles' directory and the targets' runs
  const tiling = options.tiles !== undefined ? orderByTiles(positions, options.tiles, order, listTargets) : undefined;
  if (tiling) {
    order = tiling.order;
    positions = tiling.positions;
  }
  for (let p = 0; p < n; ++p) position[order[p]] = p;

  const gather = <T extends Float32Array | Int16Array>(source: T, components: number): T => {
    const out = new (source.constructor as new (n: number) => T)(n * components);
    for (let p = 0; p < n; ++p) {
      const o = order[p] * components;
      for (let k = 0; k < components; ++k) out[p * components + k] = source[o + k];
    }
    return out;
  };

  // per target: what is stored once must be the same for its particles (by the list's order: a target's particles
  // runs where tile-ordered)
  const uniformValues = (values: ArrayLike<number> | undefined, what: string) => {
    if (!values) return undefined;
    const of = new Map<number, number>();
    for (let i = 0; i < n; ++i) {
      const t = listTargets[i];
      const v = of.get(t);
      if (v === undefined) of.set(t, values[i]);
      else if (values[i] !== v) throw new Error(`particle pack: target ${t}'s particles differ in ${what}`);
    }
    return of;
  };
  const entityOf = uniformValues(listEntities, 'entity');
  const compartmentOf = uniformValues(listCompartments, 'compartment');
  // a target's particles' radius, by definition one (the largest where they differ in rounding)
  let radiusOf: Map<number, number> | undefined;
  if (listRadii) {
    radiusOf = new Map();
    for (let i = 0; i < n; ++i) radiusOf.set(listTargets[i], Math.max(radiusOf.get(listTargets[i]) ?? 0, listRadii[i]));
  }

  const blocks: BlockData[] = [{ name: 'position', type: 'grid21', components: 2, values: positions.data }];
  if (list.rotations)
    blocks.push({
      name: 'rotation',
      type: 'snorm16',
      components: 4,
      values: gather(ParticleRotations.toSnorm16(list.rotations), 4),
    });

  // the sites of the targets that have them, in the targets' order
  let siteCount = 0;
  for (const t of targetOrder) {
    const target = list.targetMapping?.get(t);
    if (target?.kind === 'sites') siteCount += target.sites.count;
  }
  const siteData = new Int16Array(siteCount * 4);
  // the types: the distinct radii (as float32) and, where every target's sites have them, elements; an index a site
  let withElements = siteCount > 0;
  for (const t of targetOrder) {
    const target = list.targetMapping?.get(t);
    if (target?.kind === 'sites' && !ParticleSites.hasElements(target.sites)) withElements = false;
  }
  const siteRadii: number[] = [];
  const siteElements: number[] = [];
  const typeOf = new Map<number, number>();
  const sitePosition = Vec3();

  // a target's bounds from the positions as written (what a reader decodes), rounded outward
  const rangeBounds = new Float64Array(ParticleRangeBoundsValues);
  const targetBounds = (t: number): ParticlePackTarget['bounds'] => {
    const start = targetOffset.get(t)!;
    Particle.computeRangeBounds(positions, start, start + targetCounts.get(t)!, rangeBounds);
    const [x0, y0, z0, x1, y1, z1, d] = rangeBounds;
    const down = (v: number) => Math.floor(v * BoundsSteps) / BoundsSteps;
    const up = (v: number) => Math.ceil(v * BoundsSteps) / BoundsSteps;
    // the rounded box's center moved by at most half a step an axis: its distance by that more
    return {
      min: [down(x0), down(y0), down(z0)],
      max: [up(x1), up(y1), up(z1)],
      distance: up(d + 1 / BoundsSteps),
    };
  };

  const targets: ParticlePackTarget[] = [];
  const volumeOf = new Map<number, number>();
  let siteOffset = 0;
  for (const t of targetOrder) {
    const info = list.targetInfo.get(t);
    const entity = entityOf?.get(t) ?? info?.entity ?? -1;
    const compartment = compartmentOf?.get(t) ?? -1;
    const target = list.targetMapping?.get(t);
    let sites: ParticlePackSites | undefined;
    let reach: number | undefined;
    let volume: number | undefined;
    let template: ParticlePackTarget['template'];
    if (target?.kind === 'particles') {
      if (!templates) throw new Error(`particle pack: template target ${t} in a template`);
      let index = templates.indexOf(target.particles);
      if (index < 0) index = templates.push(target.particles) - 1;
      const center = ParticleTarget.getCenter(target)!;
      template = { index, center: [center[0], center[1], center[2]] };
    }
    // its volume where its sites give one: the tiles' opacities'
    if (target?.kind === 'sites') {
      const s = target.sites;
      // the box's center the origin, its larger half extent over 32767 the scale (0.01 Å for ±327 Å)
      const min = Vec3.create(Infinity, Infinity, Infinity);
      const max = Vec3.create(-Infinity, -Infinity, -Infinity);
      for (let i = 0; i < s.count; ++i) {
        ParticleSites.getPosition(s, i, sitePosition);
        Vec3.min(min, min, sitePosition);
        Vec3.max(max, max, sitePosition);
      }
      const origin = s.count > 0 ? Vec3.scale(Vec3(), Vec3.add(Vec3(), min, max), 0.5) : Vec3();
      const half = s.count > 0 ? Math.max(max[0] - min[0], max[1] - min[1], max[2] - min[2]) / 2 : 0;
      const scale = half > 0 ? half / 32767 : 1;
      for (let i = 0; i < s.count; ++i) {
        ParticleSites.getPosition(s, i, sitePosition);
        const o = (siteOffset + i) * 4;
        for (let k = 0; k < 3; ++k) {
          const q = Math.round((sitePosition[k] - origin[k]) / scale);
          siteData[o + k] = Math.max(-32767, Math.min(32767, q));
        }
        const r = Math.fround(ParticleSites.getRadius(s, i));
        const element = withElements ? ParticleSites.getElement(s, i) : 0;
        // a float32's bits and the element (a byte) as one key
        radiusBits[0] = r;
        const key = radiusBitsU32[0] * 256 + element;
        let type = typeOf.get(key);
        if (type === undefined) {
          type = siteRadii.length;
          if (type > 0xffff) throw new Error('particle pack: more than 65536 site types');
          typeOf.set(key, type);
          siteRadii.push(r);
          siteElements.push(element);
        }
        siteData[o + 3] = type;
      }
      const center = ParticleTarget.getCenter(target)!;
      let chains: ParticlePackSites['chains'];
      if (s.chains) {
        const counts: number[] = [];
        for (let c = 0; c + 1 < s.chains.offsets.length; ++c)
          counts.push(s.chains.offsets[c + 1] - s.chains.offsets[c]);
        chains = { labels: [...s.chains.labels], counts };
      }
      // each chain's extremes of the positions as written (a box decoded as a reader decodes the sites)
      const written = siteData.subarray(siteOffset * 4, (siteOffset + s.count) * 4);
      const runs = s.chains ? s.chains.offsets : Int32Array.of(0, s.count);
      const extremes = new Int16Array((runs.length - 1) * 6);
      for (let c = 0; c + 1 < runs.length; ++c) {
        ParticleSites.getQuantizedExtremes(written, runs[c], runs[c + 1], extremes, c * 6);
      }
      sites = {
        offset: siteOffset,
        count: s.count,
        origin: [origin[0], origin[1], origin[2]],
        scale,
        label: s.label,
        center: [center[0], center[1], center[2]],
        chains,
        extremes: Array.from(extremes),
      };
      // its reach and volume as a reader computes them: of its sites as written, its types' radii and elements
      const read = ParticleSites.ofQuantized(s.count, {
        data: written,
        origin,
        scale,
        typeRadii: Float32Array.from(siteRadii),
        typeElements: withElements ? Uint8Array.from(siteElements) : undefined,
      });
      const readTarget: ParticleTarget = { kind: 'sites', sites: read, center: Vec3.clone(center) };
      reach = ParticleTileCatalog.getTargetRadius(readTarget);
      volume = ParticleTileCatalog.getTargetVolume(readTarget);
      if (volume !== undefined) volumeOf.set(t, volume);
      siteOffset += s.count;
    }
    targets.push({
      id: t,
      name: info?.name,
      entity,
      compartment,
      particleOffset: tiling ? -1 : targetOffset.get(t)!,
      particleCount: targetCounts.get(t)!,
      radius: radiusOf?.get(t),
      sites,
      reach,
      volume,
      // a target a range: its bounds (tile-ordered its runs', computed where asked)
      bounds: tiling ? undefined : targetBounds(t),
      template,
    });
  }
  if (siteCount > 0) blocks.push({ name: 'site', type: 'int16', components: 4, values: siteData });
  if (tiling) {
    blocks.push({ name: 'target.runs', type: 'int32', components: 2, values: tiling.runs });
    blocks.push({ name: 'tile.ranges', type: 'int32', components: 3, values: tiling.ranges });
    blocks.push({ name: 'tile.boxes', type: 'float32', components: 6, values: tiling.boxes });
  }
  // the tiles' opacities: of their particles as written, their targets' volumes (as a tile set computes them)
  const opacityBlocks = tiling ? Math.max(0, Math.floor(options.opacityBlocks ?? 4)) : 0;
  if (tiling && opacityBlocks > 0) {
    const b3 = opacityBlocks ** 3;
    const count = tiling.ranges.length / ParticleTileDirectoryRange;
    const density = new Float32Array(count * b3);
    const grid = ParticleTileGrid.create(tiling.min, tiling.max, [
      tiling.dimensions,
      tiling.dimensions,
      tiling.dimensions,
    ]);
    const cell = Box3D();
    const targetOf = (q: number) => listTargets[order[q]];
    const volume = (id: number) => volumeOf.get(id);
    for (let k = 0; k < count; ++k) {
      const first = tiling.ranges[k * 3 + 1];
      const end = first + tiling.ranges[k * 3 + 2];
      ParticleTileGrid.getBox(grid, tiling.ranges[k * 3], cell);
      const opacity = ParticleTileOpaque.computeOpacity(
        cell,
        (sink) => ParticleTileOpaque.forEachMoleculeOf(first, end, targetOf, positions, volume, sink),
        { blocks: opacityBlocks },
      );
      if (opacity) density.set(opacity.density, k * b3);
    }
    blocks.push({ name: 'tile.opacity', type: 'float32', components: b3, values: density });
  }
  // the particles' density grid: of their positions as written, their targets' volumes; tile-ordered with the tiles'
  // opacities none (they are its), a template tile-ordered (no opacities) its own (what places its copies' occlusion,
  // `ParticleNesting.getDensityGrid`: computed where read else, the yeast cell's 5.3 s)
  const densityBlocks = tiling && opacityBlocks > 0 ? 0 : Math.max(0, Math.floor(options.densityBlocks ?? 64));
  let densityGrid: ParticlePackDocument['density'];
  if (densityBlocks > 0 && n > 0) {
    const box = ParticleTileOpaque.getCentersBox(positions, n);
    const dimensions = ParticleTileOpaque.getGridDimensions(box, densityBlocks);
    const targetOf = (q: number) => listTargets[order[q]];
    const volume = (id: number) => volumeOf.get(id);
    const density = ParticleTileOpaque.computeDensity(box, dimensions, (sink) =>
      ParticleTileOpaque.forEachMoleculeOf(0, n, targetOf, positions, volume, sink),
    );
    if (density) {
      blocks.push({ name: 'density', type: 'float32', components: 1, values: density });
      densityGrid = { dimensions, min: Array.from(box.min), max: Array.from(box.max) };
    }
  }

  // fibers: their particles at their new positions, in their order
  if (list.fibers) {
    const indices = new Int32Array(list.fibers.indices.length);
    for (let k = 0; k < indices.length; ++k) indices[k] = position[list.fibers.indices[k]];
    blocks.push({ name: 'fiber.offsets', type: 'int32', components: 1, values: Int32Array.from(list.fibers.offsets) });
    blocks.push({ name: 'fiber.indices', type: 'int32', components: 1, values: indices });
  }

  // numeric attributes
  const attributes: { name: string; label: string; min: number; max: number }[] = [];
  const attributeValues: { name: string; values: Float32Array }[] = [];
  if (list.attributes) {
    for (const [name, a] of list.attributes) {
      const type = a.column.schema.valueType;
      if (type !== 'float' && type !== 'int') continue;
      const values = new Float32Array(n);
      for (let p = 0; p < n; ++p) values[p] = a.column.value(order[p]);
      blocks.push({ name: `attribute.${name}`, type: 'float32', components: 1, values });
      attributes.push({ name, label: a.label, min: a.min, max: a.max });
      attributeValues.push({ name, values });
    }
  }

  // a template's coarse levels: copies of samples of its particles as written
  const templateLevels = options.templateLevels ?? 0;
  if (templateLevels > 0 && n > 0) {
    const rotationBlock = blocks.find((b) => b.name === 'rotation')?.values as Int16Array | undefined;
    blocks.push(...writeTemplateLevels(positions, n, (q) => listTargets[order[q]], rotationBlock, templateLevels));
  }

  // the coarse levels: their nodes' particles copies of a sample of their tiles'
  const levelCount = tiling && options.levels ? Math.log2(tiling.dimensions) : 0;
  if (tiling && levelCount > 0) {
    const rotationBlock = blocks.find((b) => b.name === 'rotation')?.values as Int16Array | undefined;
    const levels = writeLevels(tiling, levelCount, (q) => listTargets[order[q]], rotationBlock, attributeValues);
    blocks.push(...levels);
  }

  const entities: { index: number; name: string; function?: string }[] = [];
  for (const [index, e] of list.entityInfo ?? []) entities.push({ index, name: e.name, function: e.function });
  const compartments: { index: number; name: string }[] = [];
  for (const [index, c] of list.compartmentInfo ?? []) compartments.push({ index, name: c.name });

  // a template's boundary of its positions and radii as written (what a reader decodes: a radius a target), rounded
  // outward (`ParticlePackDocument.boundary`)
  let boundary: ParticlePackDocument['boundary'];
  if (sampleOrdered) {
    let radii: ParticleValues<Float32Array> | undefined;
    if (radiusOf && tiling) {
      // tile-ordered: a radius a particle (a target's runs apart)
      const values = new Float32Array(n);
      for (let q = 0; q < n; ++q) values[q] = radiusOf.get(listTargets[order[q]]) ?? 0;
      radii = values;
    } else if (radiusOf) {
      const offsets = new Int32Array(targetOrder.length + 1);
      const values = new Float32Array(targetOrder.length);
      targetOrder.forEach((t, k) => {
        offsets[k + 1] = offsets[k] + targetCounts.get(t)!;
        values[k] = radiusOf.get(t) ?? 0;
      });
      radii = ParticleValues.ofRuns(offsets, values);
    }
    const { box, sphere } = Particle.getBoundary({
      count: n,
      coordinates: positions,
      radii,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    } as unknown as ParticleList);
    const down = (v: number) => Math.floor(v * BoundsSteps) / BoundsSteps;
    const up = (v: number) => Math.ceil(v * BoundsSteps) / BoundsSteps;
    const c = sphere.center;
    boundary = {
      min: [down(box.min[0]), down(box.min[1]), down(box.min[2])],
      max: [up(box.max[0]), up(box.max[1]), up(box.max[2])],
      center: [c[0], c[1], c[2]],
      radius: up(sphere.radius + 1 / BoundsSteps),
    };
  }

  const document: Omit<ParticlePackDocument, 'blocks'> = {
    label: list.label,
    entryId: list.entryId,
    count: n,
    cell: list.cell
      ? { size: Array.from(list.cell.size), anglesInRadians: Array.from(list.cell.anglesInRadians) }
      : undefined,
    entities,
    compartments,
    targets,
    attributes,
    fiberCount: list.fibers?.count,
    siteRadii,
    siteElements: withElements ? siteElements : undefined,
    positionGrid: { origin: Array.from(positions.origin), step: positions.step },
    tiles: tiling
      ? {
          dimensions: tiling.dimensions,
          min: Array.from(tiling.min),
          max: Array.from(tiling.max),
          opacityBlocks: opacityBlocks > 0 ? opacityBlocks : undefined,
          levels: levelCount > 0 ? levelCount : undefined,
        }
      : undefined,
    density: densityGrid,
    sampleOrdered: sampleOrdered || undefined,
    templateLevels: templateLevels > 0 && n > 0 ? templateLevels : undefined,
    boundary,
  };
  return { blocks, document };
}

/**
 * A template's coarse levels' blocks (`ParticlePackDocument.templateLevels`), of its particles as written (in their
 * stored order: tiles, a target's in its sample order): level `L` of each target every `8^L`-th of its particles in
 * that order (the first kept), a run a target as first seen, its scale the square root of its particles over those
 * kept; levels 1 to `levelCount`, finest first. What `ParticleNesting.getTemplateLevels` samples where loaded
 * (`Particle.sampleByTarget` in a stored order: its picks by runs), written once
 */
function writeTemplateLevels(
  positions: GridPositions,
  n: number,
  targetOf: (q: number) => number,
  rotations: Int16Array | undefined,
  levelCount: number,
): BlockData[] {
  // the targets by id as first seen, a particle's target's slot
  const slotOf = new Map<number, number>();
  const ids: number[] = [];
  const slots = new Int32Array(n);
  for (let q = 0; q < n; ++q) {
    const t = targetOf(q);
    let s = slotOf.get(t);
    if (s === undefined) {
      s = ids.length;
      slotOf.set(t, s);
      ids.push(t);
    }
    slots[q] = s;
  }
  const targets = ids.length;
  const counts = new Int32Array(targets);
  for (let q = 0; q < n; ++q) counts[slots[q]] += 1;
  const levelSizes: number[] = [];
  let total = 0;
  for (let level = 1; level <= levelCount; ++level) {
    const stride = 8 ** level;
    let size = 0;
    for (let s = 0; s < targets; ++s) size += Math.ceil(counts[s] / stride);
    levelSizes.push(size);
    total += size;
  }
  const data = new Uint32Array(total * 2);
  const keys = new Int32Array(total);
  const rotation = rotations ? new Int16Array(total * 4) : undefined;
  const runs = new Int32Array(levelCount * targets * 2);
  const scales = new Float32Array(levelCount * targets);
  const ranges = new Int32Array(levelCount * 3);
  const seen = new Int32Array(targets);
  const cursor = new Int32Array(targets);
  let first = 0;
  for (let level = 1; level <= levelCount; ++level) {
    const stride = 8 ** level;
    // a run a target, its particles every stride-th of its own in order
    let at = first;
    for (let s = 0; s < targets; ++s) {
      cursor[s] = at;
      const kept = Math.ceil(counts[s] / stride);
      const r = (level - 1) * targets + s;
      runs[r * 2] = at;
      runs[r * 2 + 1] = ids[s];
      scales[r] = Math.sqrt(counts[s] / kept);
      at += kept;
    }
    seen.fill(0);
    for (let q = 0; q < n; ++q) {
      const s = slots[q];
      if (seen[s]++ % stride !== 0) continue;
      const k = cursor[s]++;
      data[k * 2] = positions.data[q * 2];
      data[k * 2 + 1] = positions.data[q * 2 + 1];
      keys[k] = q;
      if (rotation) rotation.set(rotations!.subarray(q * 4, q * 4 + 4), k * 4);
    }
    ranges.set([level, first, levelSizes[level - 1]], (level - 1) * 3);
    first = at;
  }
  const out: BlockData[] = [
    { name: 'level.position', type: 'grid21', components: 2, values: data },
    { name: 'level.keys', type: 'int32', components: 1, values: keys },
    { name: 'level.runs', type: 'int32', components: 2, values: runs },
    { name: 'level.scales', type: 'float32', components: 1, values: scales },
    { name: 'level.ranges', type: 'int32', components: 3, values: ranges },
  ];
  if (rotation) out.push({ name: 'level.rotation', type: 'snorm16', components: 4, values: rotation });
  return out;
}

/** a fine grid's cells an axis over the tiles' box (`writeLevels`): its Morton codes 27 bits */
const LevelSampleCells = 512;

/**
 * The coarse levels' blocks (`ParticlePackDocument.tiles.levels`) of tile-ordered particles `tiling` (`orderByTiles`;
 * a particle's index its place there, `targetOf` its target): level `L` a node a cell of the grid of
 * `dimensions >> L` cells an axis, its tiles consecutive (the tiles in Morton order), of each target every `8^L`-th
 * of its particles there (at least one) in Morton order of a fine grid over the box (`LevelSampleCells` an axis: a
 * node's sample spread over it), a run a target, its scale the cube root of its particles over those kept; the nodes
 * by level, coarsest first, a level's in Morton order. Why by target: a sample by count alone keeps a node's volume
 * only on average (its large molecules sampled more or less than its small ones), and loses a target rarer than
 * `8^L` there
 */
function writeLevels(
  tiling: ReturnType<typeof orderByTiles>,
  levelCount: number,
  targetOf: (q: number) => number,
  rotations: Int16Array | undefined,
  attributes: readonly { readonly name: string; readonly values: Float32Array }[],
): BlockData[] {
  const { positions, ranges, dimensions: d, min, max } = tiling;
  const n = positions.data.length / 2;
  // per particle its tile's Morton code and cell, by the directory's ranges
  const tileCode = new Int32Array(n);
  const tileCell = new Int32Array(n);
  for (let k = 0; k < ranges.length; k += 3) {
    const id = ranges[k];
    const code = mortonOrder3d(id % d, Math.floor(id / d) % d, Math.floor(id / (d * d)));
    tileCode.fill(code, ranges[k + 1], ranges[k + 1] + ranges[k + 2]);
    tileCell.fill(id, ranges[k + 1], ranges[k + 1] + ranges[k + 2]);
  }
  // the particles in Morton order of the fine grid, a code's by index (a stable sort: any particle count)
  const f = LevelSampleCells;
  const p = Vec3();
  const codes = new Uint32Array(n);
  const cellOf = (v: number, a: number) =>
    Math.max(0, Math.min(f - 1, Math.floor(((v - min[a]) / Math.max(1e-6, max[a] - min[a])) * f)));
  for (let q = 0; q < n; ++q) {
    GridPositions.get(positions, q, p);
    codes[q] = mortonOrder3d(cellOf(p[0], 0), cellOf(p[1], 1), cellOf(p[2], 2));
  }
  const fine = radixSortIndices(codes, f * f * f - 1);

  // the samples' values chunked, not in arrays of numbers: a few hundred million particles' samples pass an array's
  // length limit
  const ChunkSize = 1 << 20;
  const position = ChunkedArray.create(Uint32Array, 2, ChunkSize);
  const keyChunks = ChunkedArray.create(Int32Array, 1, ChunkSize);
  const nodeRanges: number[] = [];
  const nodeBoxes: number[] = [];
  const runs = ChunkedArray.create(Int32Array, 2, ChunkSize);
  const runScales = ChunkedArray.create(Float32Array, 1, ChunkSize);
  let count = 0;
  const bounds = new Float64Array(ParticleRangeBoundsValues);
  for (let level = levelCount; level >= 1; --level) {
    const stride = 8 ** level;
    const shift = 3 * level;
    // per node, per target its particles seen and those kept (every stride-th, the first kept)
    const byNode = new Map<
      number,
      { node: number; first: number; byTarget: Map<number, { seen: number; q: number[] }> }
    >();
    for (let s = 0; s < n; ++s) {
      const q = fine[s];
      const node = tileCode[q] >> shift;
      let at = byNode.get(node);
      if (!at) byNode.set(node, (at = { node, first: q, byTarget: new Map() }));
      const t = targetOf(q);
      let of = at.byTarget.get(t);
      if (!of) at.byTarget.set(t, (of = { seen: 0, q: [] }));
      if (of.seen++ % stride === 0) of.q.push(q);
    }
    const nodes = [...byNode.values()].sort((a, b) => a.node - b.node);
    const dl = d >> level;
    for (const { first: q0, byTarget } of nodes) {
      const first = count;
      for (const [t, { seen, q: list }] of byTarget) {
        ChunkedArray.add2(runs, count, t);
        ChunkedArray.add(runScales, Math.cbrt(seen / list.length));
        for (const j of list) {
          ChunkedArray.add2(position, positions.data[j * 2], positions.data[j * 2 + 1]);
          ChunkedArray.add(keyChunks, j);
          count += 1;
        }
      }
      const id = tileCell[q0];
      const x = (id % d) >> level;
      const y = (Math.floor(id / d) % d) >> level;
      const z = Math.floor(id / (d * d)) >> level;
      nodeRanges.push(level, x + y * dl + z * dl * dl, first, count - first);
    }
  }
  const data = ChunkedArray.compact(position) as Uint32Array;
  const nodeKeys = ChunkedArray.compact(keyChunks) as Int32Array;
  const written: GridPositions = { kind: 'grid', data, origin: positions.origin, step: positions.step };
  for (let k = 0; k < nodeRanges.length; k += 4) {
    Particle.computeRangeBounds(written, nodeRanges[k + 2], nodeRanges[k + 2] + nodeRanges[k + 3], bounds);
    for (let a = 0; a < 6; ++a) {
      const v = bounds[a];
      const r = Math.fround(v);
      nodeBoxes.push(a < 3 ? (r > v ? nextDown(r) : r) : r < v ? nextUp(r) : r);
    }
  }
  const out: BlockData[] = [
    { name: 'node.position', type: 'grid21', components: 2, values: data },
    { name: 'node.keys', type: 'int32', components: 1, values: nodeKeys },
    { name: 'node.runs', type: 'int32', components: 2, values: ChunkedArray.compact(runs) as Int32Array },
    { name: 'node.scales', type: 'float32', components: 1, values: ChunkedArray.compact(runScales) as Float32Array },
    { name: 'node.ranges', type: 'int32', components: 4, values: Int32Array.from(nodeRanges) },
    { name: 'node.boxes', type: 'float32', components: 6, values: Float32Array.from(nodeBoxes) },
  ];
  if (rotations) {
    const r = new Int16Array(count * 4);
    for (let i = 0; i < count; ++i) r.set(rotations.subarray(nodeKeys[i] * 4, nodeKeys[i] * 4 + 4), i * 4);
    out.push({ name: 'node.rotation', type: 'snorm16', components: 4, values: r });
  }
  for (const a of attributes) {
    out.push({
      name: `node.attribute.${a.name}`,
      type: 'float32',
      components: 1,
      values: Float32Array.from(nodeKeys, (j) => a.values[j]),
    });
  }
  return out;
}

/**
 * The order by tile of `positions` (in `order`: an old index a position, by target): a grid of `dimensions` cells an
 * axis over their box (`ParticleTileGrid`, as a tile set makes it), the tiles in Morton order of their cells, a tile's
 * particles in their order (a counting sort: by target where `order` is); the positions in that order, per run of a
 * target in a tile its first particle and target (`listTargets` by old index), per tile its cell's id, first particle,
 * count and its centers' box
 */
function orderByTiles(positions: GridPositions, dimensions: number, order: Int32Array, listTargets: Int32Array) {
  if (!(dimensions >= 1 && dimensions <= 1024) || (dimensions & (dimensions - 1)) !== 0)
    throw new Error(`particle pack: tiles ${dimensions}, a power of two up to 1024`);
  const n = order.length;
  const bounds = new Float64Array(ParticleRangeBoundsValues);
  Particle.computeRangeBounds(positions, 0, n, bounds);
  const min = Vec3.create(bounds[0], bounds[1], bounds[2]);
  const max = Vec3.create(bounds[3], bounds[4], bounds[5]);
  const grid = ParticleTileGrid.create(min, max, [dimensions, dimensions, dimensions]);
  // per particle its tile's Morton code (the cell's id once the tile is known), per code its count
  const cells = dimensions ** 3;
  const codeOf = new Int32Array(n);
  const counts = new Int32Array(cells);
  const idOf = new Int32Array(cells);
  const p = Vec3();
  for (let j = 0; j < n; ++j) {
    GridPositions.get(positions, j, p);
    const id = ParticleTileGrid.getIdOf(grid, p[0], p[1], p[2]);
    const x = id % dimensions;
    const y = Math.floor(id / dimensions) % dimensions;
    const z = Math.floor(id / (dimensions * dimensions));
    const code = mortonOrder3d(x, y, z);
    codeOf[j] = code;
    counts[code] += 1;
    idOf[code] = id;
  }
  const firsts = new Int32Array(cells);
  for (let c = 0, first = 0; c < cells; ++c) {
    firsts[c] = first;
    first += counts[c];
  }
  const next = firsts.slice();
  const tiled = new Int32Array(n);
  const data = new Uint32Array(n * 2);
  for (let j = 0; j < n; ++j) {
    const q = next[codeOf[j]]++;
    tiled[q] = order[j];
    data[q * 2] = positions.data[j * 2];
    data[q * 2 + 1] = positions.data[j * 2 + 1];
  }
  const sorted: GridPositions = { kind: 'grid', data, origin: positions.origin, step: positions.step };
  // the runs: a target in a tile
  const runs: number[] = [];
  for (let q = 0, c = 0; q < n; ++q) {
    while (q >= firsts[c] + counts[c]) ++c;
    const t = listTargets[tiled[q]];
    if (q === firsts[c] || t !== listTargets[tiled[q - 1]]) runs.push(q, t);
  }
  // the tiles with particles in Morton order, their centers' boxes
  const ranges: number[] = [];
  const boxes: number[] = [];
  for (let c = 0; c < cells; ++c) {
    if (counts[c] === 0) continue;
    ranges.push(idOf[c], firsts[c], counts[c]);
    Particle.computeRangeBounds(sorted, firsts[c], firsts[c] + counts[c], bounds);
    boxes.push(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
  }
  return {
    order: tiled,
    positions: sorted,
    runs: Int32Array.from(runs),
    ranges: Int32Array.from(ranges),
    // outward to float32: a box as read holds its particles
    boxes: Float32Array.from(boxes, (v, k) => {
      const f = Math.fround(v);
      if (k % 6 < 3) return f > v ? nextDown(f) : f;
      return f < v ? nextUp(f) : f;
    }),
    dimensions,
    min,
    max,
  };
}

const f32 = new Float32Array(1);
const u32 = new Uint32Array(f32.buffer);
/** the float32 after `f` toward +∞ (a finite one) */
function nextUp(f: number) {
  if (f === 0) return 1.401298464324817e-45;
  f32[0] = f;
  u32[0] += f > 0 ? 1 : -1;
  return f32[0];
}
/** the float32 after `f` toward -∞ */
function nextDown(f: number) {
  return -nextUp(-f);
}

function align(offset: number) {
  return Math.ceil(offset / BlockAlignment) * BlockAlignment;
}
