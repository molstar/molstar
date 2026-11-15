/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '@molstar/core/data/db';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Task } from '@molstar/core/task';
import { CustomProperties } from '../../model/custom-property.js';
import {
  Particle,
  type ParticleAttribute,
  type ParticleCompartmentInfo,
  type ParticleList,
  ParticleSites,
  type ParticleTarget,
  ParticleValues,
} from '../../model/particles/particle-list.js';
import {
  type ParticleTile,
  type ParticleTileCatalog,
  ParticleTileDirectoryBox,
  ParticleTileDirectoryRange,
  ParticleTileGrid,
  type ParticleTileSet,
  type ParticleTileTargetInfo,
  ParticleTiles,
} from '../../model/particles/particle-tiles.js';
import type { ModelFormat } from '../format.js';
import {
  createParticlePackSitesTarget,
  getParticlePackInfo,
  getParticlePackSiteTypes,
  type ParticlePackBlock,
  type ParticlePackDocument,
  type ParticlePackHeader,
  ParticlePackHeaderBytes,
  type ParticlePackTarget,
  parseParticlePackDocument,
  readParticlePackHeader,
} from './pack.js';

// A tile-ordered particle pack read by ranges: its header, its document and the blocks every
// tile needs (the directory, its opacities, the targets' runs) at once; a tile's particles (its range of the per-particle
// blocks) and a target's sites when asked (`ParticleTileProvider`). The pack itself never whole: a file larger than an
// `ArrayBuffer`, the first tiles drawn before the rest is read.

/** The bytes of a pack, read by range: a file's (`Blob.slice`), a URL's (HTTP range requests), bytes in memory */
export interface ParticlePackSource {
  /** the pack's bytes */
  readonly byteLength: number;
  /** its bytes `[offset, offset + length)`, in an array of their own (4-byte aligned) */
  read(offset: number, length: number): Promise<Uint8Array>;
}

export namespace ParticlePackSource {
  /** bytes in memory (tests, a pack downloaded whole) */
  export function fromBytes(data: Uint8Array): ParticlePackSource {
    return {
      byteLength: data.byteLength,
      read: async (offset, length) => {
        if (offset < 0 || offset + length > data.byteLength) throw new Error('particle pack read past its end');
        const bytes = data.subarray(offset, offset + length);
        return bytes.byteOffset % 4 === 0 ? bytes : bytes.slice();
      },
    };
  }

  export interface BatchedOptions {
    /** reads of `source` in flight at most */
    readonly concurrency?: number;
    /** pending ranges at most this far apart (bytes) read as one */
    readonly gap?: number;
    /** a merged read's bytes at most */
    readonly maxBytes?: number;
  }

  /**
   * `source` read at most `concurrency` ranges at a time, pending ranges near each other (`gap`) read as one (up to
   * `maxBytes`): a view's targets' sites, thousands of small ranges side by side in the sites' block, a few requests,
   * not thousands at once (a browser refuses those: `ERR_INSUFFICIENT_RESOURCES`)
   */
  export function batched(source: ParticlePackSource, options: BatchedOptions = {}): ParticlePackSource {
    const concurrency = options.concurrency ?? 6;
    const gap = options.gap ?? 256 * 1024;
    const maxBytes = options.maxBytes ?? 16 * 1024 * 1024;
    type Pending = {
      offset: number;
      length: number;
      resolve: (b: Uint8Array) => void;
      reject: (e: unknown) => void;
    };
    const pending: Pending[] = [];
    let inFlight = 0;
    const next = () => {
      while (inFlight < concurrency && pending.length > 0) {
        // the first pending by offset and those after it within the gap, as one read
        pending.sort((a, b) => a.offset - b.offset);
        const group = [pending.shift()!];
        let end = group[0].offset + group[0].length;
        while (pending.length > 0) {
          const p = pending[0];
          if (p.offset - end > gap || Math.max(end, p.offset + p.length) - group[0].offset > maxBytes) break;
          group.push(pending.shift()!);
          end = Math.max(end, p.offset + p.length);
        }
        const start = group[0].offset;
        inFlight += 1;
        source
          .read(start, end - start)
          .then(
            (bytes) => {
              for (const p of group) {
                const part = bytes.subarray(p.offset - start, p.offset - start + p.length);
                p.resolve(part.byteOffset % 4 === 0 && group.length === 1 ? part : part.slice());
              }
            },
            (e) => {
              for (const p of group) p.reject(e);
            },
          )
          .finally(() => {
            inFlight -= 1;
            next();
          });
      }
    };
    return {
      byteLength: source.byteLength,
      read: (offset, length) =>
        new Promise<Uint8Array>((resolve, reject) => {
          pending.push({ offset, length, resolve, reject });
          // the ranges asked together (a microtask's) merged before the first goes
          queueMicrotask(next);
        }),
    };
  }

  /** a file (`File`) or a blob: a range a `slice` */
  export function fromBlob(blob: Blob): ParticlePackSource {
    return {
      byteLength: blob.size,
      read: async (offset, length) => {
        if (offset < 0 || offset + length > blob.size) throw new Error('particle pack read past its end');
        return new Uint8Array(await blob.slice(offset, offset + length).arrayBuffer());
      },
    };
  }
}

/** A pack opened by range: its header and document, the blocks read when asked */
export interface ParticlePackStream {
  readonly source: ParticlePackSource;
  readonly header: ParticlePackHeader;
  readonly document: ParticlePackDocument;
}

/** The pack of `source` opened: its header and document read (two ranges) */
export async function openParticlePackStream(source: ParticlePackSource): Promise<ParticlePackStream> {
  if (source.byteLength < ParticlePackHeaderBytes) throw new Error('not a particle pack (shorter than its header)');
  const header = readParticlePackHeader(await source.read(0, ParticlePackHeaderBytes));
  const { jsonOffset, jsonLength } = header;
  if (jsonOffset + jsonLength > source.byteLength) throw new Error('particle pack truncated (its document)');
  const document = parseParticlePackDocument(await source.read(jsonOffset, jsonLength), source.byteLength);
  return { source, header, document };
}

/**
 * Whether a pack of `doc` is read by range as a tile set: tile-ordered (its directory, its targets' runs), its
 * positions on a grid, its targets' reach known (the catalog's without their sites read). Its fibers not read: their
 * particles the whole's, a tile's list none (as a whole list's tiles have none)
 */
export function isParticlePackStreamable(doc: ParticlePackDocument) {
  if (!doc.tiles || !doc.positionGrid) return false;
  const has = (name: string) => doc.blocks.some((b) => b.name === name);
  if (!has('tile.ranges') || !has('tile.boxes') || !has('target.runs') || !has('position')) return false;
  return doc.targets.every((t) => !t.sites || t.reach !== undefined);
}

export type ParticlePackStreamFormat = ModelFormat<ParticlePackStream>;

export namespace ParticlePackStreamFormat {
  export function is(x?: ModelFormat): x is ParticlePackStreamFormat {
    return x?.kind === 'particle-pack-stream';
  }

  export function create(stream: ParticlePackStream): ParticlePackStreamFormat {
    return { kind: 'particle-pack-stream', name: 'particle-pack', data: stream };
  }
}

function getBlock(doc: ParticlePackDocument, name: string) {
  return doc.blocks.find((b) => b.name === name);
}

/** items `[first, first + count)` of `block` (`components` values an item), as its type's array */
async function readItems(source: ParticlePackSource, block: ParticlePackBlock, first: number, count: number) {
  const itemBytes = block.byteLength / Math.max(1, block.count);
  const bytes = await source.read(block.byteOffset + first * itemBytes, count * itemBytes);
  const length = count * block.components;
  if (block.type === 'float32') return new Float32Array(bytes.buffer, bytes.byteOffset, length);
  if (block.type === 'int16' || block.type === 'snorm16') return new Int16Array(bytes.buffer, bytes.byteOffset, length);
  if (block.type === 'grid21') return new Uint32Array(bytes.buffer, bytes.byteOffset, length);
  return new Int32Array(bytes.buffer, bytes.byteOffset, length);
}

/** block `name` of `stream` whole, of `type`'s array (none where the pack has none) */
async function readWhole<A extends Int32Array | Float32Array>(
  stream: ParticlePackStream,
  name: string,
  type: new (...args: any[]) => A,
): Promise<A | undefined> {
  const block = getBlock(stream.document, name);
  if (!block) return undefined;
  const items = await readItems(stream.source, block, 0, block.count);
  if (!(items instanceof type)) throw new Error(`particle pack block '${name}' not of its type`);
  return items as A;
}

/**
 * Per compartment its entities (in `entityInfo`'s order) and whether some of its particles have none, by the
 * targets' (a target's particles' entity and compartment its own): `ParticleTileCatalog.getCompartmentEntities`'s, no
 * particle read
 */
function getCompartmentEntities(doc: ParticlePackDocument, compartmentInfo: ReadonlyMap<number, { name: string }>) {
  if (doc.compartments.length === 0) return undefined;
  const entities = new Map<number, Set<number>>();
  const withoutEntity = new Set<number>();
  for (const t of doc.targets) {
    if (t.particleCount === 0 || t.compartment < 0) continue;
    if (t.entity < 0) withoutEntity.add(t.compartment);
    else {
      let set = entities.get(t.compartment);
      if (!set) entities.set(t.compartment, (set = new Set()));
      set.add(t.entity);
    }
  }
  const rank = new Map<number, number>();
  for (const e of doc.entities) rank.set(e.index, rank.size);
  const rankOf = (e: number) => rank.get(e) ?? rank.size + e;
  const out = new Map<number, ParticleCompartmentInfo>();
  for (const [c, info] of compartmentInfo) {
    const ids = [...(entities.get(c) ?? [])].sort((a, b) => rankOf(a) - rankOf(b));
    out.set(c, { ...info, entities: ids, withoutEntity: withoutEntity.has(c) });
  }
  return out;
}

export interface ParticleTileSetFromPackOptions {
  /** the tile set's label, the pack's by default */
  readonly label?: string;
  /**
   * a target's sites every `siteStride`-th (`ParticleSites.coarsen`; 1, all, by default): its sites' bytes read, the
   * kept ones copied, the read dropped. Why: a view's targets the floor of its memory (the yeast's 9.5 M sites for a
   * tile of 6.8 k particles, 8 B a site read and about as many drawn), a tenth of it at 10
   */
  readonly siteStride?: number;
  /**
   * its coarse levels where written (`ParticlePackDocument.tiles.levels`; true by default): its nodes tiles of the set
   * (`ParticleTileSet.levels`), read as tiles are
   */
  readonly levels?: boolean;
}

/**
 * The tile set of tile-ordered `stream` (`isParticlePackStreamable`) on its written grid: its directory, its
 * opacities and its targets' runs read at once (the yeast's 2.3 M runs 18 MB); a tile's particles (its range of the
 * position, rotation and attribute blocks) read when its tile is asked (`provider.getTile`), a target's sites when it
 * is (`provider.getTarget`, once). A tile's list as a whole list's tile is (`ParticleTiles.fromList`): its keys its
 * particles' indices in the pack, its targets and their values runs, its bounds its centers' box padded by its radius
 */
export async function createParticleTileSetFromPack(
  stream: ParticlePackStream,
  options: ParticleTileSetFromPackOptions = {},
): Promise<ParticleTileSet> {
  const doc = stream.document;
  if (!isParticlePackStreamable(doc)) throw new Error('particle pack not read by range (not tile-ordered on a grid)');
  const tiles = doc.tiles!;
  const [ranges, boxes, runs] = await Promise.all([
    readWhole(stream, 'tile.ranges', Int32Array),
    readWhole(stream, 'tile.boxes', Float32Array),
    readWhole(stream, 'target.runs', Int32Array),
  ]);
  if (!ranges || !boxes || !runs) throw new Error('particle pack without its tile directory');
  const blocks = tiles.opacityBlocks;
  const opacity = blocks ? await readWhole(stream, 'tile.opacity', Float32Array) : undefined;
  const tileCount = ranges.length / ParticleTileDirectoryRange;
  if (boxes.length / ParticleTileDirectoryBox !== tileCount)
    throw new Error("particle pack 'tile.ranges' and 'tile.boxes' of other tile counts");

  const label = options.label ?? doc.label;
  const siteStride = Math.max(1, Math.floor(options.siteStride ?? 1));
  const info = getParticlePackInfo(doc);
  const siteTypes = getParticlePackSiteTypes(doc);
  const byId = new Map<number, ParticlePackTarget>(doc.targets.map((t) => [t.id, t]));
  const hasRadii = doc.targets.some((t) => t.radius !== undefined);
  let maxRadius = 0;
  for (const t of doc.targets) maxRadius = Math.max(maxRadius, t.radius ?? 0);

  // the catalog: a target's info from its row (its reach and volume as written, its sites not read)
  const targets = new Map<number, ParticleTileTargetInfo>();
  for (const t of doc.targets) {
    if (!t.sites) continue;
    targets.set(t.id, {
      id: t.id,
      kind: 'sites',
      sites: { quantized: true, elements: !!siteTypes.typeElements },
      radius: t.reach,
      volume: t.volume,
      byteCount: Math.ceil(t.sites.count / siteStride) * 8,
    });
  }
  const cell = doc.cell ? Cell.create(Vec3.ofArray(doc.cell.size), Vec3.ofArray(doc.cell.anglesInRadians)) : undefined;
  const catalog: ParticleTileCatalog = {
    targets,
    targetInfo: info.targetInfo,
    entityInfo: info.entityInfo.size > 0 ? info.entityInfo : undefined,
    compartmentInfo: getCompartmentEntities(doc, info.compartmentInfo),
    cell,
  };

  // the tiles, from the directory: a tile its cell's, its runs those of its range (the runs ascending, a tile's from its
  // first particle), its targets theirs
  const n = tiles.dimensions;
  const grid = ParticleTileGrid.create(Vec3.ofArray(tiles.min), Vec3.ofArray(tiles.max), [n, n, n]);
  const rotationBlock = getBlock(doc, 'rotation');
  const attributeBlocks = doc.attributes.map((a) => getBlock(doc, `attribute.${a.name}`));
  const bytesPerParticle =
    8 + (rotationBlock ? rotationBlock.byteLength / Math.max(1, rotationBlock.count) : 0) + attributeBlocks.length * 4;
  const runCount = runs.length / 2;
  const byCell = new Map<number, { readonly k: number; readonly r0: number; readonly r1: number }>();
  const tileList: ParticleTile[] = [];
  const opacityValues = blocks ? blocks ** 3 : 0;
  let r = 0;
  for (let k = 0; k < tileCount; ++k) {
    const id = ranges[k * 3];
    const first = ranges[k * 3 + 1];
    const count = ranges[k * 3 + 2];
    if (r < runCount && runs[r * 2] !== first) throw new Error("particle pack runs not from a tile's first particle");
    const r0 = r;
    while (r < runCount && runs[r * 2] < first + count) ++r;
    const ids = new Set<number>();
    for (let q = r0; q < r; ++q) ids.add(runs[q * 2 + 1]);
    byCell.set(id, { k, r0, r1: r });
    tileList.push({
      id,
      box: ParticleTileGrid.getBox(grid, id, Box3D()),
      count,
      base: first,
      targets: Int32Array.from(ids).sort(),
      opacity:
        blocks && opacity
          ? {
              dimensions: [blocks, blocks, blocks],
              density: opacity.subarray(k * opacityValues, (k + 1) * opacityValues),
            }
          : undefined,
      byteCount: count * bytesPerParticle,
    });
  }
  if (r !== runCount) throw new Error('particle pack runs past its tiles');
  tileList.sort((a, b) => a.id - b.id);
  const tileMap = new Map<number, ParticleTile>(tileList.map((t) => [t.id, t]));

  // the sphere all are within: the tiles' centers' boxes padded by the largest radius
  const helper = new BoundaryHelper('98');
  const sphere = Sphere3D();
  const tileSphere = (k: number) => {
    const b = k * ParticleTileDirectoryBox;
    Vec3.set(
      sphere.center,
      (boxes[b] + boxes[b + 3]) / 2,
      (boxes[b + 1] + boxes[b + 4]) / 2,
      (boxes[b + 2] + boxes[b + 5]) / 2,
    );
    sphere.radius =
      Math.hypot(boxes[b + 3] - boxes[b], boxes[b + 4] - boxes[b + 1], boxes[b + 5] - boxes[b + 2]) / 2 + maxRadius;
    return sphere;
  };
  helper.reset();
  for (let k = 0; k < tileCount; ++k) helper.includeSphere(tileSphere(k));
  helper.finishedIncludeStep();
  for (let k = 0; k < tileCount; ++k) helper.radiusSphere(tileSphere(k));
  const boundingSphere = tileCount > 0 ? helper.getSphere() : undefined;

  const positionBlock = getBlock(doc, 'position')!;
  const positionGrid = doc.positionGrid!;
  const gridOrigin = Vec3.ofArray(positionGrid.origin);
  const format = ParticlePackStreamFormat.create(stream);

  /** a tile's or a node's blocks (`ParticleTileLevels`): its particles' positions, rotations and attributes */
  interface RangeBlocks {
    readonly position: ParticlePackBlock;
    readonly rotation: ParticlePackBlock | undefined;
    readonly attributes: readonly (ParticlePackBlock | undefined)[];
    /** a particle's index in the pack (a node's), none: its place in the blocks */
    readonly keys?: ParticlePackBlock;
    /** its runs (first, target), ascending */
    readonly runs: Int32Array;
    /** its boxes (`ParticleTileDirectoryBox` a range) */
    readonly boxes: Float32Array;
  }
  const tileBlocks: RangeBlocks = {
    position: positionBlock,
    rotation: rotationBlock,
    attributes: attributeBlocks,
    runs,
    boxes,
  };

  /** tile `id`'s list: its ranges of the per-particle blocks read, its runs' values */
  const readTile = (id: number): Promise<ParticleList> => {
    const tile = tileMap.get(id);
    const at = byCell.get(id);
    if (!tile || !at) throw new Error(`particle tiles: no tile ${id}`);
    return readRange(`tile ${id}`, tileBlocks, tile.base, tile.count, at, 1);
  };

  /**
   * a range's list (`first`, `count` of `ranges`' blocks; its runs `at`, its box `at.k`): its ranges of the
   * per-particle blocks read, its runs' values, its particles drawn scaled by `scale`
   */
  const readRange = async (
    name: string,
    ranges: RangeBlocks,
    first: number,
    count: number,
    at: { readonly k: number; readonly r0: number; readonly r1: number },
    scale: number,
  ): Promise<ParticleList> => {
    const runs = ranges.runs;
    const boxes = ranges.boxes;
    const [positions, rotations, keyValues, ...attributeValues] = await Promise.all([
      readItems(stream.source, ranges.position, first, count),
      ranges.rotation ? readItems(stream.source, ranges.rotation, first, count) : Promise.resolve(undefined),
      ranges.keys ? readItems(stream.source, ranges.keys, first, count) : Promise.resolve(undefined),
      ...ranges.attributes.map((b) => (b ? readItems(stream.source, b, first, count) : Promise.resolve(undefined))),
    ]);
    if (!(positions instanceof Uint32Array)) throw new Error("particle pack block 'position' not on a grid");
    if (rotations instanceof Int32Array || rotations instanceof Uint32Array)
      throw new Error("particle pack block 'rotation' not snorm16 or float32");

    const runsIn = at.r1 - at.r0;
    const offsets = new Int32Array(runsIn + 1);
    const ids = new Int32Array(runsIn);
    const entityValues = new Int32Array(runsIn);
    const compartmentValues = new Int32Array(runsIn);
    const radiusValues = new Float32Array(runsIn);
    let radius = 0;
    for (let q = 0; q < runsIn; ++q) {
      const t = byId.get(runs[(at.r0 + q) * 2 + 1]);
      if (!t) throw new Error(`particle pack run of target ${runs[(at.r0 + q) * 2 + 1]}, not a target`);
      offsets[q] = runs[(at.r0 + q) * 2] - first;
      ids[q] = t.id;
      entityValues[q] = t.entity;
      compartmentValues[q] = t.compartment;
      radiusValues[q] = t.radius ?? 0;
      radius = Math.max(radius, radiusValues[q]);
    }
    offsets[runsIn] = count;
    const targetValues = ParticleValues.ofRuns(offsets, ids);
    const entities = doc.entities.length > 0 ? ParticleValues.ofRuns(offsets, entityValues) : undefined;
    let attributes: Map<string, ParticleAttribute> | undefined;
    for (let a = 0; a < doc.attributes.length; ++a) {
      const values = attributeValues[a];
      if (!(values instanceof Float32Array)) continue;
      attributes ??= new Map();
      const { name, label: attributeLabel, min, max } = doc.attributes[a];
      attributes.set(name, { column: Column.ofFloatArray(values), label: attributeLabel, min, max });
    }
    let keys: Int32Array;
    if (keyValues instanceof Int32Array) {
      keys = keyValues;
    } else {
      keys = new Int32Array(count);
      for (let i = 0; i < count; ++i) keys[i] = first + i;
    }
    const list: ParticleList = {
      entryId: doc.entryId,
      label: `${label ?? 'Particles'} ${name}`,
      count,
      cell,
      keys,
      targets: targetValues,
      targetInfo: info.targetInfo,
      entities,
      entityInfo: catalog.entityInfo,
      compartments: doc.compartments.length > 0 ? ParticleValues.ofRuns(offsets, compartmentValues) : undefined,
      compartmentInfo: catalog.compartmentInfo,
      coordinates: { kind: 'grid', data: positions, origin: gridOrigin, step: positionGrid.step },
      rotations: rotations as Float32Array | Int16Array | undefined,
      radii: hasRadii ? ParticleValues.ofRuns(offsets, radiusValues) : undefined,
      attributes,
      getParticleLabel: (index: number) => {
        const parts = [`Particle ${keys[index] + 1}`];
        const entity = entities ? info.entityInfo.get(ParticleValues.get(entities, index))?.name : undefined;
        if (entity) parts.push(entity);
        const name = info.targetName.get(ParticleValues.get(targetValues, index));
        if (name !== undefined) parts.push(`target ${name}`);
        return parts.join(' · ');
      },
      sourceData: format,
      scale: scale !== 1 ? scale : undefined,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    // its bounds: its centers' box padded by its largest radius, not a pass over it (as a whole list's tile)
    const b = at.k * ParticleTileDirectoryBox;
    const box = Box3D.create(
      Vec3.create(boxes[b] - radius, boxes[b + 1] - radius, boxes[b + 2] - radius),
      Vec3.create(boxes[b + 3] + radius, boxes[b + 4] + radius, boxes[b + 5] + radius),
    );
    const tileBounds = Sphere3D.create(
      Vec3.create((boxes[b] + boxes[b + 3]) / 2, (boxes[b + 1] + boxes[b + 4]) / 2, (boxes[b + 2] + boxes[b + 5]) / 2),
      Math.hypot(boxes[b + 3] - boxes[b], boxes[b + 4] - boxes[b + 1], boxes[b + 5] - boxes[b + 2]) / 2 + radius,
    );
    Particle.setBoundary(list, { box, sphere: tileBounds });
    return list;
  };

  // a target's sites read once (the targets' tiles share it)
  const siteBlock = getBlock(doc, 'site');
  const targetReads = new Map<number, Promise<ParticleTarget>>();
  const readTarget = (id: number): Promise<ParticleTarget> => {
    let read = targetReads.get(id);
    if (!read) {
      const s = byId.get(id)?.sites;
      if (!s || !siteBlock) return Promise.reject(new Error(`particle tiles: no target ${id}`));
      read = readItems(stream.source, siteBlock, s.offset, s.count).then((data) => {
        if (!(data instanceof Int16Array)) throw new Error("particle pack block 'site' not int16");
        const target = createParticlePackSitesTarget(s, data, siteTypes);
        return siteStride > 1 && target.kind === 'sites'
          ? { ...target, sites: ParticleSites.coarsen(target.sites, siteStride) }
          : target;
      });
      // a failed read asked again later
      read.catch(() => targetReads.delete(id));
      targetReads.set(id, read);
    }
    return read;
  };

  const base: ParticleTileSet = {
    label,
    catalog,
    grid,
    tiles: tileMap,
    count: doc.count,
    boundingSphere,
    frameCount: 1,
    provider: {
      getTile: (id, frame) =>
        Task.create(`Particle tile ${id}`, async () => {
          if (frame !== 0) throw new Error(`particle tiles: no frame ${frame}`);
          return readTile(id);
        }),
      getTarget: (id) => Task.create(`Particle target ${id}`, () => readTarget(id)),
    },
  };
  if (!(tiles.levels && tiles.levels > 0) || options.levels === false) return base;

  // the coarse levels: a node a tile of the set, its slots after the particles' (its base the pack's count and its
  // first), its id after the grid's cells
  const [nodeRanges, nodeBoxes, nodeRuns] = await Promise.all([
    readWhole(stream, 'node.ranges', Int32Array),
    readWhole(stream, 'node.boxes', Float32Array),
    readWhole(stream, 'node.runs', Int32Array),
  ]);
  const nodePosition = getBlock(doc, 'node.position');
  const nodeKeys = getBlock(doc, 'node.keys');
  if (!nodeRanges || !nodeBoxes || !nodeRuns || !nodePosition || !nodeKeys) {
    throw new Error('particle pack levels without their nodes');
  }
  const nodeBlocks: RangeBlocks = {
    position: nodePosition,
    rotation: getBlock(doc, 'node.rotation'),
    attributes: doc.attributes.map((a) => getBlock(doc, `node.attribute.${a.name}`)),
    keys: nodeKeys,
    runs: nodeRuns,
    boxes: nodeBoxes,
  };
  const cells = n * n * n;
  const nodeCount = nodeRanges.length / 4;
  const allTiles = new Map<number, ParticleTile>(tileMap);
  const nodeAt = new Map<
    number,
    { readonly k: number; readonly r0: number; readonly r1: number; readonly level: number }
  >();
  const byLevel = new Map<string, number>();
  const nodeRunCount = nodeRuns.length / 2;
  let nodeParticles = 0;
  r = 0;
  for (let k = 0; k < nodeCount; ++k) {
    const level = nodeRanges[k * 4];
    const levelCell = nodeRanges[k * 4 + 1];
    const first = nodeRanges[k * 4 + 2];
    const count = nodeRanges[k * 4 + 3];
    if (r < nodeRunCount && nodeRuns[r * 2] !== first)
      throw new Error("particle pack node runs not from a node's first");
    const r0 = r;
    while (r < nodeRunCount && nodeRuns[r * 2] < first + count) ++r;
    const ids = new Set<number>();
    for (let q = r0; q < r; ++q) ids.add(nodeRuns[q * 2 + 1]);
    const id = cells + k;
    nodeAt.set(id, { k, r0, r1: r, level });
    byLevel.set(`${level},${levelCell}`, id);
    const b = k * ParticleTileDirectoryBox;
    allTiles.set(id, {
      id,
      box: Box3D.create(
        Vec3.create(nodeBoxes[b], nodeBoxes[b + 1], nodeBoxes[b + 2]),
        Vec3.create(nodeBoxes[b + 3], nodeBoxes[b + 4], nodeBoxes[b + 5]),
      ),
      count,
      base: doc.count + first,
      targets: Int32Array.from(ids).sort(),
      byteCount: count * (bytesPerParticle + 4),
    });
    nodeParticles = Math.max(nodeParticles, first + count);
  }
  const readNode = (id: number): Promise<ParticleList> => {
    const at = nodeAt.get(id);
    const tile = allTiles.get(id);
    if (!at || !tile) throw new Error(`particle tiles: no node ${id}`);
    const first = tile.base - doc.count;
    return readRange(
      `node ${at.level}/${id - cells}`,
      nodeBlocks,
      first,
      tile.count,
      at,
      ParticleTiles.getLevelScale(at.level),
    );
  };
  return {
    ...base,
    tiles: allTiles,
    // the slots: the particles', then the nodes'
    count: doc.count + nodeParticles,
    levels: {
      base,
      count: tiles.levels,
      getNode: (level, levelCell) => byLevel.get(`${level},${levelCell}`),
    },
    provider: {
      getTile: (id, frame) =>
        Task.create(`Particle tile ${id}`, async () => {
          if (frame !== 0) throw new Error(`particle tiles: no frame ${frame}`);
          return id < cells ? readTile(id) : readNode(id);
        }),
      getTarget: base.provider.getTarget,
    },
  };
}
