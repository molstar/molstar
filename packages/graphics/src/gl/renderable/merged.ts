/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { OcclusionPyramid } from '../compute/cell-cull.js';
import { createMergedCuller, type MergedRasterParams } from './merged-culler.js';
import { isSphereRasterEligible } from '../compute/sphere-raster-pass.js';
import { packAffineRows } from '../affine-rows.js';
import { OctNormalItemSize } from '../oct-normals.js';
import { ValueCell } from '@molstar/core/util';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { createEmptyInstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { createDataArray, type DataArray } from './util.js';
import { SegmentedArray } from '../segmented-array.js';
import {
  GlobalUniformSchema,
  GlobalTextureSchema,
  GlobalDefineSchema,
  InternalSchema,
  type InternalValues,
  SegmentSchema,
  type SegmentValues,
  type GlobalDefines,
  type GlobalDefineValues,
  type RenderableSchema,
  type RenderableValues,
  type BaseValues,
} from './schema.js';
import { createSphereResolveValues, SpheresSchema } from './spheres.js';
import { MeshSchema } from './mesh.js';
import { CylindersSchema } from './cylinders.js';
import { LinesSchema } from './lines.js';
import { PointsSchema } from './points.js';
import {
  SpheresShaderCode,
  MeshShaderCode,
  CylindersShaderCode,
  LinesShaderCode,
  PointsShaderCode,
  type ShaderCode,
} from '../shader-code.js';
import type {
  DrawMode,
  GlobalBindings,
  GraphicsRenderVariant,
  MultiDrawBaseData,
  Transparency,
} from '../render-item.js';
import type { GraphicsContext } from '../context.js';

import {
  CellCullKernel,
  CellCullSchema,
  CullArgsKernel,
  CullArgsSchema,
  InstanceCopyKernel,
  InstanceCopySchema,
} from '../compute/cell-cull.js';
import { getMaxItemCount, isItemIndexing } from '../quad-indices.js';
import type { GraphicsStats } from '../stats.js';
import {
  type Renderable,
  type RenderableState,
  appendMdbDrawKey,
  createDrawVersion,
  createCullView,
  createSegmentedMdbList,
  type CullSegment,
  type LodLevelsValue,
  LodView,
  type CullValues,
  type Frame,
  createClipInfoCache,
} from '../renderable.js';

export type MergeableValues = RenderableValues & BaseValues;

export const MergedTypes = ['spheres', 'mesh', 'cylinders', 'lines', 'points'] as const;

/**
 * Vec4s per instance of merged instance data read by index (the WebGL2 `dInstanceData` texture):
 * the transform's columns, (instance, segment), the segment's sphere, its level scales
 */
export const InstanceDataVec4s = 5;

/** Packs the instance data of merged `values` into `out` (at least instances times `InstanceDataVec4s` vec4s) */
export function packInstanceData(values: RenderableValues, out: Float32Array) {
  const count: number = values.instanceCount.ref.value;
  const transform = values.aTransform.ref.value as Float32Array;
  const instance = values.aInstance.ref.value as Float32Array;
  const segment = values.aSegment.ref.value as Float32Array;
  for (let i = 0; i < count; ++i) {
    const o = i * InstanceDataVec4s * 4;
    // the transform's rows (`affine-rows.ts`), the instance and its segment's groups, the segment's index
    packAffineRows(transform.subarray(i * 16, i * 16 + 16), 1, out, o);
    out[o + 12] = instance[i];
    out[o + 13] = segment[i * 4];
    out[o + 14] = segment[i * 4 + 1];
    out[o + 15] = segment[i * 4 + 2];
    out[o + 16] = segment[i * 4 + 3];
    out[o + 17] = 0;
    out[o + 18] = 0;
    out[o + 19] = 0;
  }
}
export type MergedType = (typeof MergedTypes)[number];

export function isMergeableType(type: string): type is MergedType {
  return MergedTypes.includes(type as MergedType);
}

//

type Granularity = 'instance' | 'group' | 'groupInstance' | 'vertex' | 'item' | 'none';

type MergedTexSpec = {
  tex: string;
  itemSize: number;
  float?: boolean;
  granularity: (v: MergeableValues) => Granularity;
  /** texel count multiplier, e.g. for dual colors */
  factor?: (v: MergeableValues) => number;
};

function pickGranularity(type: string): Granularity {
  switch (type) {
    case 'instance':
      return 'instance';
    case 'group':
      return 'group';
    case 'groupInstance':
      return 'groupInstance';
    case 'vertex':
      return 'vertex';
    default:
      return 'none';
  }
}

const BaseTexSpecs: MergedTexSpec[] = [
  {
    tex: 'tColor',

    itemSize: 4,
    granularity: (v) => pickGranularity(v.colorType.ref.value),
    factor: (v) => (v.dDualColor?.ref.value ? 2 : 1),
  },
  { tex: 'tSize', itemSize: 4, granularity: (v) => pickGranularity(v.sizeType.ref.value) },
  { tex: 'tMarker', itemSize: 1, granularity: (v) => pickGranularity(v.markerType.ref.value) },
  {
    tex: 'tOverpaint',

    itemSize: 4,
    granularity: (v) => (v.dOverpaint.ref.value ? pickGranularity(v.overpaintType.ref.value) : 'none'),
  },
  {
    tex: 'tTransparency',

    itemSize: 1,
    granularity: (v) => (v.dTransparency.ref.value ? pickGranularity(v.transparencyType.ref.value) : 'none'),
  },
  {
    tex: 'tEmissive',

    itemSize: 1,
    granularity: (v) => (v.dEmissive.ref.value ? pickGranularity(v.emissiveType.ref.value) : 'none'),
  },
  {
    tex: 'tSubstance',

    itemSize: 4,
    granularity: (v) => (v.dSubstance.ref.value ? pickGranularity(v.substanceType.ref.value) : 'none'),
  },
  {
    tex: 'tClipping',

    itemSize: 1,
    granularity: (v) => (v.dClipping.ref.value ? pickGranularity(v.clippingType.ref.value) : 'none'),
  },
  {
    tex: 'tWiggle',

    itemSize: 1,
    granularity: (v) => (v.dWiggle.ref.value ? pickGranularity(v.wiggleType.ref.value) : 'none'),
  },
];

type VertexAttribute = { key: string; itemSize: number };

function getVertexAttributes(schema: RenderableSchema): VertexAttribute[] {
  const ret: VertexAttribute[] = [];
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    // `aVertex` (the emulated vertex id of WebGL1) is not member data: the render item fills it with serial ids
    if (k === 'aVertex') continue;
    if (spec.type === 'attribute' && spec.divisor === 0) {
      // the members' arrays as they are on the CPU: encoded normals three floats a vertex (`oct-normals.ts`)
      ret.push({ key: k, itemSize: spec.kind === 'oct16' ? OctNormalItemSize : spec.itemSize });
    }
  }
  return ret;
}

type MergedTypeDescriptor = {
  readonly drawMode: DrawMode;
  readonly shaderCode: ShaderCode;
  readonly schema: RenderableSchema;
  readonly hasElements: boolean;
  /** texel count of 'item'-granularity textures: the impostors' data, from their items drawn (`drawCount`) */
  itemTexels: (v: MergeableValues) => number;
  /** kind-specific textures to merge */
  readonly extraTextures: MergedTexSpec[];
  readonly vertexAttributes: VertexAttribute[];
};

/** Whether merged objects can be culled on the GPU (`gpuCulling`): compute and indirect draws by first instance */
export function canCullOnGpu(ctx: GraphicsContext | undefined) {
  const compute = ctx?.compute;
  return (
    !!compute &&
    !!ctx.capabilities.indirectInstanceRanges &&
    // the instance copy binds 8 storage buffers, the core limits' count
    compute.isKernelSupported(CellCullKernel, CellCullSchema) &&
    compute.isKernelSupported(InstanceCopyKernel, InstanceCopySchema) &&
    compute.isKernelSupported(CullArgsKernel, CullArgsSchema)
  );
}

const MergedTypeDescriptors: { [k in MergedType]: MergedTypeDescriptor } = {
  spheres: {
    drawMode: 'triangles',
    shaderCode: SpheresShaderCode,
    schema: SpheresSchema,
    hasElements: false,
    // one position texel per sphere, drawn as 6 indices
    itemTexels: (v) => v.drawCount.ref.value / 6,
    extraTextures: [{ tex: 'tPositionGroup', itemSize: 4, float: true, granularity: () => 'item' }],
    vertexAttributes: getVertexAttributes(SpheresSchema),
  },
  mesh: {
    drawMode: 'triangles',
    shaderCode: MeshShaderCode,
    schema: MeshSchema,
    hasElements: true,
    itemTexels: () => 0,
    extraTextures: [],
    vertexAttributes: getVertexAttributes(MeshSchema),
  },
  cylinders: {
    drawMode: 'triangles',
    shaderCode: CylindersShaderCode,
    schema: CylindersSchema,
    hasElements: false,
    // three data texels per cylinder, drawn as 12 indices
    itemTexels: (v) => (v.drawCount.ref.value / 12) * 3,
    extraTextures: [{ tex: 'tCylinder', itemSize: 4, float: true, granularity: () => 'item' }],
    vertexAttributes: getVertexAttributes(CylindersSchema),
  },
  lines: {
    drawMode: 'triangles',
    shaderCode: LinesShaderCode,
    schema: LinesSchema,
    hasElements: false,
    // two data texels per line, drawn as 6 indices
    itemTexels: (v) => (v.drawCount.ref.value / 6) * 2,
    extraTextures: [{ tex: 'tLine', itemSize: 4, float: true, granularity: () => 'item' }],
    vertexAttributes: getVertexAttributes(LinesSchema),
  },
  points: {
    drawMode: 'triangles',
    shaderCode: PointsShaderCode,
    schema: PointsSchema,
    hasElements: false,
    // one position texel per point, drawn as 6 indices
    itemTexels: (v) => v.drawCount.ref.value / 6,
    extraTextures: [{ tex: 'tPositionGroup', itemSize: 4, float: true, granularity: () => 'item' }],
    vertexAttributes: getVertexAttributes(PointsSchema),
  },
};

//

function isDefineKey(k: string) {
  return k.length > 1 && k[0] === 'd' && k[1] >= 'A' && k[1] <= 'Z';
}

/**
 * Signature of everything that must be equal across members sharing a merged
 * render item: define values, granularity resolution, and lod configuration.
 *
 * Note: `lodLevels[j]`'s `scale` (index 4) is intentionally excluded - it is
 * data-dependent per member (see `Spheres.getLodLevelsValue`, based on that
 * member's own sphere count), not part of the user-configured lod preset
 * (only min/max/overlap are). Per-member scales are instead carried in the
 * segment data (`tSegment`, see `SegmentSchema`) and selected in the shader
 * via `uLodLevel`, so they never need to match here.
 */
export function getMergeKey(values: MergeableValues): string {
  // checked per member on every update of a merged representation: computed again only where a cell it reads changed
  let cached = mergeKeys.get(values);
  if (!cached) {
    const keys = intern(internedDefineKeys, Object.keys(values).filter(isDefineKey).sort());
    const cells = [...keys.map((k) => values[k]), values.uInstanceGranularity, values.uLod];
    if (values.lodLevels) cells.push(values.lodLevels);
    // the index of `uLod` (compared by its components), not looked up per cell: values have too many keys for fast reads
    const lod = keys.length + 1;
    cached = { keys, cells, lod, versions: cells.map(() => -1), seen: cells.map(() => undefined), key: '' };
    mergeKeys.set(values, cached);
  }
  // updates bump versions of cells whose values stay (`ValueCell.update` of the same define): compared by value
  const { cells, versions, seen, lod } = cached;
  let changed = false;
  for (let i = 0, il = cells.length; i < il; ++i) {
    const cell = cells[i];
    if (cell.ref.version === versions[i]) continue;
    versions[i] = cell.ref.version;
    const value = cell.ref.value;
    if (i === lod) {
      const last = seen[i] as Vec4 | undefined;
      if (!last || !Vec4.exactEquals(last, value as Vec4)) {
        seen[i] = Vec4.clone(value as Vec4);
        changed = true;
      }
    } else if (value !== seen[i]) {
      seen[i] = value;
      changed = true;
    }
  }
  if (changed) cached.key = intern(internedMergeKeys, computeMergeKey(values, cached.keys));
  return cached.key;
}

/**
 * Equal merge keys as one string and equal define key lists as one array: kept per member (100k members, several
 * hundred characters each). Strings compare by value, so interning changes no result; cleared when large, the kept
 * ones then stay distinct copies
 */
const internedMergeKeys = new Map<string, string>();
const internedDefineKeys = new Map<string, string[]>();

function intern<T extends string | string[]>(table: Map<string, T>, value: T): T {
  const key = typeof value === 'string' ? value : value.join();
  const interned = table.get(key);
  if (interned !== undefined) return interned;
  if (table.size >= 1024) table.clear();
  table.set(key, value);
  return value;
}

/** `getMergeKey` per values, with the define keys and the cells it reads and their versions when it was computed */
const mergeKeys = new WeakMap<
  MergeableValues,
  { keys: string[]; cells: ValueCell<unknown>[]; lod: number; versions: number[]; seen: unknown[]; key: string }
>();

function computeMergeKey(values: MergeableValues, defineKeys: readonly string[]) {
  const parts: string[] = [];
  for (const k of defineKeys) {
    parts.push(`${k}:${values[k].ref.value}`);
  }
  parts.push(`ig:${values.uInstanceGranularity.ref.value}`);
  parts.push(`uLod:${(values.uLod.ref.value as Vec4).join(',')}`);
  const lodLevels = values.lodLevels?.ref.value as LodLevelsValue | undefined;
  parts.push(
    `lod${lodLevels?.pixelSize ? `px${lodLevels.pixelSize}` : ''}:` +
      (lodLevels ? lodLevels.map((l) => `${l[0]},${l[1]},${l[2]}`).join(';') : ''),
  );
  return parts.join('|');
}

const SupportedVertexTypes = ['uniform', 'instance', 'group', 'groupInstance', 'vertex'];
const SupportedIgTypes = ['instance', 'groupInstance'];

export function canMergeValues(type: string, members: readonly MergeableValues[]): boolean {
  if (!isMergeableType(type)) return false;
  // a single member merges as well (culled on the GPU): how many is the caller's (`canMergeRenderObjects`)
  if (members.length < 1) return false;

  const first = members[0];
  if (first.dUsePalette?.ref.value) return false;
  if (!SupportedVertexTypes.includes(first.colorType.ref.value)) return false;
  if (first.sizeType && !SupportedVertexTypes.includes(first.sizeType.ref.value)) return false;
  if (first.dOverpaint.ref.value && !SupportedIgTypes.includes(first.overpaintType.ref.value)) return false;
  if (first.dTransparency.ref.value && !SupportedIgTypes.includes(first.transparencyType.ref.value)) return false;
  if (first.dEmissive.ref.value && !SupportedIgTypes.includes(first.emissiveType.ref.value)) return false;
  if (first.dSubstance.ref.value && !SupportedIgTypes.includes(first.substanceType.ref.value)) return false;
  if (first.dClipping.ref.value && !SupportedIgTypes.includes(first.clippingType.ref.value)) return false;
  if (first.dWiggle.ref.value && !SupportedIgTypes.includes(first.wiggleType.ref.value)) return false;

  const key = getMergeKey(first);
  for (let i = 1, il = members.length; i < il; ++i) {
    if (getMergeKey(members[i]) !== key) return false;
  }
  return true;
}

export type MergedSegments = {
  count: number;
  /** first instance (into the merged instanced attributes) per segment */
  instanceBases: number[];
  /** instance count per segment */
  instanceCounts: number[];
  /** first vertex (into the merged geometry) per segment */
  vertexBases: number[];
  /** first drawn vertex per segment, for types without elements (e.g. 6 per point quad with one data vertex per point) */
  drawBases: number[];
  /** byte offset (into the merged element buffer) per segment */
  elementOffsets: number[];
};

export type Merged = {
  readonly type: MergedType;
  readonly values: MergeableValues & SegmentValues;
  readonly segments: MergedSegments;
  /** sync merged values with member values, returns true if the layout changed */
  sync: () => boolean;
  /** the members' marker cells `sync` reads, observed with the merged values' own (`Renderable.memberCells`) */
  readonly memberCells: readonly ValueCell<unknown>[];
  /**
   * advanced by `sync` where what culling reads of the members changed (layout, transforms and grids, spheres, levels):
   * a still frame's lists stand while it stays, without comparing the members' cells (six a member)
   */
  readonly cullVersion: ValueCell<number>;
};

function ensureArray<T extends Uint8Array | Uint32Array | Float32Array>(
  array: T,
  length: number,
  ctor: new (length: number) => T,
): T {
  return array.length >= length ? array : new ctor(length);
}

const boundaryHelper = new BoundaryHelper('98');

type Cell = ValueCell<any>;

/** The cells of a member `sync` reads, see `getMemberCells` */
interface MemberCells {
  readonly instanceCount: Cell;
  readonly uGroupCount: Cell;
  readonly uVertexCount: Cell;
  readonly drawCount: Cell;
  readonly aTransform: Cell;
  readonly aInstance: Cell;
  readonly boundingSphere: Cell;
  readonly invariantBoundingSphere: Cell;
  readonly markerStatus: Cell;
  readonly markerAverage: Cell;
  readonly transparencyAverage: Cell;
  readonly transparencyMin: Cell;
  readonly emissiveAverage: Cell;
  readonly wiggleAverage: Cell;
  readonly hasReflection: Cell;
  readonly lodLevels: Cell | undefined;
  readonly elements: Cell | undefined;
  /** in the order of the type's `vertexAttributes` */
  readonly attributes: readonly Cell[];
  /** in the order of the merged textures (`texSpecs`) */
  readonly textures: readonly Cell[];
}

/** The cells of `v`: they stay for the values' lifetime (updates replace their values, not the cells) */
function getMemberCells(
  v: MergeableValues,
  vertexAttributes: readonly { key: string }[],
  hasElements: boolean,
  texSpecs: readonly MergedTexSpec[],
): MemberCells {
  return {
    instanceCount: v.instanceCount,
    uGroupCount: v.uGroupCount,
    uVertexCount: v.uVertexCount,
    drawCount: v.drawCount,
    aTransform: v.aTransform,
    aInstance: v.aInstance,
    boundingSphere: v.boundingSphere,
    invariantBoundingSphere: v.invariantBoundingSphere,
    markerStatus: v.markerStatus,
    markerAverage: v.markerAverage,
    transparencyAverage: v.transparencyAverage,
    transparencyMin: v.transparencyMin,
    emissiveAverage: v.emissiveAverage,
    wiggleAverage: v.wiggleAverage,
    hasReflection: v.hasReflection,
    lodLevels: v.lodLevels,
    elements: hasElements ? v.elements : undefined,
    attributes: vertexAttributes.map((a) => v[a.key]),
    textures: texSpecs.map((s) => v[s.tex]),
  };
}

/** The versions of a member's cells `sync` saw last (-1: none) */
function createMemberVersions(attributeCount: number, textureCount: number) {
  return {
    aTransform: -1,
    aInstance: -1,
    boundingSphere: -1,
    invariantBoundingSphere: -1,
    markerStatus: -1,
    markerAverage: -1,
    transparencyAverage: -1,
    transparencyMin: -1,
    emissiveAverage: -1,
    wiggleAverage: -1,
    hasReflection: -1,
    lodLevels: -1,
    elements: -1,
    attributes: new Int32Array(attributeCount).fill(-1),
    textures: new Int32Array(textureCount).fill(-1),
  };
}

/** `gpuCulling`: culled on the GPU where available (`merged-culler.ts`), fixed for the merged object's lifetime */
export function createMergedValues(type: MergedType, members: readonly MergeableValues[], gpuCulling = false): Merged {
  const descriptor = MergedTypeDescriptors[type];
  const { hasElements, vertexAttributes } = descriptor;
  const first = members[0];
  const n = members.length;

  const texSpecs = [...descriptor.extraTextures, ...BaseTexSpecs].filter((spec) => first[spec.tex] !== undefined);

  const segments: MergedSegments = {
    count: n,
    instanceBases: new Array(n).fill(0),
    instanceCounts: new Array(n).fill(0),
    vertexBases: new Array(n).fill(0),
    drawBases: new Array(n).fill(0),
    elementOffsets: new Array(n).fill(0),
  };

  // own cells, identity must stay stable for the lifetime of the merged render item
  const own: { [k: string]: ValueCell<any> } = {
    aTransform: ValueCell.create(new Float32Array(0)),
    aInstance: ValueCell.create(new Float32Array(0)),
    aSegment: ValueCell.create(new Float32Array(0)),
    tSegment: ValueCell.create(createDataArray(n * 2, 4, Float32Array)),
    aVertex: ValueCell.create(new Float32Array(0)),
    drawCount: ValueCell.create(0),
    instanceCount: ValueCell.create(0),
    uVertexCount: ValueCell.create(0),
    uInstanceCount: ValueCell.create(0),
    uGroupCount: ValueCell.create(0),
    boundingSphere: ValueCell.create(Sphere3D()),
    invariantBoundingSphere: ValueCell.create(Sphere3D()),
    uInvariantBoundingSphere: ValueCell.create(Vec4()),
    uMarker: ValueCell.create(0),
    markerAverage: ValueCell.create(0),
    markerStatus: ValueCell.create(-1),
    transparencyAverage: ValueCell.create(0),
    transparencyMin: ValueCell.create(0),
    emissiveAverage: ValueCell.create(0),
    wiggleAverage: ValueCell.create(0),
    hasReflection: ValueCell.create(false),
    instanceGrid: ValueCell.create({ ...createEmptyInstanceGrid(), cellSize: 2 }),
    dSegmented: ValueCell.create(true),
    gpuCulling: ValueCell.create(gpuCulling),
  };
  // the members' geometry is not copied: the render item uploads their arrays as parts (`SegmentedArray`)
  for (const { key } of vertexAttributes) {
    own[key] = ValueCell.create(SegmentedArray.create(Float32Array, []));
  }
  if (hasElements) {
    own.elements = ValueCell.create(SegmentedArray.create(Uint32Array, []));
  }
  for (const spec of texSpecs) {
    const ctor = (spec.float ? Float32Array : Uint8Array) as new (length: number) => any;
    own[spec.tex] = ValueCell.create(createDataArray(1, spec.itemSize, ctor));
  }

  // the sphere raster's values set per view, the merged object's own (not its first member's)
  if ('dSphereResolve' in first) Object.assign(own, createSphereResolveValues());
  const values = {
    ...first,
    ...own,
  } as unknown as MergeableValues & SegmentValues;

  /**
   * Extent of a member's 'group'-granularity data. `uGroupCount` comes from the
   * visual's location iterator and can understate it - e.g. particle target
   * structures iterate per particle (groupCount 1) but carry per-element sizes.
   */
  function getGroupSpan(i: number): number {
    const v = members[i];
    let span = v.uGroupCount.ref.value;
    for (const spec of texSpecs) {
      if (spec.granularity(first) !== 'group') continue;
      const { count } = v[spec.tex].ref.value as DataArray<any>;
      if (count > span) span = count;
    }
    return span;
  }

  function getTexelCount(i: number, spec: MergedTexSpec): number {
    const v = members[i];
    const factor = spec.factor ? spec.factor(v) : 1;
    switch (spec.granularity(first)) {
      case 'instance':
        return v.instanceCount.ref.value * factor;
      case 'group':
        return getGroupSpan(i) * factor;
      case 'groupInstance':
        return v.instanceCount.ref.value * v.uGroupCount.ref.value * factor;
      case 'vertex':
        // the data vertices (`Geometry.getVertexCount`)
        return v.uVertexCount.ref.value * factor;
      case 'item':
        return descriptor.itemTexels(v) * factor;
      case 'none':
        return 0;
    }
  }

  /** texel base per member for each merged texture */
  const texBases = new Map<string, number[]>();
  for (const spec of texSpecs) {
    texBases.set(spec.tex, new Array(n).fill(0));
  }

  // the cells `sync` reads per member, resolved once into objects of one shape: values objects have too many keys for
  // fast named reads, and a sync reads some 30 cells of each member
  const memberCells = members.map((v) => getMemberCells(v, vertexAttributes, hasElements, texSpecs));
  const memberVersions = members.map(() => createMemberVersions(vertexAttributes.length, texSpecs.length));

  /** per member its instance, group, group span, vertex and draw counts: a change is a new layout */
  const layout = new Float64Array(n * 5).fill(-1);
  let layoutMergeKey = '';

  /** Whether the layout changed since the last call (also the first member's merge key), stores the new one */
  function updateLayout() {
    let changed = false;
    const groupTextures: number[] = [];
    for (let t = 0; t < texSpecs.length; ++t) {
      if (texSpecs[t].granularity(first) === 'group') groupTextures.push(t);
    }
    for (let i = 0; i < n; ++i) {
      const c = memberCells[i];
      let span: number = c.uGroupCount.ref.value;
      for (const t of groupTextures) {
        const { count } = c.textures[t].ref.value as DataArray<any>;
        if (count > span) span = count;
      }
      const o = i * 5;
      if (
        layout[o] !== c.instanceCount.ref.value ||
        layout[o + 1] !== c.uGroupCount.ref.value ||
        layout[o + 2] !== span ||
        layout[o + 3] !== c.uVertexCount.ref.value ||
        layout[o + 4] !== c.drawCount.ref.value
      ) {
        layout[o] = c.instanceCount.ref.value;
        layout[o + 1] = c.uGroupCount.ref.value;
        layout[o + 2] = span;
        layout[o + 3] = c.uVertexCount.ref.value;
        layout[o + 4] = c.drawCount.ref.value;
        changed = true;
      }
    }
    const mergeKey = getMergeKey(first);
    if (mergeKey !== layoutMergeKey) {
      layoutMergeKey = mergeKey;
      changed = true;
    }
    return changed;
  }

  /** the layer's array of the merged data array: the members' layers as parts (`SegmentedArray`), not copied */
  const layerArray = (spec: MergedTexSpec) => (own[spec.tex].ref.value as DataArray<any>).array as SegmentedArray;

  /** member i's part of the layer of `spec` read again from the member (its array may be another one) */
  function updateTexPart(spec: MergedTexSpec, i: number) {
    const part = layerArray(spec).parts[i] as SegmentedArray.Part;
    part.array = (members[i][spec.tex].ref.value as DataArray<any>).array;
    part.version += 1;
  }

  function mergeTextureFull(spec: MergedTexSpec) {
    const cell = own[spec.tex];
    const bases = texBases.get(spec.tex)!;
    const ctor = spec.float ? Float32Array : Uint8Array;

    if (spec.granularity(first) === 'none') {
      // no per-segment data used, the first member's layer
      const src = first[spec.tex].ref.value as DataArray<any>;
      const part = { array: src.array, offset: 0, count: src.count * spec.itemSize, bias: 0, version: 0 };
      bases.fill(0);
      ValueCell.update(cell, { array: SegmentedArray.create(ctor, [part]), count: src.count });
      return;
    }

    let total = 0;
    const parts: SegmentedArray.Part[] = [];
    for (let i = 0; i < n; ++i) {
      bases[i] = total;
      const texels = getTexelCount(i, spec);
      const { array } = members[i][spec.tex].ref.value as DataArray<any>;
      parts.push({ array, offset: total * spec.itemSize, count: texels * spec.itemSize, bias: 0, version: 0 });
      total += texels;
    }
    // whole words, as `createDataArray`
    const count = Math.max(1, total);
    const length = Math.ceil((count * spec.itemSize * ctor.BYTES_PER_ELEMENT) / 4) * (4 / ctor.BYTES_PER_ELEMENT);
    ValueCell.update(cell, { array: SegmentedArray.create(ctor, parts, length), count });
  }

  function unionSpheres(out: Sphere3D, key: 'boundingSphere' | 'invariantBoundingSphere') {
    boundaryHelper.reset();
    for (let i = 0; i < n; ++i) {
      boundaryHelper.includeSphere(members[i][key].ref.value as Sphere3D);
    }
    boundaryHelper.finishedIncludeStep();
    for (let i = 0; i < n; ++i) {
      boundaryHelper.radiusSphere(members[i][key].ref.value as Sphere3D);
    }
    return boundaryHelper.getSphere(out);
  }

  function updateBoundingSpheres() {
    ValueCell.update(own.boundingSphere, unionSpheres(own.boundingSphere.ref.value, 'boundingSphere'));
    const ibs = unionSpheres(own.invariantBoundingSphere.ref.value, 'invariantBoundingSphere');
    ValueCell.update(own.invariantBoundingSphere, ibs);
    ValueCell.update(own.uInvariantBoundingSphere, Vec4.fromSphere(own.uInvariantBoundingSphere.ref.value, ibs));
  }

  function updateMarkerAggregates() {
    let average = 0;
    let weight = 0;
    let status = members[0].markerStatus.ref.value;
    for (let i = 0; i < n; ++i) {
      const v = members[i];
      const w = Math.max(
        1,
        v.markerType.ref.value === 'instance'
          ? v.instanceCount.ref.value
          : v.instanceCount.ref.value * v.uGroupCount.ref.value,
      );
      average += v.markerAverage.ref.value * w;
      weight += w;
      if (v.markerStatus.ref.value !== status) status = -1;
    }
    ValueCell.updateIfChanged(own.markerAverage, weight > 0 ? average / weight : 0);
    ValueCell.updateIfChanged(own.markerStatus, status);
    ValueCell.updateIfChanged(own.uMarker, status);
  }

  function updateAggregates() {
    let transparencyAverage = 0;
    let transparencyMin = 1;
    let emissiveAverage = 0;
    let wiggleAverage = 0;
    let hasReflection = false;
    let weight = 0;
    for (let i = 0; i < n; ++i) {
      const v = members[i];
      const w = Math.max(1, v.instanceCount.ref.value * v.uGroupCount.ref.value);
      transparencyAverage += v.transparencyAverage.ref.value * w;
      emissiveAverage += v.emissiveAverage.ref.value * w;
      wiggleAverage += v.wiggleAverage.ref.value * w;
      transparencyMin = Math.min(transparencyMin, v.transparencyMin.ref.value);
      hasReflection = hasReflection || v.hasReflection.ref.value;
      weight += w;
    }
    ValueCell.updateIfChanged(own.transparencyAverage, weight > 0 ? transparencyAverage / weight : 0);
    ValueCell.updateIfChanged(own.transparencyMin, transparencyMin);
    ValueCell.updateIfChanged(own.emissiveAverage, weight > 0 ? emissiveAverage / weight : 0);
    ValueCell.updateIfChanged(own.wiggleAverage, weight > 0 ? wiggleAverage / weight : 0);
    ValueCell.updateIfChanged(own.hasReflection, hasReflection);
  }

  /** member i's vertex attribute as a part: its vertices at its vertex base */
  function vertexPart(attr: VertexAttribute, i: number): SegmentedArray.Part {
    return {
      array: members[i][attr.key].ref.value as Float32Array,
      offset: segments.vertexBases[i] * attr.itemSize,
      count: members[i].uVertexCount.ref.value * attr.itemSize,
      bias: 0,
      version: 0,
    };
  }

  /** member i's elements as a part: its indices at its element offset, made absolute by its vertex base */
  function elementsPart(i: number): SegmentedArray.Part {
    return {
      array: members[i].elements.ref.value as Uint32Array,
      offset: segments.elementOffsets[i] / 4,
      count: members[i].drawCount.ref.value,
      bias: segments.vertexBases[i],
      version: 0,
    };
  }

  /** member i's part of `key` read again from the member (its array may be another of the same size) */
  function updatePart(key: string, i: number) {
    const part = (own[key].ref.value as SegmentedArray).parts[i] as SegmentedArray.Part;
    part.array = members[i][key].ref.value;
    part.version += 1;
  }

  /** per-segment group counts and offsets, per instance in `aSegment` */
  const segmentGroupCounts = new Array(n).fill(0);
  const segmentGroupInstanceOffsets = new Array(n).fill(0);
  const segmentGroupBases = new Array(n).fill(0);

  /** member i's own (sphere-count-dependent) per-level scales, the second element of its segment data */
  function writeSegmentLod(i: number, array: Float32Array) {
    const member = members[i];
    // gridless members are drawn at full (undecimated) detail, see the
    // `cellSize <= 1` branch in renderable.ts's cullSegment - no per-level
    // scale compensation applies then, else spheres would render too large
    const gridless = member.instanceGrid.ref.value.cellSize <= 1;
    const lodLevels = !gridless ? (member.lodLevels?.ref.value as LodLevelsValue | undefined) : undefined;
    const s0 = lodLevels?.[0]?.[4] ?? 1;
    const s1 = lodLevels?.[1]?.[4] ?? 1;
    const s2 = lodLevels?.[2]?.[4] ?? 1;
    const s3 = lodLevels?.[3]?.[4] ?? 1;
    const o = (i * 2 + 1) * 4;
    array[o + 0] = s0;
    array[o + 1] = s1;
    array[o + 2] = s2;
    array[o + 3] = s3;
  }

  function computeSegments() {
    let instanceTotal = 0;
    let vertexTotal = 0;
    let drawTotal = 0;
    let maxGroupCount = 0;
    let groupInstanceBase = 0;
    let groupBase = 0;
    for (let i = 0; i < n; ++i) {
      const v = members[i];
      const instanceCount = v.instanceCount.ref.value;
      const groupCount = v.uGroupCount.ref.value;
      segments.instanceBases[i] = instanceTotal;
      segments.instanceCounts[i] = instanceCount;
      segments.vertexBases[i] = vertexTotal;
      segments.drawBases[i] = drawTotal;
      segments.elementOffsets[i] = drawTotal * 4;
      segmentGroupCounts[i] = groupCount;
      segmentGroupInstanceOffsets[i] = groupInstanceBase - instanceTotal * groupCount;
      segmentGroupBases[i] = groupBase;
      instanceTotal += instanceCount;
      vertexTotal += v.uVertexCount.ref.value;
      drawTotal += v.drawCount.ref.value;
      groupInstanceBase += instanceCount * groupCount;
      groupBase += getGroupSpan(i);
      if (groupCount > maxGroupCount) maxGroupCount = groupCount;
    }
    segments.count = n;
    return { instanceTotal, vertexTotal, drawTotal, maxGroupCount };
  }

  function rebuildInstances(instanceTotal: number) {
    const segmentInstanceArray = ensureArray(own.aSegment.ref.value, instanceTotal * 4, Float32Array);
    const segmentData = own.tSegment.ref.value as DataArray<Float32Array>;
    const segmentArray = segmentData.array;

    // instanced attributes in per-segment blocks
    const transformArray = ensureArray(own.aTransform.ref.value, instanceTotal * 16, Float32Array);
    const instanceArray = ensureArray(own.aInstance.ref.value, instanceTotal, Float32Array);
    for (let i = 0; i < n; ++i) {
      const v = members[i];
      const instanceCount = segments.instanceCounts[i];
      const instanceBase = segments.instanceBases[i];
      transformArray.set((v.aTransform.ref.value as Float32Array).subarray(0, instanceCount * 16), instanceBase * 16);
      const memberInstance = v.aInstance.ref.value as Float32Array;
      const ibs = v.invariantBoundingSphere.ref.value as Sphere3D;
      // the segment's data: invariant bounding sphere, level scales
      const o = i * 8;
      segmentArray[o + 0] = ibs.center[0];
      segmentArray[o + 1] = ibs.center[1];
      segmentArray[o + 2] = ibs.center[2];
      segmentArray[o + 3] = ibs.radius;
      writeSegmentLod(i, segmentArray);
      for (let k = 0; k < instanceCount; ++k) {
        instanceArray[instanceBase + k] = memberInstance[k] + instanceBase;
        const s = (instanceBase + k) * 4;
        segmentInstanceArray[s] = segmentGroupCounts[i];
        segmentInstanceArray[s + 1] = segmentGroupInstanceOffsets[i];
        segmentInstanceArray[s + 2] = segmentGroupBases[i];
        segmentInstanceArray[s + 3] = i;
      }
    }
    ValueCell.update(own.aTransform, transformArray);
    ValueCell.update(own.aInstance, instanceArray);
    // cellSize > 1 keeps the renderer on the cull path, each member's OWN instanceGrid is used instead
    ValueCell.update(own.instanceGrid, { ...createEmptyInstanceGrid(), cellSize: 2 });
    ValueCell.update(own.aSegment, segmentInstanceArray);
    ValueCell.update(own.tSegment, segmentData);
  }

  function rebuild() {
    const { instanceTotal, vertexTotal, drawTotal, maxGroupCount } = computeSegments();

    rebuildInstances(instanceTotal);

    for (const attr of vertexAttributes) {
      const parts = members.map((_, i) => vertexPart(attr, i));
      ValueCell.update(own[attr.key], SegmentedArray.create(Float32Array, parts));
    }
    if (hasElements) {
      ValueCell.update(
        own.elements,
        SegmentedArray.create(
          Uint32Array,
          members.map((_, i) => elementsPart(i)),
        ),
      );
    }

    for (const spec of texSpecs) {
      mergeTextureFull(spec);
    }

    ValueCell.updateIfChanged(own.drawCount, drawTotal);
    ValueCell.updateIfChanged(own.instanceCount, instanceTotal);
    ValueCell.updateIfChanged(own.uVertexCount, vertexTotal);
    ValueCell.updateIfChanged(own.uInstanceCount, instanceTotal);
    ValueCell.updateIfChanged(own.uGroupCount, maxGroupCount);

    updateBoundingSpheres();
    updateMarkerAggregates();
    updateAggregates();
  }

  function snapshotAll() {
    for (let i = 0; i < n; ++i) {
      const c = memberCells[i];
      const mv = memberVersions[i];
      mv.aTransform = c.aTransform.ref.version;
      mv.aInstance = c.aInstance.ref.version;
      mv.boundingSphere = c.boundingSphere.ref.version;
      mv.invariantBoundingSphere = c.invariantBoundingSphere.ref.version;
      mv.markerStatus = c.markerStatus.ref.version;
      mv.markerAverage = c.markerAverage.ref.version;
      mv.transparencyAverage = c.transparencyAverage.ref.version;
      mv.transparencyMin = c.transparencyMin.ref.version;
      mv.emissiveAverage = c.emissiveAverage.ref.version;
      mv.wiggleAverage = c.wiggleAverage.ref.version;
      mv.hasReflection = c.hasReflection.ref.version;
      mv.lodLevels = c.lodLevels ? c.lodLevels.ref.version : -1;
      for (let a = 0; a < c.attributes.length; ++a) mv.attributes[a] = c.attributes[a].ref.version;
      mv.elements = c.elements ? c.elements.ref.version : -1;
      for (let t = 0; t < c.textures.length; ++t) mv.textures[t] = c.textures[t].ref.version;
    }
  }

  const cullVersion = ValueCell.create(0);

  function sync() {
    if (updateLayout()) {
      rebuild();
      snapshotAll();
      ValueCell.update(cullVersion, cullVersion.ref.value + 1);
      return true;
    }

    let instancesChanged = false;
    let boundsChanged = false;
    let markersChanged = false;
    let aggregatesChanged = false;
    let elementsChanged = false;
    let lodChanged = false;
    const attributesChanged = new Set<string>();

    for (let i = 0; i < n; ++i) {
      const c = memberCells[i];
      const mv = memberVersions[i];
      const transform = c.aTransform.ref.version;
      const instance = c.aInstance.ref.version;
      const bounds = c.boundingSphere.ref.version;
      const invariantBounds = c.invariantBoundingSphere.ref.version;
      const markerAverage = c.markerAverage.ref.version;
      const transparencyAverage = c.transparencyAverage.ref.version;
      const transparencyMin = c.transparencyMin.ref.version;
      const emissiveAverage = c.emissiveAverage.ref.version;
      const wiggleAverage = c.wiggleAverage.ref.version;
      const hasReflection = c.hasReflection.ref.version;

      // the shared grid placement depends on transforms and invariant spheres
      if (mv.aTransform !== transform || mv.aInstance !== instance || mv.invariantBoundingSphere !== invariantBounds) {
        instancesChanged = true;
      }
      if (mv.boundingSphere !== bounds || mv.invariantBoundingSphere !== invariantBounds) {
        boundsChanged = true;
      }
      if (mv.markerStatus !== c.markerStatus.ref.version || mv.markerAverage !== markerAverage) {
        markersChanged = true;
      }
      if (
        mv.transparencyAverage !== transparencyAverage ||
        mv.transparencyMin !== transparencyMin ||
        mv.emissiveAverage !== emissiveAverage ||
        mv.wiggleAverage !== wiggleAverage ||
        mv.hasReflection !== hasReflection
      ) {
        aggregatesChanged = true;
      }
      for (let a = 0; a < c.attributes.length; ++a) {
        const version = c.attributes[a].ref.version;
        if (mv.attributes[a] !== version) {
          updatePart(vertexAttributes[a].key, i);
          attributesChanged.add(vertexAttributes[a].key);
          mv.attributes[a] = version;
        }
      }
      if (c.elements && mv.elements !== c.elements.ref.version) {
        updatePart('elements', i);
        elementsChanged = true;
        mv.elements = c.elements.ref.version;
      }

      const lodVersion = c.lodLevels ? c.lodLevels.ref.version : -1;
      if (mv.lodLevels !== lodVersion) {
        writeSegmentLod(i, (own.tSegment.ref.value as DataArray<Float32Array>).array);
        lodChanged = true;
      }

      // the versions seen, as `snapshotAll` (the textures' and `markerStatus` below: the texture loop reads them)
      mv.aTransform = transform;
      mv.aInstance = instance;
      mv.boundingSphere = bounds;
      mv.invariantBoundingSphere = invariantBounds;
      mv.markerAverage = markerAverage;
      mv.transparencyAverage = transparencyAverage;
      mv.transparencyMin = transparencyMin;
      mv.emissiveAverage = emissiveAverage;
      mv.wiggleAverage = wiggleAverage;
      mv.hasReflection = hasReflection;
      mv.lodLevels = lodVersion;
    }

    for (let t = 0; t < texSpecs.length; ++t) {
      const spec = texSpecs[t];
      let texChanged = false;
      if (spec.granularity(first) === 'none') {
        if (memberVersions[0].textures[t] !== memberCells[0].textures[t].ref.version) {
          // an array of another size is a new layout (`mergeTextureFull`)
          const src = first[spec.tex].ref.value as DataArray<any>;
          if (src.count === (own[spec.tex].ref.value as DataArray<any>).count) {
            updateTexPart(spec, 0);
            ValueCell.update(own[spec.tex], own[spec.tex].ref.value);
          } else {
            mergeTextureFull(spec);
          }
        }
        // only the first member's is read while the granularity stays (a change of it rebuilds)
        memberVersions[0].textures[t] = memberCells[0].textures[t].ref.version;
        continue;
      }
      const marker = spec.tex === 'tMarker';
      for (let i = 0; i < n; ++i) {
        const c = memberCells[i];
        const mv = memberVersions[i];
        const version = c.textures[t].ref.version;
        // uniform-status marking updates the array without bumping the texture cell
        if (mv.textures[t] !== version || (marker && mv.markerStatus !== c.markerStatus.ref.version)) {
          updateTexPart(spec, i);
          texChanged = true;
        }
        mv.textures[t] = version;
      }
      if (texChanged) {
        ValueCell.update(own[spec.tex], own[spec.tex].ref.value);
      }
    }

    if (instancesChanged) rebuildInstances(own.instanceCount.ref.value as number);
    for (const key of attributesChanged) ValueCell.update(own[key], own[key].ref.value);
    if (elementsChanged) ValueCell.update(own.elements, own.elements.ref.value);
    if (!instancesChanged && lodChanged) ValueCell.update(own.tSegment, own.tSegment.ref.value);
    if (boundsChanged) updateBoundingSpheres();
    if (markersChanged) updateMarkerAggregates();
    if (aggregatesChanged) updateAggregates();
    if (instancesChanged || boundsChanged || lodChanged) ValueCell.update(cullVersion, cullVersion.ref.value + 1);

    for (let i = 0; i < n; ++i) memberVersions[i].markerStatus = memberCells[i].markerStatus.ref.version;
    return false;
  }

  updateLayout();
  rebuild();
  snapshotAll();

  // the members' marker cells only: a mark changes members without an update of the merged object, their other changes
  // come with one; two cells a member, as there may be 100k members
  const markerTexture = texSpecs.findIndex((s) => s.tex === 'tMarker');
  const observed: ValueCell<unknown>[] = [];
  for (const c of memberCells) {
    observed.push(c.markerStatus);
    if (markerTexture >= 0) observed.push(c.textures[markerTexture]);
  }

  return { type, values, segments, sync, memberCells: observed, cullVersion };
}

//

const MergedRenderableSchemas = new Map<MergedType, RenderableSchema>();

/** The schema of merged renderables of `type` (one object per type, shared by the items) */
export function getMergedRenderableSchema(type: MergedType) {
  let schema = MergedRenderableSchemas.get(type);
  if (!schema) {
    schema = {
      ...GlobalUniformSchema,
      ...GlobalTextureSchema,
      ...GlobalDefineSchema,
      ...InternalSchema,
      ...MergedTypeDescriptors[type].schema,
      ...SegmentSchema,
    };
    MergedRenderableSchemas.set(type, schema);
  }
  return schema;
}

/** The shader code of merged renderables of `type` */
export function getMergedShaderCode(type: MergedType) {
  return MergedTypeDescriptors[type].shaderCode;
}

export function MergedRenderable(
  ctx: GraphicsContext,
  id: number,
  merged: Merged,
  members: readonly MergeableValues[],
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<MergeableValues & SegmentValues> {
  const descriptor = MergedTypeDescriptors[merged.type];
  const schema = getMergedRenderableSchema(merged.type);
  const renderValues: MergeableValues & SegmentValues & InternalValues & GlobalDefineValues = {
    ...merged.values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  /**
   * Whether the culler's records are indexed: with elements, or items drawn by the shared indices of their pattern
   * (their draws within its limit, else unindexed as the render item draws them)
   */
  const items = descriptor.shaderCode.itemIndices;
  const indexed =
    descriptor.hasElements ||
    (isItemIndexing(items) &&
      merged.values.drawCount.ref.value / items.pattern.length <= getMaxItemCount(ctx.resources.limits, items));
  // culled on the GPU where asked for and available (`merged-culler.ts`)
  const gpuCulling = !!(merged.values as { gpuCulling?: ValueCell<boolean> }).gpuCulling?.ref.value;
  const culler =
    gpuCulling && canCullOnGpu(ctx)
      ? createMergedCuller(ctx.compute!, ctx.resources, merged, members, indexed)
      : undefined;
  // drawn with the culler's compacted instanced attributes
  if (culler) Object.assign(renderValues, culler.attributes);
  const renderItem = ctx.createGraphicsRenderItem(
    descriptor.drawMode,
    descriptor.shaderCode,
    schema,
    renderValues,
    materialId,
    transparency,
  );

  const mdb = createSegmentedMdbList();
  const clipInfo = createClipInfoCache(renderValues);
  const segmentSpheres: Sphere3D[] = [];
  let mode: 'none' | 'full' | 'cull' | 'gpu' = 'none';
  let lastCullFrame: Frame | undefined;

  /** the frame culled on the GPU (its first view) and the planes it was culled with */
  let gpuFrame: Frame | undefined;
  const gpuPlanes = new Float64Array(28);
  let gpuLists: MultiDrawBaseData[] = [];
  const getPlane = (i: number, cameraPlane: Plane3D, frustum: Frustum3D) => (i < 6 ? frustum[i as 0] : cameraPlane);
  function isGpuView(cameraPlane: Plane3D, frustum: Frustum3D) {
    for (let i = 0; i < 7; ++i) {
      const { normal, constant } = getPlane(i, cameraPlane, frustum);
      const o = i * 4;
      if (
        gpuPlanes[o] !== normal[0] ||
        gpuPlanes[o + 1] !== normal[1] ||
        gpuPlanes[o + 2] !== normal[2] ||
        gpuPlanes[o + 3] !== constant
      ) {
        return false;
      }
    }
    return true;
  }
  function setGpuView(cameraPlane: Plane3D, frustum: Frustum3D) {
    for (let i = 0; i < 7; ++i) {
      const { normal, constant } = getPlane(i, cameraPlane, frustum);
      gpuPlanes.set([normal[0], normal[1], normal[2], constant], i * 4);
    }
  }

  /** the CPU culling's lists, shifted to the identity region of the compacted attributes where culled on the GPU */
  const shiftedLists: MultiDrawBaseData[] = [];
  let shiftedVersion = -1;
  function getCpuLists() {
    if (!culler) return mdb.list;
    if (shiftedVersion === drawVersion && shiftedLists.length === mdb.list.length) return shiftedLists;
    shiftedVersion = drawVersion;
    return shiftLists(mdb.list, shiftedLists);
  }
  /** `lists` shifted to the identity region of the compacted attributes (`culler.identityBase`) into `out` */
  function shiftLists(lists: MultiDrawBaseData[], out: MultiDrawBaseData[]) {
    const base = culler!.identityBase;
    out.length = lists.length;
    for (let j = 0, jl = lists.length; j < jl; ++j) {
      const l = lists[j];
      const previous = out[j];
      const baseInstances =
        previous && previous.baseInstances.length >= l.count
          ? previous.baseInstances
          : new Uint32Array(l.baseInstances.length);
      for (let i = 0; i < l.count; ++i) baseInstances[i] = l.baseInstances[i] + base;
      out[j] = { ...l, baseInstances };
    }
    return out;
  }
  /** whether the sphere raster draws the object's raster regions in the pass being drawn (`Renderable.rasterDrawn`) */
  let rasterDrawn = false;
  // the impostors draw the raster regions' instances too where the raster does not draw them
  const getLists = () => (mode === 'gpu' ? culler!.getLists(mdb.list, !rasterDrawn) : getCpuLists());

  // the members' cells culling reads every frame, resolved once: values objects have too many keys for fast named reads
  const memberCull: CullValues[] = members.map((v) => ({
    drawCount: v.drawCount,
    instanceCount: v.instanceCount,
    uLod: v.uLod,
    instanceGrid: v.instanceGrid,
    boundingSphere: v.boundingSphere,
    lodLevels: v.lodLevels,
    dClipObjects: v.dClipObjects,
    dClipping: v.dClipping,
    clipObjects: v.clipObjects,
    dClipVariant: v.dClipVariant,
    invariantBoundingSphere: v.invariantBoundingSphere,
  }));

  /**
   * what the last CPU cull read: a still frame keeps its lists (`createCullView`); the members' cells by the version
   * `sync` advances where they changed, not each (with 100k members, recording them was most of a moving frame)
   */
  const cullView = createCullView();
  const cullCells = [merged.cullVersion];

  const segment: CullSegment = { first: 0, offset: 0, instanceBase: 0 };
  function setSegment(i: number) {
    const { segments } = merged;
    segment.first = descriptor.hasElements ? 0 : segments.drawBases[i];
    segment.offset = segments.elementOffsets[i];
    segment.instanceBase = segments.instanceBases[i];
    return segment;
  }

  function getCapacity() {
    let capacity = 0;
    for (let i = 0, il = members.length; i < il; ++i) {
      capacity += Math.max(1, memberCull[i].instanceGrid.ref.value.cellCount) + 1;
    }
    return capacity;
  }

  function getMemberLod(i: number, hasLodLevels: boolean): LodLevelsValue | undefined {
    // per-member lod counts, the level ranges are equal across members
    return hasLodLevels ? (memberCull[i].lodLevels!.ref.value as LodLevelsValue) : undefined;
  }

  /** the level view of the last lists, for the caps' draws of segments (`renderSegments`) */
  let segmentLod = LodView();

  function buildFullList(lod: LodView) {
    segmentLod = lod;
    const hasLodLevels = !!mdb.prepare(memberCull[0].lodLevels, getCapacity(), true, lod);
    for (let i = 0, il = members.length; i < il; ++i) {
      mdb.fullSegment(memberCull[i], setSegment(i), getMemberLod(i, hasLodLevels));
    }
    mode = 'full';
  }

  function cull(
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    isOccluded: ((s: Sphere3D) => boolean) | null,
    stats: GraphicsStats,
    frame: Frame,
    lod: LodView,
  ) {
    segmentLod = lod;
    // skip recomputation if this renderable was already culled for the current frame
    if (frame === lastCullFrame) {
      stats.cacheHits.cull++;
      return;
    }

    if (culler && frame === gpuFrame && isGpuView(cameraPlane, frustum)) {
      // culled on the GPU before the passes (`prepareCull`): the level lists for their uniforms only
      mdb.prepare(memberCull[0].lodLevels, 0, true, lod);
      gpuLists = culler.getLists(mdb.list);
      mode = 'gpu';
      lastCullFrame = frame;
      drawVersion = getDrawVersion();
      return;
    }

    // the same view and members' data as the last CPU cull (a still frame): its lists stand
    if (mode === 'cull' && cullView.isSame(lod, isOccluded, cullCells)) {
      lastCullFrame = frame;
      stats.cacheHits.cull++;
      return;
    }
    cullView.set(lod, cullCells);
    const lodLevels = mdb.prepare(memberCull[0].lodLevels, getCapacity(), true, lod);
    mode = 'cull';

    for (let i = 0, il = members.length; i < il; ++i) {
      mdb.cullSegment(
        memberCull[i],
        setSegment(i),
        getMemberLod(i, !!lodLevels),
        cameraPlane,
        frustum,
        isOccluded,
        stats,
        true,
      );
    }

    lastCullFrame = frame;
    drawVersion = getDrawVersion();
  }

  function cullSimple(d: number, radius: number, scale: number, lod: LodView) {
    segmentLod = lod;
    cullView.invalidate();
    const lodLevels = mdb.prepare(memberCull[0].lodLevels, getCapacity(), true, lod);
    if (!lodLevels) {
      buildFullList(lod);
    } else {
      mode = 'cull';
      for (let i = 0, il = members.length; i < il; ++i) {
        mdb.cullSimpleSegment(memberCull[i], setSegment(i), getMemberLod(i, true)!, d, radius, scale);
      }
    }
    drawVersion = getDrawVersion();
  }

  buildFullList(LodView());

  /**
   * The draws of the listed segments (ascending), each whole at its levels, as lists of their own (the main lists stay):
   * the solid-interior caps of the segments a cap cuts (`Renderable.renderSegments`), a few a frame
   */
  const segmentMdb = createSegmentedMdbList();
  const shiftedSegmentLists: MultiDrawBaseData[] = [];
  function renderSegments(variant: GraphicsRenderVariant, globals: GlobalBindings, segments: readonly number[]) {
    const hasLodLevels = !!segmentMdb.prepare(memberCull[0].lodLevels, segments.length, true, segmentLod);
    for (let k = 0, kl = segments.length; k < kl; ++k) {
      const i = segments[k];
      segmentMdb.fullSegment(memberCull[i], setSegment(i), getMemberLod(i, hasLodLevels));
    }
    renderItem.render(variant, globals, culler ? shiftLists(segmentMdb.list, shiftedSegmentLists) : segmentMdb.list);
  }

  let programVersion = 0;
  const getDrawVersion = createDrawVersion((out) => appendMdbDrawKey(mode === 'gpu' ? gpuLists : mdb.list, out));
  let drawVersion = getDrawVersion();
  return {
    id,
    materialId,
    renderItem,
    values: renderValues,
    state,

    cull,
    uncull: (lod: LodView) => {
      lastCullFrame = undefined;
      // levels by size on screen convert their ranges with the view (`prepare`)
      const pixel = ((members[0].lodLevels?.ref.value as LodLevelsValue | undefined)?.pixelSize ?? 0) > 0;
      if (mode === 'full' && !pixel) return;
      buildFullList(lod);
      drawVersion = getDrawVersion();
    },
    cullSimple,
    prepareCull: culler
      ? (
          cameraPlane: Plane3D,
          frustum: Frustum3D,
          frame: Frame,
          lod: LodView,
          pyramid?: OcclusionPyramid,
          raster?: MergedRasterParams,
        ) => {
          // once per frame, for its first view; other views of the frame draw the CPU culling's lists
          if (frame === gpuFrame) return true;
          // the sphere raster's split per instance for objects it can draw
          const split = raster && 'dSphereResolve' in renderValues && isSphereRasterEligible(renderValues as never);
          gpuFrame = culler.dispatch(cameraPlane, frustum, lod, pyramid, split ? raster : undefined)
            ? frame
            : undefined;
          if (gpuFrame) setGpuView(cameraPlane, frustum);
          return !!gpuFrame;
        }
      : undefined,
    needsCullSecond: culler ? (frame: Frame) => frame === gpuFrame && culler.secondPending : undefined,
    prepareCullSecond: culler
      ? (frame: Frame, pyramid: OcclusionPyramid | undefined) => frame === gpuFrame && culler.dispatchSecond(pyramid)
      : undefined,
    // where the view draws the GPU culling's lists, as `cull` decides: a frame culled already keeps its lists (views of
    // a reused frame, as the samples of temporal multi-sampling), else the view culled on the GPU
    getRasterInput: culler
      ? (frame: Frame, cameraPlane: Plane3D, frustum: Frustum3D) => {
          const gpu = frame === lastCullFrame ? mode === 'gpu' : frame === gpuFrame && isGpuView(cameraPlane, frustum);
          return gpu ? culler.getRasterInput() : undefined;
        }
      : undefined,
    render: (variant: GraphicsRenderVariant, globals: GlobalBindings) => {
      renderItem.render(variant, globals, getLists());
    },
    renderSecond: culler
      ? (variant: GraphicsRenderVariant, globals: GlobalBindings) => {
          if (mode !== 'gpu' || !culler.secondDispatched) return;
          renderItem.render(variant, globals, culler.getSecondLists(mdb.list, !rasterDrawn));
        }
      : undefined,
    renderSegments,
    getClipInfo: clipInfo.get,
    memberCells: merged.memberCells,
    getSegmentSpheres: () => {
      // the members' cells, read again: a member's sphere may be replaced
      for (let i = 0, il = members.length; i < il; ++i) segmentSpheres[i] = members[i].boundingSphere.ref.value;
      return segmentSpheres;
    },
    getByteCount: () => renderItem.getByteCount(),
    getProgram: (variant: GraphicsRenderVariant) => renderItem.getProgram(variant),
    setTransparency: (value: Transparency) => {
      renderItem.setTransparency(value);
      programVersion += 1;
    },
    get programVersion() {
      return programVersion;
    },
    get rasterDrawn() {
      return rasterDrawn;
    },
    set rasterDrawn(value: boolean) {
      rasterDrawn = value;
    },
    get recordsCompact() {
      return mode === 'gpu' ? culler!.recordsCompact : undefined;
    },
    get recordsVersion() {
      return culler?.recordsVersion ?? 0;
    },
    get drawVersion() {
      // the lists without the raster regions are other draws (recordings of the pass replay while it is the same)
      return drawVersion * 2 + (rasterDrawn ? 1 : 0);
    },
    update: () => {
      merged.sync();
      culler?.sync();
      gpuFrame = undefined;
      mode = 'none';
      renderItem.update();
      lastCullFrame = undefined;
    },
    dispose: () => {
      renderItem.destroy();
      culler?.destroy();
    },
  };
}
