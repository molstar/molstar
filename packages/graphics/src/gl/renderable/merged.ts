/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createInstanceStorage, usesInstanceStorage, withInstanceStorage } from '../instance-storage.js';
import { getCulledSegmentData } from '../culled-data.js';
import type { PagedCullChanges, PagedCullSource } from '../compute/paged-cull.js';
import { createPagedMergedCuller } from './merged-paged-culler.js';
import { fillSerial, getSharedSerial, isSharedSerial } from '@molstar/core/util/array';
import { isDebugMode } from '@molstar/core/util/debug';
import type { OcclusionPyramid } from '../compute/cell-cull.js';
import { createMergedCuller, type MergedRasterParams } from './merged-culler.js';
import { isSphereRasterEligible } from '../compute/sphere-raster-pass.js';
import {
  applyPendingGrid,
  buildPendingGrid,
  getObjectTransforms,
  hasInstanceCells,
  type InstanceGridValues,
  type TransformPart,
  type TransformValues,
  writeAffineRowsRange,
} from '../transform-parts.js';
import { OctNormalItemSize } from '../oct-normals.js';
import { ValueCell } from '@molstar/core/util';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import type { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { createEmptyInstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { createDataArray, type DataArray, getDrawnBoundingSphere } from './util.js';
import { SegmentedArray, type SegmentedArrayType } from '../segmented-array.js';
import {
  GlobalUniformSchema,
  GlobalTextureSchema,
  GlobalDefineSchema,
  InternalSchema,
  type InternalValues,
  SegmentSchema,
  SegmentTexels,
  type SegmentValues,
  type GlobalDefines,
  type GlobalDefineValues,
  type RenderableSchema,
  type RenderableValues,
  type BaseValues,
  DefineSpec,
  type Values,
} from './schema.js';
import { StorageBufferSpec } from '../compute.js';
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
  createDrawVersionCounter,
  createCullView,
  type DrawKey,
  createSegmentedMdbList,
  type CullSegment,
  type LodLevelsValue,
  LodView,
  type CullValues,
  type CullBatch,
  type Frame,
  createClipInfoCache,
  createInstancePositions,
  createMarkedKey,
  getMarkedDrawVersion,
  isMarkedFew,
} from '../renderable.js';

export type MergeableValues = RenderableValues & BaseValues;

/** The cells of a member a merged object reads (`MergeRecord`), named; with them its geometry and layers by key */
export const MergeRecordKeys = [
  'instanceCount',
  'uGroupCount',
  'uVertexCount',
  'drawCount',
  'aTransform',
  'aInstance',
  'boundingSphere',
  'invariantBoundingSphere',
  'instanceGrid',
  'instanceCellSize',
  'markerStatus',
  'markerAverage',
  'markerType',
  'transparencyAverage',
  'transparencyMin',
  'emissiveAverage',
  'wiggleAverage',
  'hasReflection',
] as const;

/** The cells of a member a merged object reads where its values have them (`MergeRecord`) */
export const MergeRecordOptionalKeys = [
  'sizeType',
  'tSize',
  'tSiteColor',
  'instancePoses',
  'lodLevels',
  'padding',
  'uSiteOrigin',
  'uSiteScale',
  'markerInstanceCounts',
  'markedInstanceCount',
] as const;

/**
 * What a merged object reads of a member (`createMergedValues`, its renderable and cullers): its counts, instances,
 * bounds, grid, levels, markers' and layers' aggregates, and by key its geometry (the type's vertex attributes,
 * `elements`) and layers (the merged textures). A member's values are its record (a record view,
 * `createMergedRenderObject`); what all members share (defines, uniforms, props, the clip) is the object's shared values, not a member's
 */
export type MergeRecord = Pick<BaseValues, (typeof MergeRecordKeys)[number]> & {
  readonly [key: string]: ValueCell<any> | undefined;
};

export const MergedTypes = ['spheres', 'mesh', 'cylinders', 'lines', 'points'] as const;

/**
 * Vec4s per instance of merged instance data read by index (the WebGL2 `dInstanceData` texture): the transform's rows,
 * the instance and its segment's index (the integers as their bits, `floatBitsToUint`, `floatBitsToInt`: exact past
 * 2^24). The segments' groups (`aSegment.xyz`, the same for each instance of a segment) in a table
 * after the instances, a slot each (its first vec4): 64 bytes an instance (80 with the groups per instance)
 */
export const InstanceDataVec4s = 4;

/** The segments of merged `values`' instances: the largest segment index + 1 */
export function getInstanceDataSegmentCount(values: RenderableValues) {
  // an object not merged (an expanded one, `dInstanceExpansion`): none
  if (!values.aSegment) return 0;
  const count: number = values.instanceCount.ref.value;
  const segment = values.aSegment.ref.value as Int32Array;
  let max = -1;
  for (let i = 0; i < count; ++i) if (segment[i * 4 + 3] > max) max = segment[i * 4 + 3];
  return max + 1;
}

/**
 * Packs the instance data of the instances `first` to `first + count` of merged `values` into `out` (slots of
 * `InstanceDataVec4s` vec4s): per instance its transform's rows and (instance, segment index); per segment of these
 * instances, at slot `segmentSlot` + its index, its groups (`aSegment.xyz`)
 */
export function packInstanceData(
  values: RenderableValues,
  out: Float32Array,
  first: number,
  count: number,
  segmentSlot: number,
) {
  const instance = values.aInstance.ref.value as Uint32Array;
  // an expanded object's: per instance its template particle (`expansionLocals`), the third word
  const locals = values.expansionLocals?.ref.value as Uint32Array | undefined;
  // an object not merged (an expanded one, `dInstanceExpansion`): its instances' alone
  const segment = values.aSegment?.ref.value as Int32Array | undefined;
  const bits = new Uint32Array(out.buffer, out.byteOffset, out.length);
  const signed = new Int32Array(out.buffer, out.byteOffset, out.length);
  const stride = InstanceDataVec4s * 4;
  // the transform's rows (`affine-rows.ts`, by the encoding), the instance and its segment's index
  writeAffineRowsRange(values as unknown as TransformValues, first, count, out, 0, stride);
  for (let k = 0; k < count; ++k) {
    const i = first + k;
    const o = k * stride;
    bits[o + 12] = instance[i];
    if (locals) bits[o + 14] = locals[i];
    else out[o + 14] = 0;
    out[o + 15] = 0;
    if (!segment) {
      signed[o + 13] = 0;
      continue;
    }
    const s = segment[i * 4 + 3];
    signed[o + 13] = s;
    // its segment's groups, the same from each of its instances
    const t = (segmentSlot + s) * stride;
    signed[t] = segment[i * 4];
    signed[t + 1] = segment[i * 4 + 1];
    signed[t + 2] = segment[i * 4 + 2];
    signed[t + 3] = 0;
  }
}
export type MergedType = (typeof MergedTypes)[number];

export function isMergeableType(type: string): type is MergedType {
  return MergedTypes.includes(type as MergedType);
}

//

/** `table`: a member's own entries packed after the members' before (spheres' `siteType` sizes, a table by type) */
type Granularity = 'instance' | 'group' | 'groupInstance' | 'vertex' | 'item' | 'table' | 'none';

type MergedTexSpec = {
  tex: string;
  itemSize: number;
  float?: boolean;
  /** int16 components (`rgba16i`: spheres' quantized centers) */
  int16?: boolean;
  granularity: (v: MergeableValues) => Granularity;
  /** texel count multiplier, e.g. for dual colors */
  factor?: (v: MergeableValues) => number;
  /**
   * whether its table is read by the segment's index (a member's texel at its index: `ColorTypeIndex` 'object'), not
   * by a base in the segment's data: merged by member, never by allocation (`isAllocTexture`)
   */
  bySegment?: (v: MergeableValues) => boolean;
};

/** The array type of the layer of `spec` */
function getTexArrayCtor(spec: MergedTexSpec) {
  return spec.float ? Float32Array : spec.int16 ? Int16Array : Uint8Array;
}

/** The granularity of spheres' second color (`SiteColorSchema`) */
function getSiteColorGranularity(v: MergeableValues): Granularity {
  const type = (v as { dSiteColor?: ValueCell<string> }).dSiteColor?.ref.value;
  return type === 'siteType' ? 'table' : type === 'group' ? 'group' : 'none';
}

/** Whether spheres' values draw quantized centers (`dSites`): their `tSite` merged, else `tPositionGroup` */
/** a page store's `aInstance`: none (its ids its culled poses') */
const EmptyInstances = new Uint32Array(0);

const hasSites = (v: MergeableValues) => !!(v as { dSites?: ValueCell<boolean> }).dSites?.ref.value;

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
    case 'siteType':
    // a texel a member: packed in order, member i's at i (its segment, `ColorTypeIndex` 'object')
    case 'object':
      return 'table';
    default:
      return 'none';
  }
}

/** The layers a source may hold while off (`MergedInstanceSource.layers`): on later by their define alone */
const SourceOptionalLayers = new Set([
  'tOverpaint',
  'tTransparency',
  'tClipping',
  'tEmissive',
  'tSubstance',
  'tWiggle',
]);

const BaseTexSpecs: MergedTexSpec[] = [
  {
    tex: 'tColor',

    itemSize: 4,
    granularity: (v) => pickGranularity(v.colorType.ref.value),
    factor: (v) => (v.dDualColor?.ref.value ? 2 : 1),
    bySegment: (v) => v.colorType.ref.value === 'object',
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

type VertexAttribute = {
  key: string;
  itemSize: number;
  ctor: Float32ArrayConstructor | Uint32ArrayConstructor | Int32ArrayConstructor;
};

function getVertexAttributes(schema: RenderableSchema): VertexAttribute[] {
  const ret: VertexAttribute[] = [];
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    // `aVertex` (the emulated vertex id of WebGL1) is not member data: the render item fills it with serial ids
    if (k === 'aVertex') continue;
    if (spec.type === 'attribute' && spec.divisor === 0) {
      // the members' arrays as they are on the CPU: encoded normals three floats a vertex (`oct-normals.ts`)
      // the arrays by the kind: integers as integers (`aGroup`)
      const ctor = spec.kind === 'uint32' ? Uint32Array : spec.kind === 'int32' ? Int32Array : Float32Array;
      ret.push({ key: k, itemSize: spec.kind === 'oct16' ? OctNormalItemSize : spec.itemSize, ctor });
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
  itemTexels: (v: MergeRecord) => number;
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
    // the instance copy binds 6 storage buffers, within the core limits' 8
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
    // one position texel per sphere, drawn as 6 indices; quantized centers (`dSites`) their sites', the other the first
    // member's placeholder (`none`)
    itemTexels: (v) => v.drawCount.ref.value / 6,
    extraTextures: [
      { tex: 'tPositionGroup', itemSize: 4, float: true, granularity: (v) => (hasSites(v) ? 'none' : 'item') },
      { tex: 'tSite', itemSize: 4, int16: true, granularity: (v) => (hasSites(v) ? 'item' : 'none') },
      // the second color: a table by site type (a member's packed after the others', its first entry in its segment
      // data, `siteColorTableBases`) or a texel a group
      { tex: 'tSiteColor', itemSize: 4, granularity: (v) => getSiteColorGranularity(v) },
    ],
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

/** The keys of the records of a merged object of `type` with the shared values `shared` (`MergeRecord`) */
export function getMergeRecordKeys(type: MergedType, shared: MergeableValues): string[] {
  const { vertexAttributes, hasElements, extraTextures } = MergedTypeDescriptors[type];
  const keys: string[] = [...MergeRecordKeys];
  for (const k of MergeRecordOptionalKeys) if (shared[k] !== undefined) keys.push(k);
  for (const { key } of vertexAttributes) keys.push(key);
  if (hasElements) keys.push('elements');
  for (const { tex } of [...extraTextures, ...BaseTexSpecs]) {
    if (shared[tex] !== undefined && !keys.includes(tex)) keys.push(tex);
  }
  return keys;
}

/** The record of `values` (of `type`, `getMergeRecordKeys`): its cells of the record's keys alone, not copied */
export function pickMergeRecord(type: MergedType, values: MergeableValues, keys = getMergeRecordKeys(type, values)) {
  const record: { [k: string]: ValueCell<any> | undefined } = {};
  for (const k of keys) record[k] = values[k];
  return record as MergeRecord;
}

//

function isDefineKey(k: string) {
  return k.length > 1 && k[0] === 'd' && k[1] >= 'A' && k[1] <= 'Z';
}

/**
 * Signature of everything that must be equal across members sharing a merged
 * render item: define values, granularity resolution, and lod configuration.
 *
 * Note: `lodLevels[j]`'s `scale` (index 3) is intentionally excluded - it is
 * data-dependent per member (see `Spheres.getLodLevelsValue`, based on that
 * member's own sphere count), not part of the user-configured lod preset
 * (only min/max are). Per-member scales are instead carried in the
 * segment data (`tSegment`, see `SegmentSchema`) and selected in the shader
 * via `uLodLevel`, so they never need to match here.
 */
export function getMergeKey(values: MergeableValues): string {
  // checked per member on every update of a merged representation: computed again only where a cell it reads changed
  let cached = mergeKeys.get(values);
  if (!cached) {
    const keys = intern(internedDefineKeys, Object.keys(values).filter(isDefineKey).sort());
    // the color type: types of one define read their layer by another granularity (`object` and `instance`)
    const cells = [...keys.map((k) => values[k]), values.uInstanceGranularity, values.colorType, values.uLod];
    if (values.lodLevels) cells.push(values.lodLevels);
    // the index of `uLod` (compared by its components), not looked up per cell: values have too many keys for fast reads
    const lod = keys.length + 2;
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
  parts.push(`ct:${values.colorType.ref.value}`);
  parts.push(`uLod:${(values.uLod.ref.value as Vec4).join(',')}`);
  const lodLevels = values.lodLevels?.ref.value as LodLevelsValue | undefined;
  parts.push(
    `lod${lodLevels?.pixelSize ? `px${lodLevels.pixelSize}` : ''}:` +
      (lodLevels ? lodLevels.map((l) => `${l[0]},${l[1]}`).join(';') : ''),
  );
  return parts.join('|');
}

const SupportedVertexTypes = ['uniform', 'object', 'instance', 'group', 'groupInstance', 'vertex', 'siteType'];
const SupportedIgTypes = ['instance', 'groupInstance'];

/**
 * Whether `members` merge into one object of `type`; `layers` false: their layers' types not checked (values whose
 * layers the merged object sets again its own way, e.g. per-group visuals' smoothed layers, `UnitsRepresentation`)
 */
export function canMergeValues(type: string, members: readonly MergeableValues[], layers = true): boolean {
  if (!isMergeableType(type)) return false;
  // a single member merges as well (culled on the GPU): how many is the caller's (`canMergeRenderObjects`)
  if (members.length < 1) return false;

  // a palette the theme's (its texture, domain and default shared values), looked up after `tColor`, read by segment
  const first = members[0];
  if (!SupportedVertexTypes.includes(first.colorType.ref.value)) return false;
  if (first.sizeType && !SupportedVertexTypes.includes(first.sizeType.ref.value)) return false;
  if (layers && !canMergeLayers(first)) return false;

  const key = getMergeKey(first);
  for (let i = 1, il = members.length; i < il; ++i) {
    if (getMergeKey(members[i]) !== key) return false;
  }
  return true;
}

function canMergeLayers(first: MergeableValues) {
  if (first.dOverpaint.ref.value && !SupportedIgTypes.includes(first.overpaintType.ref.value)) return false;
  if (first.dTransparency.ref.value && !SupportedIgTypes.includes(first.transparencyType.ref.value)) return false;
  if (first.dEmissive.ref.value && !SupportedIgTypes.includes(first.emissiveType.ref.value)) return false;
  if (first.dSubstance.ref.value && !SupportedIgTypes.includes(first.substanceType.ref.value)) return false;
  if (first.dClipping.ref.value && !SupportedIgTypes.includes(first.clippingType.ref.value)) return false;
  if (first.dWiggle.ref.value && !SupportedIgTypes.includes(first.wiggleType.ref.value)) return false;
  return true;
}

export type MergedSegments = {
  count: number;
  /** first instance (into the merged instanced attributes) per segment */
  instanceBases: number[];
  /**
   * per segment its first slot in each culling region (`k * instanceCount + regionBases[s]`, `MergedCuller`): where its
   * culled instances go, apart from where its data is (a store's blocks: a member's instances in many of them, its
   * region a range of its live count); `instanceBases` without (a member's instances its own range)
   */
  regionBases?: number[];
  /** instance count per segment */
  instanceCounts: number[];
  /** first vertex (into the merged geometry) per segment */
  vertexBases: number[];
  /** first drawn vertex per segment, for types without elements (e.g. 6 per point quad with one data vertex per point) */
  drawBases: number[];
  /** byte offset (into the merged element buffer) per segment */
  elementOffsets: number[];
  /**
   * per segment its groups as every instance of it reads them (`aSegment.xyz`): the group count, the group-instance
   * offset (its first group-instance less its first instance times the group count) and the group offset
   */
  groupCounts: number[];
  groupInstanceOffsets: number[];
  groupBases: number[];
};

export type Merged = {
  readonly type: MergedType;
  readonly values: MergeableValues & SegmentValues;
  /** the members' records, in the segments' order (`MergeRecord`) */
  readonly records: readonly MergeRecord[];
  readonly segments: MergedSegments;
  /**
   * sync merged values with member values, returns true if the layout changed; `observedOnly`: only observed cells
   * changed (`memberCells`, a mark: `Renderable.update`), the members read those of which changed
   */
  sync: (observedOnly?: boolean) => boolean;
  /** the members' marker cells `sync` reads, observed with the merged values' own (`Renderable.memberCells`) */
  readonly memberCells: readonly ValueCell<unknown>[];
  /**
   * the members with a marker average other than 0 as the last sync saw them, ascending; undefined for a source's
   * markers (`MergedInstanceSource.markers`)
   */
  readonly getMarkedMembers?: () => readonly number[];
  /**
   * whether the records' owner wrote none of their cells since the last sync but observed ones (its `recordsVersion`
   * the same: e.g. a clip moved, the shared values alone); undefined without one
   */
  readonly recordsUnchanged?: () => boolean;
  /**
   * advanced by `sync` where what culling reads of the members changed (layout, transforms and grids, spheres, levels):
   * a still frame's lists stand while it stays, without comparing the members' cells (six a member)
   */
  readonly cullVersion: ValueCell<number>;
  /** advanced by `sync` where a member's markers changed (or the layout): what the marked lists are built from */
  readonly markerVersion: ValueCell<number>;
  /** the instances kept elsewhere (`MergedInstanceSource`): the culler's data theirs */
  readonly source?: MergedInstanceSource;
  /**
   * per member a count `sync` advances where it read the member's cells as changed (all at a whole layout): what the
   * culler reads again of a member (with a source: those its `takeMembers` gave), the rest kept
   */
  readonly memberChanges: Uint32Array;
  /**
   * With a source: member `i`'s geometry is another (its cells filled from another target's): allocated again at the
   * arrays' ends at the next `sync`, the others' kept and not written again (not a new layout)
   */
  replaceGeometry(i: number): void;
};

/**
 * With a source, the room the geometry's arrays are laid out with (a share more, at least so many elements more: a small
 * object's too): members replaced fill it before a new layout
 */
const SourceHeadroom = 1.5;
const SourceMinRoom = 1 << 14;
const sourceRoom = (used: number) => Math.max(Math.ceil(used * SourceHeadroom), used + SourceMinRoom);

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
  /** spheres' quantized centers' origin and scale (`dSites`), in the segment data */
  readonly uSiteOrigin: Cell | undefined;
  readonly uSiteScale: Cell | undefined;
  readonly elements: Cell | undefined;
}

/** The cells of `v`: they stay for the values' lifetime (updates replace their values, not the cells) */
function getMemberCells(v: MergeRecord, hasElements: boolean): MemberCells {
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
    uSiteOrigin: (v as { uSiteOrigin?: Cell }).uSiteOrigin,
    uSiteScale: (v as { uSiteScale?: Cell }).uSiteScale,
    elements: hasElements ? v.elements : undefined,
  };
}

/** The versions of a member's sites' origin and scale (`MemberCells`), one number: -1 without */
function getSitesVersion(c: MemberCells) {
  return c.uSiteOrigin && c.uSiteScale ? c.uSiteOrigin.ref.version * 0x100000 + c.uSiteScale.ref.version : -1;
}

/** The versions of a member's cells `sync` saw last (-1: none); of its attributes and textures in flat arrays */
function createMemberVersions() {
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
    sites: -1,
    elements: -1,
  };
}

/**
 * The instances of a merged object kept elsewhere (a page store, `page-store.ts`), culled and composed on the GPU by a
 * culler of its own (`merged-paged-culler.ts`): a member's instances its `instanceCount` (the store's live count of
 * it), so the segments come out as the store's members; its per-instance layers the source's cells (written and logged
 * by range, `TextureChanges`). An instance count changed is not a new layout of the geometry: the segments and the
 * segment table made again alone
 */
export interface MergedInstanceSource {
  /** the layers of instance granularity by texture key (`tColor`, `tMarker`, ...): the merged object's own cells */
  readonly layers: ReadonlyMap<string, ValueCell<DataArray<any>>>;
  /**
   * the markers' aggregates of its marker layer (`Visual.mark` on the layer sets them): the merged object's own cells,
   * not the members' (their instances have no markers of their own); the marked lists not used (no counts a member)
   */
  readonly markers?: {
    readonly markerAverage: ValueCell<number>;
    readonly markerStatus: ValueCell<number>;
    readonly uMarker: ValueCell<number>;
  };
  /**
   * the members whose cells it changed since the last call (their counts, bounds, a fill, `fillMergeRecord`): what
   * `sync` reads of the members, the others' cells not read (undefined: all). A values object's reads are slow and
   * scattered, thousands of members an update with a few of them changed
   */
  takeMembers?(): Iterable<number> | undefined;
  /**
   * per member its instance count and placed bounds as its cells hold them, the source's typed copies: what `sync` and
   * the culler read of the members `take` tells of (their count or bounds changed alone), not their cells
   */
  readonly counts?: MergedSourceCounts;
  /**
   * a page store's instances (`ParticlePageStore`, `paged-cull.ts`): culled and composed on the GPU by a culler of its
   * own (`merged-paged-culler.ts`); its changes taken at each sync, its members' reaches (a member's
   * geometry's from its center: a particle's sphere its scale times it)
   */
  readonly paged: {
    readonly source: PagedCullSource;
    takeChanges(): PagedCullChanges;
    readonly reaches: ArrayLike<number>;
    /**
     * its copies' markers (a nested list's: a template particle drawn once a copy from one store index, its store
     * marker every copy's), uploaded where its version changed: the words (two a copy), then the blocks' bytes; none
     * without copies
     */
    readonly copyMarkers?: {
      readonly version: number;
      readonly words: Uint32Array;
      readonly bytes: Uint8Array;
      readonly byteLength: number;
      /** the bytes changed since the last call (`[start, end)` pairs in `bytes`), undefined all of them */
      takeChanged(): readonly number[] | undefined;
    };
    /**
     * its copies' own layers (a nested list's, `CulledSchema.sCopyLayers`): uploaded whole where its version changed;
     * none without a copy's
     */
    readonly copyLayers?: { readonly version: number; readonly buffer: Uint32Array };
    /** the copies' instances a frame draws at most, by their levels (`PagedCullView.maxInstances`; 0: no cap) */
    readonly maxInstances?: number;
  };
}

/** A source's members' instance counts and placed bounds (`MergedInstanceSource.counts`) */
export interface MergedSourceCounts {
  /** per member its instance count, as its `instanceCount` */
  readonly instanceCounts: ArrayLike<number>;
  /** per member its placed bounds (center, radius: four a member), as its `boundingSphere` */
  readonly spheres: ArrayLike<number>;
  /** the members whose count or bounds it changed since the last call (their cells with them, not read) */
  take(): Iterable<number>;
}

/**
 * Values like `like` (of `type`) without geometry: no vertices nor draws, no instances, every cell its own (`like`'s
 * values in cells of its own: a cell written later changes these alone). A source's shared values and its empty record
 * (`pickMergeRecord` of them: a member not filled yet, a segment that draws nothing)
 */
export function createEmptyMergeMember(type: MergedType, like: MergeableValues): MergeableValues {
  const { vertexAttributes, hasElements } = MergedTypeDescriptors[type];
  const values: { [k: string]: ValueCell<any> } = {};
  for (const [k, cell] of Object.entries(like as unknown as { [k: string]: ValueCell<any> })) {
    values[k] = ValueCell.create(cell.ref.value);
  }
  for (const { key, ctor } of vertexAttributes) values[key] = ValueCell.create(new ctor(0));
  if (hasElements) values.elements = ValueCell.create(new Uint32Array(0));
  values.drawCount = ValueCell.create(0);
  values.uVertexCount = ValueCell.create(0);
  values.instanceCount = ValueCell.create(0);
  // no grid of its own (`like`'s may be pending): a cell of all its instances
  values.instanceGrid = ValueCell.create(createEmptyInstanceGrid());
  values.markerAverage = ValueCell.create(0);
  values.markerStatus = ValueCell.create(0);
  return values as unknown as MergeableValues;
}

/** What a source holds of its members (`MergedInstanceSource`): their instances, grids and markers, not a record's */
const SourceHeldKeys = new Set<string>([
  'instanceCount',
  'aTransform',
  'aInstance',
  'instancePoses',
  'instanceGrid',
  'instanceCellSize',
  'markerStatus',
  'markerAverage',
  'markerType',
  'markerInstanceCounts',
  'markedInstanceCount',
]);

/**
 * The keys of a record (`MergeRecord`) a member of a source fills, of a merged object of `type` with the shared values
 * `shared` and the source's layers `sourceLayers`: its geometry, counts, bounds, levels, sites, aggregates and layers
 * of its own; not what the source holds (`SourceHeldKeys`, its layers of instances) nor the layers of no granularity
 * (the shared values'). The other keys' cells may be shared by the source's records
 */
export function getSourceRecordFillKeys(
  type: MergedType,
  shared: MergeableValues,
  sourceLayers: ReadonlySet<string>,
): string[] {
  const { extraTextures } = MergedTypeDescriptors[type];
  const layers = new Map<string, MergedTexSpec>();
  for (const spec of [...extraTextures, ...BaseTexSpecs]) layers.set(spec.tex, spec);
  return getMergeRecordKeys(type, shared).filter((k) => {
    if (SourceHeldKeys.has(k)) return false;
    const spec = layers.get(k);
    if (!spec) return true;
    const granularity = spec.granularity(shared);
    return granularity !== 'none' && !(granularity === 'instance' && sourceLayers.has(k));
  });
}

/**
 * `record` filled with `values` (of the same type and defines): its cells of `keys` given `values`' value where another
 * (`getSourceRecordFillKeys`: its own cells); with a source, its geometry then another (`Merged.replaceGeometry`)
 */
export function fillMergeRecord(record: MergeRecord, values: MergeableValues, keys: readonly string[]) {
  for (const k of keys) {
    const cell = values[k];
    const target = record[k];
    if (cell && target && target.ref.value !== cell.ref.value) ValueCell.update(target, cell.ref.value);
  }
}

/**
 * A merged object of `members` (their records) with the values all of them share `shared` (of their type and merge
 * key: its defines, uniforms, props and clip the object's); `gpuCulling`: culled on the GPU where available
 * (`merged-culler.ts`), fixed for the merged object's lifetime;
 * `rasterOnly`: for the sphere raster only, culled on the CPU where the raster does not take it (`Renderable.rasterSpheres`);
 * `source`: the instances kept elsewhere (`MergedInstanceSource`)
 */
export function createMergedValues(
  type: MergedType,
  shared: MergeableValues,
  members: readonly MergeRecord[],
  gpuCulling = false,
  rasterOnly = false,
  source?: MergedInstanceSource,
  recordsVersion?: ValueCell<number>,
): Merged {
  const descriptor = MergedTypeDescriptors[type];
  const { hasElements, vertexAttributes } = descriptor;
  const n = members.length;

  const texSpecs = [...descriptor.extraTextures, ...BaseTexSpecs].filter((spec) => shared[spec.tex] !== undefined);
  /** the second color's type (`dSiteColor`): a define, the shared values' (read again: a define changed rebuilds) */
  const siteColorTypeOf = () => (shared as { dSiteColor?: ValueCell<string> }).dSiteColor?.ref.value;

  const segments: MergedSegments = {
    count: n,
    instanceBases: new Array(n).fill(0),
    instanceCounts: new Array(n).fill(0),
    vertexBases: new Array(n).fill(0),
    drawBases: new Array(n).fill(0),
    elementOffsets: new Array(n).fill(0),
    groupCounts: new Array(n).fill(0),
    groupInstanceOffsets: new Array(n).fill(0),
    groupBases: new Array(n).fill(0),
  };

  // own cells, identity must stay stable for the lifetime of the merged render item
  const own: { [k: string]: ValueCell<any> } = {
    aTransform: ValueCell.create(new Float32Array(0)),
    // the members' parts where one of them has poses (`transform-parts.ts`), none where all have matrices
    instancePoses: ValueCell.create(undefined),
    aInstance: ValueCell.create(new Uint32Array(0)),
    aSegment: ValueCell.create(new Int32Array(0)),
    tSegment: ValueCell.create(createDataArray(n * SegmentTexels, 4, Float32Array)),
    aVertex: ValueCell.create(new Float32Array(0)),
    drawCount: ValueCell.create(0),
    drawCountMax: ValueCell.create(0),
    instanceCount: ValueCell.create(0),
    uVertexCount: ValueCell.create(0),
    uInstanceCount: ValueCell.create(0),
    uGroupCount: ValueCell.create(0),
    boundingSphere: ValueCell.create(Sphere3D()),
    invariantBoundingSphere: ValueCell.create(Sphere3D()),
    uInvariantBoundingSphere: ValueCell.create(Vec4()),
    // the transient transform of the whole object, not a member's: a merged object moves as one (`Merged`)
    uObjectTransform: ValueCell.create(Mat4.identity()),
    uMarker: ValueCell.create(0),
    markerAverage: ValueCell.create(0),
    markerStatus: ValueCell.create(-1),
    // the members' are read (`Renderable.renderMarked`), not the first's
    markerInstanceCounts: ValueCell.create(new Uint32Array(0)),
    markedInstanceCount: ValueCell.create(0),
    transparencyAverage: ValueCell.create(0),
    transparencyMin: ValueCell.create(0),
    emissiveAverage: ValueCell.create(0),
    wiggleAverage: ValueCell.create(0),
    hasReflection: ValueCell.create(false),
    instanceGrid: ValueCell.create({ ...createEmptyInstanceGrid(), cellSize: 2 }),
    instanceCellSize: ValueCell.create(2),
    instanceBatchSize: ValueCell.create(0),
    dSegmented: ValueCell.create(true),
    // drawn with instance indices where culled on the GPU: set by the renderable's culler
    dCulled: ValueCell.create(false),
    dCulledPaged: ValueCell.create(false),
    dCopyMarkers: ValueCell.create(false),
    dCopyLayers: ValueCell.create(false),
    gpuCulling: ValueCell.create(gpuCulling),
    gpuCullingRasterOnly: ValueCell.create(rasterOnly),
  };
  // the members' geometry is not copied: the render item uploads their arrays as parts (`SegmentedArray`)
  for (const { key, ctor } of vertexAttributes) {
    own[key] = ValueCell.create(SegmentedArray.create<SegmentedArrayType>(ctor, []));
  }
  if (hasElements) {
    own.elements = ValueCell.create(SegmentedArray.create(Uint32Array, []));
  }
  /** the layers the source holds (`MergedInstanceSource.layers`): its cells, not merged from the members */
  const sourced = new Set<string>();
  for (const spec of texSpecs) {
    const layer = source?.layers.get(spec.tex);
    // a layer off now too (`SourceOptionalLayers`): turned on later (its define, its type the instance) it reads them
    const granularity = spec.granularity(shared);
    if (layer && (granularity === 'instance' || (granularity === 'none' && SourceOptionalLayers.has(spec.tex)))) {
      own[spec.tex] = layer;
      sourced.add(spec.tex);
      continue;
    }
    const ctor = getTexArrayCtor(spec) as new (length: number) => any;
    own[spec.tex] = ValueCell.create(createDataArray(1, spec.itemSize, ctor));
  }

  // a source's markers' aggregates its cells (`MergedInstanceSource.markers`)
  if (source?.markers) Object.assign(own, source.markers);

  // the sphere raster's values set per view, the merged object's own (not the shared values')
  if ('dSphereResolve' in shared) Object.assign(own, createSphereResolveValues());
  const values = {
    ...shared,
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
      if (spec.granularity(shared) !== 'group') continue;
      const { count } = v[spec.tex]!.ref.value as DataArray<any>;
      if (count > span) span = count;
    }
    return span;
  }

  function getTexelCount(i: number, spec: MergedTexSpec): number {
    const v = members[i];
    const factor = spec.factor ? spec.factor(shared) : 1;
    switch (spec.granularity(shared)) {
      case 'instance':
        return memberInstanceCount(i) * factor;
      case 'group':
        return getGroupSpan(i) * factor;
      case 'groupInstance':
        return memberInstanceCount(i) * v.uGroupCount.ref.value * factor;
      case 'vertex':
        // the data vertices (`Geometry.getVertexCount`)
        return v.uVertexCount.ref.value * factor;
      case 'item':
        return descriptor.itemTexels(v) * factor;
      case 'table':
        return (v[spec.tex]!.ref.value as DataArray<any>).count * factor;
      case 'none':
        return 0;
    }
  }

  /** texel base per member for each merged texture */
  const texBases = new Map<string, number[]>();
  for (const spec of texSpecs) {
    texBases.set(spec.tex, new Array(n).fill(0));
  }

  /**
   * With a source: the members' geometry by allocation, not by member (a pool's member filled with another target's
   * geometry, `replaceGeometry`): an allocation a member's counts at the time, at the end of each array (its parts
   * added, `SegmentedArray.extend`: the render item writes those alone), a member's last its own, those before it
   * unused until the arrays are laid out again (whole, with room, `SourceHeadroom`: where the room is used up or the
   * unused draw more than the used)
   */
  type GeometryAlloc = {
    readonly member: number;
    readonly vertexCount: number;
    readonly drawCount: number;
    readonly groupSpan: number;
    readonly sizeTable: number;
    readonly siteColorTable: number;
    /** per texture spec its texels (none for a layer of the source's or of no granularity) */
    readonly texels: number[];
    vertexBase: number;
    drawBase: number;
    groupBase: number;
    sizeBase: number;
    siteColorBase: number;
    readonly texBases: number[];
  };
  let allocs: GeometryAlloc[] = [];
  const liveAlloc = new Int32Array(n);
  /** the arrays' room (their lengths): vertices, draws, per texture spec texels */
  let roomVertices = 0;
  let roomDraws = 0;
  let roomTexels: number[] = [];
  /** the allocations' ends */
  let allocEnd = { vertex: 0, draw: 0 };
  /** every array's end after the allocations laid out (`layoutAllocs`): those appended laid out from there alone */
  const layoutEnd = { vertex: 0, draw: 0, group: 0, size: 0, siteColor: 0, tex: [] as number[] };
  /**
   * the members whose allocation is another since `computeSourceSegments` (their bases to write): all after a whole
   * layout, else those appended
   */
  let segmentAllocs: Set<number> | 'all' = 'all';
  /** per member its live allocation's draw count (`drawCountMax`: a max over a typed array, not the allocations) */
  const memberDraws = new Int32Array(n);
  /** members whose geometry is another since the last sync (`replaceGeometry`) */
  const replaced = new Set<number>();
  /** the segment table's rows to write at the next `rebuildSourceInstances`: all after a whole layout, else those replaced */
  let segmentRows: Set<number> | 'all' = 'all';

  // the cells `sync` reads per member, resolved once into objects of one shape: values objects have too many keys for
  // fast named reads, and a sync reads some 30 cells of each member
  const memberCells = members.map((v) => getMemberCells(v, hasElements));
  const memberVersions = members.map(() => createMemberVersions());
  // a member's vertex attributes' and merged textures' cells and versions, flat (member `i`'s attribute `a` at
  // `i * attrCount + a`, texture `t` at `i * texCount + t`): no arrays a member (some 0.4 KB each with 100k members)
  const attrCount = vertexAttributes.length;
  const texCount = texSpecs.length;
  const attributeCells: Cell[] = new Array(n * attrCount);
  const textureCells: Cell[] = new Array(n * texCount);
  for (let i = 0; i < n; ++i) {
    const v = members[i];
    for (let a = 0; a < attrCount; ++a) attributeCells[i * attrCount + a] = v[vertexAttributes[a].key]!;
    for (let t = 0; t < texCount; ++t) textureCells[i * texCount + t] = v[texSpecs[t].tex]!;
  }
  const attributeVersions = new Int32Array(n * attrCount).fill(-1);
  const textureVersions = new Int32Array(n * texCount).fill(-1);

  /** per member its instance, group, group span, vertex and draw counts: a change is a new layout */
  const layout = new Float64Array(n * 5).fill(-1);
  let layoutMergeKey = '';
  /** with a source: an instance count changed alone at the last `updateLayout` (not a new layout of the geometry) */
  let instanceCountsChanged = false;
  /** per member a count advanced where `sync` read its cells as changed (`Merged.memberChanges`) */
  const memberChanges = new Uint32Array(n);
  const allMembers = fillSerial(new Int32Array(n));
  /**
   * the members a sync reads: all, or with a source those it changed (`MergedInstanceSource.takeMembers`) after the
   * first sync; the others' cells as the last sync saw them
   */
  let readList: ArrayLike<number> = allMembers;
  /** with a source's counts (`MergedInstanceSource.counts`): the members a sync reads the count and bounds of alone */
  let countList: ArrayLike<number> = [];
  const counts = source?.counts;
  /**
   * member `i`'s instance count: a paged source's copy (`MergedSourceCounts.instanceCounts`, its members' cells not
   * written: thousands a step), else its cell
   */
  const pagedCounts = source?.paged ? counts : undefined;
  function memberInstanceCount(i: number): number {
    return pagedCounts ? pagedCounts.instanceCounts[i] : members[i].instanceCount.ref.value;
  }
  let synced = false;

  const markerTexture = texSpecs.findIndex((s) => s.tex === 'tMarker');
  const observedChanged: number[] = [];
  /**
   * the members whose observed cells changed (`memberCells`): two versions a member compared, not all a sync reads (a
   * hover over 100k members: most of its frame)
   */
  function getObservedChanged() {
    observedChanged.length = 0;
    const t = markerTexture >= 0 && texSpecs[markerTexture].granularity(shared) !== 'none' ? markerTexture : -1;
    for (let i = 0; i < n; ++i) {
      if (
        memberVersions[i].markerStatus !== memberCells[i].markerStatus.ref.version ||
        (t >= 0 && textureVersions[i * texCount + t] !== textureCells[i * texCount + t].ref.version)
      ) {
        observedChanged.push(i);
      }
    }
    return observedChanged;
  }

  /** Whether the layout changed since the last call (also the shared values' merge key), stores the new one */
  function updateLayout() {
    instanceCountsChanged = false;
    let changed = false;
    const groupTextures: number[] = [];
    for (let t = 0; t < texSpecs.length; ++t) {
      if (texSpecs[t].granularity(shared) === 'group') groupTextures.push(t);
    }
    for (let k = 0; k < readList.length; ++k) {
      const i = readList[k];
      const c = memberCells[i];
      setAggregatesDefault(i);
      let span: number = c.uGroupCount.ref.value;
      for (const t of groupTextures) {
        const { count } = textureCells[i * texCount + t].ref.value as DataArray<any>;
        if (count > span) span = count;
      }
      const o = i * 5;
      // a source's instance counts their own change: the segments again, not the geometry
      const instanceCount = memberInstanceCount(i);
      if (source && layout[o] !== instanceCount && layout[o] !== -1) {
        layout[o] = instanceCount;
        instanceCountsChanged = true;
      }
      // a source's member of another geometry: allocated again (`appendSourceGeometry`), not a new layout
      if (replaced.has(i) && layout[o] !== -1) {
        layout[o + 1] = c.uGroupCount.ref.value;
        layout[o + 2] = span;
        layout[o + 3] = c.uVertexCount.ref.value;
        layout[o + 4] = c.drawCount.ref.value;
        continue;
      }
      if (
        layout[o] !== instanceCount ||
        layout[o + 1] !== c.uGroupCount.ref.value ||
        layout[o + 2] !== span ||
        layout[o + 3] !== c.uVertexCount.ref.value ||
        layout[o + 4] !== c.drawCount.ref.value
      ) {
        layout[o] = instanceCount;
        layout[o + 1] = c.uGroupCount.ref.value;
        layout[o + 2] = span;
        layout[o + 3] = c.uVertexCount.ref.value;
        layout[o + 4] = c.drawCount.ref.value;
        changed = true;
      }
    }
    // the count changed alone: the source's copy, not the member's cell
    for (let k = 0; k < countList.length; ++k) {
      const o = countList[k] * 5;
      const count = counts!.instanceCounts[countList[k]];
      if (layout[o] !== -1 && layout[o] !== count) {
        layout[o] = count;
        instanceCountsChanged = true;
      }
    }
    const mergeKey = getMergeKey(shared);
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
    const part = layerArray(spec).parts[source && isAllocTexture(spec) ? liveAlloc[i] : i] as SegmentedArray.Part;
    part.array = (members[i][spec.tex]!.ref.value as DataArray<any>).array;
    part.version += 1;
  }

  function mergeTextureFull(spec: MergedTexSpec) {
    // the source's layer: its cell as it is
    if (sourced.has(spec.tex)) return;
    const cell = own[spec.tex];
    const bases = texBases.get(spec.tex)!;
    const ctor = getTexArrayCtor(spec);

    if (spec.granularity(shared) === 'none') {
      // no per-segment data used, the shared values' layer
      const src = shared[spec.tex].ref.value as DataArray<any>;
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
      const { array } = members[i][spec.tex]!.ref.value as DataArray<any>;
      parts.push({ array, offset: total * spec.itemSize, count: texels * spec.itemSize, bias: 0, version: 0 });
      total += texels;
    }
    // whole words, as `createDataArray`
    const count = Math.max(1, total);
    const length = Math.ceil((count * spec.itemSize * ctor.BYTES_PER_ELEMENT) / 4) * (4 / ctor.BYTES_PER_ELEMENT);
    ValueCell.update(cell, { array: SegmentedArray.create(ctor, parts, length), count });
  }

  /** the members' spheres, those drawing nothing left out (as the scene leaves out such objects, `Scene`) */
  function unionSpheres(out: Sphere3D, sphere: (i: number) => Sphere3D) {
    boundaryHelper.reset();
    for (let i = 0; i < n; ++i) {
      if (!memberCells[i].drawCount.ref.value) continue;
      const s = sphere(i);
      if (s.radius) boundaryHelper.includeSphere(s);
    }
    boundaryHelper.finishedIncludeStep();
    for (let i = 0; i < n; ++i) {
      if (!memberCells[i].drawCount.ref.value) continue;
      const s = sphere(i);
      if (s.radius) boundaryHelper.radiusSphere(s);
    }
    return boundaryHelper.getSphere(out);
  }

  // by the cells resolved once (`memberCells`): a member's values object's named reads slow
  const placedSphere = (i: number) => memberCells[i].boundingSphere.ref.value as Sphere3D;
  const invariantSphere = (i: number) => memberCells[i].invariantBoundingSphere.ref.value as Sphere3D;

  /** with a source: per member its placed and invariant spheres (center, radius), of `readList` read again a sync */
  const placedSpheres = new Float64Array(source ? n * 4 : 0);
  const invariantSpheres = new Float64Array(source ? n * 4 : 0);
  /**
   * a source's: whether a member's sphere read changed since its object's spheres were last unioned
   * (`updateBoundingSpheres`: a pass over every member, not again for a step's counts alone)
   */
  let spheresDirty = true;
  function readSpheres() {
    for (let k = 0; k < readList.length; ++k) {
      const i = readList[k];
      const p = placedSphere(i);
      const v = invariantSphere(i);
      const o = i * 4;
      // NaN as itself (a member of no bounds)
      if (
        !Object.is(placedSpheres[o], p.center[0]) ||
        !Object.is(placedSpheres[o + 1], p.center[1]) ||
        !Object.is(placedSpheres[o + 2], p.center[2]) ||
        !Object.is(placedSpheres[o + 3], p.radius) ||
        !Object.is(invariantSpheres[o], v.center[0]) ||
        !Object.is(invariantSpheres[o + 1], v.center[1]) ||
        !Object.is(invariantSpheres[o + 2], v.center[2]) ||
        !Object.is(invariantSpheres[o + 3], v.radius)
      ) {
        spheresDirty = true;
      }
      placedSpheres[o] = p.center[0];
      placedSpheres[o + 1] = p.center[1];
      placedSpheres[o + 2] = p.center[2];
      placedSpheres[o + 3] = p.radius;
      invariantSpheres[o] = v.center[0];
      invariantSpheres[o + 1] = v.center[1];
      invariantSpheres[o + 2] = v.center[2];
      invariantSpheres[o + 3] = v.radius;
    }
    readCountSpheres();
  }

  /** the placed spheres of `countList` from the source's copies (`MergedInstanceSource.counts`), whether one changed */
  function readCountSpheres() {
    let changed = false;
    if (!counts) return changed;
    const { spheres } = counts;
    for (let k = 0; k < countList.length; ++k) {
      const o = countList[k] * 4;
      for (let j = 0; j < 4; ++j) {
        // NaN as itself (a member of no bounds of its own: its copy's NaN)
        if (!Object.is(placedSpheres[o + j], spheres[o + j])) {
          placedSpheres[o + j] = spheres[o + j];
          changed = true;
        }
      }
    }
    if (changed) spheresDirty = true;
    return changed;
  }

  /**
   * The box of the spheres' boxes as a sphere: a source's (its members' spheres boxes of their particles, thousands of
   * members at each change of their counts; the boundary helper's tighter sphere costs a pass more and more a sphere)
   */
  function unionSpheresBox(out: Sphere3D, spheres: Float64Array) {
    let x0 = Infinity;
    let y0 = Infinity;
    let z0 = Infinity;
    let x1 = -Infinity;
    let y1 = -Infinity;
    let z1 = -Infinity;
    for (let o = 0; o < spheres.length; o += 4) {
      const cx = spheres[o];
      const cy = spheres[o + 1];
      const cz = spheres[o + 2];
      const r = spheres[o + 3];
      if (cx - r < x0) x0 = cx - r;
      if (cy - r < y0) y0 = cy - r;
      if (cz - r < z0) z0 = cz - r;
      if (cx + r > x1) x1 = cx + r;
      if (cy + r > y1) y1 = cy + r;
      if (cz + r > z1) z1 = cz + r;
    }
    const c = out.center;
    if (!(x1 >= x0)) {
      c[0] = c[1] = c[2] = 0;
      out.radius = 0;
      return out;
    }
    c[0] = (x0 + x1) / 2;
    c[1] = (y0 + y1) / 2;
    c[2] = (z0 + z1) / 2;
    // about the box's center, as far as the farthest sphere reaches: not the box's own sphere (a sphere's box's corner
    // a third again as far, √3 its radius: a page store's object, its sphere one member's, so far too large)
    let radius = 0;
    for (let o = 0; o < spheres.length; o += 4) {
      const r = spheres[o + 3];
      if (!(r >= 0)) continue;
      const d = Math.hypot(spheres[o] - c[0], spheres[o + 1] - c[1], spheres[o + 2] - c[2]) + r;
      if (d > radius) radius = d;
    }
    out.radius = radius;
    return out;
  }

  /** as placed: the object's transient transform (`uObjectTransform`) moves what is drawn, not these */
  function updateBoundingSpheres() {
    if (source) {
      readSpheres();
      if (!spheresDirty) return;
      spheresDirty = false;
    }
    const ibs = source
      ? unionSpheresBox(own.invariantBoundingSphere.ref.value, invariantSpheres)
      : unionSpheres(own.invariantBoundingSphere.ref.value, invariantSphere);
    ValueCell.update(
      own.boundingSphere,
      source
        ? unionSpheresBox(own.boundingSphere.ref.value, placedSpheres)
        : unionSpheres(own.boundingSphere.ref.value, placedSphere),
    );
    ValueCell.update(own.invariantBoundingSphere, ibs);
    ValueCell.update(own.uInvariantBoundingSphere, Vec4.fromSphere(own.uInvariantBoundingSphere.ref.value, ibs));
  }

  /**
   * per member its marker average by its weight, its weight and its status, their sums and the members by status: a
   * sync of a few members (a mark) updates theirs (`updateMarkerAggregates`), not a pass over every member (a hover over
   * 100k members: a tenth of its frame)
   */
  const markerSums = new Float64Array(n);
  const markerWeights = new Float64Array(n);
  /** `NoStatus`: not counted yet */
  const markerStatuses = new Int8Array(n);
  const NoStatus = 127;
  const markerStatusCounts = new Map<number, number>();
  let markerSum = 0;
  let markerWeight = 0;
  /** members with a marker average other than 0: none, the object's average 0 whatever the sums' rounding */
  let markedMembers = 0;
  let markerAggregatesInit = false;
  /** the members with a marker average other than 0 (`Merged.getMarkedMembers`), their list made where read */
  const markedSet = new Set<number>();
  let markedList: number[] | undefined;
  function getMarkedMembers() {
    return (markedList ??= Array.from(markedSet).sort((a, b) => a - b));
  }
  function setMarkerAggregate(i: number) {
    const v = members[i];
    const count = memberInstanceCount(i);
    const w = Math.max(1, v.markerType.ref.value === 'instance' ? count : count * v.uGroupCount.ref.value);
    const sum = v.markerAverage.ref.value * w;
    const status = v.markerStatus.ref.value;
    const previous = markerStatuses[i];
    if (previous !== NoStatus) markerStatusCounts.set(previous, markerStatusCounts.get(previous)! - 1);
    if (markerStatusCounts.get(previous) === 0) markerStatusCounts.delete(previous);
    markerStatusCounts.set(status, (markerStatusCounts.get(status) ?? 0) + 1);
    markerStatuses[i] = status;
    if ((sum !== 0) !== (markerSums[i] !== 0)) {
      markedMembers += sum !== 0 ? 1 : -1;
      if (sum !== 0) markedSet.add(i);
      else markedSet.delete(i);
      markedList = undefined;
    }
    markerSum += sum - markerSums[i];
    markerWeight += w - markerWeights[i];
    markerSums[i] = sum;
    markerWeights[i] = w;
  }

  /** the object's marker aggregates by its members': those of `list` changed alone (all where none given) */
  function updateMarkerAggregates(list?: ArrayLike<number>) {
    // a source's: set where its layer is marked
    if (source?.markers) return;
    // all of them, summed again: no rounding carried over from the updates of a few
    if (!list || list.length * 8 > n || !markerAggregatesInit) {
      markerAggregatesInit = true;
      markerSums.fill(0);
      markerWeights.fill(0);
      markerStatuses.fill(NoStatus);
      markerStatusCounts.clear();
      markedSet.clear();
      markedList = undefined;
      markerSum = markerWeight = markedMembers = 0;
      for (let i = 0; i < n; ++i) setMarkerAggregate(i);
    } else {
      for (let k = 0, kl = list.length; k < kl; ++k) setMarkerAggregate(list[k]);
    }
    let status = -1;
    if (markerStatusCounts.size === 1) for (const s of markerStatusCounts.keys()) status = s;
    ValueCell.updateIfChanged(own.markerAverage, markedMembers > 0 && markerWeight > 0 ? markerSum / markerWeight : 0);
    ValueCell.updateIfChanged(own.markerStatus, status);
    ValueCell.updateIfChanged(own.uMarker, status);
  }

  /**
   * per member whether its aggregates' values are the defaults (no transparency, emission, wiggle nor reflection), and
   * how many are not: none, the object's aggregates the defaults whatever the weights (`updateAggregates`, a pass over
   * every member, not a member's geometry filled: a tile set's new targets as its first tiles load)
   */
  const aggregatesOther = new Uint8Array(n);
  let aggregatesOtherCount = 0;
  /** every member's read once (the first `updateAggregates`, before a layout read them all) */
  let aggregatesInit = false;
  function setAggregatesDefault(i: number) {
    const v = members[i];
    const other =
      v.transparencyAverage.ref.value !== 0 ||
      v.emissiveAverage.ref.value !== 0 ||
      v.wiggleAverage.ref.value !== 0 ||
      v.transparencyMin.ref.value !== 1 ||
      !!v.hasReflection.ref.value
        ? 1
        : 0;
    aggregatesOtherCount += other - aggregatesOther[i];
    aggregatesOther[i] = other;
  }

  function updateAggregates() {
    if (!aggregatesInit) {
      aggregatesInit = true;
      for (let i = 0; i < n; ++i) setAggregatesDefault(i);
    }
    if (aggregatesOtherCount === 0) {
      ValueCell.updateIfChanged(own.transparencyAverage, 0);
      ValueCell.updateIfChanged(own.transparencyMin, 1);
      ValueCell.updateIfChanged(own.emissiveAverage, 0);
      ValueCell.updateIfChanged(own.wiggleAverage, 0);
      ValueCell.updateIfChanged(own.hasReflection, false);
      return;
    }
    let transparencyAverage = 0;
    let transparencyMin = 1;
    let emissiveAverage = 0;
    let wiggleAverage = 0;
    let hasReflection = false;
    let weight = 0;
    for (let i = 0; i < n; ++i) {
      const v = members[i];
      const w = Math.max(1, memberInstanceCount(i) * v.uGroupCount.ref.value);
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
      array: members[i][attr.key]!.ref.value as Float32Array | Uint32Array | Int32Array,
      offset: segments.vertexBases[i] * attr.itemSize,
      count: members[i].uVertexCount.ref.value * attr.itemSize,
      bias: 0,
      version: 0,
    };
  }

  /** member i's elements as a part: its indices at its element offset, made absolute by its vertex base */
  function elementsPart(i: number): SegmentedArray.Part {
    return {
      array: members[i].elements!.ref.value as Uint32Array,
      offset: segments.elementOffsets[i] / 4,
      count: members[i].drawCount.ref.value,
      bias: segments.vertexBases[i],
      version: 0,
    };
  }

  /** member i's part of `key` read again from the member (its array may be another of the same size) */
  function updatePart(key: string, i: number) {
    const part = (own[key].ref.value as SegmentedArray).parts[source ? liveAlloc[i] : i] as SegmentedArray.Part;
    part.array = members[i][key]!.ref.value;
    part.version += 1;
  }

  /**
   * member i's own (sphere-count-dependent) per-level scales and strides, its first item in the merged items and its
   * quantized centers' origin and scale (spheres' `dSites`; elements 1 to 4 of its segment data): quad `p` of level `j`
   * draws the member's item `p * stride_j` (the impostors' quads)
   */
  function writeSegmentLod(i: number, array: Float32Array) {
    const member = members[i];
    // gridless members are drawn at full (undecimated) detail, see the
    // `hasInstanceCells` branch in renderable.ts's cullSegment - no per-level
    // scale compensation or stride applies then, else spheres would render too large, or items twice. A source's
    // members by their levels whatever their grids (culled by level on the GPU, a CPU cull building their grids): the
    // culler's records count a level's items alike (`layoutRegions`), else a stride over the whole count reads past
    // the member's items
    const gridless = !source && !hasInstanceCells(member);
    const lodLevels = !gridless ? (member.lodLevels?.ref.value as LodLevelsValue | undefined) : undefined;
    const o = (i * SegmentTexels + 1) * 4;
    for (let j = 0; j < 4; ++j) {
      array[o + j] = lodLevels?.[j]?.[3] ?? 1;
      array[o + 4 + j] = lodLevels?.[j]?.[4] ?? 1;
    }
    // high and low 16 bits: exact as floats
    const first = segments.drawBases[i] / 6;
    array[o + 8] = Math.floor(first / 65536);
    array[o + 9] = first % 65536;
    // the member's first entry of the merged size table (`siteType` sizes: the members' tables packed)
    array[o + 10] = sizeTableBases[i];
    // the member's first entry of the merged second color table (`dSiteColor` 'siteType')
    array[o + 11] = siteColorTableBases[i];
    const origin = (member as { uSiteOrigin?: ValueCell<Vec3> }).uSiteOrigin?.ref.value;
    array[o + 12] = origin?.[0] ?? 0;
    array[o + 13] = origin?.[1] ?? 0;
    array[o + 14] = origin?.[2] ?? 0;
    array[o + 15] = (member as { uSiteScale?: ValueCell<number> }).uSiteScale?.ref.value ?? 1;
  }

  /** per member its first entry of the merged size table (`siteType` sizes, `table` granularity), 0 otherwise */
  const sizeTableBases = new Array<number>(n).fill(0);
  /** per member its first entry of the merged second color table (`dSiteColor` 'siteType'), 0 otherwise */
  const siteColorTableBases = new Array<number>(n).fill(0);

  /**
   * whether a texture of `spec` is merged by allocation: the source's layers, those of no granularity and those read by
   * the segment's index (`bySegment`: a member's texel its index, an allocation's end not) not
   */
  const isAllocTexture = (spec: MergedTexSpec) =>
    !sourced.has(spec.tex) && spec.granularity(shared) !== 'none' && !spec.bySegment?.(shared);
  /** whether a texture of `spec` is a source's table by member (`bySegment`): its member's part read again by itself */
  const isMemberTable = (spec: MergedTexSpec) =>
    !!source && !sourced.has(spec.tex) && spec.granularity(shared) !== 'none' && !!spec.bySegment?.(shared);

  /**
   * The members' parts of the tables by member (`isMemberTable`) of `list` read again where their textures changed:
   * a member allocated again (`replaceGeometry`, `relayoutSourceInstances`) its table's texel as its index's, not
   * an allocation's (the texture loop of `sync` not reached then)
   */
  function updateMemberTables(list: ArrayLike<number>) {
    for (let t = 0; t < texSpecs.length; ++t) {
      const spec = texSpecs[t];
      if (!isMemberTable(spec)) continue;
      let changed = false;
      for (let k = 0; k < list.length; ++k) {
        const i = list[k];
        const version = textureCells[i * texCount + t].ref.version;
        if (textureVersions[i * texCount + t] === version) continue;
        updateTexPart(spec, i);
        textureVersions[i * texCount + t] = version;
        changed = true;
      }
      if (changed) ValueCell.update(own[spec.tex], own[spec.tex].ref.value);
    }
  }

  function createAlloc(i: number): GeometryAlloc {
    const v = members[i] as MergeRecord & {
      sizeType?: ValueCell<string>;
      tSize?: ValueCell<DataArray<Uint8Array>>;
      tSiteColor?: ValueCell<DataArray<Uint8Array>>;
    };
    return {
      member: i,
      vertexCount: v.uVertexCount.ref.value,
      drawCount: v.drawCount.ref.value,
      groupSpan: getGroupSpan(i),
      sizeTable: v.sizeType?.ref.value === 'siteType' ? v.tSize!.ref.value.count : 0,
      siteColorTable: siteColorTypeOf() === 'siteType' ? v.tSiteColor!.ref.value.count : 0,
      texels: texSpecs.map((spec) => (isAllocTexture(spec) ? getTexelCount(i, spec) : 0)),
      vertexBase: 0,
      drawBase: 0,
      groupBase: 0,
      sizeBase: 0,
      siteColorBase: 0,
      texBases: texSpecs.map(() => 0),
    };
  }

  /**
   * the allocations' bases (each array's in their order) from `from` on, those before as they were (appended alone:
   * thousands laid out a member filled otherwise), and the arrays' ends
   */
  function layoutAllocs(from = 0) {
    const e = layoutEnd;
    let vertex = from > 0 ? e.vertex : 0;
    let draw = from > 0 ? e.draw : 0;
    let group = from > 0 ? e.group : 0;
    let size = from > 0 ? e.size : 0;
    let siteColor = from > 0 ? e.siteColor : 0;
    const tex = texSpecs.map((_, t) => (from > 0 ? e.tex[t] : 0));
    for (let k = from; k < allocs.length; ++k) {
      const a = allocs[k];
      a.vertexBase = vertex;
      a.drawBase = draw;
      a.groupBase = group;
      a.sizeBase = size;
      a.siteColorBase = siteColor;
      vertex += a.vertexCount;
      draw += a.drawCount;
      group += a.groupSpan;
      size += a.sizeTable;
      siteColor += a.siteColorTable;
      for (let t = 0; t < tex.length; ++t) {
        a.texBases[t] = tex[t];
        tex[t] += a.texels[t];
      }
    }
    e.vertex = vertex;
    e.draw = draw;
    e.group = group;
    e.size = size;
    e.siteColor = siteColor;
    e.tex = tex.slice();
    return { vertex, draw, tex };
  }

  /** an allocation's part of a vertex attribute: its member's vertices at its base */
  const allocVertexPart = (attr: VertexAttribute, a: GeometryAlloc): SegmentedArray.Part => ({
    array: members[a.member][attr.key]!.ref.value as Float32Array | Uint32Array | Int32Array,
    offset: a.vertexBase * attr.itemSize,
    count: a.vertexCount * attr.itemSize,
    bias: 0,
    version: 0,
  });
  /** an allocation's part of the elements: its member's indices at its draw base, made absolute by its vertex base */
  const allocElementsPart = (a: GeometryAlloc): SegmentedArray.Part => ({
    array: members[a.member].elements!.ref.value as Uint32Array,
    offset: a.drawBase,
    count: a.drawCount,
    bias: a.vertexBase,
    version: 0,
  });
  /** an allocation's part of the texture of `spec` (the `t`-th) */
  const allocTexPart = (spec: MergedTexSpec, t: number, a: GeometryAlloc): SegmentedArray.Part => ({
    array: (members[a.member][spec.tex]!.ref.value as DataArray<any>).array,
    offset: a.texBases[t] * spec.itemSize,
    count: a.texels[t] * spec.itemSize,
    bias: 0,
    version: 0,
  });

  /** A source's geometry laid out whole: an allocation a member, in order, the arrays with room */
  function buildSourceGeometry() {
    allocs = members.map((_, i) => createAlloc(i));
    for (let i = 0; i < n; ++i) liveAlloc[i] = i;
    replaced.clear();
    segmentRows = 'all';
    segmentAllocs = 'all';
    const end = layoutAllocs();
    allocEnd = end;
    roomVertices = sourceRoom(end.vertex);
    roomDraws = sourceRoom(end.draw);
    roomTexels = end.tex.map(sourceRoom);
    for (const attr of vertexAttributes) {
      const parts = allocs.map((a) => allocVertexPart(attr, a));
      const array = SegmentedArray.create<SegmentedArrayType>(attr.ctor, parts, roomVertices * attr.itemSize);
      ValueCell.update(own[attr.key], array);
    }
    if (hasElements) {
      ValueCell.update(own.elements, SegmentedArray.create(Uint32Array, allocs.map(allocElementsPart), roomDraws));
    }
    for (let t = 0; t < texSpecs.length; ++t) {
      const spec = texSpecs[t];
      if (!isAllocTexture(spec)) {
        mergeTextureFull(spec);
        continue;
      }
      const ctor = getTexArrayCtor(spec);
      const count = Math.max(1, roomTexels[t]);
      const length = Math.ceil((count * spec.itemSize * ctor.BYTES_PER_ELEMENT) / 4) * (4 / ctor.BYTES_PER_ELEMENT);
      const parts = allocs.map((a) => allocTexPart(spec, t, a));
      ValueCell.update(own[spec.tex], { array: SegmentedArray.create(ctor, parts, length), count });
    }
  }

  /**
   * A source's members of `replaced` allocated again at the arrays' ends, the arrays extended by their parts alone
   * (longer where the room is used up: room again, the parts where they were); laid out whole where the unused
   * allocations draw more than the used
   */
  function appendSourceGeometry() {
    const firstNew = allocs.length;
    for (const i of replaced) {
      liveAlloc[i] = allocs.length;
      allocs.push(createAlloc(i));
      if (segmentRows !== 'all') segmentRows.add(i);
      if (segmentAllocs !== 'all') segmentAllocs.add(i);
    }
    const replacedFirst = replaced.has(0);
    replaced.clear();
    const end = layoutAllocs(firstNew);
    let used = 0;
    for (let i = 0; i < n; ++i) used += allocs[liveAlloc[i]].drawCount;
    if (end.draw - used > used) {
      buildSourceGeometry();
      return;
    }
    // the room used up: more (by half again), the arrays extended longer, not laid out whole (every layer written
    // again: a source's members filled a few at a time, its sites hundreds of MB, at each growth)
    if (end.vertex > roomVertices) roomVertices = sourceRoom(end.vertex);
    if (end.draw > roomDraws) roomDraws = sourceRoom(end.draw);
    for (let t = 0; t < roomTexels.length; ++t) if (end.tex[t] > roomTexels[t]) roomTexels[t] = sourceRoom(end.tex[t]);
    allocEnd = end;
    const added = allocs.slice(firstNew);
    for (const attr of vertexAttributes) {
      const array = own[attr.key].ref.value as SegmentedArray;
      ValueCell.update(
        own[attr.key],
        SegmentedArray.extend(
          array,
          added.map((a) => allocVertexPart(attr, a)),
          Math.max(array.length, roomVertices * attr.itemSize),
        ),
      );
    }
    if (hasElements) {
      const array = own.elements.ref.value as SegmentedArray;
      ValueCell.update(
        own.elements,
        SegmentedArray.extend(array, added.map(allocElementsPart), Math.max(array.length, roomDraws)),
      );
    }
    for (let t = 0; t < texSpecs.length; ++t) {
      const spec = texSpecs[t];
      if (!isAllocTexture(spec)) {
        // the shared values' layer where it has no granularity: read again where it is another
        if (replacedFirst && !sourced.has(spec.tex)) mergeTextureFull(spec);
        continue;
      }
      const data = own[spec.tex].ref.value as DataArray<any>;
      const ctor = getTexArrayCtor(spec);
      const count = Math.max(data.count, roomTexels[t]);
      const length = Math.ceil((count * spec.itemSize * ctor.BYTES_PER_ELEMENT) / 4) * (4 / ctor.BYTES_PER_ELEMENT);
      const segmented = data.array as SegmentedArray;
      const array = SegmentedArray.extend(
        segmented,
        added.map((a) => allocTexPart(spec, t, a)),
        Math.max(segmented.length, length),
      );
      ValueCell.update(own[spec.tex], { array, count });
    }
  }

  /**
   * A source's segments: the instances by member, the geometry by the members' allocations (of those allocated again
   * alone, `segmentAllocs`: an allocation's bases fixed once laid out, the others' as they were)
   */
  function computeSourceSegments() {
    const bases = texSpecs.map((spec) => texBases.get(spec.tex)!);
    const writeAlloc = (i: number) => {
      const a = allocs[liveAlloc[i]];
      segments.vertexBases[i] = a.vertexBase;
      segments.drawBases[i] = a.drawBase;
      segments.elementOffsets[i] = a.drawBase * 4;
      segments.groupBases[i] = a.groupBase;
      sizeTableBases[i] = a.sizeBase;
      siteColorTableBases[i] = a.siteColorBase;
      for (let t = 0; t < bases.length; ++t) bases[t][i] = a.texBases[t];
      memberDraws[i] = a.drawCount;
    };
    if (segmentAllocs === 'all') for (let i = 0; i < n; ++i) writeAlloc(i);
    else for (const i of segmentAllocs) writeAlloc(i);
    segmentAllocs = new Set();

    let instanceTotal = 0;
    let drawMax = 0;
    let maxGroupCount = 0;
    let groupInstanceBase = 0;
    for (let i = 0; i < n; ++i) {
      // the counts as `updateLayout` read them (of `readList`): not the members' cells, thousands of them
      const instanceCount = layout[i * 5];
      const groupCount = layout[i * 5 + 1];
      segments.instanceBases[i] = instanceTotal;
      segments.instanceCounts[i] = instanceCount;
      segments.groupCounts[i] = groupCount;
      segments.groupInstanceOffsets[i] = groupInstanceBase - instanceTotal * groupCount;
      instanceTotal += instanceCount;
      groupInstanceBase += instanceCount * groupCount;
      if (memberDraws[i] > drawMax) drawMax = memberDraws[i];
      if (groupCount > maxGroupCount) maxGroupCount = groupCount;
    }
    segments.count = n;
    sourceTotals = { instanceTotal, vertexTotal: allocEnd.vertex, drawTotal: allocEnd.draw, drawMax, maxGroupCount };
    return sourceTotals;
  }

  /** `computeSourceSegments`' totals as last computed: what a count change alone keeps */
  let sourceTotals:
    | { instanceTotal: number; vertexTotal: number; drawTotal: number; drawMax: number; maxGroupCount: number }
    | undefined;

  /**
   * A source's instance counts changed alone (no member allocated again): the segments' instance bases, counts and
   * group offsets, from `layout`; the allocations' (bases, draws, tables) as they were (a member's allocation records
   * scattered and slow to read, thousands a change)
   */
  function computeSourceCounts() {
    if (!sourceTotals) return computeSourceSegments();
    let instanceTotal = 0;
    let groupInstanceBase = 0;
    for (let i = 0; i < n; ++i) {
      const instanceCount = layout[i * 5];
      const groupCount = layout[i * 5 + 1];
      segments.instanceBases[i] = instanceTotal;
      segments.instanceCounts[i] = instanceCount;
      segments.groupInstanceOffsets[i] = groupInstanceBase - instanceTotal * groupCount;
      instanceTotal += instanceCount;
      groupInstanceBase += instanceCount * groupCount;
    }
    sourceTotals = { ...sourceTotals, instanceTotal };
    return sourceTotals;
  }

  function computeSegments() {
    if (source) return computeSourceSegments();
    let sizeTableTotal = 0;
    let siteColorTableTotal = 0;
    for (let i = 0; i < n; ++i) {
      sizeTableBases[i] = sizeTableTotal;
      siteColorTableBases[i] = siteColorTableTotal;
      const v = members[i] as {
        sizeType?: ValueCell<string>;
        tSize?: ValueCell<DataArray<Uint8Array>>;
        tSiteColor?: ValueCell<DataArray<Uint8Array>>;
      };
      if (v.sizeType?.ref.value === 'siteType') sizeTableTotal += v.tSize!.ref.value.count;
      if (siteColorTypeOf() === 'siteType') siteColorTableTotal += v.tSiteColor!.ref.value.count;
    }
    let instanceTotal = 0;
    let vertexTotal = 0;
    let drawTotal = 0;
    let drawMax = 0;
    let maxGroupCount = 0;
    let groupInstanceBase = 0;
    let groupBase = 0;
    for (let i = 0; i < n; ++i) {
      const v = members[i];
      const instanceCount = memberInstanceCount(i);
      const groupCount = v.uGroupCount.ref.value;
      segments.instanceBases[i] = instanceTotal;
      segments.instanceCounts[i] = instanceCount;
      segments.vertexBases[i] = vertexTotal;
      segments.drawBases[i] = drawTotal;
      segments.elementOffsets[i] = drawTotal * 4;
      segments.groupCounts[i] = groupCount;
      segments.groupInstanceOffsets[i] = groupInstanceBase - instanceTotal * groupCount;
      segments.groupBases[i] = groupBase;
      instanceTotal += instanceCount;
      vertexTotal += v.uVertexCount.ref.value;
      drawTotal += v.drawCount.ref.value;
      if (v.drawCount.ref.value > drawMax) drawMax = v.drawCount.ref.value;
      groupInstanceBase += instanceCount * groupCount;
      groupBase += getGroupSpan(i);
      if (groupCount > maxGroupCount) maxGroupCount = groupCount;
    }
    segments.count = n;
    return { instanceTotal, vertexTotal, drawTotal, drawMax, maxGroupCount };
  }

  /** `aSegment` of `instanceTotal` instances: per instance its segment's groups and its segment's index */
  function getSegmentInstanceArray(instanceTotal: number) {
    const out = new Int32Array(instanceTotal * 4);
    for (let i = 0; i < n; ++i) {
      const base = segments.instanceBases[i];
      const g0 = segments.groupCounts[i];
      const g1 = segments.groupInstanceOffsets[i];
      const g2 = segments.groupBases[i];
      for (let k = 0, kl = segments.instanceCounts[i]; k < kl; ++k) {
        const s = (base + k) * 4;
        out[s] = g0;
        out[s + 1] = g1;
        out[s + 2] = g2;
        out[s + 3] = i;
      }
    }
    return out;
  }

  function rebuildInstances(instanceTotal: number) {
    if (source) return rebuildSourceInstances(instanceTotal);
    const segmentData = own.tSegment.ref.value as DataArray<Float32Array>;
    const segmentArray = segmentData.array;

    // instanced attributes in per-segment blocks; the transforms the members' as parts (`SegmentedArray`), not copied
    // (millions of instances: their rows packed where read, `packAffineRowsOf`)
    const transformParts: SegmentedArray.Part<Float32Array>[] = [];
    // each member's by its encoding, kept where one of them has poses (its `aTransform` empty)
    const poseParts: TransformPart[] = [];
    let hasPoses = false;
    // every member's instances in their own order (`getSharedSerial`): the merged ones too, a view as well, no copy;
    // else an array of its own (never a shared view, read only)
    const allSerial = members.every((v) => isSharedSerial(v.aInstance.ref.value as Uint32Array));
    const ownInstance = own.aInstance.ref.value as Uint32Array;
    const instanceArray = allSerial
      ? getSharedSerial(instanceTotal)
      : isSharedSerial(ownInstance)
        ? new Uint32Array(instanceTotal)
        : ensureArray(ownInstance, instanceTotal, Uint32Array);
    for (let i = 0; i < n; ++i) {
      const v = members[i];
      const instanceCount = segments.instanceCounts[i];
      const instanceBase = segments.instanceBases[i];
      const memberTransforms = getObjectTransforms(v as unknown as TransformValues);
      if (memberTransforms.kind === 'poses') hasPoses = true;
      poseParts.push({ transforms: memberTransforms, first: instanceBase });
      transformParts.push({
        array: v.aTransform.ref.value as Float32Array,
        offset: instanceBase * 16,
        // a member with poses: its `aTransform` empty, its instances in `poseParts`
        count: memberTransforms.kind === 'poses' ? 0 : instanceCount * 16,
        bias: 0,
        version: 0,
      });
      const memberInstance = v.aInstance.ref.value as Uint32Array;
      const ibs = v.invariantBoundingSphere.ref.value as Sphere3D;
      // the segment's data: invariant bounding sphere, level scales and order
      const o = i * SegmentTexels * 4;
      segmentArray[o + 0] = ibs.center[0];
      segmentArray[o + 1] = ibs.center[1];
      segmentArray[o + 2] = ibs.center[2];
      segmentArray[o + 3] = ibs.radius;
      writeSegmentLod(i, segmentArray);
      if (!allSerial) {
        for (let k = 0; k < instanceCount; ++k) {
          instanceArray[instanceBase + k] = memberInstance[k] + instanceBase;
        }
      }
    }
    ValueCell.update(
      own.instancePoses,
      hasPoses ? { kind: 'segmented', count: instanceTotal, parts: poseParts } : undefined,
    );
    ValueCell.update(own.aTransform, SegmentedArray.create(Float32Array, transformParts, instanceTotal * 16));
    ValueCell.update(own.aInstance, instanceArray);
    // cells (`hasInstanceCells`) keep the renderer on the cull path, each member's OWN instanceGrid is used instead
    ValueCell.update(own.instanceGrid, { ...createEmptyInstanceGrid(), cellSize: 2 });
    ValueCell.updateIfChanged(own.instanceCellSize, 2);
    // made where read (16 B an instance, the same for every instance of a segment): the instanced attribute of objects
    // culled on the CPU, WebGL2's instance data; GPU culling reads the segments' groups (`CulledData`), not this
    let segmentInstanceArray: Int32Array | undefined;
    const box = own.aSegment.ref;
    ValueCell.set(own.aSegment, {
      id: box.id,
      version: box.version + 1,
      metadata: box.metadata,
      get value() {
        return (segmentInstanceArray ??= getSegmentInstanceArray(instanceTotal));
      },
    });
    ValueCell.update(own.tSegment, segmentData);
  }

  /**
   * A source's instances (`MergedInstanceSource`): the segment table alone (a member's invariant sphere, levels and
   * order); no transforms, poses nor `aInstance` of the members (the source's data: the culler's, written by range), an
   * instance's `aInstance` its index (the source's layers' texel)
   */
  function rebuildSourceInstances(instanceTotal: number) {
    const segmentData = own.tSegment.ref.value as DataArray<Float32Array>;
    const segmentArray = segmentData.array;
    // the rows of what changed (`segmentRows`): an instance count changed alone changes none
    const writeRow = (i: number) => {
      const ibs = members[i].invariantBoundingSphere.ref.value as Sphere3D;
      const o = i * SegmentTexels * 4;
      segmentArray[o + 0] = ibs.center[0];
      segmentArray[o + 1] = ibs.center[1];
      segmentArray[o + 2] = ibs.center[2];
      segmentArray[o + 3] = ibs.radius;
      writeSegmentLod(i, segmentArray);
    };
    if (segmentRows === 'all') for (let i = 0; i < n; ++i) writeRow(i);
    else for (const i of segmentRows) writeRow(i);
    segmentRows = new Set();
    ValueCell.update(own.instancePoses, undefined);
    ValueCell.update(own.aTransform, SegmentedArray.create(Float32Array, [], 0));
    // a page store's: its instances' ids the culled poses' (`culledInstance`), no array of them (the shared serial of the
    // store's instances: 4 B each, a list with templates' every copy's)
    ValueCell.update(own.aInstance, source?.paged ? EmptyInstances : getSharedSerial(instanceTotal));
    ValueCell.update(own.instanceGrid, { ...createEmptyInstanceGrid(), cellSize: 2 });
    ValueCell.updateIfChanged(own.instanceCellSize, 2);
    let segmentInstanceArray: Int32Array | undefined;
    const box = own.aSegment.ref;
    ValueCell.set(own.aSegment, {
      id: box.id,
      version: box.version + 1,
      metadata: box.metadata,
      get value() {
        return (segmentInstanceArray ??= getSegmentInstanceArray(instanceTotal));
      },
    });
    ValueCell.update(own.tSegment, segmentData);
  }

  /** A source's instance counts changed alone: the segments, the segment table and the counts again, not the geometry */
  function relayoutSourceInstances() {
    // the members replaced since: allocated again
    const geometryChanged = replaced.size > 0;
    if (geometryChanged) appendSourceGeometry();
    const { instanceTotal, vertexTotal, drawTotal, drawMax, maxGroupCount } = geometryChanged
      ? computeSegments()
      : computeSourceCounts();
    rebuildSourceInstances(instanceTotal);
    ValueCell.updateIfChanged(own.instanceCount, instanceTotal);
    ValueCell.updateIfChanged(own.uInstanceCount, instanceTotal);
    ValueCell.updateIfChanged(own.drawCount, drawTotal);
    ValueCell.updateIfChanged(own.drawCountMax, drawMax);
    ValueCell.updateIfChanged(own.uVertexCount, vertexTotal);
    ValueCell.updateIfChanged(own.uGroupCount, maxGroupCount);
    updateBoundingSpheres();
    updateMarkerAggregates();
    // their averages by the members' (an instance count changed alone moves their weights alone: kept till a member's
    // geometry or a whole layout)
    if (geometryChanged) updateAggregates();
  }

  function rebuild() {
    if (source) buildSourceGeometry();
    const { instanceTotal, vertexTotal, drawTotal, drawMax, maxGroupCount } = computeSegments();

    rebuildInstances(instanceTotal);

    // a source's arrays by allocation (`buildSourceGeometry`)
    if (!source) rebuildGeometry();

    ValueCell.updateIfChanged(own.drawCount, drawTotal);
    ValueCell.updateIfChanged(own.drawCountMax, drawMax);
    ValueCell.updateIfChanged(own.instanceCount, instanceTotal);
    ValueCell.updateIfChanged(own.uVertexCount, vertexTotal);
    ValueCell.updateIfChanged(own.uInstanceCount, instanceTotal);
    ValueCell.updateIfChanged(own.uGroupCount, maxGroupCount);

    updateBoundingSpheres();
    updateMarkerAggregates();
    updateAggregates();
  }

  /** the members' vertex attributes, elements and layers as parts, by member */
  function rebuildGeometry() {
    for (const attr of vertexAttributes) {
      const parts = members.map((_, i) => vertexPart(attr, i));
      ValueCell.update(own[attr.key], SegmentedArray.create<SegmentedArrayType>(attr.ctor, parts));
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
  }

  /** the versions of the members read (`readList`) as seen, their changes counted (`memberChanges`) */
  function snapshotAll() {
    for (let k = 0; k < readList.length; ++k) {
      const i = readList[k];
      memberChanges[i] += 1;
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
      mv.sites = getSitesVersion(c);
      for (let a = 0, o = i * attrCount; a < attrCount; ++a)
        attributeVersions[o + a] = attributeCells[o + a].ref.version;
      mv.elements = c.elements ? c.elements.ref.version : -1;
      for (let t = 0, o = i * texCount; t < texCount; ++t) textureVersions[o + t] = textureCells[o + t].ref.version;
    }
  }

  const cullVersion = ValueCell.create(0);
  const markerVersion = ValueCell.create(0);

  /** the owner's `recordsVersion` at the last sync */
  let syncedRecordsVersion = -1;
  function recordsUnchanged() {
    return synced && !source && recordsVersion!.ref.value === syncedRecordsVersion;
  }

  function sync(observedOnly = false) {
    const version = recordsVersion?.ref.value;
    // with a source: the members it changed alone after the first sync, the others' cells as the last saw them
    const taken = source?.takeMembers?.();
    // and those it changed the count or bounds of alone: from its copies
    const counted = counts?.take();
    countList = synced && counted ? Int32Array.from(counted) : [];
    if (synced && taken) {
      // and those of another geometry (`replaceGeometry`)
      const read = new Set(taken);
      for (const i of replaced) read.add(i);
      readList = Int32Array.from(read);
    } else if (synced && observedOnly && !source) {
      readList = getObservedChanged();
    } else {
      readList = allMembers;
    }
    synced = true;
    // the records' writes seen by a read of all members alone: one of the observed cells (a mark between a write and
    // its update) leaves them to the next sync
    if (version !== undefined && readList === allMembers) syncedRecordsVersion = version;
    if (updateLayout()) {
      // a whole layout: all of them seen, and counted as changed
      readList = allMembers;
      if (version !== undefined) syncedRecordsVersion = version;
      rebuild();
      snapshotAll();
      ValueCell.update(cullVersion, cullVersion.ref.value + 1);
      ValueCell.update(markerVersion, markerVersion.ref.value + 1);
      return true;
    }
    if (instanceCountsChanged || replaced.size > 0) {
      relayoutSourceInstances();
      updateMemberTables(readList);
      snapshotAll();
      ValueCell.update(cullVersion, cullVersion.ref.value + 1);
      ValueCell.update(markerVersion, markerVersion.ref.value + 1);
      return true;
    }

    let instancesChanged = false;
    // the bounds of those the source changed the count or bounds of alone (`counts`)
    let boundsChanged = readCountSpheres();
    let markersChanged = false;
    let aggregatesChanged = false;
    let elementsChanged = false;
    let lodChanged = false;
    const attributesChanged = new Set<string>();

    for (let k = 0; k < readList.length; ++k) {
      const i = readList[k];
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
      for (let a = 0, o = i * attrCount; a < attrCount; ++a) {
        const version = attributeCells[o + a].ref.version;
        if (attributeVersions[o + a] !== version) {
          updatePart(vertexAttributes[a].key, i);
          attributesChanged.add(vertexAttributes[a].key);
          attributeVersions[o + a] = version;
        }
      }
      if (c.elements && mv.elements !== c.elements.ref.version) {
        updatePart('elements', i);
        elementsChanged = true;
        mv.elements = c.elements.ref.version;
      }

      const lodVersion = c.lodLevels ? c.lodLevels.ref.version : -1;
      const sitesVersion = getSitesVersion(c);
      if (mv.lodLevels !== lodVersion || mv.sites !== sitesVersion) {
        writeSegmentLod(i, (own.tSegment.ref.value as DataArray<Float32Array>).array);
        lodChanged = true;
      }
      // what the culler reads of it (`Merged.memberChanges`)
      if (
        mv.aTransform !== transform ||
        mv.aInstance !== instance ||
        mv.boundingSphere !== bounds ||
        mv.invariantBoundingSphere !== invariantBounds ||
        mv.lodLevels !== lodVersion ||
        mv.sites !== sitesVersion
      ) {
        memberChanges[i] += 1;
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
      mv.sites = sitesVersion;
    }

    for (let t = 0; t < texSpecs.length; ++t) {
      const spec = texSpecs[t];
      // the source's layers its own cells: updated (and their ranges logged) by it
      if (sourced.has(spec.tex)) continue;
      let texChanged = false;
      if (spec.granularity(shared) === 'none') {
        if (textureVersions[t] !== textureCells[t].ref.version) {
          // an array of another size is a new layout (`mergeTextureFull`)
          const src = shared[spec.tex].ref.value as DataArray<any>;
          if (src.count === (own[spec.tex].ref.value as DataArray<any>).count) {
            updateTexPart(spec, 0);
            ValueCell.update(own[spec.tex], own[spec.tex].ref.value);
          } else {
            mergeTextureFull(spec);
          }
        }
        // only the shared values' is read while the granularity stays (a change of it rebuilds)
        textureVersions[t] = textureCells[t].ref.version;
        continue;
      }
      const marker = spec.tex === 'tMarker';
      for (let k = 0; k < readList.length; ++k) {
        const i = readList[k];
        const c = memberCells[i];
        const mv = memberVersions[i];
        const version = textureCells[i * texCount + t].ref.version;
        // uniform-status marking updates the array without bumping the texture cell
        if (
          textureVersions[i * texCount + t] !== version ||
          (marker && mv.markerStatus !== c.markerStatus.ref.version)
        ) {
          updateTexPart(spec, i);
          texChanged = true;
        }
        textureVersions[i * texCount + t] = version;
      }
      if (texChanged) {
        ValueCell.update(own[spec.tex], own[spec.tex].ref.value);
        if (marker) markersChanged = true;
      }
    }

    // a source's members have no instance data of their own: their transforms' versions not read
    if (instancesChanged && !source) rebuildInstances(own.instanceCount.ref.value as number);
    for (const key of attributesChanged) ValueCell.update(own[key], own[key].ref.value);
    if (elementsChanged) ValueCell.update(own.elements, own.elements.ref.value);
    if (!instancesChanged && lodChanged) ValueCell.update(own.tSegment, own.tSegment.ref.value);
    if (boundsChanged) updateBoundingSpheres();
    if (markersChanged) {
      updateMarkerAggregates(readList);
      ValueCell.update(markerVersion, markerVersion.ref.value + 1);
    }
    if (aggregatesChanged) updateAggregates();
    if (instancesChanged || boundsChanged || lodChanged) ValueCell.update(cullVersion, cullVersion.ref.value + 1);

    for (let k = 0; k < readList.length; ++k) {
      const i = readList[k];
      memberVersions[i].markerStatus = memberCells[i].markerStatus.ref.version;
    }
    return false;
  }

  updateLayout();
  rebuild();
  snapshotAll();

  // the members' marker cells only: a mark changes members without an update of the merged object, their other changes
  // come with one; two cells a member, as there may be 100k members
  const observed: ValueCell<unknown>[] = [];
  for (let i = 0; i < n; ++i) {
    observed.push(memberCells[i].markerStatus);
    if (markerTexture >= 0) observed.push(textureCells[i * texCount + markerTexture]);
  }

  /** a source's member of another geometry: allocated again at the next sync */
  function replaceGeometry(i: number) {
    if (!source) throw new Error('merged values: a geometry replaced without a source');
    replaced.add(i);
  }

  return {
    type,
    values,
    records: members,
    segments,
    sync,
    memberCells: observed,
    getMarkedMembers: source?.markers ? undefined : getMarkedMembers,
    recordsUnchanged: recordsVersion ? recordsUnchanged : undefined,
    cullVersion,
    markerVersion,
    source,
    memberChanges,
    replaceGeometry,
  };
}

//

/**
 * Merged objects culled on the GPU (`dCulled`, WebGPU, `merged-culler.ts`): drawn from indirect records whose first
 * instances index the culler's regions of instance indices (`sCulledIndices`), the instances' data read from storage by
 * those (`sCulledData`, `CulledData`) in place of the instanced attributes, which can't follow an index (the program's
 * vertex inputs filled by a generated reader, `webgpu/culled.ts`)
 */
export const CulledSchema = {
  sCulledIndices: StorageBufferSpec('u32', 'read', ['dCulled']),
  sCulledData: StorageBufferSpec('vec4<f32>', 'read', ['dCulled']),
  dCulled: DefineSpec('boolean'),
  /**
   * a page store's instances (`paged-cull.ts`): `sCulledIndices` the regions' pose slots (an instance's at
   * `instance_index`), `sCulledData` the segments and the composed poses (a slot each, its store index its `aInstance`)
   */
  dCulledPaged: DefineSpec('boolean'),
  /**
   * a page store's copies' markers (`MergedInstanceSource.paged.copyMarkers`): two words a copy, its marker or its
   * block's offset, then the blocks' bytes; a copy's instance its marker by them, not by its store index's
   */
  sCopyMarkers: StorageBufferSpec('u32', 'read', ['dCopyMarkers']),
  dCopyMarkers: DefineSpec('boolean'),
  /**
   * a page store's copies' own layers (`MergedInstanceSource.paged.copyLayers`): two words a copy's kind (overpaint,
   * transparency, clipping, emissive, substance, wiggle), a value's or a block's by store index, then their bytes; a copy's
   * instance its layer by them where its copy has one, not by its store index's
   */
  sCopyLayers: StorageBufferSpec('u32', 'read', ['dCopyLayers']),
  dCopyLayers: DefineSpec('boolean'),
} as const;
export type CulledValues = Values<typeof CulledSchema>;

/** The instanced attributes culled programs read from `sCulledData` (`CulledData`), not as vertex inputs */
export const CulledAttributes = ['aTransform', 'aInstance', 'aSegment'] as const;

const MergedRenderableSchemas = new Map<MergedType, RenderableSchema>();

/** The segments of `segments` and their groups (`aSegment.xyz`) as the instance data holds them (`CulledData`) */
function getStorageSegments(segments: MergedSegments, tSegment: DataArray<Float32Array>) {
  return { segments, groups: getCulledSegmentData(segments, tSegment) };
}

/** The schema of merged renderables of `type` (one object per type, shared by the items) */
export function getMergedRenderableSchema(type: MergedType) {
  let schema = MergedRenderableSchemas.get(type);
  if (!schema) {
    const merged: { [k: string]: RenderableSchema[string] } = {
      ...GlobalUniformSchema,
      ...GlobalTextureSchema,
      ...GlobalDefineSchema,
      ...InternalSchema,
      ...MergedTypeDescriptors[type].schema,
      ...SegmentSchema,
      ...CulledSchema,
    };
    // the instanced attributes: vertex inputs but in culled programs (merged objects are always segmented) and those of
    // poses read from storage (`instance-storage.ts`)
    schema = withInstanceStorage(merged, CulledAttributes, 'dCulled');
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
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<MergeableValues & SegmentValues> {
  const descriptor = MergedTypeDescriptors[merged.type];
  const schema = getMergedRenderableSchema(merged.type);
  const members = merged.records;
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
  // the largest segment's items: draws of items address theirs by their first vertex (`ItemIndexPattern`), the shared
  // indices then cover one segment's, not the merged geometry's
  const indexed =
    descriptor.hasElements ||
    (isItemIndexing(items) &&
      merged.values.drawCountMax.ref.value / items.pattern.length <= getMaxItemCount(ctx.resources.limits, items));
  // culled on the GPU where asked for and available (`merged-culler.ts`)
  const gpuCulling = !!(merged.values as { gpuCulling?: ValueCell<boolean> }).gpuCulling?.ref.value;
  const itemPattern = !descriptor.hasElements && isItemIndexing(items) ? items : undefined;
  const culler =
    gpuCulling && canCullOnGpu(ctx)
      ? merged.source?.paged
        ? createPagedMergedCuller(ctx.compute!, ctx.resources, merged, members, indexed, itemPattern)
        : createMergedCuller(ctx.compute!, ctx.resources, merged, members, indexed, itemPattern)
      : undefined;
  // members with a pending grid (`InstanceGrid.pending`) drawn as stored where nothing reads them by index: their grids
  // built and applied (the data in their order), as they are for objects culled on the CPU
  const applyPendingGrids = () => {
    let applied = false;
    for (const m of members) applied = applyPendingGrid(m as unknown as InstanceGridValues) || applied;
    return applied;
  };
  if (!culler && applyPendingGrids()) merged.sync();
  // drawn through the culler's instance indices (`CulledSchema`)
  if (culler) Object.assign(renderValues, culler.values);
  // else poses read by index from storage where the context does (`instance-storage.ts`), not packed as rows
  const storage =
    !culler && usesInstanceStorage(ctx, merged.values)
      ? createInstanceStorage(ctx.compute!, merged.values, () =>
          getStorageSegments(merged.segments, merged.values.tSegment.ref.value as DataArray<Float32Array>),
        )
      : undefined;
  Object.assign(renderValues, storage?.values ?? { dInstanceStorage: ValueCell.create(false) });
  const rasterOnly =
    !!culler && !!(merged.values as { gpuCullingRasterOnly?: ValueCell<boolean> }).gpuCullingRasterOnly?.ref.value;
  // the members' spheres as `resolveGpuCulling` counts them, again where they changed (`cullVersion`)
  let rasterSpheres = 0;
  let rasterSpheresVersion = -1;
  function getRasterSpheres() {
    const version = merged.cullVersion.ref.value;
    if (version !== rasterSpheresVersion) {
      rasterSpheresVersion = version;
      rasterSpheres = 0;
      // the members' counts as the segments hold them (a paged source's from its copies: its members' cells not written)
      const { instanceCounts } = merged.segments;
      for (let i = 0; i < members.length; ++i) {
        const v = members[i];
        rasterSpheres += instanceCounts[i] * ((v as { uVertexCount?: ValueCell<number> }).uVertexCount?.ref.value ?? 0);
      }
    }
    return rasterSpheres;
  }
  const renderItem = ctx.createGraphicsRenderItem(
    descriptor.drawMode,
    descriptor.shaderCode,
    schema,
    renderValues,
    materialId,
    transparency,
  );

  const clipInfo = createClipInfoCache(renderValues);
  const segmentSpheres: Sphere3D[] = [];
  /**
   * the lists of the visible passes and of the ghost pass (`clipGhost`): each with its own levels' uniforms and caches,
   * so that neither rewrites the other's every frame; drawn from the set of the last cull. The ghost's culled on the
   * GPU where its region was dispatched (`culler.ghostDispatched`), on the CPU otherwise
   */
  const createSet = () => ({
    mdb: createSegmentedMdbList(),
    mode: 'none' as 'none' | 'full' | 'cull' | 'gpu',
    lastCullFrame: undefined as Frame | undefined,
    cullView: createCullView(),
    drawVersion: 0,
    /** the CPU culling's lists, shifted to the identity region of the compacted attributes where culled on the GPU */
    shiftedLists: [] as MultiDrawBaseData[],
    shiftedVersion: -1,
    /** the identity region they were shifted to (`culler.identityBase`: moves with the ghost's and raster regions) */
    shiftedBase: -1,
  });
  const visibleSet = createSet();
  const ghostSet = createSet();
  let set = visibleSet;
  const getSet = (clipGhost: boolean | undefined) => (clipGhost ? ghostSet : visibleSet);

  /** the frame culled on the GPU (its first view) and the planes it was culled with */
  let gpuFrame: Frame | undefined;
  /** told once (debug mode) that a view not culled on the GPU for it drew none of it (`prepareCullView` not called) */
  let warnedUnculledView = false;
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
  /**
   * whether the planes culled on the GPU are of a camera at `cameraPlane`: a page store's view (its buffers that view's,
   * whatever the frame; a jittered sample of temporal multi-sampling the same)
   */
  function isGpuCamera(cameraPlane: Plane3D) {
    const { normal, constant } = cameraPlane;
    return (
      gpuPlanes[24] === normal[0] &&
      gpuPlanes[25] === normal[1] &&
      gpuPlanes[26] === normal[2] &&
      gpuPlanes[27] === constant
    );
  }
  function setGpuView(cameraPlane: Plane3D, frustum: Frustum3D) {
    for (let i = 0; i < 7; ++i) {
      const { normal, constant } = getPlane(i, cameraPlane, frustum);
      gpuPlanes.set([normal[0], normal[1], normal[2], constant], i * 4);
    }
  }

  function getCpuLists() {
    if (!culler) return set.mdb.list;
    if (
      set.shiftedVersion === set.drawVersion &&
      set.shiftedBase === culler.identityBase &&
      set.shiftedLists.length === set.mdb.list.length
    ) {
      return set.shiftedLists;
    }
    set.shiftedVersion = set.drawVersion;
    set.shiftedBase = culler.identityBase;
    return shiftLists(set.mdb.list, set.shiftedLists);
  }
  /**
   * per list its shifted one (`shiftLists`): the same object each call, its fields the list's now; a level's draw
   * block is made once a level object (WebGPU's), a new one each call would leave one behind each call
   */
  const shiftedListOf = new WeakMap<MultiDrawBaseData, MultiDrawBaseData>();
  /** `lists` shifted to the identity region of the compacted attributes (`culler.identityBase`) into `out` */
  function shiftLists(lists: MultiDrawBaseData[], out: MultiDrawBaseData[]) {
    const base = culler!.identityBase;
    out.length = lists.length;
    for (let j = 0, jl = lists.length; j < jl; ++j) {
      const l = lists[j];
      let shifted = shiftedListOf.get(l);
      const baseInstances =
        shifted && shifted.baseInstances.length >= l.count
          ? shifted.baseInstances
          : new Uint32Array(l.baseInstances.length);
      for (let i = 0; i < l.count; ++i) baseInstances[i] = l.baseInstances[i] + base;
      if (shifted) {
        Object.assign(shifted, l);
        shifted.baseInstances = baseInstances;
      } else {
        shiftedListOf.set(l, (shifted = { ...l, baseInstances }));
      }
      out[j] = shifted;
    }
    return out;
  }
  /** whether the sphere raster draws the object's raster regions in the pass being drawn (`Renderable.rasterDrawn`) */
  let rasterDrawn = false;
  // the impostors draw the raster regions' instances too where the raster does not draw them
  const getGpuLists = () =>
    set === ghostSet ? culler!.getGhostLists(set.mdb.list) : culler!.getLists(set.mdb.list, !rasterDrawn);
  const getLists = () => (set.mode === 'gpu' ? getGpuLists() : getCpuLists());

  /**
   * The members culled on the GPU with a pending grid (`InstanceGrid.pending`): their data
   * in their own order, drawn by the culler's indices; their grids built where a CPU cull needs their cells
   * (`ensureIndexes`) into cells of their own, their orders the identity region's (`MergedCuller.setIdentityOrder`)
   */
  const pending = members.map((v) => !!culler && !!v.instanceGrid.ref.value.pending);
  /** a source's instances (`MergedInstanceSource.getMember`): the identity region its indices, by member */
  const hasPending = pending.some((p) => p);
  const indexGrids = members.map((v, i) => (pending[i] ? ValueCell.create(v.instanceGrid.ref.value) : v.instanceGrid));
  /**
   * the order of a member's built grid, the draw slots its marked lists read (`createInstancePositions`): made for the
   * pending members as their grids are built (not a cell a member: 100k members)
   */
  const indexOrders: (ValueCell<Uint32Array> | undefined)[] = new Array(members.length);
  /** `cullVersion` of the built grids, -1 for none built */
  let indexVersion = -1;
  const isIndexCurrent = () => indexVersion === merged.cullVersion.ref.value;
  function ensureIndexes() {
    if (!hasPending) return;
    if (indexVersion >= 0 && isIndexCurrent()) return;
    indexVersion = merged.cullVersion.ref.value;
    // the members' counts as the segments hold them (a paged source's from its copies: its members' cells not written)
    const { instanceBases, instanceCounts } = merged.segments;
    const order = new Uint32Array(merged.values.instanceCount.ref.value as number);
    for (let i = 0, il = members.length; i < il; ++i) {
      const base = instanceBases[i];
      const n = instanceCounts[i];
      const grid = pending[i] ? members[i].instanceGrid.ref.value : undefined;
      if (grid?.pending) {
        const built = buildPendingGrid(members[i] as unknown as InstanceGridValues);
        ValueCell.update(indexGrids[i], built);
        if (indexOrders[i]) ValueCell.update(indexOrders[i]!, built.cellInstance);
        else indexOrders[i] = ValueCell.create(built.cellInstance);
        const cells = built.cellInstance;
        for (let k = 0; k < n; ++k) order[base + k] = base + cells[k];
      } else {
        for (let k = 0; k < n; ++k) order[base + k] = base + k;
      }
    }
    culler!.setIdentityOrder(order);
  }
  /** the built grids dropped where the members changed: the identity region in the instances' order again */
  function syncIndexes() {
    if (indexVersion < 0 || isIndexCurrent()) return;
    indexVersion = -1;
    culler!.setIdentityOrder(undefined);
    for (let i = 0, il = members.length; i < il; ++i) {
      if (pending[i]) ValueCell.update(indexGrids[i], members[i].instanceGrid.ref.value);
    }
  }
  const isIndexBuilt = (i: number) => pending[i] && indexVersion >= 0;

  /**
   * a paged source's members' counts as the CPU cull reads them: cells of their own, set from the segments' counts as a
   * CPU cull begins (`syncCullCounts`), the members' cells not written (a step: thousands of them)
   */
  const pagedCountCells = merged.source?.paged ? members.map(() => ValueCell.create(0)) : undefined;
  /** the merged object's cull version the counts were synced at: a count changed advances it */
  let cullCountsVersion = -1;
  function syncCullCounts() {
    if (!pagedCountCells || cullCountsVersion === merged.cullVersion.ref.value) return;
    cullCountsVersion = merged.cullVersion.ref.value;
    const { instanceCounts } = merged.segments;
    for (let i = 0, il = pagedCountCells.length; i < il; ++i) {
      ValueCell.updateIfChanged(pagedCountCells[i], instanceCounts[i] ?? 0);
    }
  }
  // the members' cells culling reads every frame, resolved once: values objects have too many keys for fast named reads.
  // The level of detail and the clip the object's (its shared values', of the merge key), not a member's
  const { uLod, dClipObjects, dClipping, clipObjects, dClipVariant } = merged.values;
  const memberCull: CullValues[] = members.map((v, i) => ({
    drawCount: v.drawCount,
    instanceCount: pagedCountCells ? pagedCountCells[i] : v.instanceCount,
    uLod,
    instanceGrid: indexGrids[i],
    instanceCellSize: v.instanceCellSize,
    boundingSphere: v.boundingSphere,
    lodLevels: v.lodLevels,
    dClipObjects,
    dClipping,
    clipObjects,
    dClipVariant,
    invariantBoundingSphere: v.invariantBoundingSphere,
    // the object's transient transform, not the member's: its grid moved by it (`getCullGrid`)
    uObjectTransform: merged.values.uObjectTransform,
  }));

  /**
   * what the last CPU cull read: a still frame keeps its lists (`createCullView`); the members' cells by the version
   * `sync` advances where they changed, not each (with 100k members, recording them was most of a moving frame); the
   * object's transient transform, which moves every member's cells; its clip objects (the shared values', a clip moved
   * alone syncs as a mark: neither of the others moves, nor are the modes reset where a culler culls on the GPU)
   */
  const cullCells = [
    merged.cullVersion,
    merged.values.uObjectTransform,
    (merged.values as { clipObjects?: ValueCell<unknown> }).clipObjects,
  ];

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
  /** the view of the last cull for the marked lists (`renderMarked`), undefined where not culled by it */
  let markedLod: LodView | undefined;
  const markedPlane = Plane3D();
  const markedFrustum = Frustum3D();
  function setMarkedView(cameraPlane: Plane3D, frustum: Frustum3D, lod: LodView) {
    markedLod = lod;
    Plane3D.copy(markedPlane, cameraPlane);
    Frustum3D.copy(markedFrustum, frustum);
  }

  function buildFullList(lod: LodView) {
    syncCullCounts();
    if (set === visibleSet) {
      segmentLod = lod;
      markedLod = undefined;
    }
    const hasLodLevels = !!set.mdb.prepare(memberCull[0].lodLevels, getCapacity(), true, lod);
    for (let i = 0, il = members.length; i < il; ++i) {
      set.mdb.fullSegment(memberCull[i], setSegment(i), getMemberLod(i, hasLodLevels));
    }
    set.mode = 'full';
  }

  function cull(
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    isOccluded: ((s: Sphere3D) => boolean) | null,
    stats: GraphicsStats,
    frame: Frame,
    lod: LodView,
    clipGhost = false,
  ) {
    syncCullCounts();
    set = getSet(clipGhost);
    if (!clipGhost) {
      segmentLod = lod;
      setMarkedView(cameraPlane, frustum, lod);
    }
    // skip recomputation if this renderable was already culled for the current frame
    if (frame === set.lastCullFrame) {
      stats.cacheHits.cull++;
      return;
    }

    // a view of its own culled on the GPU (`prepareCullView`: the thickness pass) replaces the frame's culling: the
    // frame's passes after it draw that (without the occlusion a superset) where of the same view; a page store's where
    // of the same camera (whatever the frame, a jittered sample too), the others' of the same planes
    const gpuView =
      frame === gpuFrame
        ? isGpuView(cameraPlane, frustum)
        : gpuFrame !== undefined && (culler?.paged ? isGpuCamera(cameraPlane) : isGpuView(cameraPlane, frustum));
    if (culler && gpuView && (!clipGhost || culler.ghostDispatched)) {
      // culled on the GPU before the passes (`prepareCull`): the level lists for their uniforms only
      set.mdb.prepare(memberCull[0].lodLevels, 0, true, lod);
      gpuLists = clipGhost ? culler.getGhostLists(set.mdb.list) : culler.getLists(set.mdb.list);
      set.mode = 'gpu';
      set.lastCullFrame = frame;
      set.drawVersion = getDrawVersion();
      return;
    }

    // culled on the GPU alone (a page store's: nothing of its instances on the CPU; the others but for the sphere raster
    // alone, `rasterOnly`): a view not culled on the GPU for it draws none, as while the kernels compile (`'none'`; no
    // CPU lists, no grids built for them; `prepareCullView` culls a view of its own)
    if (culler && (culler.paged || !rasterOnly)) {
      // its kernels ready: a pass of a view of its own (a shadow map, a custom camera) that did not cull it on the GPU
      // (`Renderer.prepareCullView`), no CPU fallback: said once in debug mode, the cause hard to see otherwise
      if (isDebugMode && !warnedUnculledView && gpuFrame !== undefined && culler.areKernelsReady()) {
        warnedUnculledView = true;
        console.warn(
          'merged object culled on the GPU: drawn in a view not culled for it (Renderer.prepareCullView), nothing drawn',
        );
      }
      set.mdb.prepare(memberCull[0].lodLevels, 0, true, lod);
      set.mode = 'none';
      set.lastCullFrame = frame;
      set.drawVersion = getDrawVersion();
      return;
    }

    // for the sphere raster alone: members' grids pending while the culling's kernels compile (an app's first frames):
    // nothing drawn this frame, not every grid built for it; the kernels are pending work, the frames after draw (mode
    // 'none': no CPU lists kept as a still frame's)
    if (hasPending && indexVersion < 0 && !culler!.areKernelsReady()) {
      set.mdb.prepare(memberCull[0].lodLevels, 0, true, lod);
      set.mode = 'none';
      set.lastCullFrame = frame;
      set.drawVersion = getDrawVersion();
      return;
    }

    // the same view and members' data as the last CPU cull (a still frame): its lists stand
    if (set.mode === 'cull' && set.cullView.isSame(lod, isOccluded, cullCells)) {
      set.lastCullFrame = frame;
      stats.cacheHits.cull++;
      return;
    }
    set.cullView.set(lod, cullCells);
    // the members' pending grids built: their cells for the lists, their order the identity region's
    ensureIndexes();
    const lodLevels = set.mdb.prepare(memberCull[0].lodLevels, getCapacity(), true, lod);
    set.mode = 'cull';

    for (let i = 0, il = members.length; i < il; ++i) {
      set.mdb.cullSegment(
        memberCull[i],
        setSegment(i),
        getMemberLod(i, !!lodLevels),
        cameraPlane,
        frustum,
        isOccluded,
        stats,
        true,
        clipGhost,
      );
    }

    set.lastCullFrame = frame;
    set.drawVersion = getDrawVersion();
  }

  // the marked instances alone (`renderMarked`), at the levels of the view of the last cull
  const markedMdb = createSegmentedMdbList();
  const shiftedMarked: MultiDrawBaseData[] = [];
  let markedLists: MultiDrawBaseData[] = [];
  let markedUsed = false;
  const markedKey = createMarkedKey();
  /** per member the draw slots of its instances (`createInstancePositions`), made as a member is first marked */
  const markedPositions: (ReturnType<typeof createInstancePositions> | undefined)[] = new Array(members.length);
  const getMarkedVersion = createDrawVersion((out) => appendMdbDrawKey(markedLists, out));
  let markedVersion = 0;
  /** the marked lists for the members' markers and the view, where they are few (`isMarkedFew`) */
  function updateMarked() {
    syncCullCounts();
    const lod = markedLod;
    if (!lod) {
      markedUsed = false;
      markedKey.reset();
      return;
    }
    markedKey.begin(lod);
    markedKey.add(merged.markerVersion.ref.value);
    markedKey.add(merged.cullVersion.ref.value);
    markedKey.add(merged.values.uObjectTransform.ref.version);
    // shifted to the identity region (`shiftLists`): moved where the ghost's or raster regions are added
    markedKey.add(culler ? culler.identityBase : -1);
    if (markedKey.isSame()) return;

    const levelCount = (memberCull[0].lodLevels?.ref.value as LodLevelsValue | undefined)?.length ?? 0;
    let marked = 0;
    // a source's markers have no counts a member: not known where
    markedUsed = !merged.source;
    // the marked members, not a pass over all (100k members: half a hover's frame)
    const markedMembers = merged.getMarkedMembers?.();
    if (!markedUsed || !markedMembers) {
      markedUsed = false;
      return;
    }
    for (let k = 0, kl = markedMembers.length; k < kl; ++k) {
      const v = members[markedMembers[k]];
      // markers without counts: not known where
      if (!v.markerInstanceCounts) {
        markedUsed = false;
        break;
      }
      marked += v.markedInstanceCount!.ref.value;
    }
    if (markedUsed) markedUsed = isMarkedFew(marked, levelCount, merged.values.instanceCount.ref.value as number);
    if (!markedUsed) return;

    const hasLodLevels = !!markedMdb.prepare(memberCull[0].lodLevels, Math.max(1, marked), true, lod);
    for (let k = 0, kl = markedMembers.length; k < kl; ++k) {
      const i = markedMembers[k];
      const v = members[i];
      markedMdb.markedSegment(
        memberCull[i],
        setSegment(i),
        getMemberLod(i, hasLodLevels),
        v.markerInstanceCounts!.ref.value,
        // the draw slots of the identity region: a built grid's order where a pending member's was built
        (markedPositions[i] ??= createInstancePositions())(
          isIndexBuilt(i) ? indexOrders[i]! : v.aInstance,
          v.instanceCount.ref.value,
        ),
        getObjectTransforms(v as unknown as TransformValues),
        v.invariantBoundingSphere.ref.value,
        markedPlane,
        markedFrustum,
      );
    }
    markedLists = culler ? shiftLists(markedMdb.list, shiftedMarked) : markedMdb.list;
    markedVersion = getMarkedVersion();
  }

  buildFullList(LodView());

  /**
   * The draws of the listed segments (ascending), each whole at its levels, as lists of their own (the main lists stay):
   * the solid-interior caps of the segments a cap cuts (`Renderable.renderSegments`), a few a frame
   */
  const segmentMdb = createSegmentedMdbList();
  const shiftedSegmentLists: MultiDrawBaseData[] = [];
  function renderSegments(variant: GraphicsRenderVariant, globals: GlobalBindings, segments: readonly number[]) {
    syncCullCounts();
    const hasLodLevels = !!segmentMdb.prepare(memberCull[0].lodLevels, segments.length, true, segmentLod);
    for (let k = 0, kl = segments.length; k < kl; ++k) {
      const i = segments[k];
      segmentMdb.fullSegment(memberCull[i], setSegment(i), getMemberLod(i, hasLodLevels));
    }
    renderItem.render(variant, globals, culler ? shiftLists(segmentMdb.list, shiftedSegmentLists) : segmentMdb.list);
  }

  let programVersion = 0;
  /** clip objects, not opted out of the ghost (as `DrawFact.ClipGhost`) */
  const hasClipGhost = () => {
    const v = renderValues as { dClipObjects?: ValueCell<boolean>; clipGhost?: ValueCell<boolean> };
    return !!v.dClipObjects?.ref.value && v.clipGhost?.ref.value !== false;
  };
  // a version a set (one key each: the visible and the ghost set culled in one frame would each advance a shared one
  // every frame), from one counter (a set's version never another's with other lists)
  const drawVersions = createDrawVersionCounter();
  const appendDrawKey = (out: DrawKey) => appendMdbDrawKey(set.mode === 'gpu' ? gpuLists : set.mdb.list, out);
  const getVisibleDrawVersion = createDrawVersion(appendDrawKey, drawVersions);
  const getGhostDrawVersion = createDrawVersion(appendDrawKey, drawVersions);
  const getDrawVersion = () => (set === ghostSet ? getGhostDrawVersion() : getVisibleDrawVersion());
  visibleSet.drawVersion = ghostSet.drawVersion = getDrawVersion();
  return {
    id,
    materialId,
    renderItem,
    values: renderValues,
    state,

    cull,
    uncull: (lod: LodView, clipGhost = false) => {
      set = getSet(clipGhost);
      set.lastCullFrame = undefined;
      // levels by size on screen convert their ranges with the view (`prepare`)
      const pixel = ((members[0].lodLevels?.ref.value as LodLevelsValue | undefined)?.pixelSize ?? 0) > 0;
      if (set.mode === 'full' && !pixel) return;
      buildFullList(lod);
      set.drawVersion = getDrawVersion();
    },
    prepareCull: culler
      ? (
          cameraPlane: Plane3D,
          frustum: Frustum3D,
          frame: Frame,
          lod: LodView,
          pyramid?: OcclusionPyramid,
          raster?: MergedRasterParams,
          batch?: CullBatch,
          clipGhost = false,
          twoPhase = false,
        ) => {
          // once per view of a frame: another view of it (a stereo or XR eye) culled for it too; a page store's by its
          // camera (a jittered sample the same), the others' by their planes (jittered samples cull with the camera's)
          if (frame === gpuFrame && (culler.paged ? isGpuCamera(cameraPlane) : isGpuView(cameraPlane, frustum))) {
            return true;
          }
          // the sphere raster's split per instance for objects it can draw
          const split = raster && 'dSphereResolve' in renderValues && isSphereRasterEligible(renderValues as never);
          // the ghost's region where the object has a ghost (`DrawFact.ClipGhost`)
          const ghost = clipGhost && hasClipGhost();
          gpuFrame = culler.dispatch(
            cameraPlane,
            frustum,
            lod,
            pyramid,
            split ? raster : undefined,
            batch,
            ghost,
            twoPhase,
          )
            ? frame
            : undefined;
          if (gpuFrame) setGpuView(cameraPlane, frustum);
          return !!gpuFrame;
        }
      : undefined,
    prepareCullView: culler
      ? (cameraPlane: Plane3D, frustum: Frustum3D, frame: Frame, lod: LodView) => {
          // a still view of its own culls nothing again (`MergedCuller.dispatch`'s key); the ghost kept (drawn after, the
          // illumination's; not tested against the pyramid: the same)
          const ghost = culler.ghostDispatched;
          gpuFrame = culler.dispatch(cameraPlane, frustum, lod, undefined, undefined, undefined, ghost)
            ? frame
            : undefined;
          if (gpuFrame) setGpuView(cameraPlane, frustum);
          return !!gpuFrame;
        }
      : undefined,
    isGpuCullStale: culler
      ? (frame: Frame) => visibleSet.lastCullFrame === frame && visibleSet.mode === 'gpu' && gpuFrame !== frame
      : undefined,
    get rasterSpheres() {
      return rasterOnly ? getRasterSpheres() : undefined;
    },
    needsCullSecond: culler ? (frame: Frame) => frame === gpuFrame && culler.secondPending : undefined,
    prepareCullSecond: culler
      ? (frame: Frame, pyramid: OcclusionPyramid | undefined, batch?: CullBatch) =>
          frame === gpuFrame && culler.dispatchSecond(pyramid, batch)
      : undefined,
    // where the view draws the GPU culling's lists, as `cull` decides: a frame culled already keeps its lists (views of
    // a reused frame, as the samples of temporal multi-sampling), else the view culled on the GPU
    getRasterInput: culler
      ? (frame: Frame, cameraPlane: Plane3D, frustum: Frustum3D) => {
          const v = visibleSet;
          const gpu =
            frame === v.lastCullFrame ? v.mode === 'gpu' : frame === gpuFrame && isGpuView(cameraPlane, frustum);
          return gpu ? culler.getRasterInput() : undefined;
        }
      : undefined,
    render: (variant: GraphicsRenderVariant, globals: GlobalBindings) => {
      renderItem.render(variant, globals, getLists());
    },
    renderMarked: (variant: GraphicsRenderVariant, globals: GlobalBindings) => {
      updateMarked();
      renderItem.render(variant, globals, markedUsed ? markedLists : getLists());
    },
    prepareMarked: (cameraPlane: Plane3D, frustum: Frustum3D, lod: LodView) => {
      setMarkedView(cameraPlane, frustum, lod);
      updateMarked();
      return markedUsed;
    },
    get markedDrawVersion() {
      updateMarked();
      return getMarkedDrawVersion(markedUsed, markedVersion, visibleSet.drawVersion * 2 + (rasterDrawn ? 1 : 0));
    },
    renderSecond: culler
      ? (variant: GraphicsRenderVariant, globals: GlobalBindings) => {
          if (visibleSet.mode !== 'gpu' || !culler.secondDispatched) return;
          renderItem.render(variant, globals, culler.getSecondLists(visibleSet.mdb.list, !rasterDrawn));
        }
      : undefined,
    renderSegments,
    getClipInfo: clipInfo.get,
    memberCells: merged.memberCells,
    getSegmentSpheres: () => {
      // the members' cells, read again: a member's sphere may be replaced; as drawn (moved by the object's transform)
      for (let i = 0, il = members.length; i < il; ++i) segmentSpheres[i] = getDrawnBoundingSphere(memberCull[i]);
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
      return visibleSet.mode === 'gpu' ? culler!.recordsCompact : undefined;
    },
    get recordsVersion() {
      return culler?.recordsVersion ?? 0;
    },
    get drawVersion() {
      // the lists without the raster regions are other draws (recordings of the pass replay while it is the same); of
      // the visible passes' lists: draw lists record those
      return visibleSet.drawVersion * 2 + (rasterDrawn ? 1 : 0);
    },
    update: (observedOnly?: boolean) => {
      // the records unwritten but observed cells (`Merged.recordsUnchanged`: a clip moved): synced as a mark, but the
      // CPU culling's lists culled again (their clip objects are the object's)
      const membersObserved = !!observedOnly || !!merged.recordsUnchanged?.();
      // none pending where only observed cells changed: a grid's change comes with a full update, which applied it
      if (!culler && !membersObserved) applyPendingGrids();
      const layoutChanged = merged.sync(membersObserved);
      const dataVersion = culler?.dataVersion;
      culler?.sync(membersObserved);
      if (hasPending) syncIndexes();
      storage?.sync();
      // the frame's GPU culling stands where what it culled is the same (e.g. markers changed: marking redrawn over the
      // frame with its camera); the CPU culling's lists where only observed cells changed (a mark: what they read moves
      // `cullVersion`, which their view's key holds)
      if (culler ? layoutChanged || culler.dataVersion !== dataVersion : layoutChanged || !observedOnly) {
        gpuFrame = undefined;
        visibleSet.mode = ghostSet.mode = 'none';
      }
      renderItem.update();
      visibleSet.lastCullFrame = ghostSet.lastCullFrame = undefined;
    },
    dispose: () => {
      renderItem.destroy();
      culler?.destroy();
      storage?.destroy();
    },
  };
}
