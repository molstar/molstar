/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ValueCell } from '@molstar/core/util';
import type { AttributeItemSize, ElementsKind, AttributeValues, AttributeKind, DataTypeArrayType } from '../buffer.js';
import type { UniformKind, UniformValues, UniformKindValue } from '../uniform.js';
import type { DefineKind, DefineValues } from '../shader-code.js';
import type { Mat4 } from '@molstar/core/math/linear-algebra';
import type {
  TextureValues,
  TextureType,
  TextureFormat,
  TextureFilter,
  TextureKind,
  TextureKindValue,
} from '../texture.js';
import type { Sphere3D } from '@molstar/core/math/geometry';
import type { InstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { Clip } from '@molstar/core/util/clip';
import type { DrawIndirect, StorageBuffer, StorageBufferSpec } from '../compute.js';
import { convertDepthRange, type DepthRangeConversion, type DepthRangeTag } from '../depth-range.js';
import { flipWindowRects, type WindowRectTag } from '../window-rect.js';
import type { DataArray } from './util.js';

export type ValueKindType = {
  number: number;
  string: string;
  boolean: boolean;
  any: any;
  unknown: unknown;

  m4: Mat4;
  float32: Float32Array;
  uint32: Uint32Array;
  sphere: Sphere3D;
  instanceGrid: InstanceGrid;
  drawIndirect: DrawIndirect | undefined;
  clipObjects: Clip.Objects;
  /** storage buffers of render schemas (`StorageBufferSpec`), unset on backends without them */
  'storage-buffer': StorageBuffer | undefined;
};
export type ValueKind = keyof ValueKindType;

//

export type DataKindValue = {
  'data-uint8': DataArray<Uint8Array>;
  'data-float32': DataArray<Float32Array>;
};
export type DataKind = keyof DataKindValue;
export type DataValues = { [k: string]: ValueCell<DataArray<any>> };

export type KindValue = UniformKindValue & DataTypeArrayType & TextureKindValue & ValueKindType & DataKindValue;

export type Values<S extends RenderableSchema> = { readonly [k in keyof S]: ValueCell<KindValue[S[k]['kind']]> };
export type UnboxedValues<S extends RenderableSchema> = { readonly [k in keyof S]: KindValue[S[k]['kind']] };

/** How a backend maps tagged uniforms where it binds them */
export interface UniformConversions {
  /** the projections (`UniformSpec.depthRange`), see `DepthRangeTag` */
  readonly depthRange?: DepthRangeConversion;
  /** the height of the render pass target, to flip window rectangles (`UniformSpec.windowRect`), see `WindowRectTag` */
  readonly targetHeight?: () => number;
}

/** The values by kind, the tagged uniforms mapped with `conversions` */
export function splitValues(schema: RenderableSchema, values: RenderableValues, conversions: UniformConversions = {}) {
  const attributeValues: AttributeValues = {};
  const defineValues: DefineValues = {};
  const textureValues: TextureValues = {};
  const materialTextureValues: TextureValues = {};
  const dataValues: DataValues = {};
  const uniformValues: UniformValues = {};
  const materialUniformValues: UniformValues = {};
  const bufferedUniformValues: UniformValues = {};
  const storageBufferValues: { [k: string]: ValueCell<StorageBuffer | undefined> } = {};
  Object.keys(schema).forEach((k) => {
    if (schema[k].type === 'storage-buffer' && values[k] !== undefined) storageBufferValues[k] = values[k];
    const spec = schema[k];
    if (spec.type === 'attribute') attributeValues[k] = values[k];
    if (spec.type === 'define') defineValues[k] = values[k];
    if (spec.type === 'data' && values[k] !== undefined) dataValues[k] = values[k];
    // check if k exists in values to exclude global textures
    if (spec.type === 'texture' && values[k] !== undefined) {
      if (spec.variant === 'material') materialTextureValues[k] = values[k];
      else textureValues[k] = values[k];
    }
    // check if k exists in values to exclude global uniforms
    if (spec.type === 'uniform' && values[k] !== undefined) {
      if (spec.variant === 'material') materialUniformValues[k] = values[k];
      else if (spec.variant === 'buffered') bufferedUniformValues[k] = values[k];
      else uniformValues[k] = values[k];
    }
  });
  return {
    attributeValues,
    defineValues,
    textureValues,
    materialTextureValues,
    dataValues,
    uniformValues: convertUniforms(schema, uniformValues, conversions),
    materialUniformValues: convertUniforms(schema, materialUniformValues, conversions),
    bufferedUniformValues: convertUniforms(schema, bufferedUniformValues, conversions),
    storageBufferValues,
  };
}

function convertUniforms(schema: RenderableSchema, values: UniformValues, conversions: UniformConversions) {
  const converted = convertDepthRange(schema, values, conversions.depthRange);
  return flipWindowRects(schema, converted, conversions.targetHeight);
}

//

/**
 * Condition for a texture or attribute to be used by a program: any of the given flattened
 * define names (as tested in shader code, e.g. 'dUsePalette', 'dColorType_grid',
 * 'emulateVertexId') is defined, a '!' prefix negates a name.
 * Must cover all uses in the shader code, a superset is fine.
 */
export type UsedIf = readonly string[];

/** Whether a flattened define name is defined, `undefined` if not known (yet) */
export type DefineLookup = (name: string) => boolean | undefined;

export function getDefineLookup(names: Iterable<string>): DefineLookup {
  const set = new Set(names);
  return (name) => set.has(name);
}

/** Flattened names of the defined schema defines, like `#define`d by the shader header (string `dX: 'v'` → `dX_v`) */
export function getDefineValueNames(defineValues: DefineValues): string[] {
  const names: string[] = [];
  for (const k of Object.keys(defineValues)) {
    const v = defineValues[k].ref.value;
    if (typeof v === 'string') names.push(`${k}_${v}`);
    else if (v !== false && v !== undefined) names.push(k);
  }
  return names;
}

/** Lookup of the flattened names of the schema defines, see `getDefineValueNames` */
export function getDefineValuesLookup(defineValues: DefineValues): DefineLookup {
  return getDefineLookup(getDefineValueNames(defineValues));
}

/** `false` only if `spec` is known to be unused, i.e. no `usedIf` condition can hold */
export function isUsed(spec: { readonly usedIf?: UsedIf }, defines: DefineLookup): boolean {
  if (!spec.usedIf) return true;
  for (const condition of spec.usedIf) {
    const negate = condition.startsWith('!');
    const defined = defines(negate ? condition.slice(1) : condition);
    if (defined === undefined || defined !== negate) return true;
  }
  return false;
}

/** `kind`, `itemSize` and `divisor` are the complete logical description, backends decide the transport */
export type AttributeSpec<K extends AttributeKind> = {
  type: 'attribute';
  kind: K;
  itemSize: AttributeItemSize;
  divisor: number;
  usedIf?: UsedIf;
};
export function AttributeSpec<K extends AttributeKind>(
  kind: K,
  itemSize: AttributeItemSize,
  divisor: number,
  usedIf?: UsedIf,
): AttributeSpec<K> {
  return { type: 'attribute', kind, itemSize, divisor, usedIf };
}

type UniformVariant = 'material' | 'buffered';
/**
 * On WebGL2, uniforms in a block are backed by uniform buffer objects:
 * - 'global': one buffer per renderer, updated by the renderer
 * - 'material': one buffer per materialId, updated on material change
 * - 'object': one buffer per render item, updated on value change
 *
 * On WebGL1 they are set as classic uniforms.
 */
export type UniformBlock = 'global' | 'material' | 'object';
/**
 * Where a uniform lives; `undefined` and 'draw' are classic uniforms on all backends:
 * 'draw' for values that change between the draws of a render item (see `MultiDrawBaseData`),
 * `undefined` for textures and arrays without a fixed `length`
 */
export type UniformScope = UniformBlock | 'draw';
export type UniformSpec<K extends UniformKind> = {
  type: 'uniform';
  kind: K;
  variant?: UniformVariant;
  block?: UniformScope;
  length?: number;
  /** number define with the classic (WebGL1) array size, `length` is the max. (block) size */
  lengthDefine?: string;
  /** a loop bound: a number define on WebGL1, whose loops need constant bounds */
  loopBound?: boolean;
  /** a matrix holding the projection, mapped by each backend to its reversed clip depth, see `DepthRangeTag` */
  depthRange?: DepthRangeTag;
  /** a window rectangle (top-left origin), flipped by backends with a bottom-left origin, see `WindowRectTag` */
  windowRect?: WindowRectTag;
};
/**
 * Unless given, `block` defaults to 'material' for the 'material' variant and to 'object' otherwise.
 * Array uniforms need a (max) `length` to be in a block.
 */
export function UniformSpec<K extends UniformKind>(
  kind: K,
  variant?: UniformVariant,
  block?: UniformScope,
  length?: number,
  lengthDefine?: string,
): UniformSpec<K> {
  const isArray = kind.endsWith('[]');
  if ((length !== undefined || lengthDefine !== undefined) && !isArray) {
    throw new Error(`'length' requires an array uniform kind, got '${kind}'`);
  }
  if (kind === 't' || (isArray && length === undefined)) {
    if (block !== undefined && block !== 'draw') {
      throw new Error(
        `uniform kind '${kind}' cannot be assigned to block '${block}'${isArray ? ' without a length' : ''}`,
      );
    }
  } else if (block === undefined) {
    block = variant === 'material' ? 'material' : 'object';
  }
  if (variant === 'buffered' && block !== undefined && block !== 'draw' && block !== 'object') {
    throw new Error(`'buffered' uniforms cannot be assigned to block '${block}'`);
  }
  if (block === 'material') {
    variant = 'material';
  } else if (variant === 'material' && block !== undefined && block !== 'draw') {
    throw new Error(`'material' uniforms cannot be assigned to block '${block}'`);
  }
  return { type: 'uniform', kind, variant, block, length, lengthDefine };
}

function GlobalUniformSpec<K extends UniformKind>(kind: K, length?: number, lengthDefine?: string): UniformSpec<K> {
  return UniformSpec(kind, undefined, 'global', length, lengthDefine);
}
function MaterialUniformSpec<K extends UniformKind>(kind: K, length?: number, lengthDefine?: string): UniformSpec<K> {
  return UniformSpec(kind, undefined, 'material', length, lengthDefine);
}
function ObjectUniformSpec<K extends UniformKind>(kind: K): UniformSpec<K> {
  return UniformSpec(kind, undefined, 'object');
}
function DrawUniformSpec<K extends UniformKind>(kind: K): UniformSpec<K> {
  return UniformSpec(kind, undefined, 'draw');
}

/** max. array lengths in uniform blocks */
export const MaxLightCount = 16;

type TextureVariant = 'material';
export type TextureDimension = '2d' | '3d' | 'cube';
export const TextureDimensions: readonly TextureDimension[] = ['2d', '3d', 'cube'];
/** as in WebGPU; integer types map to float samplers on WebGL1 */
export type TextureSampleType = 'float' | 'unfilterable-float' | 'sint' | 'uint' | 'depth';
export type TextureSpec<K extends TextureKind> = {
  type: 'texture';
  kind: K;
  format: TextureFormat;
  dataType: TextureType;
  filter: TextureFilter;
  variant?: TextureVariant;
  /** `{ define }`: given by the value of a string define, e.g. 'dGridTexType' with values '2d' and '3d' */
  dimension: TextureDimension | { readonly define: string };
  sampleType: TextureSampleType;
  usedIf?: UsedIf;
  /** see `TextureSpecOptions.shared` */
  shared?: boolean;
  /** a backend's texture of a data array (`DataSpec`) of this format, read by its accessor (WebGL data textures) */
  data?: DataFormat;
};
export type TextureSpecOptions = {
  /** defaults to '3d' for 'volume-*' kinds, '2d' otherwise */
  dimension?: TextureSpec<TextureKind>['dimension'];
  /** defaults from `format`, `dataType` and `filter`, see `getDefaultSampleType` */
  sampleType?: TextureSampleType;
  usedIf?: UsedIf;
  /**
   * Render items of a context share one texture per image object (`TextureImage`, `TextureVolume`), for images shared
   * by many objects (the font atlas); not for kind 'texture' (owned outside the items), see
   * `.v6/plans/gl/shared-textures.md`
   */
  shared?: boolean;
};
export function TextureSpec<K extends TextureKind>(
  kind: K,
  format: TextureFormat,
  dataType: TextureType,
  filter: TextureFilter,
  variant?: TextureVariant,
  options?: TextureSpecOptions,
): TextureSpec<K> {
  const dimension = options?.dimension ?? (kind.startsWith('volume-') ? '3d' : '2d');
  const sampleType = options?.sampleType ?? getDefaultSampleType(format, dataType, filter);
  const spec: TextureSpec<K> = {
    type: 'texture',
    kind,
    format,
    dataType,
    filter,
    variant,
    dimension,
    sampleType,
    usedIf: options?.usedIf,
  };
  if (options?.shared) {
    if (kind === 'texture') throw new Error("textures of kind 'texture' cannot be shared");
    spec.shared = true;
  }
  return spec;
}

function getDefaultSampleType(format: TextureFormat, dataType: TextureType, filter: TextureFilter): TextureSampleType {
  if (format === 'depth' || format === 'depth-stencil') return 'depth';
  if (dataType === 'int') return 'sint';
  // 32-bit float textures are not filterable everywhere
  if (dataType === 'float' && filter === 'nearest') return 'unfilterable-float';
  return 'float';
}

export function getTextureDimension(spec: TextureSpec<TextureKind>, defines?: DefineLookup): TextureDimension {
  const { dimension } = spec;
  if (typeof dimension === 'string') return dimension;
  const resolved = defines && TextureDimensions.find((d) => defines(`${dimension.define}_${d}`));
  if (!resolved) throw new Error(`texture dimension depends on define '${dimension.define}', which is not set`);
  return resolved;
}

/**
 * Elements of a data array by index (`DataSpec`): `rgba8` four bytes, `r8` one, `rgba32f` four floats; read as `vec4`
 * (`r8` in `.r`) by the generated accessors `read_<name>(i)`
 */
export type DataFormat = 'rgba8' | 'r8' | 'rgba32f';
export const DataFormatItemSize: { readonly [F in DataFormat]: number } = { rgba8: 4, r8: 1, rgba32f: 4 };
/**
 * Data read by element index (group, instance or vertex), a flat array without a 2D layout (`DataArray`): the backends
 * lay it out, WebGL as a data texture (its size uniform its own), WebGPU in storage pages; shaders read it with the
 * accessor `read_<name>(i)`
 */
export type DataSpec<K extends DataKind = DataKind> = { type: 'data'; kind: K; format: DataFormat; usedIf?: UsedIf };
export function DataSpec<F extends DataFormat>(
  format: F,
  usedIf?: UsedIf,
): DataSpec<F extends 'rgba32f' ? 'data-float32' : 'data-uint8'> {
  const kind = (format === 'rgba32f' ? 'data-float32' : 'data-uint8') as F extends 'rgba32f'
    ? 'data-float32'
    : 'data-uint8';
  return { type: 'data', kind, format, usedIf };
}

export type ElementsSpec<K extends ElementsKind> = { type: 'elements'; kind: K };
export function ElementsSpec<K extends ElementsKind>(kind: K): ElementsSpec<K> {
  return { type: 'elements', kind };
}

export type DefineSpec<K extends DefineKind> = { type: 'define'; kind: K; options?: string[] };
export function DefineSpec<K extends DefineKind>(kind: K, options?: string[]): DefineSpec<K> {
  return { type: 'define', kind, options };
}

export type ValueSpec<K extends ValueKind> = { type: 'value'; kind: K };
export function ValueSpec<K extends ValueKind>(kind: K): ValueSpec<K> {
  return { type: 'value', kind };
}

//

export type RenderableSchema = {
  readonly [k: string]:
    | AttributeSpec<AttributeKind>
    | UniformSpec<UniformKind>
    | TextureSpec<TextureKind>
    | DataSpec
    | ValueSpec<ValueKind>
    | DefineSpec<DefineKind>
    | ElementsSpec<ElementsKind>
    | StorageBufferSpec;
};
export type RenderableValues = { readonly [k: string]: ValueCell<any> };

//

export const GlobalUniformSchema = {
  uDrawId: DrawUniformSpec('i'),

  uModel: GlobalUniformSpec('m4'),
  uView: GlobalUniformSpec('m4'),
  uInvView: GlobalUniformSpec('m4'),
  uModelView: GlobalUniformSpec('m4'),
  uInvModelView: GlobalUniformSpec('m4'),
  uProjection: { ...GlobalUniformSpec('m4'), depthRange: 'projection' as const },
  uInvProjection: { ...GlobalUniformSpec('m4'), depthRange: 'inverse' as const },
  uModelViewProjection: { ...GlobalUniformSpec('m4'), depthRange: 'projection' as const },
  uInvModelViewProjection: { ...GlobalUniformSpec('m4'), depthRange: 'inverse' as const },
  uHasHeadRotation: GlobalUniformSpec('b'),
  uInvHeadRotation: GlobalUniformSpec('m4'),
  uIsAsymmetricProjection: GlobalUniformSpec('b'),
  uHasEyeCamera: GlobalUniformSpec('b'),
  uModelViewEye: GlobalUniformSpec('m4'),
  uInvModelViewEye: GlobalUniformSpec('m4'),

  uIsOrtho: GlobalUniformSpec('f'),
  uPixelRatio: GlobalUniformSpec('f'),
  uViewport: GlobalUniformSpec('v4'),
  uViewOffset: GlobalUniformSpec('v2'),
  uModelScale: GlobalUniformSpec('f'),
  uDrawingBufferSize: GlobalUniformSpec('v2'),

  uCameraPosition: GlobalUniformSpec('v3'),
  uCameraDir: GlobalUniformSpec('v3'),
  uCameraPlane: GlobalUniformSpec('v4'),
  uNear: GlobalUniformSpec('f'),
  uFar: GlobalUniformSpec('f'),
  uFog: GlobalUniformSpec('b'),
  uFogNear: GlobalUniformSpec('f'),
  uFogFar: GlobalUniformSpec('f'),
  uFogColor: GlobalUniformSpec('v3'),

  uTransparentBackground: GlobalUniformSpec('b'),

  uLightDirection: GlobalUniformSpec('v3[]', MaxLightCount, 'dLightCount'),
  uLightColor: GlobalUniformSpec('v3[]', MaxLightCount, 'dLightCount'),
  uAmbientColor: GlobalUniformSpec('v3'),

  uPickingAlphaThreshold: GlobalUniformSpec('f'),

  uHighlightColor: GlobalUniformSpec('v3'),
  uSelectColor: GlobalUniformSpec('v3'),
  uDimColor: GlobalUniformSpec('v3'),
  uHighlightStrength: GlobalUniformSpec('f'),
  uSelectStrength: GlobalUniformSpec('f'),
  uDimStrength: GlobalUniformSpec('f'),
  uMarkerPriority: GlobalUniformSpec('i'),
  uMarkerAverage: GlobalUniformSpec('f'),

  uXrayEdgeFalloff: GlobalUniformSpec('f'),
  uCelSteps: GlobalUniformSpec('f'),
  uExposure: GlobalUniformSpec('f'),

  uRenderMask: GlobalUniformSpec('i'),
  uMarkingDepthTest: GlobalUniformSpec('b'),
  uDepthBack: GlobalUniformSpec('b'),
  uMarkingType: GlobalUniformSpec('i'),
  uPickType: GlobalUniformSpec('i'),
  uSolidInteriorPass: DrawUniformSpec('i'),
  uSolidInteriorClip: DrawUniformSpec('i'),

  uTime: GlobalUniformSpec('f'),
  uEnableAnimation: GlobalUniformSpec('b'),
} as const;
export type GlobalUniformSchema = typeof GlobalUniformSchema;
export type GlobalUniformValues = Values<GlobalUniformSchema>;

export const GlobalTextureSchema = {
  /** opaque scene depth: a depth texture where supported (`depthTextureSupport`), packed `rgba` otherwise */
  tDepth: TextureSpec('texture', 'depth', 'ushort', 'nearest'),
  /** depth with alpha packed into `rgba` (`packDepthWithAlphaToRGBA`): front transparent depth, marking depth */
  tDepthAlpha: TextureSpec('texture', 'rgba', 'ubyte', 'nearest', undefined, {
    usedIf: ['dRenderVariant_emissive', 'dRenderVariant_marking'],
  }),

  // dpoit
  tDpoitDepth: TextureSpec('texture', 'rg', 'float', 'nearest', undefined, { usedIf: ['dRenderVariant_colorDpoit'] }),
  tDpoitFrontColor: TextureSpec('texture', 'rgba', 'float', 'nearest', undefined, {
    usedIf: ['dRenderVariant_colorDpoit'],
  }),
} as const;
export type GlobalTextureSchema = typeof GlobalTextureSchema;
export type GlobalTextureValues = Values<GlobalTextureSchema>;

export const GlobalDefineSchema = {
  dLightCount: DefineSpec('number'),
  dColorMarker: DefineSpec('boolean'),
} as const;
export type GlobalDefineSchema = typeof GlobalDefineSchema;
export type GlobalDefineValues = Values<GlobalDefineSchema>;
export type GlobalDefines = UnboxedValues<GlobalDefineSchema>;

export const InternalSchema = {
  uObjectId: ObjectUniformSpec('i'),
} as const;
export type InternalSchema = typeof InternalSchema;
export type InternalValues = Values<InternalSchema>;

/**
 * For render items merged from multiple sets of member values (segments),
 * see `createMergedRenderObject`.
 */
export const SegmentSchema = {
  /**
   * per-instance segment data:
   * x: segment groupCount,
   * y: groupInstance texel offset (base - instanceBase * groupCount),
   * z: group texel offset
   */
  aSegment: AttributeSpec('float32', 3, 1, ['dSegmented']),
  /** per-instance invariant bounding sphere (xyz: center, w: radius) of the segment */
  aSegmentSphere: AttributeSpec('float32', 4, 1, ['dSegmented']),
  /** per-instance, per-lod-level radius scale (see spheres' `lodLevels`), duplicated across a segment's instances */
  aSegmentLod: AttributeSpec('float32', 4, 1, ['dSegmented']),
  /** index of the lod level currently being drawn, selects a component of `aSegmentLod` */
  uLodLevel: DrawUniformSpec('i'),
  dSegmented: DefineSpec('boolean'),
} as const;
export type SegmentSchema = typeof SegmentSchema;
export type SegmentValues = Values<SegmentSchema>;

export const ColorSchema = {
  uColor: MaterialUniformSpec('v3'),
  /** the index of `tColor` the 'texture' color type reads, see `ColorTypeIndex` */
  uColorIndexType: ObjectUniformSpec('i'),
  uColorGridDim: ObjectUniformSpec('v3'),
  /** the size of the 2D layout of the grid (`tColorGrid`) */
  uColorGridTexDim: ObjectUniformSpec('v2'),
  uColorGridTransform: ObjectUniformSpec('v4'),
  uPaletteDomain: ObjectUniformSpec('v2'),
  uPaletteDefault: ObjectUniformSpec('v3'),
  tColor: DataSpec('rgba8'),
  tPalette: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest', undefined, { usedIf: ['dUsePalette'] }),
  tColorGrid: TextureSpec('texture', 'rgba', 'ubyte', 'linear', undefined, {
    usedIf: ['dColorType_grid'],
  }),
  /**
   * the color type's define (`getColorTypeDefine`): the texture types share 'texture', the grid types 'grid',
   * `uColorIndexType` selects
   */
  dColorType: DefineSpec('string', ['uniform', 'texture', 'grid', 'direct']),
  /** the color type (`ColorType`) */
  colorType: ValueSpec('string'),
  dUsePalette: DefineSpec('boolean'),
} as const;
export type ColorSchema = typeof ColorSchema;
export type ColorValues = Values<ColorSchema>;

export const SizeSchema = {
  uSize: MaterialUniformSpec('f'),
  tSize: DataSpec('rgba8'),
  /** the index of `tSize` by the size type (`SizeTypeIndex`) */
  uSizeIndexType: ObjectUniformSpec('i'),
  /** the size type's define (`getSizeTypeDefine`): the texture types share 'texture', `uSizeIndexType` selects */
  dSizeType: DefineSpec('string', ['uniform', 'texture']),
  /** the size type (`SizeType`) */
  sizeType: ValueSpec('string'),
  uSizeFactor: MaterialUniformSpec('f'),
} as const;
export type SizeSchema = typeof SizeSchema;
export type SizeValues = Values<SizeSchema>;

export const MarkerSchema = {
  uMarker: ObjectUniformSpec('f'),
  tMarker: DataSpec('r8'),
  markerAverage: ValueSpec('number'),
  markerStatus: ValueSpec('number'),
  /** the marker type (`MarkerType`); `uInstanceGranularity` selects the index */
  markerType: ValueSpec('string'),
} as const;
export type MarkerSchema = typeof MarkerSchema;
export type MarkerValues = Values<MarkerSchema>;

export const OverpaintSchema = {
  tOverpaint: DataSpec('rgba8'),
  dOverpaint: DefineSpec('boolean'),

  uOverpaintGridDim: ObjectUniformSpec('v3'),
  uOverpaintGridTexDim: ObjectUniformSpec('v2'),
  uOverpaintGridTransform: ObjectUniformSpec('v4'),
  tOverpaintGrid: TextureSpec('texture', 'rgba', 'ubyte', 'linear', undefined, {
    usedIf: ['dOverpaintType_volumeInstance'],
  }),
  /** the layer type (`LayerType`) */
  overpaintType: ValueSpec('string'),
  /** the layer type's define (`getLayerTypeDefine`): instance and groupInstance share 'texture' */
  dOverpaintType: DefineSpec('string', ['texture', 'vertexInstance', 'volumeInstance']),
  uOverpaintStrength: MaterialUniformSpec('f'),
} as const;
export type OverpaintSchema = typeof OverpaintSchema;
export type OverpaintValues = Values<OverpaintSchema>;

export const TransparencySchema = {
  tTransparency: DataSpec('r8'),
  dTransparency: DefineSpec('boolean'),
  transparencyAverage: ValueSpec('number'),
  transparencyMin: ValueSpec('number'),

  uTransparencyGridDim: ObjectUniformSpec('v3'),
  uTransparencyGridTexDim: ObjectUniformSpec('v2'),
  uTransparencyGridTransform: ObjectUniformSpec('v4'),
  tTransparencyGrid: TextureSpec('texture', 'r', 'ubyte', 'linear', undefined, {
    usedIf: ['dTransparencyType_volumeInstance'],
  }),
  /** the layer type (`LayerType`) */
  transparencyType: ValueSpec('string'),
  /** the layer type's define (`getLayerTypeDefine`): instance and groupInstance share 'texture' */
  dTransparencyType: DefineSpec('string', ['texture', 'vertexInstance', 'volumeInstance']),
  uTransparencyStrength: MaterialUniformSpec('f'),
} as const;
export type TransparencySchema = typeof TransparencySchema;
export type TransparencyValues = Values<TransparencySchema>;

export const EmissiveSchema = {
  tEmissive: DataSpec('r8'),
  dEmissive: DefineSpec('boolean'),
  emissiveAverage: ValueSpec('number'),

  uEmissiveGridDim: ObjectUniformSpec('v3'),
  uEmissiveGridTexDim: ObjectUniformSpec('v2'),
  uEmissiveGridTransform: ObjectUniformSpec('v4'),
  tEmissiveGrid: TextureSpec('texture', 'r', 'ubyte', 'linear', undefined, {
    usedIf: ['dEmissiveType_volumeInstance'],
  }),
  /** the layer type (`LayerType`) */
  emissiveType: ValueSpec('string'),
  /** the layer type's define (`getLayerTypeDefine`): instance and groupInstance share 'texture' */
  dEmissiveType: DefineSpec('string', ['texture', 'vertexInstance', 'volumeInstance']),
  uEmissiveStrength: MaterialUniformSpec('f'),
} as const;
export type EmissiveSchema = typeof EmissiveSchema;
export type EmissiveValues = Values<EmissiveSchema>;

export const SubstanceSchema = {
  tSubstance: DataSpec('rgba8'),
  dSubstance: DefineSpec('boolean'),

  uSubstanceGridDim: ObjectUniformSpec('v3'),
  uSubstanceGridTexDim: ObjectUniformSpec('v2'),
  uSubstanceGridTransform: ObjectUniformSpec('v4'),
  tSubstanceGrid: TextureSpec('texture', 'rgba', 'ubyte', 'linear', undefined, {
    usedIf: ['dSubstanceType_volumeInstance'],
  }),
  /** the layer type (`LayerType`) */
  substanceType: ValueSpec('string'),
  /** the layer type's define (`getLayerTypeDefine`): instance and groupInstance share 'texture' */
  dSubstanceType: DefineSpec('string', ['texture', 'vertexInstance', 'volumeInstance']),
  uSubstanceStrength: MaterialUniformSpec('f'),
} as const;
export type SubstanceSchema = typeof SubstanceSchema;
export type SubstanceValues = Values<SubstanceSchema>;

export const ClippingSchema = {
  tClipping: DataSpec('r8'),
  dClipping: DefineSpec('boolean'),
  /** the clipping type (`ClippingType`); `uInstanceGranularity` selects the index */
  clippingType: ValueSpec('string'),
} as const;
export type ClippingSchema = typeof ClippingSchema;
export type ClippingValues = Values<ClippingSchema>;

export const WiggleSchema = {
  tWiggle: DataSpec('r8'),
  dWiggle: DefineSpec('boolean'),
  wiggleAverage: ValueSpec('number'),
  /** the wiggle type (`WiggleType`); `uInstanceGranularity` selects the index */
  wiggleType: ValueSpec('string'),
  uWiggleStrength: MaterialUniformSpec('f'),
} as const;
export type WiggleSchema = typeof WiggleSchema;
export type WiggleValues = Values<WiggleSchema>;

export const BaseSchema = {
  dGeometryType: DefineSpec('string', [
    'cylinders',
    'directVolume',
    'image',
    'lines',
    'mesh',
    'points',
    'spheres',
    'text',
    'textureMesh',
  ]),

  ...ColorSchema,
  ...MarkerSchema,
  ...OverpaintSchema,
  ...TransparencySchema,
  ...EmissiveSchema,
  ...SubstanceSchema,
  ...ClippingSchema,
  ...WiggleSchema,

  /** whether there are clip objects */
  dClipObjects: DefineSpec('boolean'),
  /**
   * the clip object count, the bound of the clip loops: on WebGPU a pipeline constant (`oClipObjectCount`, the number
   * define), unrolled (a uniform bound cost 6-23% of clipped fragment-bound frames there); on WebGL the uniform
   * (`dClipObjectCount` ignored: one program for all counts; a define on WebGL1, `loopBound`), its loops faster so
   */
  dClipObjectCount: DefineSpec('number'),
  uClipObjectCount: { ...MaterialUniformSpec('i'), loopBound: true },
  dClipVariant: DefineSpec('string', ['instance', 'pixel']),
  /** the clip objects as given (`Clip.Objects`), for the CPU (cap culling); the shaders read frame and shape */
  clipObjects: ValueSpec('clipObjects'),
  /** per clip object the transform of a point into its frame (`Clip.Objects.frame`) */
  uClipObjectFrame: MaterialUniformSpec('m4[]', Clip.MaxObjectCount, 'uClipObjectCount'),
  /** per clip object its half size, the type in `w`, negative where inverted (`Clip.Objects.shape`) */
  uClipObjectShape: MaterialUniformSpec('v4[]', Clip.MaxObjectCount, 'uClipObjectCount'),

  aInstance: AttributeSpec('float32', 1, 1),
  /**
   * final per-instance transform calculated for instance `i` as
   * `aTransform[i] = matrix * transform[i] * extraTransform[i]`
   */
  aTransform: AttributeSpec('float32', 16, 1),
  /** vertex index, emulates `gl_VertexID`; the render item fills it when used */
  aVertex: AttributeSpec('float32', 1, 0, ['emulateVertexId']),

  /** final alpha, the base alpha clamped to [0, 1] (set with it) */
  uAlpha: MaterialUniformSpec('f'),
  uMetalness: MaterialUniformSpec('f'),
  uRoughness: MaterialUniformSpec('f'),
  uBumpiness: MaterialUniformSpec('f'),
  uEmissive: MaterialUniformSpec('f'),

  /** density value to estimate object thickness */
  uDensity: MaterialUniformSpec('f'),

  uVertexCount: ObjectUniformSpec('i'),
  uInstanceCount: ObjectUniformSpec('i'),
  uGroupCount: ObjectUniformSpec('i'),
  uInvariantBoundingSphere: ObjectUniformSpec('v4'),
  /** overridden per lod level via `MultiDrawBaseData.uniforms` */
  uLod: DrawUniformSpec('v4'),

  drawCount: ValueSpec('number'),
  instanceCount: ValueSpec('number'),
  /**
   * Vertex count from a GPU buffer (WebGPU), drawn instead of `drawCount` without multi-draws; `drawCount` is then an
   * upper bound (multi-draws draw it, beyond the vertex count the geometry must not be drawn, e.g. clipped)
   */
  drawIndirect: ValueSpec('drawIndirect'),

  /** base alpha, see uAlpha  */
  alpha: ValueSpec('number'),

  /** global transform, see aTransform */
  matrix: ValueSpec('m4'),
  /** base per-instance transform, see aTransform */
  transform: ValueSpec('float32'),
  /** additional per-instance transform, see aTransform */
  extraTransform: ValueSpec('float32'),
  /** whether extraTransform is used, see aTransform */
  hasExtraTransform: ValueSpec('boolean'),
  /** denotes reflection in transform */
  hasReflection: ValueSpec('boolean'),
  /**
   * instance granularity: the marker and layer data per instance, else per group instance; selects their index at
   * runtime, one program for both
   */
  uInstanceGranularity: ObjectUniformSpec('b'),

  /** bounding sphere taking aTransform into account and encompases all instances */
  boundingSphere: ValueSpec('sphere'),
  /** bounding sphere NOT taking aTransform into account */
  invariantBoundingSphere: ValueSpec('sphere'),

  instanceGrid: ValueSpec('instanceGrid'),
} as const;
export type BaseSchema = typeof BaseSchema;
export type BaseValues = Values<BaseSchema>;

//

export const InteriorSchema = {
  uInteriorColor: ObjectUniformSpec('v4'),
  uInteriorSubstance: ObjectUniformSpec('v4'),
} as const;
export type InteriorSchema = typeof InteriorSchema;
export type InteriorValues = Values<InteriorSchema>;

export const AnimationSchema = {
  uWiggleSpeed: MaterialUniformSpec('f'),
  uWiggleAmplitude: MaterialUniformSpec('f'),
  uWiggleFrequency: MaterialUniformSpec('f'),
  uWiggleMode: MaterialUniformSpec('i'),
  uTumbleSpeed: MaterialUniformSpec('f'),
  uTumbleAmplitude: MaterialUniformSpec('f'),
  uTumbleFrequency: MaterialUniformSpec('f'),
} as const;
export type AnimationSchema = typeof AnimationSchema;
export type AnimationValues = Values<AnimationSchema>;
