/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Geometry, GeometryUtils } from '@molstar/graphics/geo/geometry/geometry';
import { type ColorType, getColorTypeDefine } from '@molstar/graphics/geo/geometry/color-data';
import { getSizeTypeDefine, type SizeType } from '@molstar/graphics/geo/geometry/size-data';
import { getRenderableProgramSource, type RenderObjectType } from '@molstar/graphics/gl/render-object';
import { isComputeIsosurfaceSupported } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import { isMergeableType } from '@molstar/graphics/gl/renderable/merged';
import {
  ContextDefines,
  getSchemaLinearTextures,
  getSignatureKey,
  mergeSignatures,
  type ObjectSignature,
} from '@molstar/graphics/gl/precompile';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import { getQualityProps } from './util.js';
import type { VisualQuality } from '@molstar/graphics/geo/geometry/base';

type DefineValue = string | number | boolean;

/** Values a define may take, each combination with the other defines' values gives a signature */
export type SignatureDefines = { readonly [define: string]: readonly DefineValue[] };

/**
 * What the signatures of a representation's objects depend on besides its props, without data: the context the objects
 * are drawn in (capabilities choose visuals, e.g. impostors or GPU surfaces), the themes and the layers set on the
 * representation (overpaint, transparency, ...)
 */
export interface SignatureContext {
  readonly graphics?: GraphicsContext;
  /** color types the color theme may give (its granularity, see `getColorThemeSignature`) */
  readonly color: {
    readonly types: readonly ColorType[];
    readonly palette: boolean;
    readonly preferSmoothing: boolean;
  };
  /** size types the size theme may give */
  readonly size: { readonly types: readonly SizeType[] };
  readonly layers: {
    readonly overpaint: boolean;
    readonly transparency: boolean;
    readonly emissive: boolean;
    readonly substance: boolean;
    readonly clipping: boolean;
    readonly wiggle: boolean;
  };
}

export const NoSignatureLayers: SignatureContext['layers'] = {
  overpaint: false,
  transparency: false,
  emissive: false,
  substance: false,
  clipping: false,
  wiggle: false,
};

/** The color types of a granularity: 'instance' also as 'group' (objects that can't be instanced, see `createColors`) */
export function getColorTypes(granularity: ColorType): ColorType[] {
  return granularity === 'instance' ? ['instance', 'group'] : [granularity];
}

/** The color types a theme may give, see `getColorTypes` */
export function getColorThemeSignature(theme: ColorTheme<any, any>): SignatureContext['color'] {
  return {
    types: getColorTypes(theme.granularity),
    palette: !!theme.palette,
    preferSmoothing: !!theme.preferSmoothing,
  };
}

export function getSizeThemeSignature(theme: SizeTheme<any>): SignatureContext['size'] {
  return { types: [theme.granularity] };
}

/**
 * The color type smoothing gives for `type`: with a context a grid, else colors interpolated per vertex (see
 * `calcMeshColorSmoothing`)
 */
function getSmoothedColorType(type: ColorType, graphics?: GraphicsContext): ColorType | undefined {
  if (type === 'group') return graphics ? 'volume' : 'vertex';
  if (type === 'groupInstance') return graphics ? 'volumeInstance' : 'vertexInstance';
}

/** Whether colors of `props` may be smoothed (the resolution, which also decides, is data) */
function maySmooth(props: PD.Values<any>, ctx: SignatureContext) {
  const smoothColors = props.smoothColors as { name: string } | undefined;
  if (!smoothColors) return false;
  return smoothColors.name === 'on' || (smoothColors.name === 'auto' && ctx.color.preferSmoothing);
}

/** Every combination of the alternatives' values with `base` */
function expand(base: { [k: string]: DefineValue }, alternatives: [string, readonly DefineValue[]][]) {
  let result = [base];
  for (const [k, values] of alternatives) {
    const next: { [k: string]: DefineValue }[] = [];
    for (const r of result) for (const v of values) next.push({ ...r, [k]: v });
    result = next;
  }
  return result;
}

/**
 * The signatures of the objects of geometry `utils`' kind that `props` give in `ctx`: the define values of an empty
 * geometry with the props, and for what data decides (color, size and marker types, the layers' presence and types,
 * smoothing) every value it may take; `visualDefines` adds what a visual decides (e.g. varying groups)
 */
export function getGeometrySignatures<G extends Geometry>(
  utils: GeometryUtils<G>,
  props: PD.Values<any>,
  ctx: SignatureContext,
  visualDefines?: SignatureDefines,
): ObjectSignature[] {
  const geometry = utils.createEmpty();
  const kind: RenderObjectType = geometry.kind;
  const source = getRenderableProgramSource({ type: kind, merged: false });
  if (!source) return [];
  const { schema } = source;
  const values = utils.createValuesSimple(geometry, props, Color(0), 1) as unknown as {
    [k: string]: { ref: { value: unknown } } | undefined;
  };
  const base: { [k: string]: DefineValue } = {};
  for (const k of Object.keys(schema)) {
    if (schema[k].type !== 'define' || k === 'dRenderVariant' || ContextDefines.has(k)) continue;
    const value = values[k]?.ref.value;
    if (typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean') base[k] = value;
  }

  const smooth = maySmooth(props, ctx) && (kind === 'mesh' || kind === 'texture-mesh');
  // as define values: the types read by index share one (`getColorTypeDefine`)
  const colorTypes = new Set<DefineValue>(ctx.color.types.map(getColorTypeDefine));
  if (smooth) {
    for (const t of ctx.color.types) {
      const s = getSmoothedColorType(t, ctx.graphics);
      if (s) colorTypes.add(getColorTypeDefine(s));
    }
  }
  // instance and groupInstance share a define (`getLayerTypeDefine`), the instance granularity is a uniform
  const layerTypes = smooth ? ['texture', ctx.graphics ? 'volumeInstance' : 'vertexInstance'] : ['texture'];
  const { layers } = ctx;
  const alternatives: [string, readonly DefineValue[]][] = [
    ['dColorType', [...colorTypes]],
    ['dUsePalette', [ctx.color.palette]],
    ['dSizeType', [...new Set(ctx.size.types.map(getSizeTypeDefine))]],
    ['dOverpaint', [layers.overpaint]],
    ['dTransparency', [layers.transparency]],
    ['dEmissive', [layers.emissive]],
    ['dSubstance', [layers.substance]],
    ['dClipping', [layers.clipping]],
    ['dWiggle', [layers.wiggle]],
  ];
  // a layer's type is that of its default (empty) data unless it is set
  if (layers.overpaint) alternatives.push(['dOverpaintType', layerTypes]);
  if (layers.transparency) alternatives.push(['dTransparencyType', layerTypes]);
  if (layers.emissive) alternatives.push(['dEmissiveType', layerTypes]);
  if (layers.substance) alternatives.push(['dSubstanceType', layerTypes]);
  // grids are 3d textures where supported (`createDirectVolume`, the Gaussian density)
  if (kind === 'direct-volume')
    alternatives.push(['dGridTexType', [ctx.graphics?.capabilities.texture3d ? '3d' : '2d']]);
  // texture meshes are extracted into chunks where compute isosurfaces are supported (`TextureMesh.createChunked`)
  if (kind === 'texture-mesh') {
    alternatives.push(['dIsoChunks', [!!ctx.graphics && isComputeIsosurfaceSupported(ctx.graphics)]]);
  }
  if (visualDefines) for (const k of Object.keys(visualDefines)) alternatives.push([k, visualDefines[k]]);

  const doubleSided = !!values.uDoubleSided?.ref.value;
  // how the objects may draw: by the props' alpha and emission, any way the layers allow
  const alpha = (values.alpha?.ref.value as number | undefined) ?? 1;
  const opaque = alpha === 1;
  const transparent = alpha < 1 || layers.transparency;
  const emissive = ((values.uEmissive?.ref.value as number | undefined) ?? 0) > 0 || layers.emissive;
  const signatures: ObjectSignature[] = [];
  const present = alternatives.filter(([k]) => schema[k]?.type === 'define');
  for (const defines of expand(base, present)) {
    signatures.push({
      type: kind,
      merged: false,
      defines,
      linearTextures: getSchemaLinearTextures(schema, defines),
      doubleSided,
      hasReflection: false,
      opaque,
      transparent,
      emissive,
    });
  }
  return signatures;
}

/**
 * `signatures` and their merged variants where `props` merge the objects: `renderMerged` on, or `gpuCulling` on (an
 * object alone merged to be culled on the GPU, `gpuCulled`); not with a palette (`canMergeValues`). 'auto' (either) is
 * taken as unmerged and culled on the CPU: by the objects' counts, large objects that compile while they build
 */
export function withMergedSignatures(
  signatures: readonly ObjectSignature[],
  props: { readonly renderMerged?: unknown; readonly gpuCulling?: unknown },
): ObjectSignature[] {
  const gpuCulled = props.gpuCulling === true;
  if (props.renderMerged !== true && !gpuCulled) return [...signatures];
  const out = [...signatures];
  for (const s of signatures) {
    if (!isMergeableType(s.type) || s.defines.dUsePalette) continue;
    const source = getRenderableProgramSource({ type: s.type, merged: true });
    if (!source) continue;
    const defines = { ...s.defines, dSegmented: true };
    const linearTextures = getSchemaLinearTextures(source.schema, defines);
    out.push({ ...s, merged: true, defines, linearTextures, gpuCulled });
  }
  return out;
}

/** `signatures` without duplicates, the draw flags of those of the same key merged (`mergeSignatures`) */
export function uniqueSignatures(signatures: Iterable<ObjectSignature>): ObjectSignature[] {
  const map = new Map<string, ObjectSignature>();
  for (const s of signatures) {
    const key = getSignatureKey(s);
    const known = map.get(key);
    map.set(key, known ? mergeSignatures(known, s) : s);
  }
  return [...map.values()];
}

/**
 * Surfaces keep groups varying across triangles where flat interpolation is supported, else they make them uniform
 * (`Mesh.uniformTriangleGroup`)
 */
export function getSurfaceSignatureDefines(_props: unknown, graphics?: GraphicsContext): SignatureDefines {
  return { dVaryingGroup: [!graphics || graphics.capabilities.flatInterpolation] };
}

/**
 * The signatures `get` gives for the quality props `props` resolve to (`getQualityProps`, as representations apply
 * them): 'auto' quality, which data decides, as medium (low only for large structures, which compile while they build;
 * the qualities differ in double-sidedness)
 */
export function getQualitySignatures<P extends { quality?: VisualQuality }>(
  props: P,
  get: (props: P) => ObjectSignature[],
): ObjectSignature[] {
  const quality: VisualQuality = !props.quality || props.quality === 'auto' ? 'medium' : props.quality;
  return get({ ...props, ...getQualityProps({ ...props, quality }) });
}
