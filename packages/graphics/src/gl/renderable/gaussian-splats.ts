/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { GraphicsContext } from '../context.js';
import type { Transparency } from '../render-item.js';
import {
  GlobalUniformSchema,
  BaseSchema,
  UniformSpec,
  DataSpec,
  DefineSpec,
  type Values,
  InternalSchema,
  type InternalValues,
  GlobalTextureSchema,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  TextureSpec,
  ValueSpec,
} from './schema.js';
import { GaussianSplatsShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';
import { StorageBufferSpec } from '../compute.js';
import { disposeSplatOrder } from './gaussian-splats-order.js';

/** A SOG image read in place (`sog` encoding): RGBA bytes, a pixel a splat (`GaussianSplats.Sog`) */
const SogImageSpec = (usedIf: string[] = ['dSplatEncoding_sog']) =>
  TextureSpec('texture', 'rgba', 'ubyte', 'nearest', undefined, { usedIf });

export const GaussianSplatsSchema = {
  ...BaseSchema,
  /**
   * `packed`: the splats in the layouts below (`tPositionGroup`, `tSplat`, `tSplatSh`); `sog`: SOG's images read in
   * place (`tSog*`), the packed data a placeholder each
   */
  dSplatEncoding: DefineSpec('string', ['packed', 'sog']),
  /** the centers on the CPU, for the sorts (the GPU's `tPositionGroup` a placeholder in the `sog` encoding) */
  splatCenters: ValueSpec('any'),
  tPositionGroup: DataSpec('rgb32f-i32', ['!dSplatEncoding_sog']),
  /** rotation, then log-scales and opacity: two elements a splat (`geo/geometry/gaussian-splats/splat-data.ts`) */
  tSplat: DataSpec('rgba16i', ['!dSplatEncoding_sog']),
  /** the splat of each draw slot, per instance a block, back to front (`gaussian-splats-order.ts`) */
  tSplatOrder: DataSpec('rgba8'),
  /** the view-dependent color: per splat a texel per basis function past the zeroth (`SplatSh`) */
  tSplatSh: DataSpec('rgba8', ['!dSplatEncoding_sog']),
  /** its degree (0: none) and the magnitude of its largest coefficient */
  uSplatShDegree: UniformSpec('i'),
  uSplatShScale: UniformSpec('f'),
  /** SOG's images (`sog`): the means' low and high bytes, scales, rotations, base color and opacity, all `uSogImageDim` */
  tSogMeansLow: SogImageSpec(),
  tSogMeansHigh: SogImageSpec(),
  tSogScales: SogImageSpec(),
  tSogQuats: SogImageSpec(),
  tSogSh0: SogImageSpec(),
  uSogImageDim: UniformSpec('v2'),
  /** the means' ranges, in the log domain */
  uSogMeansMin: UniformSpec('v3'),
  uSogMeansMax: UniformSpec('v3'),
  /** the codebooks of the scales (r), the base color (g) and the view-dependent color (b): 256 × 1 */
  tSogCodebook: TextureSpec('image-float32', 'rgba', 'float', 'nearest', undefined, { usedIf: ['dSplatEncoding_sog'] }),
  /** version 3's rotations (an exponential map) and opacities (`(a / 204)²`), else version 2's */
  dSogVersion3: DefineSpec('boolean'),
  /** the palette of view-dependent color (`shN`): a label a splat, the entries of `uSogShCoefficients` pixels */
  dSogSh: DefineSpec('boolean'),
  tSogShLabels: SogImageSpec(['dSplatEncoding_sog&dSogSh']),
  tSogShCentroids: SogImageSpec(['dSplatEncoding_sog&dSogSh']),
  uSogShCoefficients: UniformSpec('i'),
  uSogShCount: UniformSpec('i'),
  /** the centroids image's size (64 entries a row) */
  uSogShDim: UniformSpec('v2'),
  /** the order sorted on the GPU (WebGPU: `gl/compute/splat-sort.ts`), read in place of `tSplatOrder` */
  dSplatGpuOrder: DefineSpec('boolean'),
  sSplatOrder: StorageBufferSpec('u32', 'read', ['dSplatGpuOrder']),
  /** scales every splat's standard deviations */
  uSizeFactor: UniformSpec('f'),
  /** the quad's extent in standard deviations */
  uKernelSize: UniformSpec('f'),
  /** fragments of a lower coverage (the Gaussian's falloff times the opacity) discarded */
  uSplatAlphaCutoff: UniformSpec('f'),
};
export type GaussianSplatsSchema = typeof GaussianSplatsSchema;
export type GaussianSplatsValues = Values<GaussianSplatsSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const GaussianSplatsRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...GaussianSplatsSchema,
};

export function GaussianSplatsRenderable(
  ctx: GraphicsContext,
  id: number,
  values: GaussianSplatsValues,
  state: RenderableState,
  materialId: number,
  _transparency: Transparency,
  globals: GlobalDefines,
): Renderable<GaussianSplatsValues> {
  const schema = GaussianSplatsRenderableSchema;
  const renderValues: GaussianSplatsValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  // each splat is a quad of two triangles, see the gaussian splats vertex shader; blended whatever the scene's
  // transparency (`_transparency` unused): their own pass blends them back to front over the other transparency
  const renderItem = ctx.createGraphicsRenderItem(
    'triangles',
    GaussianSplatsShaderCode,
    schema,
    renderValues,
    materialId,
    BlendedTransparency,
  );
  const renderable = createRenderable(renderItem, renderValues, state);
  renderable.setTransparency = () => {};
  // the GPU sort's buffers are the order's (`updateSplatOrder`), not the render item's
  const dispose = renderable.dispose;
  renderable.dispose = () => {
    dispose();
    disposeSplatOrder(values);
  };
  return renderable;
}

const BlendedTransparency: Transparency = 'blended';
