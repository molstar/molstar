/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export type GraphicsBackend = 'webgl' | 'webgl2' | 'webgpu';

/**
 * Backend-neutral optional features, named by what they enable (not by extension).
 * Features required by Mol* (e.g. instancing, float textures, min/max blending) are not listed.
 */
export const GraphicsFeatures = [
  /** sampling 3D textures and rendering into their slices */
  'texture3d',
  /** `flat` varyings */
  'flatInterpolation',
  /** `int` textures: sampling, rendering and readback */
  'integerTextures',
  /** color attachments of one framebuffer with different formats or types, and `rg` formats */
  'mixedAttachmentFormats',
  'floatTextureLinear',
  'halfFloatTexture',
  'halfFloatTextureLinear',
  /** rendering into float32 color targets, without blending */
  'floatRenderTarget',
  /** rendering into float32 color targets with blending (e.g. accumulation, min/max) */
  'floatBlend',
  /** rendering into float16 color targets, with blending */
  'halfFloatRenderTarget',
  'depthTexture',
  'fragDepth',
  'multipleRenderTargets',
  'textureLod',
  'multiDraw',
  'clipControl',
  'depthClamp',
  /** compute shaders (`GraphicsContext.compute`): kernels, storage buffers and storage textures */
  'compute',
  /**
   * indirect draws of instance ranges: a first instance other than 0 in the arguments (`MultiDrawBaseData.indirect`,
   * GPU culling)
   */
  'indirectInstanceRanges',
] as const;
export type GraphicsFeature = (typeof GraphicsFeatures)[number];

export type GraphicsCapabilities = { readonly [K in GraphicsFeature]: boolean } & {
  /** for reports only, never check features with it */
  readonly backend: GraphicsBackend;
};
