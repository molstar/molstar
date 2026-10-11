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
    'floatRenderTarget',
    'halfFloatRenderTarget',
    'depthTexture',
    'fragDepth',
    'multipleRenderTargets',
    'textureLod',
    'multiDraw',
    'clipControl',
    'depthClamp',
] as const;
export type GraphicsFeature = (typeof GraphicsFeatures)[number];

export type GraphicsCapabilities = { readonly [K in GraphicsFeature]: boolean } & {
    /** for reports only, never check features with it */
    readonly backend: GraphicsBackend;
};
