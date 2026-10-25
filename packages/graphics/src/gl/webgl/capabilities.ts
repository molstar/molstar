/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsCapabilities } from '../capabilities.js';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import type { WebGLExtensions } from './extensions.js';

export function getWebGLCapabilities(gl: GLRenderingContext, extensions: WebGLExtensions): GraphicsCapabilities {
  const webgl2 = isWebGL2(gl);
  const e = extensions;
  return {
    backend: webgl2 ? 'webgl2' : 'webgl',
    texture3d: webgl2,
    flatInterpolation: webgl2,
    integerTextures: webgl2,
    mixedAttachmentFormats: webgl2,
    floatTextureLinear: !!e.textureFloatLinear,
    halfFloatTexture: !!e.textureHalfFloat,
    halfFloatTextureLinear: !!e.textureHalfFloat && !!e.textureHalfFloatLinear,
    floatRenderTarget: !!e.colorBufferFloat,
    floatBlend: !!e.colorBufferFloat && !!e.floatBlend,
    halfFloatRenderTarget: !!e.textureHalfFloat && !!e.colorBufferHalfFloat,
    depthTexture: !!e.depthTexture,
    fragDepth: !!e.fragDepth,
    multipleRenderTargets: !!e.drawBuffers,
    textureLod: !!e.shaderTextureLod,
    multiDraw: !!e.multiDraw,
    clipControl: !!e.clipControl,
    depthClamp: !!e.depthClamp,
    indirectInstanceRanges: false,
    compute: false,
  };
}
