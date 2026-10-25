/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { GraphicsContext } from '../context.js';
import type { GlobalBindings } from '../render-item.js';
import { GlobalUniformSchema } from '../renderable/schema.js';
import { createNullTexture, type Texture } from '../texture.js';
import { isWebGLContext } from '../webgl/context.js';
import { asWebGLGlobalBindings, isWebGLGlobalBindings } from '../webgl/global-bindings.js';
import { asGLTexture, isGLTexture } from '../webgl/texture.js';
import { tryGetGLContext } from './gl.js';

/** stand-ins for objects of another backend */
const foreignTexture = { id: -1 } as Texture;
const foreignGlobals: GlobalBindings = {
  textures: [],
  invalidate: () => {},
  setSolidInteriorCap: () => {},
  destroy: () => {},
};

/** whether the optional `gl` package can create a context (`pnpm native:install`) */
const hasGl = (() => {
  const ctx = tryGetGLContext(1, 1);
  ctx?.destroy();
  return !!ctx;
})();

describe('graphics context', () => {
  it('isWebGLContext rejects other contexts', () => {
    expect(isWebGLContext({} as GraphicsContext)).toBe(false);
  });

  it.skipIf(!hasGl)('isWebGLContext accepts WebGL contexts', { tags: ['medium'] }, () => {
    const ctx = tryGetGLContext(32, 32)!;
    expect(isWebGLContext(ctx)).toBe(true);
    ctx.destroy();
  });

  it('isGLTexture rejects textures not of the WebGL resources', () => {
    expect(isGLTexture(foreignTexture)).toBe(false);
    expect(isGLTexture(createNullTexture())).toBe(false);
    expect(() => asGLTexture(foreignTexture)).toThrow(/another backend/);
  });

  it.skipIf(!hasGl)('isGLTexture accepts textures of the WebGL resources', { tags: ['medium'] }, () => {
    const ctx = tryGetGLContext(32, 32)!;
    const texture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    expect(isGLTexture(texture)).toBe(true);
    expect(asGLTexture(texture)).toBe(texture);
    ctx.destroy();
  });

  it('isWebGLGlobalBindings rejects global bindings of other contexts', () => {
    expect(isWebGLGlobalBindings(foreignGlobals)).toBe(false);
    expect(() => asWebGLGlobalBindings(foreignGlobals)).toThrow(/another backend/);
  });

  it.skipIf(!hasGl)('isWebGLGlobalBindings accepts global bindings of the WebGL context', { tags: ['medium'] }, () => {
    const ctx = tryGetGLContext(32, 32)!;
    const globals = ctx.createGlobalBindings(GlobalUniformSchema, [], []);
    expect(isWebGLGlobalBindings(globals)).toBe(true);
    globals.destroy();
    ctx.destroy();
  });
});
