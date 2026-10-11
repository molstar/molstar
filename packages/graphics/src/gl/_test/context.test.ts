/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

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
const foreignGlobals: GlobalBindings = { textures: [], invalidate: () => {}, destroy: () => {} };

describe('graphics context', () => {
  it('isWebGLContext narrows only WebGL contexts', () => {
    expect(isWebGLContext({} as GraphicsContext)).toBe(false);

    const ctx = tryGetGLContext(32, 32);
    if (!ctx) return;
    expect(isWebGLContext(ctx)).toBe(true);
    ctx.destroy();
  });

  it('isGLTexture narrows only textures of the WebGL resources', () => {
    expect(isGLTexture(foreignTexture)).toBe(false);
    expect(isGLTexture(createNullTexture())).toBe(false);
    expect(() => asGLTexture(foreignTexture)).toThrow(/another backend/);

    const ctx = tryGetGLContext(32, 32);
    if (!ctx) return;
    const texture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    expect(isGLTexture(texture)).toBe(true);
    expect(asGLTexture(texture)).toBe(texture);
    expect(() => ctx.resources.getValueChannel(foreignTexture)).toThrow(/another backend/);
    ctx.destroy();
  });

  it('isWebGLGlobalBindings narrows only global bindings of the WebGL context', () => {
    expect(isWebGLGlobalBindings(foreignGlobals)).toBe(false);
    expect(() => asWebGLGlobalBindings(foreignGlobals)).toThrow(/another backend/);

    const ctx = tryGetGLContext(32, 32);
    if (!ctx) return;
    const globals = ctx.createGlobalBindings(GlobalUniformSchema, [], []);
    expect(isWebGLGlobalBindings(globals)).toBe(true);
    globals.destroy();
    ctx.destroy();
  });
});
