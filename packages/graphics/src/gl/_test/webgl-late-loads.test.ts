/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// @vitest-environment jsdom

import { describe, expect, it, vi } from 'vitest';
import type { CubeFaces } from '../texture.js';
import { createTextureImage } from '../renderable/util.js';
import { createContext } from '../webgl/context.js';
import { createCubeTexture } from '../webgl/texture.js';
import { createGl } from './gl.shim.js';

// What arrives after a WebGL texture was destroyed (an image decoded, a skybox face fetched): uploaded nowhere, as the
// WebGPU textures already do; binding a deleted texture fails, so the upload would go to whatever texture is bound.

/** The texture uploads `gl` makes from now on */
function countUploads(gl: WebGLRenderingContext) {
  let count = 0;
  for (const name of ['texImage2D', 'texSubImage2D', 'texImage3D', 'texSubImage3D'] as const) {
    const f = (gl as any)[name] as ((...args: any[]) => any) | undefined;
    if (!f) continue;
    (gl as any)[name] = (...args: any[]) => {
      count += 1;
      return f.apply(gl, args);
    };
  }
  return () => count;
}

describe('WebGL textures, a load after destroy', () => {
  it('uploads nothing: an image, data, rows or mipmaps', () => {
    const gl = createGl(32, 32, {});
    const ctx = createContext(gl);
    const texture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    const image = createTextureImage(4, 4, Uint8Array);
    texture.load(image);
    texture.destroy();

    const uploads = countUploads(gl);
    texture.load(image);
    texture.loadRows?.(image, [0, 2]);
    texture.mipmap();
    expect(uploads()).toBe(0);
    ctx.destroy();
  });

  it('settles a skybox whose face fails to load, once, as an error', async () => {
    const gl = createGl(32, 32, {});
    const onload = vi.fn();
    let fail: (e: Error) => void = () => {};
    const failing = new Promise<Blob>((_, reject) => (fail = reject));
    const faces: CubeFaces = {
      nx: failing,
      ny: failing,
      nz: new Promise<Blob>(() => {}),
      px: new Promise<Blob>(() => {}),
      py: new Promise<Blob>(() => {}),
      pz: new Promise<Blob>(() => {}),
    };
    const texture = createCubeTexture(gl, faces, false, onload);
    fail(new Error('404'));
    await Promise.resolve();
    await Promise.resolve();
    // two faces failed: reported once, as an error
    expect(onload).toHaveBeenCalledTimes(1);
    expect(onload).toHaveBeenCalledWith(true);
    texture.destroy();
  });
});
