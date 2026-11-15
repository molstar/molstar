/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import type { WebGLContext } from '../webgl/context.js';
import { tryGetGLContext } from './gl.js';

/** A WebGL context of the native backend (the attachments are tracked by GL object: the shim's are numbers) */
function createCase(ctx: WebGLContext) {
  const textures = vi.spyOn(ctx.gl, 'framebufferTexture2D');
  const renderbuffers = vi.spyOn(ctx.gl, 'framebufferRenderbuffer');
  const color = () => {
    const t = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    t.define(8, 8);
    return t;
  };
  return { textures, renderbuffers, color };
}

// passes attach their targets every frame: the same attachment again is no GL call (no end of the open pass, no
// framebuffer check by the driver); any change attaches
describe('WebGL framebuffer attachments', () => {
  const ctx = tryGetGLContext(32, 32);

  (ctx ? it : it.skip)('skips attaching what an attachment point holds already', { tags: ['medium'] }, () => {
    const { textures, color } = createCase(ctx!);
    const framebuffer = ctx!.resources.framebuffer();
    const a = color();
    const b = color();

    a.attachFramebuffer(framebuffer, 0);
    expect(textures).toHaveBeenCalledTimes(1);
    a.attachFramebuffer(framebuffer, 0);
    a.attachFramebuffer(framebuffer, 0);
    expect(textures).toHaveBeenCalledTimes(1);

    // another texture, another point, or the same after a detach: attached
    b.attachFramebuffer(framebuffer, 0);
    expect(textures).toHaveBeenCalledTimes(2);
    a.attachFramebuffer(framebuffer, 1);
    expect(textures).toHaveBeenCalledTimes(3);
    a.detachFramebuffer(framebuffer, 1);
    a.attachFramebuffer(framebuffer, 1);
    expect(textures).toHaveBeenCalledTimes(5);
    // per framebuffer
    a.attachFramebuffer(ctx!.resources.framebuffer(), 1);
    expect(textures).toHaveBeenCalledTimes(6);
    vi.restoreAllMocks();
  });

  (ctx ? it : it.skip)(
    'attaches a texture again where a renderbuffer took its point meanwhile',
    { tags: ['medium'] },
    () => {
      const { textures, renderbuffers } = createCase(ctx!);
      const framebuffer = ctx!.resources.framebuffer();
      const depthTexture = ctx!.resources.texture('image-depth', 'depth', 'ushort', 'nearest');
      depthTexture.define(8, 8);
      const depthBuffer = ctx!.resources.renderbuffer('depth16', 'depth', 8, 8);

      depthTexture.attachFramebuffer(framebuffer, 'depth');
      depthBuffer.attachFramebuffer(framebuffer);
      depthBuffer.attachFramebuffer(framebuffer);
      expect(renderbuffers).toHaveBeenCalledTimes(1);
      depthTexture.attachFramebuffer(framebuffer, 'depth');
      expect(textures).toHaveBeenCalledTimes(2);
      vi.restoreAllMocks();
    },
  );
});
