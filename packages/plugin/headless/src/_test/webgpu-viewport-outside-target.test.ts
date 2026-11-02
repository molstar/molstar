/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { CopySchema, CopyShaderCode, QuadValues } from '@molstar/graphics/gl/compute/util';
import { createComputeRenderable } from '@molstar/graphics/gl/renderable';
import type { Values } from '@molstar/graphics/gl/renderable/schema';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { testViewportOutsideTarget } from '@molstar/graphics/gl/webgpu/viewport-outside-target';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 32;

/**
 * Draws a white quad over the bottom-left quarter of a viewport partly outside the target, forcing whether the context
 * sets such a viewport as given; the red channel of the target, rows top to bottom
 */
async function draw(viewportOutsideTarget: boolean, viewport?: (graphics: WebGPUContext) => number[]) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu', webgpu: { viewportOutsideTarget } },
  );
  try {
    await plugin.init();
    const graphics = plugin.canvas3d!.graphics as WebGPUContext;
    const white = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    white.load({ array: new Uint8Array(Size * Size * 4).fill(255), width: Size, height: Size });
    const target = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    target.define(Size, Size);
    const framebuffer = graphics.resources.framebuffer();
    target.attachFramebuffer(framebuffer, 'color0');
    const values: Values<typeof CopySchema> = {
      ...QuadValues,
      uQuadScale: ValueCell.create(Vec2.create(0.5, 0.5)),
      tColor: ValueCell.create(white),
      uTexSize: ValueCell.create(Vec2.create(Size, Size)),
      dPackDepth: ValueCell.create(false),
    };
    const copy = createComputeRenderable(
      graphics.createComputeRenderItem('triangles', CopyShaderCode, { ...CopySchema }, values),
      values,
    );
    const pass = () => {
      graphics.disableScissor();
      graphics.beginRenderPass(framebuffer, { label: 'viewport outside the target', clear: { color: [0, 0, 0, 1] } });
      graphics.setRenderState(RenderState.Default);
      // a viewport twice the target's size, starting left of and above it
      const [x, y, width, height] = viewport?.(graphics) ?? [-Size / 2, -Size, 2 * Size, 2 * Size];
      graphics.setViewport(x, y, width, height);
      copy.render();
      graphics.endRenderPass();
    };
    // the first draw requests the pipeline, which compiles in the background (draws skipped meanwhile)
    pass();
    await graphics.resources.pending.whenReady();
    pass();
    const array = await graphics.readTextureAsync(target, 0, 0, Size, Size);
    const errors = graphics.takeErrors();
    copy.dispose();
    framebuffer.destroy();
    target.destroy();
    white.destroy();
    return { red: (x: number, y: number) => array[(y * Size + x) * 4], errors };
  } finally {
    plugin.dispose();
  }
}

// Safari's WebGPU rejects a viewport outside the attachments (the pass encoder invalidated, the pass lost): the probe
// finds it, and the context clips such a viewport to the target there
describe.skipIf(!WebGPU.module)('viewports outside the target on WebGPU', () => {
  it('the probe finds them taken on Dawn', async () => {
    // `gpu` kept alive while its device is in use (collecting it crashes Dawn's Node.js bindings)
    const gpu = WebGPU.module!.create([]);
    const adapter = await gpu.requestAdapter();
    const device = await adapter!.requestDevice();
    try {
      expect(await testViewportOutsideTarget(device)).toBe(true);
    } finally {
      device.destroy();
      expect(gpu).toBeDefined();
    }
  });

  it('are set as given where taken (cropped), clipped to the target where not (squashed, the pass drawn)', async () => {
    // as given, the quad over the viewport's bottom-left quarter is the target's left half
    const taken = await draw(true);
    expect(taken.errors).toEqual([]);
    expect(taken.red(4, 4)).toBe(255);
    expect(taken.red(4, Size - 4)).toBe(255);
    expect(taken.red(Size - 4, 4)).toBe(0);
    // clipped to the target, the quad is over the target's bottom-left quarter
    const clipped = await draw(false);
    expect(clipped.errors).toEqual([]);
    expect(clipped.red(4, 4)).toBe(0);
    expect(clipped.red(4, Size - 4)).toBe(255);
    expect(clipped.red(Size - 4, Size - 4)).toBe(0);
  }, 60_000);

  it('are clipped beyond the spec bounds where taken (larger than maxTextureDimension2D)', async () => {
    // as given, the size would invalidate the pass on any device
    const oversized = await draw(true, (graphics) => {
      const max = graphics.device.limits.maxTextureDimension2D;
      return [-max / 2, -Size, max + 2, 2 * Size];
    });
    expect(oversized.errors).toEqual([]);
    expect(oversized.red(4, Size - 4)).toBe(255);
    expect(oversized.red(4, 4)).toBe(0);
  }, 60_000);
});
