/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TextureUsage } from './constants.js';

/**
 * Whether the device takes a viewport partly outside the pass's attachments, as the spec allows (`setViewport`: down to
 * `-2 * maxTextureDimension2D`). Safari 26 (iOS) validates the older rule, the viewport inside the attachments: the
 * `setViewport` invalidates the pass encoder, so the whole pass is lost. An encoder setting a viewport starting at -1
 * on a 1x1 target, its validation error caught; false where the probe fails.
 */
export async function testViewportOutsideTarget(device: GPUDevice): Promise<boolean> {
  let color: GPUTexture | undefined;
  try {
    color = device.createTexture({ size: [1, 1], format: 'rgba8unorm', usage: TextureUsage.RENDER_ATTACHMENT });
    device.pushErrorScope('validation');
    const encoder = device.createCommandEncoder();
    const pass = encoder.beginRenderPass({
      colorAttachments: [{ view: color.createView(), loadOp: 'clear', storeOp: 'store' }],
    });
    pass.setViewport(-1, 0, 2, 1, 0, 1);
    pass.end();
    // not submitted: the invalid encoder errors on `finish`
    encoder.finish();
    return (await device.popErrorScope()) === null;
  } catch {
    return false;
  } finally {
    color?.destroy();
  }
}
