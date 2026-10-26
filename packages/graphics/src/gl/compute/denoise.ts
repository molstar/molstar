/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import type { GraphicsContext } from '../context.js';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type KernelRequest,
  StorageTextureSpec,
} from '../compute.js';
import { TextureSpec, UniformSpec } from '../renderable/schema.js';
import type { Texture } from '../texture.js';

/**
 * The illumination denoiser as a kernel (WebGPU): tiles of 16 x 8 pixels, each loading its tile and a border of 6
 * texels once into workgroup memory (28 x 20 texels of color and normal, 15.7 KiB of the 16 KiB of the core limits)
 */
export const DenoiseKernel = KernelCode('illumination-denoise', [16, 8]);

export const DenoiseSchema = {
  uDenoiseThreshold: UniformSpec('f'),
  /** the accumulated color */
  tColor: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  /** octahedral-encoded normals, `packNormalToRGBA` */
  tNormal: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  outColor: StorageTextureSpec('rgba', 'fp16', '2d', 'write'),
};

/** The kernel of `Denoiser`, see `ComputeContext.prepareKernels` */
export function getDenoiseKernels(): KernelRequest[] {
  return [{ code: DenoiseKernel, schema: DenoiseSchema }];
}

export function isComputeDenoiseSupported(
  graphics: GraphicsContext,
): graphics is GraphicsContext & { compute: ComputeContext } {
  return !!graphics.compute?.isKernelSupported(DenoiseKernel, DenoiseSchema);
}

/** Denoises `color` (by `normal`) into `output`, an rgba fp16 texture of the same size */
export class Denoiser {
  private readonly kernel: ComputeKernel;
  private readonly values: {
    uDenoiseThreshold: ValueCell<number>;
    tColor: ValueCell<Texture>;
    tNormal: ValueCell<Texture>;
    outColor: ValueCell<Texture>;
  };

  constructor(
    private readonly compute: ComputeContext,
    color: Texture,
    normal: Texture,
    readonly output: Texture,
  ) {
    this.values = {
      uDenoiseThreshold: ValueCell.create(0.1),
      tColor: ValueCell.create(color),
      tNormal: ValueCell.create(normal),
      outColor: ValueCell.create(output),
    };
    this.kernel = compute.createKernel(DenoiseKernel, DenoiseSchema, this.values);
  }

  run(threshold: number) {
    ValueCell.updateIfChanged(this.values.uDenoiseThreshold, threshold);
    const [x, y] = DenoiseKernel.workgroupSize;
    this.compute.beginComputePass({ label: 'denoise' });
    this.kernel.dispatch(Math.ceil(this.output.getWidth() / x), Math.ceil(this.output.getHeight() / y));
    this.compute.endComputePass();
  }

  destroy() {
    this.kernel.destroy();
  }
}
