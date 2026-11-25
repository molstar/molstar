/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
  StorageTextureSpec,
} from '../compute.js';
import { TextureSpec, UniformSpec } from '../renderable/schema.js';
import type { Texture } from '../texture.js';

// The identity of a page store's picked instances (`paged-cull.ts`), gathered on the GPU in the pick's own submission:
// a page store's object picks its instance's composed pose slot (the frame's: a copy's instances composed each cull),
// what its record holds read before the pick's buffers are: its store index (the pick's instance) and, a copy's
// instance, the copy plus one (the cull's, above the record's segment: `PagedCopyShift`) into a texture of its own.
// Read from the slots the pick drew, not later (the next cull rewrites a frame's slots), no second readback.

export const PagedPickKernel = KernelCode('paged-pick', [8, 8]);

export const PagedPickSchema = {
  /** the pick's texels resolved: from `uOffset`, `uSize` of them */
  uOffset: UniformSpec('iv2'),
  uSize: UniformSpec('iv2'),
  /** the page store's object: its texels' instance its slot */
  uObjectId: UniformSpec('i'),
  /** 1 the first object of a pick: every other texel's instance as picked, its copy none */
  uFirst: UniformSpec('i'),
  tObject: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tInstance: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  /** the object's `CulledData` (`CulledPagedFlag`): its header, its segments, its composed poses */
  poses: StorageBufferSpec('vec4<u32>', 'read'),
  outInstance: StorageTextureSpec('rgba', 'ubyte', '2d', 'write'),
  outCopy: StorageTextureSpec('rgba', 'ubyte', '2d', 'write'),
};

export const PagedPickKernels: readonly KernelRequest[] = [{ code: PagedPickKernel, schema: PagedPickSchema }];

/** A page store's object a pick resolves: its id and its `CulledData` */
export interface PagedPickObject {
  readonly objectId: number;
  readonly poses: StorageBuffer;
}

export interface PagedPickResolver {
  /**
   * The pick's texels `[x, x + width) x [y, y + height)` resolved into `outInstance` (a page store's object's its
   * store index, others' as picked) and `outCopy` (a copy's instance its copy plus one, else none), in compute passes of
   * their own after the pick's; false where the kernel is not ready (compiling)
   */
  resolve(
    textures: { object: Texture; instance: Texture; outInstance: Texture; outCopy: Texture },
    viewport: { x: number; y: number; width: number; height: number },
    objects: readonly PagedPickObject[],
  ): boolean;
  destroy(): void;
}

export function createPagedPickResolver(compute: ComputeContext): PagedPickResolver {
  let kernel: ComputeKernel | undefined;
  const values = {
    uOffset: ValueCell.create(Vec2()),
    uSize: ValueCell.create(Vec2()),
    uObjectId: ValueCell.create(0),
    uFirst: ValueCell.create(1),
    tObject: ValueCell.create(undefined as unknown as Texture),
    tInstance: ValueCell.create(undefined as unknown as Texture),
    poses: ValueCell.create(undefined as unknown as StorageBuffer),
    outInstance: ValueCell.create(undefined as unknown as Texture),
    outCopy: ValueCell.create(undefined as unknown as Texture),
  };
  return {
    resolve(textures, viewport, objects) {
      if (objects.length === 0) return true;
      if (!compute.hasKernels(PagedPickKernels)) {
        void compute.prepareKernels(PagedPickKernels);
        return false;
      }
      ValueCell.update(values.uOffset, Vec2.set(values.uOffset.ref.value, viewport.x, viewport.y));
      ValueCell.update(values.uSize, Vec2.set(values.uSize.ref.value, viewport.width, viewport.height));
      ValueCell.updateIfChanged(values.tObject, textures.object);
      ValueCell.updateIfChanged(values.tInstance, textures.instance);
      ValueCell.updateIfChanged(values.outInstance, textures.outInstance);
      ValueCell.updateIfChanged(values.outCopy, textures.outCopy);
      kernel ??= compute.createKernel(PagedPickKernel, PagedPickSchema, values);
      const x = Math.ceil(viewport.width / 8);
      const y = Math.ceil(viewport.height / 8);
      for (let i = 0; i < objects.length; ++i) {
        ValueCell.updateIfChanged(values.uObjectId, objects[i].objectId);
        ValueCell.updateIfChanged(values.uFirst, i === 0 ? 1 : 0);
        ValueCell.updateIfChanged(values.poses, objects[i].poses);
        compute.beginComputePass({ label: 'paged pick' });
        kernel.dispatch(x, y);
        compute.endComputePass();
      }
      return true;
    },
    destroy() {
      kernel?.destroy();
      kernel = undefined;
    },
  };
}
