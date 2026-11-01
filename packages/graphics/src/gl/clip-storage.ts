/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util/value-cell';

/** Vec4s per clip object in storage (`clip_storage_wgsl`): its frame's columns, its shape */
export const ClipStorageVec4s = 5;

/** The clip values of an object a kernel reads (`BaseSchema`'s, `dClipPrimitive` the spheres') */
export interface ClipStorageValues {
  readonly dClipObjects: ValueCell<boolean>;
  readonly dClipVariant: ValueCell<string>;
  readonly dClipping: ValueCell<boolean>;
  readonly uClipObjectCount: ValueCell<number>;
  readonly uClipObjectFrame: ValueCell<ArrayLike<number>>;
  readonly uClipObjectShape: ValueCell<ArrayLike<number>>;
}

/** The clip objects of an object, for the kernels (`packStorageClipObjects`), 0 without */
export function getStorageClipCount(values: ClipStorageValues) {
  return values.dClipObjects.ref.value ? values.uClipObjectCount.ref.value : 0;
}

/**
 * Those of instance clipping (the vertex stage's, `clip_instance_wgsl`): an instance's center tested; none with a
 * clipping layer (its flags per group exempt the spheres of an instance apart, not the instance)
 */
export function getInstanceClipCount(values: ClipStorageValues) {
  if (values.dClipVariant.ref.value !== 'instance' || values.dClipping.ref.value) return 0;
  return getStorageClipCount(values);
}

/**
 * Those of pixel clipping (the fragment stage's, `clip_pixel_wgsl`, or a primitive's center with `dClipPrimitive`), a
 * clipping layer's flags read per sphere
 */
export function getPixelClipCount(values: ClipStorageValues) {
  return values.dClipVariant.ref.value === 'pixel' ? getStorageClipCount(values) : 0;
}

/**
 * Those culling spheres they clip whole (`clipSphereClass`, `clip-culling.md`), either variant: none with a clipping
 * layer (its flags exempt groups)
 */
export function getCullClipCount(values: ClipStorageValues) {
  return values.dClipping.ref.value ? 0 : getStorageClipCount(values);
}

/** Whether the kernels clip `values` as the programs do: not instance clipping with a clipping layer */
export function isStorageClipSupported(values: ClipStorageValues) {
  return !values.dClipObjects.ref.value || values.dClipVariant.ref.value !== 'instance' || !values.dClipping.ref.value;
}

/** The versions of what `packStorageClipObjects` reads: its result is the same while they are */
export function getStorageClipVersion(values: ClipStorageValues) {
  return `${getStorageClipCount(values)},${values.uClipObjectFrame.ref.version},${values.uClipObjectShape.ref.version}`;
}

/** The clip objects of `values` (`getStorageClipCount` of them) into `out` from vec4 `base`, `ClipStorageVec4s` each */
export function packStorageClipObjects(values: ClipStorageValues, out: Float32Array, base = 0) {
  const count = getStorageClipCount(values);
  const frame = values.uClipObjectFrame.ref.value;
  const shape = values.uClipObjectShape.ref.value;
  for (let i = 0; i < count; ++i) {
    const o = (base + i * ClipStorageVec4s) * 4;
    for (let k = 0; k < 16; ++k) out[o + k] = frame[i * 16 + k];
    for (let k = 0; k < 4; ++k) out[o + 16 + k] = shape[i * 4 + k];
  }
  return count;
}
