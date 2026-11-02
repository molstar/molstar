/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { createNullTexture, type Texture } from '../../gl/texture.js';

/**
 * The types of the layers' data (overpaint, transparency, emissive, substance): per instance or group instance (the
 * object's instance granularity), smoothed per vertex (without a graphics context) or into a grid
 */
export type LayerType = 'instance' | 'groupInstance' | 'vertexInstance' | 'volumeInstance';

/**
 * The define value of the layer type `type`: instance and groupInstance share 'texture', the object uniform
 * `uInstanceGranularity` selects the index; the smoothed types keep theirs (the grid types read their own textures)
 */
export function getLayerTypeDefine(type: LayerType | string): string {
  return type === 'instance' || type === 'groupInstance' ? 'texture' : type;
}

let nullGridTexture: Texture | undefined;

/**
 * The values of a layer's grid where it has none (not smoothed into a grid), shared by every object: a null texture has
 * no state, and the vectors are frozen, their writers replace them (`ValueCell.update`, as the smoothing does). A copy
 * an object (about 2 KB) adds up with 100k objects
 */
export const NoLayerGrid = {
  get texture(): Texture {
    return (nullGridTexture ??= createNullTexture());
  },
  dim: Object.freeze(Vec3.create(1, 1, 1)) as Vec3,
  texDim: Object.freeze(Vec2.create(1, 1)) as Vec2,
  transform: Object.freeze(Vec4.create(0, 0, 0, 1)) as Vec4,
};
