/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

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
