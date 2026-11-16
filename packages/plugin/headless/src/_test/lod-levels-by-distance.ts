/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { getLodLevelsValue, type LodLevels } from '@molstar/graphics/geo/geometry/lod-levels';
import type { SpheresValues } from '@molstar/graphics/gl/renderable/spheres';

/**
 * Levels of detail by distance, set on `values` as the renderable takes them (the props take a pixel size only): a level
 * every stride-th sphere, the strides adjusted for `sizeFactor`
 */
export function setLodLevelsByDistance(values: SpheresValues, levels: LodLevels, sizeFactor: number) {
  const count = values.uVertexCount.ref.value;
  ValueCell.update(values.lodLevels, getLodLevelsValue(levels, sizeFactor, count, 0, 0));
}
