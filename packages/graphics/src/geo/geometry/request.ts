/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Theme } from '@molstar/graphics/theme/theme';
import { isPositionColorType } from './color-data.js';
import { isPositionSizeType } from './size-data.js';

/** What creating or updating values needs of a geometry, see `GeometryUtils.prepare` */
export interface GeometryRequest {
  /** the position iterator is used (vertex-granularity color or size themes) */
  readonly positions: boolean;
}

export namespace GeometryRequest {
  export const None: GeometryRequest = { positions: false };

  /** The request of `theme`: positions for a vertex-granularity color or size theme */
  export function fromTheme(theme: Theme): GeometryRequest {
    return {
      positions: isPositionColorType(theme.color.granularity) || isPositionSizeType(theme.size.granularity),
    };
  }
}
