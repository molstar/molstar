/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { namedCatalog } from '@molstar/graphics/util/named-catalog';
import { DirectVolumeRepresentationProvider } from './direct-volume.js';
import { DotRepresentationProvider } from './dot.js';
import { IsosurfaceRepresentationProvider } from './isosurface.js';
import { SegmentRepresentationProvider } from './segment.js';
import { SliceRepresentationProvider } from './slice.js';

export const BuiltInVolumeRepresentations = namedCatalog({
  'direct-volume': DirectVolumeRepresentationProvider,
  dot: DotRepresentationProvider,
  isosurface: IsosurfaceRepresentationProvider,
  segment: SegmentRepresentationProvider,
  slice: SliceRepresentationProvider,
});
