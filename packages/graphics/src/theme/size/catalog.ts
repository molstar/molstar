/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { namedCatalog } from '@molstar/graphics/util/named-catalog';
import { ParticleSizeThemeProvider } from './particle-size.js';
import { PhysicalSizeThemeProvider } from './physical.js';
import { ShapeGroupSizeThemeProvider } from './shape-group.js';
import { UncertaintySizeThemeProvider } from './uncertainty.js';
import { UniformSizeThemeProvider } from './uniform.js';
import { VolumeValueSizeThemeProvider } from './volume-value.js';

export const BuiltInSizeThemes = namedCatalog({
  'particle-size': ParticleSizeThemeProvider,
  physical: PhysicalSizeThemeProvider,
  'shape-group': ShapeGroupSizeThemeProvider,
  uncertainty: UncertaintySizeThemeProvider,
  uniform: UniformSizeThemeProvider,
  'volume-value': VolumeValueSizeThemeProvider,
});
