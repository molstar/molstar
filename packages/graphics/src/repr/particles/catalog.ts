/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { namedCatalog } from '@molstar/graphics/util/named-catalog';
import { OrientationParticlesRepresentationProvider } from './representation/orientation.js';
import { SpacefillParticlesRepresentationProvider } from './representation/spacefill.js';
import { FibersRepresentationProvider } from './representation/fibers.js';
import { ParticleTargetRepresentationProvider } from './representation/target/representation.js';

export const BuiltInParticleRepresentations = namedCatalog({
  spacefill: SpacefillParticlesRepresentationProvider,
  orientation: OrientationParticlesRepresentationProvider,
  fibers: FibersRepresentationProvider,
  target: ParticleTargetRepresentationProvider,
});
