/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { objectForEach } from '@molstar/core/util/object';
import { type Representation, type RepresentationProvider, RepresentationRegistry } from '../representation.js';
import { OrientationParticlesRepresentationProvider } from './representation/orientation.js';
import { SpacefillParticlesRepresentationProvider } from './representation/spacefill.js';
import { FibersRepresentationProvider } from './representation/fibers.js';
import { ParticleTargetRepresentationProvider } from './representation/target/representation.js';

export class ParticleRepresentationRegistry extends RepresentationRegistry<ParticleList, Representation.State> {
  constructor() {
    super();
    objectForEach(ParticleRepresentationRegistry.BuiltIn, (p, k) => {
      if (p.name !== k) throw new Error(`Fix BuiltInParticleRepresentations to have matching names. ${p.name} ${k}`);
      this.add(p as any);
    });
  }
}

export namespace ParticleRepresentationRegistry {
  export const BuiltIn = {
    spacefill: SpacefillParticlesRepresentationProvider,
    orientation: OrientationParticlesRepresentationProvider,
    fibers: FibersRepresentationProvider,
    target: ParticleTargetRepresentationProvider,
  };

  type _BuiltIn = typeof BuiltIn;
  export type BuiltIn = keyof _BuiltIn;
  export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>;
}
