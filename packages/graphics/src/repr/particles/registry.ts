/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { type Representation, type RepresentationProvider, RepresentationRegistry } from '../representation.js';
import type { BuiltInParticleRepresentations } from './catalog.js';

export class ParticleRepresentationRegistry extends RepresentationRegistry<ParticleList, Representation.State> {}

export namespace ParticleRepresentationRegistry {
  type _BuiltIn = typeof BuiltInParticleRepresentations;
  export type BuiltIn = keyof _BuiltIn;
  export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>;
}
