/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { objectForEach } from '@molstar/core/util/object';
import { RepresentationRegistry, type Representation, type RepresentationProvider } from '../representation.js';
import { BuiltInVolumeRepresentations } from './catalog.js';
import type { Volume } from '@molstar/model/model/volume';

export class VolumeRepresentationRegistry extends RepresentationRegistry<Volume, Representation.State> {
  constructor() {
    super();
    objectForEach(BuiltInVolumeRepresentations, (p) => this.add(p as any));
  }
}

export namespace VolumeRepresentationRegistry {
  type _BuiltIn = typeof BuiltInVolumeRepresentations;
  export type BuiltIn = keyof _BuiltIn;
  export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>;
}
