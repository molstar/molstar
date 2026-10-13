/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RepresentationRegistry, type Representation, type RepresentationProvider } from '../representation.js';
import type { BuiltInVolumeRepresentations } from './catalog.js';
import type { Volume } from '@molstar/model/model/volume';

export class VolumeRepresentationRegistry extends RepresentationRegistry<Volume, Representation.State> {}

export namespace VolumeRepresentationRegistry {
  type _BuiltIn = typeof BuiltInVolumeRepresentations;
  export type BuiltIn = keyof _BuiltIn;
  export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>;
}
