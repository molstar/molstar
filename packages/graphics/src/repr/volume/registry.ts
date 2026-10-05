/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RepresentationRegistry, type Representation, type RepresentationProvider } from '../representation.js';
import type { Volume } from '@molstar/model/model/volume';
import { IsosurfaceRepresentationProvider } from './isosurface.js';
import { objectForEach } from '@molstar/core/util/object';
import { SliceRepresentationProvider } from './slice.js';
import { DirectVolumeRepresentationProvider } from './direct-volume.js';
import { SegmentRepresentationProvider } from './segment.js';
import { DotRepresentationProvider } from './dot.js';

export class VolumeRepresentationRegistry extends RepresentationRegistry<Volume, Representation.State> {
  constructor() {
    super();
    objectForEach(VolumeRepresentationRegistry.BuiltIn, (p, k) => {
      if (p.name !== k) throw new Error(`Fix BuiltInVolumeRepresentations to have matching names. ${p.name} ${k}`);
      this.add(p as any);
    });
  }
}

export namespace VolumeRepresentationRegistry {
  export const BuiltIn = {
    'direct-volume': DirectVolumeRepresentationProvider,
    dot: DotRepresentationProvider,
    isosurface: IsosurfaceRepresentationProvider,
    segment: SegmentRepresentationProvider,
    slice: SliceRepresentationProvider,
  };

  type _BuiltIn = typeof BuiltIn;
  export type BuiltIn = keyof _BuiltIn;
  export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>;
}
