/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RepresentationRegistry, Representation, RepresentationProvider } from '../representation';
import { Volume } from '../../mol-model/volume';
import { IsosurfaceRepresentationProvider } from './isosurface';
import { objectForEach } from '../../mol-util/object';
import { SliceRepresentationProvider } from './slice';

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
        'isosurface': IsosurfaceRepresentationProvider,
        'slice': SliceRepresentationProvider,
        // 'direct-volume': DirectVolumeRepresentationProvider, // TODO disabled for now, needs more work
    };

    type _BuiltIn = typeof BuiltIn
    export type BuiltIn = keyof _BuiltIn
    export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>
}