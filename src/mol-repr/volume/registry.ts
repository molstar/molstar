/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RepresentationProvider, RepresentationRegistry, Representation } from '../representation';
import { VolumeData } from 'mol-model/volume';
import { IsosurfaceRepresentationProvider } from './isosurface';
import { DirectVolumeRepresentationProvider } from './direct-volume';

export class VolumeRepresentationRegistry extends RepresentationRegistry<VolumeData, Representation.State> {
    constructor() {
        super()
        Object.keys(BuiltInVolumeRepresentations).forEach(name => {
            const p = (BuiltInVolumeRepresentations as { [k: string]: RepresentationProvider<VolumeData, any, Representation.State> })[name]
            this.add(name, p)
        })
    }
}

export const BuiltInVolumeRepresentations = {
    'isosurface': IsosurfaceRepresentationProvider,
    'direct-volume': DirectVolumeRepresentationProvider,
}
export type BuiltInVolumeRepresentationsName = keyof typeof BuiltInVolumeRepresentations
export const BuiltInVolumeRepresentationsNames = Object.keys(BuiltInVolumeRepresentations)
export const BuiltInVolumeRepresentationsOptions = BuiltInVolumeRepresentationsNames.map(n => [n, n] as [BuiltInVolumeRepresentationsName, string])