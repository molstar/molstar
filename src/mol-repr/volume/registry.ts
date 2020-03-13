/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RepresentationRegistry, Representation } from '../representation';
import { VolumeData } from '../../mol-model/volume';
import { IsosurfaceRepresentationProvider } from './isosurface';
import { objectForEach } from '../../mol-util/object';

export class VolumeRepresentationRegistry extends RepresentationRegistry<VolumeData, Representation.State> {
    constructor() {
        super()
        Object.keys(BuiltInVolumeRepresentations).forEach(name => {
            objectForEach(BuiltInVolumeRepresentations, (p, k) => {
                if (p.name !== k) throw new Error(`Fix BuiltInVolumeRepresentations to have matching names. ${p.name} ${k}`);
                this.add(p as any)
            })
        })
    }
}

export const BuiltInVolumeRepresentations = {
    'isosurface': IsosurfaceRepresentationProvider,
    // 'direct-volume': DirectVolumeRepresentationProvider, // TODO disabled for now, needs more work
}
export type BuiltInVolumeRepresentationsName = keyof typeof BuiltInVolumeRepresentations
export const BuiltInVolumeRepresentationsNames = Object.keys(BuiltInVolumeRepresentations)
export const BuiltInVolumeRepresentationsOptions = BuiltInVolumeRepresentationsNames.map(n => [n, n] as [BuiltInVolumeRepresentationsName, string])