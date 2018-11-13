/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RepresentationProvider, RepresentationRegistry } from '../representation';
import { VolumeData } from 'mol-model/volume';

export class VolumeRepresentationRegistry extends RepresentationRegistry<VolumeData> {
    constructor() {
        super()
        Object.keys(BuiltInVolumeRepresentations).forEach(name => {
            const p = (BuiltInVolumeRepresentations as { [k: string]: RepresentationProvider<VolumeData, any> })[name]
            this.add(name, p.factory, p.getParams)
        })
    }
}

export const BuiltInVolumeRepresentations = {
    // TODO
}
export type BuiltInVolumeRepresentationsName = keyof typeof BuiltInVolumeRepresentations
export const BuiltInVolumeRepresentationsNames = Object.keys(BuiltInVolumeRepresentations)
export const BuiltInVolumeRepresentationsOptions = BuiltInVolumeRepresentationsNames.map(n => [n, n] as [BuiltInVolumeRepresentationsName, string])