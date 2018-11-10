/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { RepresentationProvider, RepresentationRegistry } from '../representation';
import { CartoonRepresentationProvider } from './representation/cartoon';
import { BallAndStickRepresentationProvider } from './representation/ball-and-stick';

export class StructureRepresentationRegistry extends RepresentationRegistry<Structure> {
    constructor() {
        super()
        Object.keys(BuiltInStructureRepresentations).forEach(name => {
            const p = (BuiltInStructureRepresentations as { [k: string]: RepresentationProvider<Structure, any> })[name]
            this.add(name, p.factory, p.params)
        })
    }
}

export const BuiltInStructureRepresentations = {
    'cartoon': CartoonRepresentationProvider,
    'ball-and-stick': BallAndStickRepresentationProvider,
}
export type BuiltInStructureRepresentationsName = keyof typeof BuiltInStructureRepresentations
export const BuiltInStructureRepresentationsNames = Object.keys(BuiltInStructureRepresentations)
export const BuiltInStructureRepresentationsOptions = BuiltInStructureRepresentationsNames.map(n => [n, n] as [BuiltInStructureRepresentationsName, string])