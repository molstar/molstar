/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../mol-model/structure';
import { RepresentationProvider, RepresentationRegistry } from '../representation';
import { CartoonRepresentationProvider } from './representation/cartoon';
import { BallAndStickRepresentationProvider } from './representation/ball-and-stick';
import { GaussianSurfaceRepresentationProvider } from './representation/gaussian-surface';
import { CarbohydrateRepresentationProvider } from './representation/carbohydrate';
import { SpacefillRepresentationProvider } from './representation/spacefill';
import { PointRepresentationProvider } from './representation/point';
import { StructureRepresentationState } from './representation';
import { PuttyRepresentationProvider } from './representation/putty';
import { MolecularSurfaceRepresentationProvider } from './representation/molecular-surface';
import { EllipsoidRepresentationProvider } from './representation/ellipsoid';
import { OrientationRepresentationProvider } from './representation/orientation';
import { LabelRepresentationProvider } from './representation/label';

export class StructureRepresentationRegistry extends RepresentationRegistry<Structure, StructureRepresentationState> {
    constructor() {
        super()
        Object.keys(BuiltInStructureRepresentations).forEach(name => {
            const p = (BuiltInStructureRepresentations as { [k: string]: RepresentationProvider<Structure, any, StructureRepresentationState> })[name]
            this.add(name, p)
        })
    }
}

export const BuiltInStructureRepresentations = {
    'cartoon': CartoonRepresentationProvider,
    'ball-and-stick': BallAndStickRepresentationProvider,
    'carbohydrate': CarbohydrateRepresentationProvider,
    'ellipsoid': EllipsoidRepresentationProvider,
    'gaussian-surface': GaussianSurfaceRepresentationProvider,
    // 'gaussian-volume': GaussianVolumeRepresentationProvider, // TODO disabled for now, needs more work
    'label': LabelRepresentationProvider,
    'molecular-surface': MolecularSurfaceRepresentationProvider,
    'orientation': OrientationRepresentationProvider,
    'point': PointRepresentationProvider,
    'putty': PuttyRepresentationProvider,
    'spacefill': SpacefillRepresentationProvider,
}
export type BuiltInStructureRepresentations = typeof BuiltInStructureRepresentations
export type BuiltInStructureRepresentationsName = keyof typeof BuiltInStructureRepresentations
export const BuiltInStructureRepresentationsNames = Object.keys(BuiltInStructureRepresentations)
export const BuiltInStructureRepresentationsOptions = BuiltInStructureRepresentationsNames.map(n => [n, n] as [BuiltInStructureRepresentationsName, string])