/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../mol-model/structure';
import { objectForEach } from '../../mol-util/object';
import { RepresentationRegistry, RepresentationProvider } from '../representation';
import { StructureRepresentationState } from './representation';
import { BallAndStickRepresentationProvider } from './representation/ball-and-stick';
import { CarbohydrateRepresentationProvider } from './representation/carbohydrate';
import { CartoonRepresentationProvider } from './representation/cartoon';
import { EllipsoidRepresentationProvider } from './representation/ellipsoid';
import { GaussianSurfaceRepresentationProvider } from './representation/gaussian-surface';
import { LabelRepresentationProvider } from './representation/label';
import { MolecularSurfaceRepresentationProvider } from './representation/molecular-surface';
import { OrientationRepresentationProvider } from './representation/orientation';
import { PointRepresentationProvider } from './representation/point';
import { PuttyRepresentationProvider } from './representation/putty';
import { SpacefillRepresentationProvider } from './representation/spacefill';

export class StructureRepresentationRegistry extends RepresentationRegistry<Structure, StructureRepresentationState> {
    constructor() {
        super();
        objectForEach(StructureRepresentationRegistry.BuiltIn, (p, k) => {
            if (p.name !== k) throw new Error(`Fix BuiltInStructureRepresentations to have matching names. ${p.name} ${k}`);
            this.add(p as any);
        });
    }
}

export namespace StructureRepresentationRegistry {
    export const BuiltIn = {
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
    };

    type _BuiltIn = typeof BuiltIn
    export type BuiltIn = keyof _BuiltIn
    export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>
}