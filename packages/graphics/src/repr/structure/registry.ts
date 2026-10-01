/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Structure } from '@molstar/model/model/structure';
import { objectForEach } from '@molstar/core/util/object';
import { RepresentationRegistry, type RepresentationProvider } from '../representation.js';
import type { StructureRepresentationState } from './representation.js';
import { BallAndStickRepresentationProvider } from './representation/ball-and-stick.js';
import { BlobSurfaceRepresentationProvider } from './representation/blob-surface.js';
import { CarbohydrateRepresentationProvider } from './representation/carbohydrate.js';
import { CartoonRepresentationProvider } from './representation/cartoon.js';
import { EllipsoidRepresentationProvider } from './representation/ellipsoid.js';
import { GaussianSurfaceRepresentationProvider } from './representation/gaussian-surface.js';
import { LabelRepresentationProvider } from './representation/label.js';
import { MolecularSurfaceRepresentationProvider } from './representation/molecular-surface.js';
import { OrientationRepresentationProvider } from './representation/orientation.js';
import { PointRepresentationProvider } from './representation/point.js';
import { PuttyRepresentationProvider } from './representation/putty.js';
import { SpacefillRepresentationProvider } from './representation/spacefill.js';
import { LineRepresentationProvider } from './representation/line.js';
import { GaussianVolumeRepresentationProvider } from './representation/gaussian-volume.js';
import { BackboneRepresentationProvider } from './representation/backbone.js';
import { PolyhedronRepresentationProvider } from './representation/polyhedron.js';
import { PlaneRepresentationProvider } from './representation/plane.js';

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
        'backbone': BackboneRepresentationProvider,
        'ball-and-stick': BallAndStickRepresentationProvider,
        'blob-surface': BlobSurfaceRepresentationProvider,
        'carbohydrate': CarbohydrateRepresentationProvider,
        'ellipsoid': EllipsoidRepresentationProvider,
        'gaussian-surface': GaussianSurfaceRepresentationProvider,
        'gaussian-volume': GaussianVolumeRepresentationProvider,
        'label': LabelRepresentationProvider,
        'line': LineRepresentationProvider,
        'molecular-surface': MolecularSurfaceRepresentationProvider,
        'orientation': OrientationRepresentationProvider,
        'plane': PlaneRepresentationProvider,
        'point': PointRepresentationProvider,
        'putty': PuttyRepresentationProvider,
        'spacefill': SpacefillRepresentationProvider,
        'polyhedron': PolyhedronRepresentationProvider,
    };

    type _BuiltIn = typeof BuiltIn
    export type BuiltIn = keyof _BuiltIn
    export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>
}