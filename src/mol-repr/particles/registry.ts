/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParticleList } from '../../mol-model/particles/particle-list';
import { objectForEach } from '../../mol-util/object';
import { Representation, RepresentationProvider, RepresentationRegistry } from '../representation';
import { OrientationParticlesRepresentationProvider } from './representation/orientation';
import { SpacefillParticlesRepresentationProvider } from './representation/spacefill';
import { ParticleTargetRepresentationProvider } from './representation/target/representation';
import { FibersRepresentationProvider } from './representation/fibers';

export class ParticleRepresentationRegistry extends RepresentationRegistry<ParticleList, Representation.State> {
    constructor() {
        super();
        objectForEach(ParticleRepresentationRegistry.BuiltIn, (p, k) => {
            if (p.name !== k) throw new Error(`Fix BuiltInParticleRepresentations to have matching names. ${p.name} ${k}`);
            this.add(p as any);
        });
    }
}

export namespace ParticleRepresentationRegistry {
    export const BuiltIn = {
        'spacefill': SpacefillParticlesRepresentationProvider,
        'orientation': OrientationParticlesRepresentationProvider,
        'target': ParticleTargetRepresentationProvider,
        'fibers': FibersRepresentationProvider,
    };

    type _BuiltIn = typeof BuiltIn
    export type BuiltIn = keyof _BuiltIn
    export type BuiltInParams<T extends BuiltIn> = Partial<RepresentationProvider.ParamValues<_BuiltIn[T]>>
}
