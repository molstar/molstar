/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { SphericalHarmonicSurfaceRepresentationProvider } from './representation';

export const SphericalHarmonicSurface = PluginBehavior.create<{}>({
    name: 'spherical-harmonic-surface',
    category: 'representation',
    display: {
        name: 'Spherical Harmonic Surface',
        description: 'A smooth shape envelope approximated by a spherical harmonic expansion of the atom positions.',
    },
    ctor: class extends PluginBehavior.Handler<{}> {
        register(): void {
            this.ctx.representation.structure.registry.add(SphericalHarmonicSurfaceRepresentationProvider);
        }

        unregister() {
            this.ctx.representation.structure.registry.remove(SphericalHarmonicSurfaceRepresentationProvider);
        }
    },
});
