/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateAnimation } from '../model';

export const AnimateTime = PluginStateAnimation.create({
    name: 'built-in.animate-time',
    display: { name: 'Animate Time', description: 'Animate the passage of time in the 3D scene' },
    isExportable: true,
    params: () => ({
        durationInMs: PD.Numeric(4000, { min: 100, max: 20000, step: 100 }),
    }),
    initialState: () => ({ }),
    getDuration: p => ({ kind: 'fixed', durationMs: p.durationInMs }),

    async apply(animState, t, ctx) {
        return t.current < ctx.params.durationInMs
            ? { kind: 'next', state: animState }
            : { kind: 'finished' };
    }
});
