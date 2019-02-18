/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginContext } from 'mol-plugin/context';

export { PluginStateAnimation }

interface PluginStateAnimation<P extends PD.Params, S> {
    id: string,
    params: (ctx: PluginContext) => P,
    initialState(params: PD.Values<P>, ctx: PluginContext): S,

    /**
     * Apply the current frame and modify the state.
     * @param t Current absolute time since the animation started.
     */
    apply(state: S, t: number, ctx: PluginStateAnimation.Context<P>): Promise<PluginStateAnimation.ApplyResult<S>>,

    /**
     * The state must be serializable to JSON. If JSON.stringify is not enough,
     * custom serializer can be provided.
     */
    stateSerialization?: { toJson?(state: S): any, fromJson?(data: any): S }
}

namespace PluginStateAnimation {
    export type ApplyResult<S> = { kind: 'finished' } | { kind: 'next', state: S }
    export interface Context<P extends PD.Params> {
        params: PD.Values<P>,
        plugin: PluginContext
    }
}

