/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginContext } from 'mol-plugin/context';

export { PluginStateAnimation }

// TODO: helpers for building animations (once more animations are added)
//       for example "composite animation"

interface PluginStateAnimation<P extends PD.Params = any, S = any> {
    name: string,
    readonly display: { readonly name: string, readonly description?: string },
    params: (ctx: PluginContext) => P,
    initialState(params: PD.Values<P>, ctx: PluginContext): S,

    /**
     * Apply the current frame and modify the state.
     * @param t Current absolute time since the animation started.
     */
    apply(state: S, t: PluginStateAnimation.Time, ctx: PluginStateAnimation.Context<P>): Promise<PluginStateAnimation.ApplyResult<S>>,

    /**
     * The state must be serializable to JSON. If JSON.stringify is not enough,
     * custom converted to an object that works with JSON.stringify can be provided.
     */
    stateSerialization?: { toJSON(state: S): any, fromJSON(data: any): S }
}

namespace PluginStateAnimation {
    export interface Time {
        lastApplied: number,
        current: number
    }

    export type ApplyResult<S> = { kind: 'finished' } | { kind: 'skip' } | { kind: 'next', state: S }
    export interface Context<P extends PD.Params> {
        params: PD.Values<P>,
        plugin: PluginContext
    }

    export function create<P extends PD.Params, S>(params: PluginStateAnimation<P, S>) {
        return params;
    }
}