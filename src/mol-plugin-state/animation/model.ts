/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginContext } from '../../mol-plugin/context';

export { PluginStateAnimation };

// TODO: helpers for building animations (once more animations are added)
//       for example "composite animation"

interface PluginStateAnimation<P = any, S = any> {
    name: string,
    readonly display: { readonly name: string, readonly description?: string },
    readonly isExportable?: boolean,

    params(ctx: PluginContext): PD.For<P>,
    canApply?(ctx: PluginContext): PluginStateAnimation.CanApply,
    initialState(params: P, ctx: PluginContext): S,
    getDuration?(params: P, ctx: PluginContext): PluginStateAnimation.Duration,

    // Optionally returns new state
    setup?(params: P, state: S, ctx: PluginContext): void | Promise<S> | Promise<void>,
    teardown?(params: P, state: S, ctx: PluginContext): void | Promise<void>,

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
    export type CanApply = { canApply: true } | { canApply: false, reason?: string }
    export type Duration = { kind: 'unknown' } | { kind: 'infinite' } | { kind: 'fixed', durationMs: number }

    export interface Instance<A extends PluginStateAnimation> {
        definition: PluginStateAnimation,
        params: Params<A>,
        customDurationMs?: number
    }

    export interface Time {
        lastApplied: number,
        current: number,
        animation?: { currentFrame: number, frameCount: number }
    }

    export type ApplyResult<S> = { kind: 'finished' } | { kind: 'skip' } | { kind: 'next', state: S }
    export interface Context<P> {
        params: P,
        plugin: PluginContext
    }

    export type Params<A extends PluginStateAnimation> = A extends PluginStateAnimation<infer P> ? P : never

    export function create<P, S>(params: PluginStateAnimation<P, S>) {
        return params;
    }

    export function getDuration<A extends PluginStateAnimation>(ctx: PluginContext, instance: Instance<A>) {
        if (instance.customDurationMs) return instance.customDurationMs;
        const d = instance.definition.getDuration?.(instance.params, ctx);
        if (d?.kind === 'fixed') return d.durationMs;
    }
}