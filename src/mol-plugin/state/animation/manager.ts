/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginComponent } from 'mol-plugin/component';
import { PluginContext } from 'mol-plugin/context';
import { PluginStateAnimation } from './model';
import { ParamDefinition as PD } from 'mol-util/param-definition';

export { PluginAnimationManager }

class PluginAnimationManager extends PluginComponent<PluginAnimationManager.State> {
    private map = new Map<string, PluginStateAnimation>();
    private animations: PluginStateAnimation[] = [];
    private _current: PluginAnimationManager.Current;
    private _params?: PD.For<PluginAnimationManager.State['params']> = void 0;

    get isEmpty() { return this.animations.length === 0; }
    get current() { return this._current!; }

    getParams(): PD.Params {
        if (!this._params) {
            this._params = {
                current: PD.Select(this.animations[0] && this.animations[0].name,
                    this.animations.map(a => [a.name, a.display.name] as [string, string]),
                    { label: 'Animation' })
            };
        }
        return this._params as any as PD.Params;
    }

    updateParams(newParams: Partial<PluginAnimationManager.State['params']>) {
        this.updateState({ params: { ...this.latestState.params, ...newParams } });
        const anim = this.map.get(this.latestState.params.current)!;
        const params = anim.params(this.context);
        this._current = {
            anim,
            params,
            paramValues: PD.getDefaultValues(params),
            state: {},
            startedTime: -1,
            lastTime: 0
        }
        this.triggerUpdate();
    }

    updateCurrentParams(values: any) {
        this._current.paramValues = { ...this._current.paramValues, ...values };
        this.triggerUpdate();
    }

    register(animation: PluginStateAnimation) {
        if (this.map.has(animation.name)) {
            this.context.log.error(`Animation '${animation.name}' is already registered.`);
            return;
        }
        this._params = void 0;
        this.map.set(animation.name, animation);
        this.animations.push(animation);
        if (this.animations.length === 1) {
            this.updateParams({ current: animation.name });
        } else {
            this.triggerUpdate();
        }
    }

    start() {
        this.updateState({ animationState: 'playing' });
        this.triggerUpdate();

        this._current.lastTime = 0;
        this._current.startedTime = -1;
        this._current.state = this._current.anim.initialState(this._current.paramValues, this.context);

        requestAnimationFrame(this.animate);
    }

    stop() {
        this.updateState({ animationState: 'stopped' });
        this.triggerUpdate();
    }

    private animate = async (t: number) => {
        if (this._current.startedTime < 0) this._current.startedTime = t;
        const newState = await this._current.anim.apply(
            this._current.state,
            { lastApplied: this._current.lastTime, current: t - this._current.startedTime },
            { params: this._current.paramValues, plugin: this.context });

        if (newState.kind === 'finished') {
            this.stop();
        } else if (newState.kind === 'next') {
            this._current.state = newState.state;
            this._current.lastTime = t - this._current.startedTime;
            if (this.latestState.animationState === 'playing') requestAnimationFrame(this.animate);
        } else if (newState.kind === 'skip') {
            if (this.latestState.animationState === 'playing') requestAnimationFrame(this.animate);
        }
    }

    getSnapshot(): PluginAnimationManager.Snapshot {
        if (!this.current) return { state: this.latestState };

        return {
            state: this.latestState,
            current: {
                paramValues: this._current.paramValues,
                state: this._current.anim.stateSerialization ? this._current.anim.stateSerialization.toJSON(this._current.state) : this._current.state
            }
        };
    }

    setSnapshot(snapshot: PluginAnimationManager.Snapshot) {
        this.updateState({ animationState: snapshot.state.animationState });
        this.updateParams(snapshot.state.params);

        if (snapshot.current) {
            this.current.paramValues = snapshot.current.paramValues;
            this.current.state = this._current.anim.stateSerialization
                ? this._current.anim.stateSerialization.fromJSON(snapshot.current.state)
                : snapshot.current.state;
            this.triggerUpdate();
            if (this.latestState.animationState === 'playing') this.resume();
        }
    }

    private resume() {
        this._current.lastTime = 0;
        this._current.startedTime = -1;
        requestAnimationFrame(this.animate);
    }

    constructor(ctx: PluginContext) {
        super(ctx, { params: { current: '' }, animationState: 'stopped' });
    }
}

namespace PluginAnimationManager {
    export interface Current {
        anim: PluginStateAnimation
        params: PD.Params,
        paramValues: any,
        state: any,
        startedTime: number,
        lastTime: number
    }

    export interface State {
        params: { current: string },
        animationState: 'stopped' | 'playing'
    }

    export interface Snapshot {
        state: State,
        current?: {
            paramValues: any,
            state: any
        }
    }
}