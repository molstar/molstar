/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StatefulPluginComponent } from '../component.js';
import type { PluginContext } from '@molstar/plugin/context';
import type { PluginStateAnimation } from '../animation/model.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { isDebugMode } from '@molstar/core/util/debug';

export { PluginAnimationManager };

// TODO: pause functionality (this needs to reset if the state tree changes)
// TODO: handle unregistered animations on state restore
// TODO: better API

class PluginAnimationManager extends StatefulPluginComponent<PluginAnimationManager.State> {
  private map = new Map<string, PluginStateAnimation>();
  private _animations: PluginStateAnimation[] = [];
  /** Registration counts by name; 0 for an animation `play` added on demand */
  private counts = new Map<string, number>();
  private currentTime: number = 0;

  private _current: PluginAnimationManager.Current | undefined = void 0;
  private _params?: PD.For<PluginAnimationManager.State['params']> = void 0;

  readonly events = {
    updated: this.ev(),
    applied: this.ev(),
  };

  get isEmpty() {
    return this._animations.length === 0;
  }
  /** The selected animation, or `undefined` when none is registered. */
  get current(): PluginAnimationManager.Current | undefined {
    return this._current;
  }

  get animations() {
    return this._animations;
  }

  get isAnimatingStateTransition() {
    return !!this._current && this._current.anim.name === 'built-in.animate-state-snapshot-transition';
  }

  private triggerUpdate() {
    this.events.updated.next(void 0);
  }

  private triggerApply() {
    this.events.applied.next(void 0);
  }

  getParams(): PD.Params {
    if (!this._params) {
      this._params = {
        current: PD.Select(
          this._animations[0] && this._animations[0].name,
          this._animations.map((a) => [a.name, a.display.name] as [string, string]),
          { label: 'Animation' },
        ),
      };
    }
    return this._params as any as PD.Params;
  }

  updateParams(newParams: Partial<PluginAnimationManager.State['params']>) {
    if (this.isEmpty) return;
    this.updateState({ params: { ...this.state.params, ...newParams } });
    const anim = this.map.get(this.state.params.current);
    if (!anim) {
      if (isDebugMode) {
        console.warn(
          `Animation '${this.state.params.current}' not found. State might be from a different plugin instance.`,
        );
      }
      return;
    }

    const params = anim.params(this.context) as PD.Params;
    this._current = {
      anim,
      params,
      paramValues: PD.getDefaultValues(params),
      state: {},
      startedTime: -1,
      lastTime: 0,
    };
    this.triggerUpdate();
  }

  updateCurrentParams(values: any) {
    if (!this._current) return;
    this._current.paramValues = { ...this._current.paramValues, ...values };
    this.triggerUpdate();
  }

  /** Returns the message `register` would throw for `animation`, without changing anything. */
  findConflict(animation: PluginStateAnimation): string | undefined {
    const existing = this.map.get(animation.name);
    if (existing && existing !== animation) {
      return `PluginAnimationManager: a different animation is already registered under the name '${animation.name}'.`;
    }
    return undefined;
  }

  /**
   * Animations are keyed by `name`. Registering the same object again increments a count,
   * a different object under an existing name throws.
   */
  register(animation: PluginStateAnimation) {
    const conflict = this.findConflict(animation);
    if (conflict) throw new Error(conflict);

    if (this.map.has(animation.name)) {
      // also adopts an animation added on demand by `play`
      this.counts.set(animation.name, (this.counts.get(animation.name) ?? 0) + 1);
      return;
    }
    this.counts.set(animation.name, 1);
    this.add(animation);
  }

  /**
   * Decrements the count and removes the animation at zero. No-op for an unknown animation
   * and for one that `play` added on demand. Removing the current animation selects the first remaining one,
   * or none, and resets the cached params.
   */
  unregister(animation: PluginStateAnimation) {
    if (this.map.get(animation.name) !== animation) return;
    const count = this.counts.get(animation.name) ?? 0;
    if (count <= 0) return;
    if (count > 1) {
      this.counts.set(animation.name, count - 1);
      return;
    }

    this.counts.delete(animation.name);
    this.remove(animation);
  }

  private add(animation: PluginStateAnimation) {
    this._params = void 0;
    this.map.set(animation.name, animation);
    this._animations.push(animation);
    if (this._animations.length === 1) {
      this.updateParams({ current: animation.name });
    } else {
      this.triggerUpdate();
    }
  }

  private remove(animation: PluginStateAnimation) {
    this.map.delete(animation.name);
    const idx = this._animations.indexOf(animation);
    if (idx >= 0) this._animations.splice(idx, 1);
    this._params = void 0;

    const current = this._current;
    if (!current || current.anim !== animation) {
      this.triggerUpdate();
      return;
    }

    // the removed animation is the current one
    this.isStopped = true;
    if (this.state.animationState !== 'stopped') {
      this.updateState({ animationState: 'stopped' });
      if (this.context.behaviors.state.isAnimating.value) {
        this.context.behaviors.state.isAnimating.next(false);
      }
      if (animation.teardown) {
        Promise.resolve(animation.teardown(current.paramValues, current.state, this.context)).catch((e) =>
          console.error(`Failed to tear down animation '${animation.name}'`, e),
        );
      }
    }

    const next = this._animations[0];
    if (next) {
      this.updateParams({ current: next.name });
    } else {
      this._current = void 0;
      this.updateState({ params: { ...this.state.params, current: '' } });
      this.triggerUpdate();
    }
  }

  async play<P>(animation: PluginStateAnimation<P>, params: P) {
    await this.stop();
    if (!this.map.has(animation.name)) {
      // registered on demand, outside reference counting
      this.counts.set(animation.name, 0);
      this.add(animation);
    }
    this.updateParams({ current: animation.name });
    this.updateCurrentParams(params);
    await this.start();
  }

  async tick(t: number, isSynchronous?: boolean, animation?: PluginAnimationManager.AnimationInfo) {
    this.currentTime = t;
    if (this.isStopped) return;

    if (isSynchronous || animation) {
      await this.applyFrame(animation);
    } else {
      this.applyAsync();
    }
  }

  private isStopped = true;
  private isApplying = false;

  async start() {
    if (!this._current) return;

    this.updateState({ animationState: 'playing' });
    if (!this.context.behaviors.state.isAnimating.value) {
      this.context.behaviors.state.isAnimating.next(true);
    }
    this.triggerUpdate();

    const current = this._current;
    if (!current) return;

    const anim = current.anim;
    let initialState = anim.initialState(current.paramValues, this.context);
    if (anim.setup) {
      const state = await anim.setup(current.paramValues, initialState, this.context);
      if (state) initialState = state;
    }

    current.lastTime = 0;
    current.startedTime = -1;
    current.state = initialState;
    this.isStopped = false;
  }

  async stop() {
    this.isStopped = true;
    if (this.state.animationState !== 'stopped') {
      const current = this._current;
      if (current?.anim.teardown) {
        await current.anim.teardown(current.paramValues, current.state, this.context);
      }

      this.updateState({ animationState: 'stopped' });
      this.triggerUpdate();
    }

    if (this.context.behaviors.state.isAnimating.value) {
      this.context.behaviors.state.isAnimating.next(false);
    }
  }

  stopStateTransitionAnimation() {
    if (!this.isAnimatingStateTransition) return;
    return this.stop();
  }

  get isAnimating() {
    return this.state.animationState === 'playing';
  }

  private async applyAsync() {
    if (this.isApplying) return;

    this.isApplying = true;
    try {
      await this.applyFrame();
    } finally {
      this.isApplying = false;
    }
  }

  private async applyFrame(animation?: PluginAnimationManager.AnimationInfo) {
    const current = this._current;
    if (!current) return;

    const t = this.currentTime;
    if (current.startedTime < 0) current.startedTime = t;
    const newState = await current.anim.apply(
      current.state,
      { lastApplied: current.lastTime, current: t - current.startedTime, animation },
      { params: current.paramValues, plugin: this.context },
    );

    if (newState.kind === 'finished') {
      this.stop();
    } else if (newState.kind === 'next') {
      current.state = newState.state;
      current.lastTime = t - current.startedTime;
    }
    this.triggerApply();
  }

  getSnapshot(): PluginAnimationManager.Snapshot {
    const current = this._current;
    if (!current) return { state: this.state };

    return {
      state: this.state,
      current: {
        paramValues: current.paramValues,
        state: current.anim.stateSerialization ? current.anim.stateSerialization.toJSON(current.state) : current.state,
      },
    };
  }

  setSnapshot(snapshot: PluginAnimationManager.Snapshot) {
    if (this.isEmpty) return;
    this.updateState({ animationState: snapshot.state.animationState });
    this.updateParams(snapshot.state.params);

    if (snapshot.current) {
      const name = snapshot.state.params.current;
      const current = this._current;
      if (!current || current.anim.name !== name) {
        // `updateParams` ignores an unregistered name, so `current` would still be the previous animation
        this.context.log.warn(`Animation '${name}' of the snapshot is not registered; skipping its state.`);
        return;
      }

      current.paramValues = snapshot.current.paramValues;
      current.state = current.anim.stateSerialization
        ? current.anim.stateSerialization.fromJSON(snapshot.current.state)
        : snapshot.current.state;
      this.triggerUpdate();
      if (this.state.animationState === 'playing') this.resume();
    }
  }

  private async resume() {
    const current = this._current;
    if (!current) return;

    current.lastTime = 0;
    current.startedTime = -1;
    const anim = current.anim;
    if (!this.context.behaviors.state.isAnimating.value) {
      this.context.behaviors.state.isAnimating.next(true);
    }
    if (anim.setup) {
      await anim.setup(current.paramValues, current.state, this.context);
    }
    this.isStopped = false;
  }

  constructor(private context: PluginContext) {
    super({ params: { current: '' }, animationState: 'stopped' });
  }
}

namespace PluginAnimationManager {
  export interface AnimationInfo {
    currentFrame: number;
    frameCount: number;
  }

  export interface Current {
    anim: PluginStateAnimation;
    params: PD.Params;
    paramValues: any;
    state: any;
    startedTime: number;
    lastTime: number;
  }

  export interface State {
    params: { current: string };
    animationState: 'stopped' | 'playing';
  }

  export interface Snapshot {
    state: State;
    current?: {
      paramValues: any;
      state: any;
    };
  }
}
