/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { now } from '@molstar/core/util/now';
import type { PluginAnimationManager } from '@molstar/plugin/state/manager/animation';
import { isTimingMode } from '@molstar/core/util/debug';
import { printTimerResults } from '@molstar/graphics/gl/timer';
import { PluginConfig } from '@molstar/plugin/config';
import type { Subscription } from 'rxjs';

const MaxProperFrameDelta = 1000 / 30;
/**
 * With `PluginConfig.General.PauseLoopWhenIdle`, the canvas idle (`Canvas3D.isIdle`) this long before the loop stops
 * asking for frames: slack for work that is signalled a frame late.
 */
const IdleGraceMs = 250;

export class PluginAnimationLoop {
  private lastTickT: number = 0;
  // Proper time is used to prevent animations from skipping
  // if there is a blocking operation, e.g., shader compilation
  // The drawback of this is that sometimes the animation will take
  // longer than intended, but hopefully that's a reasonable tradeoff
  private properTimeT: number = 0;

  private currentFrame: number | undefined = undefined;
  private _isAnimating = false;

  private pausesWhenIdle = false;
  /** Started, but asking for no frames while the canvas is idle: the next `wake` asks again. */
  private sleeping = false;
  private idleSince = -1;
  private wakeSubs: Subscription[] | undefined = undefined;

  /** Started and not stopped: also while it sleeps (`isSleeping`). */
  get isAnimating() {
    return this._isAnimating;
  }

  /** Started, but asking for no frames: the canvas is idle (`PluginConfig.General.PauseLoopWhenIdle`). */
  get isSleeping() {
    return this.sleeping;
  }

  async tick(
    t: number,
    options?: {
      isSynchronous?: boolean;
      manualDraw?: boolean;
      animation?: PluginAnimationManager.AnimationInfo;
      updateControls?: boolean;
      xrFrame?: XRFrame;
    },
  ) {
    await this.plugin.managers.animation.tick(t, options?.isSynchronous, options?.animation);
    this.plugin.canvas3d?.tick(t as now.Timestamp, options);

    if (isTimingMode) {
      const timerResults = this.plugin.canvas3d?.graphics.timer.resolve();
      if (timerResults) {
        for (const result of timerResults) {
          printTimerResults([result]);
        }
      }
    }
  }

  private frame = (_timestamp?: number, xrFrame?: XRFrame) => {
    const t = now();
    const dt = t - this.lastTickT;
    this.lastTickT = t;
    this.properTimeT += Math.min(dt, MaxProperFrameDelta);
    this.currentFrame = undefined;
    // A frame that runs is awake, whoever asked for it.
    this.sleeping = false;
    this.tick(this.properTimeT, { xrFrame });
    if (!this._isAnimating) return;
    if (this.idleFor(t) >= IdleGraceMs) {
      this.sleeping = true;
      // Canvas3D keeps the last frame callback and asks for it again on `resume()` (`start`, an XR session change):
      // cleared, so a sleeping loop is not run twice a frame by that and a `wake` together.
      this.plugin.canvas3d?.cancelAnimationFrame(0);
      return;
    }
    this.currentFrame = this.plugin.canvas3d?.requestAnimationFrame(this.frame);
  };

  /** How long the canvas has been idle with nothing animating (ms), when the loop may sleep; else 0. */
  private idleFor(t: number) {
    const canvas = this.plugin.canvas3d;
    if (
      !this.pausesWhenIdle ||
      !canvas?.isIdle ||
      canvas.xr.isPresenting.value ||
      this.plugin.managers.animation.isAnimating
    ) {
      this.idleSince = -1;
      return 0;
    }
    if (this.idleSince < 0) this.idleSince = t;
    return t - this.idleSince;
  }

  /** Something gave the canvas work (or a state animation started): a sleeping loop asks for frames again. */
  private wake = () => {
    this.idleSince = -1;
    if (!this.sleeping || !this._isAnimating) return;
    this.sleeping = false;
    this.lastTickT = now();
    this.currentFrame = this.plugin.canvas3d?.requestAnimationFrame(this.frame);
  };

  private followWake() {
    const canvas = this.plugin.canvas3d;
    if (this.wakeSubs || !canvas) return;
    this.wakeSubs = [
      canvas.wake.subscribe(this.wake),
      this.plugin.behaviors.state.isAnimating.subscribe((animating) => {
        if (animating) this.wake();
      }),
    ];
  }

  resetTime(t: number) {
    this.plugin.canvas3d?.resetTime(t);
  }

  start(options?: { immediate?: boolean }) {
    this.plugin.canvas3d?.resume();
    this.pausesWhenIdle = !!this.plugin.config.get(PluginConfig.General.PauseLoopWhenIdle);
    if (this.pausesWhenIdle) this.followWake();
    this.sleeping = false;
    this.idleSince = -1;
    this._isAnimating = true;
    this.resetTime(0);
    this.properTimeT = 0;
    this.lastTickT = now();
    if (options?.immediate) this.frame();
    else this.currentFrame = this.plugin.canvas3d?.requestAnimationFrame(this.frame);
  }

  stop(options?: { noDraw?: boolean }) {
    this._isAnimating = false;
    this.sleeping = false;
    if (this.currentFrame !== undefined) {
      this.plugin.canvas3d?.cancelAnimationFrame(this.currentFrame);
      this.currentFrame = undefined;
    }
    if (options?.noDraw) {
      this.plugin.canvas3d?.pause(options?.noDraw);
    }
  }

  dispose() {
    this.stop();
    if (this.wakeSubs) for (const s of this.wakeSubs) s.unsubscribe();
    this.wakeSubs = undefined;
  }

  constructor(private plugin: PluginContext) {}
}
