/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from './context';
import { now } from '../mol-util/now';
import { PluginAnimationManager } from '../mol-plugin-state/manager/animation';
import { isTimingMode } from '../mol-util/debug';
import { printTimerResults } from '../mol-gl/webgl/timer';

const MaxProperFrameDelta = 1000 / 30;

export class PluginAnimationLoop {
    private lastTickT: number = 0;
    // Proper time is used to prevent animations from skipping
    // if there is a blocking operation, e.g., shader compilation
    // The drawback of this is that sometimes the animation will take
    // longer than intended, but hopefully that's a reasonable tradeoff
    private properTimeT: number = 0;

    private currentFrame: any = void 0;
    private _isAnimating = false;

    get isAnimating() {
        return this._isAnimating;
    }

    async tick(t: number, options?: { isSynchronous?: boolean, manualDraw?: boolean, animation?: PluginAnimationManager.AnimationInfo, updateControls?: boolean}) {
        await this.plugin.managers.animation.tick(t, options?.isSynchronous, options?.animation);
        this.plugin.canvas3d?.tick(t as now.Timestamp, options);

        if (isTimingMode) {
            const timerResults = this.plugin.canvas3d?.webgl.timer.resolve();
            if (timerResults) {
                for (const result of timerResults) {
                    printTimerResults([result]);
                }
            }
        }
    }

    private frame = () => {
        const t = now();
        const dt = t - this.lastTickT;
        this.lastTickT = t;
        this.properTimeT += Math.min(dt, MaxProperFrameDelta);
        this.tick(this.properTimeT);
        if (this._isAnimating) {
            this.currentFrame = requestAnimationFrame(this.frame);
        }
    };

    resetTime(t: number) {
        this.plugin.canvas3d?.resetTime(t);
    }

    start(options?: { immediate?: boolean }) {
        this.plugin.canvas3d?.resume();
        this._isAnimating = true;
        this.resetTime(0);
        this.properTimeT = 0;
        this.lastTickT = now();
        if (options?.immediate) this.frame();
        else this.currentFrame = requestAnimationFrame(this.frame);
    }

    stop(options?: { noDraw?: boolean }) {
        this._isAnimating = false;
        if (this.currentFrame !== void 0) {
            cancelAnimationFrame(this.currentFrame);
            this.currentFrame = void 0;
        }
        if (options?.noDraw) {
            this.plugin.canvas3d?.pause(options?.noDraw);
        }
    }

    constructor(private plugin: PluginContext) {

    }
}