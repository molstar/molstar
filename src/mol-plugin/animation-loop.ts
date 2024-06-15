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

export class PluginAnimationLoop {
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
        this.tick(now());
        if (this._isAnimating) {
            this.currentFrame = requestAnimationFrame(this.frame);
        }
    };

    resetTime(t: number = now()) {
        this.plugin.canvas3d?.resetTime(t);
    }

    start(options?: { immediate?: boolean }) {
        this.plugin.canvas3d?.resume();
        this._isAnimating = true;
        this.resetTime();
        // TODO: should immediate be the default mode?
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