/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from './context';
import { now } from '../mol-util/now';

export class PluginAnimationLoop {
    private currentFrame: any = void 0;
    private _isAnimating = false;

    tick(t: number, options?: { isSynchronous?: boolean, manualDraw?: boolean }) {
        this.plugin.canvas3d?.tick(t as now.Timestamp, options);
        return this.plugin.managers.animation.tick(t, options?.isSynchronous);
    }

    private frame = () => {
        this.tick(now());
        if (this._isAnimating) {
            this.currentFrame = requestAnimationFrame(this.frame);
        }
    }

    resetTime(t: number = now()) {
        this.plugin.canvas3d?.resetTime(t);
    }

    start() {
        this._isAnimating = true;
        this.resetTime();
        this.currentFrame = requestAnimationFrame(this.frame);
    }

    stop() {
        this._isAnimating = false;
        if (this.currentFrame !== void 0) {
            cancelAnimationFrame(this.currentFrame);
            this.currentFrame = void 0;
        }
    }

    constructor(private plugin: PluginContext) {

    }
}