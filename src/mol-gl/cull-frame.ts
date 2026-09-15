/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * 
 * @author Taylor Hoffmann <taylor@hoffmann.io>
 */

import { WebGLStats } from './webgl/context';

let cullFrameId = 0;

export function beginCullFrame(stats?: WebGLStats) {
    cullFrameId++;
    if (stats) {
        stats.cull.cached = 0;
        stats.cull.computed = 0;
        stats.uniforms.uploaded = 0;
        stats.uniforms.skipped = 0;
    }
}

export function getCullFrameId() {
    return cullFrameId;
}
