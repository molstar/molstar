/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createContext } from '../webgl/context';

export function getGLContext(width: number, height: number) {
    const gl = require('gl')(width, height, {
        alpha: true,
        depth: true,
        premultipliedAlpha: true,
        preserveDrawingBuffer: true,
        antialias: true,
    });
    return createContext(gl);
}

export function tryGetGLContext(width: number, height: number, requiredExtensions?: { fragDepth?: boolean, textureFloat?: boolean }) {
    try {
        const ctx = getGLContext(width, height);
        if (requiredExtensions?.fragDepth && !ctx.extensions.fragDepth) return;
        if (requiredExtensions?.textureFloat && !ctx.extensions.textureFloat) return;
        return ctx;
    } catch (e) {
        return;
    }
}