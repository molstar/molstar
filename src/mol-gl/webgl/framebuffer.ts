/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLStats } from './context'
import { idFactory } from 'mol-util/id-factory';
import { ReferenceCache, createReferenceCache } from 'mol-util/reference-cache';
import { GLRenderingContext } from './compat';

const getNextFramebufferId = idFactory()

export interface Framebuffer {
    readonly id: number

    bind: () => void
    destroy: () => void
}

export function createFramebuffer (gl: GLRenderingContext, stats: WebGLStats): Framebuffer {
    const _framebuffer = gl.createFramebuffer()
    if (_framebuffer === null) {
        throw new Error('Could not create WebGL framebuffer')
    }

    let destroyed = false
    stats.framebufferCount += 1

    return {
        id: getNextFramebufferId(),

        bind: () => gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer),
        destroy: () => {
            if (destroyed) return
            gl.deleteFramebuffer(_framebuffer)
            destroyed = true
            stats.framebufferCount -= 1
        }
    }
}

export type FramebufferCache = ReferenceCache<Framebuffer, string>

export function createFramebufferCache(gl: GLRenderingContext, stats: WebGLStats): FramebufferCache {
    return createReferenceCache(
        (name: string) => name,
        () => createFramebuffer(gl, stats),
        (framebuffer: Framebuffer) => { framebuffer.destroy() }
    )
}