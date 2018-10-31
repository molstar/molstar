/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'
import { idFactory } from 'mol-util/id-factory';
import { ReferenceCache, createReferenceCache } from 'mol-util/reference-cache';

const getNextFramebufferId = idFactory()

export interface Framebuffer {
    readonly id: number

    bind: () => void
    destroy: () => void
}

export function createFramebuffer (ctx: Context): Framebuffer {
    const { gl } = ctx
    const _framebuffer = gl.createFramebuffer()
    if (_framebuffer === null) {
        throw new Error('Could not create WebGL framebuffer')
    }

    let destroyed = false
    ctx.framebufferCount += 1

    return {
        id: getNextFramebufferId(),

        bind: () => gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer),
        destroy: () => {
            if (destroyed) return
            gl.deleteFramebuffer(_framebuffer)
            destroyed = true
            ctx.framebufferCount -= 1
        }
    }
}

export type FramebufferCache = ReferenceCache<Framebuffer, string, Context>

export function createFramebufferCache(): FramebufferCache {
    return createReferenceCache(
        (name: string) => name,
        (ctx: Context) => createFramebuffer(ctx),
        (framebuffer: Framebuffer) => { framebuffer.destroy() }
    )
}