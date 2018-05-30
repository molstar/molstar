/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'
import { idFactory } from 'mol-util/id-factory';
import { createTexture } from './texture';
import { createFramebuffer } from './framebuffer';
// import { createRenderbuffer } from './renderbuffer';

const getNextRenderTargetId = idFactory()

export interface RenderTarget {
    readonly id: number

    bind: () => void
    setSize: (width: number, height: number) => void
    readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => void
    destroy: () => void
}

export function createRenderTarget (ctx: Context, _width: number, _height: number): RenderTarget {
    const { gl } = ctx

    const targetTexture = createTexture(ctx, 'rgba', 'ubyte')
    targetTexture.setSize(_width, _height)

    const framebuffer = createFramebuffer(ctx)
    framebuffer.bind()

    // attach the texture as the first color attachment
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, targetTexture.texture, 0);

    // const depthRenderbuffer = createRenderbuffer(ctx)
    // depthRenderbuffer.bind()
    // // make a depth buffer and the same size as the targetTexture
    // gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, targetTexture.width, targetTexture.height);
    // gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, depthRenderbuffer);

    let destroyed = false

    return {
        id: getNextRenderTargetId(),

        bind: () => {
            framebuffer.bind()
            gl.viewport(0, 0, _width, _height);
        },
        setSize: (width: number, height: number) => {
            _width = width
            _height = height
            targetTexture.setSize(_width, _height)
        },
        readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => {
            framebuffer.bind()
            ctx.readPixels(x, y, width, height, buffer)
            ctx.unbindFramebuffer()
        },
        destroy: () => {
            if (destroyed) return
            targetTexture.destroy()
            framebuffer.destroy()
            // depthRenderbuffer.destroy()
            destroyed = true
        }
    }
}