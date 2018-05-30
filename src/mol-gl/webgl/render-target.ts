/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context, createImageData } from './context'
import { idFactory } from 'mol-util/id-factory';
import { createTexture } from './texture';
import { createFramebuffer } from './framebuffer';
// import { createRenderbuffer } from './renderbuffer';

const getNextRenderTargetId = idFactory()

export interface RenderTarget {
    readonly id: number

    bind: () => void
    setSize: (width: number, height: number) => void
    getImageData: () => ImageData
    destroy: () => void
}

export function createRenderTarget (ctx: Context, _width: number, _height: number): RenderTarget {
    const { gl } = ctx

    const image = {
        array: new Uint8Array(_width * _height * 4),
        width: _width,
        height: _height
    }

    const targetTexture = createTexture(ctx, 'rgba', 'ubyte')
    targetTexture.load(image)

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
            image.array = new Uint8Array(_width * _height * 4)
            image.width = _width
            image.height = _height
            targetTexture.load(image)
        },
        getImageData: () => {
            framebuffer.bind()
            ctx.readPixels(0, 0, _width, _height, image.array)
            return createImageData(image.array, _width, _height)
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