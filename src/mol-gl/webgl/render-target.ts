/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context, createImageData } from './context'
import { idFactory } from 'mol-util/id-factory';
import { createTexture } from './texture';
import { createFramebuffer } from './framebuffer';
import { createRenderbuffer } from './renderbuffer';
import { TextureImage } from '../renderable/util';

const getNextRenderTargetId = idFactory()

export interface RenderTarget {
    readonly id: number
    readonly width: number
    readonly height: number
    readonly image: Readonly<TextureImage>

    bind: () => void
    setSize: (width: number, height: number) => void
    readBuffer: (x: number, y: number, width: number, height: number, dst: Uint8Array) => void
    getBuffer: () => Uint8Array
    getImageData: () => ImageData
    destroy: () => void
}

export function createRenderTarget (ctx: Context, _width: number, _height: number): RenderTarget {
    const { gl } = ctx

    const image: TextureImage = {
        array: new Uint8Array(_width * _height * 4),
        width: _width,
        height: _height
    }

    const targetTexture = createTexture(ctx, 'rgba', 'ubyte')
    targetTexture.load(image)

    const framebuffer = createFramebuffer(ctx)

    // attach the texture as the first color attachment
    targetTexture.attachFramebuffer(framebuffer, 'color0')

    // make a depth renderbuffer of the same size as the targetTexture
    const depthRenderbuffer = createRenderbuffer(ctx, 'depth16', 'depth', _width, _height)

    let destroyed = false

    function readBuffer(x: number, y: number, width: number, height: number, dst: Uint8Array) {
        framebuffer.bind()
        ctx.readPixels(x, y, width, height, dst)
    }

    function getBuffer() {
        readBuffer(0, 0, _width, _height, image.array)
        return image.array
    }

    return {
        id: getNextRenderTargetId(),
        get width () { return _width },
        get height () { return _height },
        image,

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

            depthRenderbuffer.setSize(_width, _height)
        },
        readBuffer,
        getBuffer,
        getImageData: () => createImageData(getBuffer(), _width, _height),
        destroy: () => {
            if (destroyed) return
            targetTexture.destroy()
            framebuffer.destroy()
            depthRenderbuffer.destroy()
            destroyed = true
        }
    }
}