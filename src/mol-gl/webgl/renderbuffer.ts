/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from './context'
import { idFactory } from 'mol-util/id-factory';

const getNextRenderbufferId = idFactory()

export type RenderbufferFormat = 'depth16' | 'stencil8' | 'rgba4' | 'depth-stencil'
export type RenderbufferAttachment = 'depth' | 'stencil' | 'depth-stencil' | 'color0'

export function getFormat(ctx: WebGLContext, format: RenderbufferFormat) {
    const { gl } = ctx
    switch (format) {
        case 'depth16': return gl.DEPTH_COMPONENT16
        case 'stencil8': return gl.STENCIL_INDEX8
        case 'rgba4': return gl.RGBA4
        case 'depth-stencil': return gl.DEPTH_STENCIL
    }
}

export function getAttachment(ctx: WebGLContext, attachment: RenderbufferAttachment) {
    const { gl } = ctx
    switch (attachment) {
        case 'depth': return gl.DEPTH_ATTACHMENT
        case 'stencil': return gl.STENCIL_ATTACHMENT
        case 'depth-stencil': return gl.DEPTH_STENCIL_ATTACHMENT
        case 'color0': return gl.COLOR_ATTACHMENT0
    }
}

export interface Renderbuffer {
    readonly id: number

    bind: () => void
    setSize: (width: number, height: number) => void

    destroy: () => void
}

export function createRenderbuffer (ctx: WebGLContext, format: RenderbufferFormat, attachment: RenderbufferAttachment, _width: number, _height: number): Renderbuffer {
    const { gl } = ctx
    const _renderbuffer = gl.createRenderbuffer()
    if (_renderbuffer === null) {
        throw new Error('Could not create WebGL renderbuffer')
    }

    const bind = () => gl.bindRenderbuffer(gl.RENDERBUFFER, _renderbuffer)
    const _format = getFormat(ctx, format)
    const _attachment = getAttachment(ctx, attachment)

    bind()
    gl.renderbufferStorage(gl.RENDERBUFFER, _format, _width, _height)
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, _attachment, gl.RENDERBUFFER, _renderbuffer)

    let destroyed = false
    ctx.renderbufferCount += 1

    return {
        id: getNextRenderbufferId(),

        bind,
        setSize: (_width: number, _height: number) => {
            bind()
            gl.renderbufferStorage(gl.RENDERBUFFER, _format, _width, _height)
        },

        destroy: () => {
            if (destroyed) return
            gl.deleteRenderbuffer(_renderbuffer)
            destroyed = true
            ctx.framebufferCount -= 1
        }
    }
}