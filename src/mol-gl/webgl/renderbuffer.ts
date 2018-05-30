/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'
import { idFactory } from 'mol-util/id-factory';

const getNextRenderbufferId = idFactory()

export interface Renderbuffer {
    readonly id: number

    bind: () => void
    destroy: () => void
}

export function createRenderbuffer (ctx: Context): Renderbuffer {
    const { gl } = ctx
    const _renderbuffer = gl.createRenderbuffer()
    if (_renderbuffer === null) {
        throw new Error('Could not create WebGL renderbuffer')
    }

    let destroyed = false
    ctx.renderbufferCount += 1

    return {
        id: getNextRenderbufferId(),

        bind: () => {
            gl.bindRenderbuffer(gl.RENDERBUFFER, _renderbuffer)
        },

        destroy: () => {
            if (destroyed) return
            gl.deleteRenderbuffer(_renderbuffer)
            destroyed = true
            ctx.framebufferCount -= 1
        }
    }
}