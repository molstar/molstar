/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '../../mol-util/id-factory';
import { GLRenderingContext } from './compat';
import { Framebuffer, checkFramebufferStatus } from './framebuffer';
import { isDebugMode } from '../../mol-util/debug';

const getNextRenderbufferId = idFactory();

export type RenderbufferFormat = 'depth16' | 'stencil8' | 'rgba4' | 'depth-stencil'
export type RenderbufferAttachment = 'depth' | 'stencil' | 'depth-stencil' | 'color0'

export function getFormat(gl: GLRenderingContext, format: RenderbufferFormat) {
    switch (format) {
        case 'depth16': return gl.DEPTH_COMPONENT16;
        case 'stencil8': return gl.STENCIL_INDEX8;
        case 'rgba4': return gl.RGBA4;
        case 'depth-stencil': return gl.DEPTH_STENCIL;
    }
}

export function getAttachment(gl: GLRenderingContext, attachment: RenderbufferAttachment) {
    switch (attachment) {
        case 'depth': return gl.DEPTH_ATTACHMENT;
        case 'stencil': return gl.STENCIL_ATTACHMENT;
        case 'depth-stencil': return gl.DEPTH_STENCIL_ATTACHMENT;
        case 'color0': return gl.COLOR_ATTACHMENT0;
    }
}

export interface Renderbuffer {
    readonly id: number

    bind: () => void
    attachFramebuffer: (framebuffer: Framebuffer) => void
    setSize: (width: number, height: number) => void
    reset: () => void
    destroy: () => void
}

function getRenderbuffer(gl: GLRenderingContext) {
    const renderbuffer = gl.createRenderbuffer();
    if (renderbuffer === null) {
        throw new Error('Could not create WebGL renderbuffer');
    }
    return renderbuffer;
}

export function createRenderbuffer (gl: GLRenderingContext, format: RenderbufferFormat, attachment: RenderbufferAttachment, _width: number, _height: number): Renderbuffer {
    let _renderbuffer = getRenderbuffer(gl);

    const bind = () => gl.bindRenderbuffer(gl.RENDERBUFFER, _renderbuffer);
    const _format = getFormat(gl, format);
    const _attachment = getAttachment(gl, attachment);

    function init() {
        bind();
        gl.renderbufferStorage(gl.RENDERBUFFER, _format, _width, _height);
    }
    init();

    let destroyed = false;

    return {
        id: getNextRenderbufferId(),

        bind,
        attachFramebuffer: (framebuffer: Framebuffer) => {
            framebuffer.bind();
            bind();
            gl.framebufferRenderbuffer(gl.FRAMEBUFFER, _attachment, gl.RENDERBUFFER, _renderbuffer);
            if (isDebugMode) checkFramebufferStatus(gl);
        },
        setSize: (width: number, height: number) => {
            _width = width;
            _height = height;
            init();
        },
        reset: () => {
            _renderbuffer = getRenderbuffer(gl);
            init();
        },
        destroy: () => {
            if (destroyed) return;
            gl.deleteRenderbuffer(_renderbuffer);
            destroyed = true;
        }
    };
}