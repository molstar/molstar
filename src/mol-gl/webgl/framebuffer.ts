/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '../../mol-util/id-factory';
import { GLRenderingContext, isWebGL2 } from './compat';

const getNextFramebufferId = idFactory();

function getFramebufferStatusDescription(gl: GLRenderingContext, status: number) {
    switch (status) {
        case gl.FRAMEBUFFER_COMPLETE: return 'complete';
        case gl.FRAMEBUFFER_INCOMPLETE_ATTACHMENT: return 'incomplete attachment';
        case gl.FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: return 'incomplete missing attachment';
        case gl.FRAMEBUFFER_INCOMPLETE_DIMENSIONS: return 'incomplete dimensions';
        case gl.FRAMEBUFFER_UNSUPPORTED: return 'unsupported';
    }
    if (isWebGL2(gl)) {
        switch (status) {
            case gl.FRAMEBUFFER_INCOMPLETE_MULTISAMPLE: return 'incomplete multisample';
            case gl.RENDERBUFFER_SAMPLES: return 'renderbuffer samples';
        }
    }
    return 'unknown error';
}

export function checkFramebufferStatus(gl: GLRenderingContext) {
    const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (status !== gl.FRAMEBUFFER_COMPLETE) {
        const description = getFramebufferStatusDescription(gl, status);
        throw new Error(`Framebuffer status: ${description}`);
    }
}

export interface Framebuffer {
    readonly id: number

    bind: () => void
    reset: () => void
    destroy: () => void
}

function getFramebuffer(gl: GLRenderingContext) {
    const framebuffer = gl.createFramebuffer();
    if (framebuffer === null) {
        throw new Error('Could not create WebGL framebuffer');
    }
    return framebuffer;
}

export function createFramebuffer (gl: GLRenderingContext): Framebuffer {
    let _framebuffer = getFramebuffer(gl);

    let destroyed = false;

    return {
        id: getNextFramebufferId(),
        bind: () => gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer),

        reset: () => {
            _framebuffer = getFramebuffer(gl);
        },
        destroy: () => {
            if (destroyed) return;
            gl.deleteFramebuffer(_framebuffer);
            destroyed = true;
        }
    };
}