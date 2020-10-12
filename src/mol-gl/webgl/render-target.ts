/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '../../mol-util/id-factory';
import { Texture, TextureFilter } from './texture';
import { Framebuffer } from './framebuffer';
import { WebGLResources } from './resources';
import { GLRenderingContext } from './compat';

const getNextRenderTargetId = idFactory();

export interface RenderTarget {
    readonly id: number
    readonly texture: Texture
    readonly framebuffer: Framebuffer

    getWidth: () => number
    getHeight: () => number
    /** binds framebuffer and sets viewport to rendertarget's width and height */
    bind: () => void
    setSize: (width: number, height: number) => void
    reset: () => void
    destroy: () => void
}

export function createRenderTarget(gl: GLRenderingContext, resources: WebGLResources, _width: number, _height: number, depth = true, type: 'uint8' | 'float32' = 'uint8', filter: TextureFilter = 'nearest'): RenderTarget {

    const framebuffer = resources.framebuffer();
    const targetTexture = type === 'float32'
        ? resources.texture('image-float32', 'rgba', 'float', filter)
        : resources.texture('image-uint8', 'rgba', 'ubyte', filter);
    // make a depth renderbuffer of the same size as the targetTexture
    const depthRenderbuffer = depth
        ? resources.renderbuffer('depth16', 'depth', _width, _height)
        : null;

    function init() {
        targetTexture.define(_width, _height);
        targetTexture.attachFramebuffer(framebuffer, 'color0');
        if (depthRenderbuffer) depthRenderbuffer.attachFramebuffer(framebuffer);
    }
    init();

    let destroyed = false;

    return {
        id: getNextRenderTargetId(),
        texture: targetTexture,
        framebuffer,

        getWidth: () => _width,
        getHeight: () => _height,
        bind: () => {
            framebuffer.bind();
            gl.viewport(0, 0, _width, _height);
        },
        setSize: (width: number, height: number) => {
            _width = width;
            _height = height;
            targetTexture.define(_width, _height);
            if (depthRenderbuffer) depthRenderbuffer.setSize(_width, _height);
        },
        reset: () => {
            init();
        },
        destroy: () => {
            if (destroyed) return;
            targetTexture.destroy();
            framebuffer.destroy();
            if (depthRenderbuffer) depthRenderbuffer.destroy();
            destroyed = true;
        }
    };
}
