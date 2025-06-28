/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '../../mol-util/id-factory';
import { Texture, TextureFilter } from './texture';
import { Framebuffer } from './framebuffer';
import { WebGLResources } from './resources';
import { GLRenderingContext, isWebGL2 } from './compat';
import { Renderbuffer } from './renderbuffer';

const getNextRenderTargetId = idFactory();

export interface RenderTarget {
    readonly id: number
    readonly texture: Texture
    readonly framebuffer: Framebuffer
    readonly depthRenderbuffer: Renderbuffer | null

    getWidth: () => number
    getHeight: () => number
    /** binds framebuffer */
    bind: () => void
    setSize: (width: number, height: number) => void
    reset: () => void
    destroy: () => void
}

export function createRenderTarget(gl: GLRenderingContext, resources: WebGLResources, _width: number, _height: number, depth = true, type: 'uint8' | 'float32' | 'fp16' = 'uint8', filter: TextureFilter = 'nearest', format: 'rgba' | 'alpha' = 'rgba'): RenderTarget {

    if (format === 'alpha' && !isWebGL2(gl)) {
        throw new Error('cannot render to alpha format in webgl1');
    }

    const framebuffer = resources.framebuffer();
    const targetTexture = type === 'fp16'
        ? resources.texture('image-float16', format, 'fp16', filter)
        : type === 'float32'
            ? resources.texture('image-float32', format, 'float', filter)
            : resources.texture('image-uint8', format, 'ubyte', filter);
    // make a depth renderbuffer of the same size as the targetTexture
    const depthRenderbuffer = !depth
        ? null
        : isWebGL2(gl)
            ? resources.renderbuffer('depth32f', 'depth', _width, _height)
            : resources.renderbuffer('depth16', 'depth', _width, _height);

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
        depthRenderbuffer,

        getWidth: () => _width,
        getHeight: () => _height,
        bind: () => {
            framebuffer.bind();
        },
        setSize: (width: number, height: number) => {
            if (_width === width && _height === height) {
                return;
            }

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
