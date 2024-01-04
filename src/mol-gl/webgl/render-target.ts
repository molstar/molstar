/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '../../mol-util/id-factory';
import { createNullTexture, Texture, TextureFilter } from './texture';
import { createNullFramebuffer, Framebuffer } from './framebuffer';
import { WebGLResources } from './resources';
import { GLRenderingContext, isWebGL2 } from './compat';

const getNextRenderTargetId = idFactory();

export interface RenderTarget {
    readonly id: number
    readonly framebuffer: Framebuffer
    readonly texture: Texture
    readonly depthTexture: Texture | null

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
    const colorTexture = type === 'fp16'
        ? resources.texture('image-float16', format, 'fp16', filter)
        : type === 'float32'
            ? resources.texture('image-float32', format, 'float', filter)
            : resources.texture('image-uint8', format, 'ubyte', filter);
    const depthTexture = !depth
        ? null
        : isWebGL2(gl)
            ? resources.texture('image-depth', 'depth', 'float', 'nearest')
            : resources.texture('image-depth', 'depth', 'ushort', 'nearest');

    function init() {
        colorTexture.define(_width, _height);
        colorTexture.attachFramebuffer(framebuffer, 'color0');
        if (depthTexture) {
            depthTexture.define(_width, _height);
            depthTexture.attachFramebuffer(framebuffer, 'depth');
        }
    }
    init();

    let destroyed = false;

    return {
        id: getNextRenderTargetId(),
        framebuffer,
        texture: colorTexture,
        depthTexture,

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
            colorTexture.define(_width, _height);
            if (depthTexture) depthTexture.define(_width, _height);
        },
        reset: () => {
            init();
        },
        destroy: () => {
            if (destroyed) return;
            framebuffer.destroy();
            colorTexture.destroy();
            if (depthTexture) depthTexture.destroy();
            destroyed = true;
        }
    };
}

//

export function createNullRenderTarget(gl: GLRenderingContext): RenderTarget {
    return {
        id: getNextRenderTargetId(),
        framebuffer: createNullFramebuffer(),
        texture: createNullTexture(gl),
        depthTexture: null,

        getWidth: () => 0,
        getHeight: () => 0,
        bind: () => {
            gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        },
        setSize: () => {},
        reset: () => {},
        destroy: () => {}
    };
}