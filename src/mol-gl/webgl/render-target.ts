/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readPixels } from './context';
import { idFactory } from '../../mol-util/id-factory';
import { Texture } from './texture';
import { Framebuffer } from './framebuffer';
import { TextureImage } from '../renderable/util';
import { Mutable } from '../../mol-util/type-helpers';
import { PixelData } from '../../mol-util/image';
import { WebGLResources } from './resources';
import { GLRenderingContext } from './compat';

const getNextRenderTargetId = idFactory();

export interface RenderTarget {
    readonly id: number
    readonly image: TextureImage<any>
    readonly texture: Texture
    readonly framebuffer: Framebuffer

    getWidth: () => number
    getHeight: () => number
    /** binds framebuffer and sets viewport to rendertarget's width and height */
    bind: () => void
    setSize: (width: number, height: number) => void
    readBuffer: (x: number, y: number, width: number, height: number, dst: Uint8Array) => void
    getBuffer: () => Uint8Array
    getPixelData: () => PixelData
    reset: () => void
    destroy: () => void
}

export function createRenderTarget(gl: GLRenderingContext, resources: WebGLResources, _width: number, _height: number): RenderTarget {
    const image: Mutable<TextureImage<Uint8Array>> = {
        array: new Uint8Array(_width * _height * 4),
        width: _width,
        height: _height
    };

    const framebuffer = resources.framebuffer();
    const targetTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    // make a depth renderbuffer of the same size as the targetTexture
    const depthRenderbuffer = resources.renderbuffer('depth16', 'depth', _width, _height);

    function init() {
        targetTexture.load(image);
        targetTexture.attachFramebuffer(framebuffer, 'color0');
        depthRenderbuffer.attachFramebuffer(framebuffer);
    }
    init();

    let destroyed = false;

    function readBuffer(x: number, y: number, width: number, height: number, dst: Uint8Array) {
        framebuffer.bind();
        gl.viewport(0, 0, _width, _height);
        readPixels(gl, x, y, width, height, dst);
    }

    function getBuffer() {
        readBuffer(0, 0, _width, _height, image.array);
        return image.array;
    }

    return {
        id: getNextRenderTargetId(),
        image,
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
            image.array = new Uint8Array(_width * _height * 4);
            image.width = _width;
            image.height = _height;
            targetTexture.load(image);
            depthRenderbuffer.setSize(_width, _height);
        },
        readBuffer,
        getBuffer,
        getPixelData: () => PixelData.flipY(PixelData.create(getBuffer(), _width, _height)),
        reset: () => {
            init();
        },
        destroy: () => {
            if (destroyed) return;
            targetTexture.destroy();
            framebuffer.destroy();
            depthRenderbuffer.destroy();
            destroyed = true;
        }
    };
}