/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from './context';
import { TextureImage, TextureVolume } from '../renderable/util';
import { ValueCell } from '../../mol-util';
import { RenderableSchema } from '../renderable/schema';
import { idFactory } from '../../mol-util/id-factory';
import { Framebuffer } from './framebuffer';
import { isWebGL2, GLRenderingContext } from './compat';
import { ValueOf } from '../../mol-util/type-helpers';
import { WebGLExtensions } from './extensions';

const getNextTextureId = idFactory();

export type TextureKindValue = {
    'image-uint8': TextureImage<Uint8Array>
    'image-float32': TextureImage<Float32Array>
    'image-float16': TextureImage<Float32Array>
    'image-int32': TextureImage<Int32Array>
    'image-depth': TextureImage<Uint8Array> // TODO should be Uint32Array
    'volume-uint8': TextureVolume<Uint8Array>
    'volume-float32': TextureVolume<Float32Array>
    'volume-float16': TextureVolume<Float32Array>
    'texture': Texture
}
export type TextureValueType = ValueOf<TextureKindValue>
export type TextureKind = keyof TextureKindValue
export type TextureType = 'ubyte' | 'ushort' | 'float' | 'fp16' | 'int'
export type TextureFormat = 'alpha' | 'rgb' | 'rgba' | 'depth'
/** Numbers are shortcuts for color attachment */
export type TextureAttachment = 'depth' | 'stencil' | 'color0' | 'color1' | 'color2' | 'color3' | 'color4' | 'color5' | 'color6' | 'color7' | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7
export type TextureFilter = 'nearest' | 'linear'

export function getTarget(gl: GLRenderingContext, kind: TextureKind): number {
    switch (kind) {
        case 'image-uint8': return gl.TEXTURE_2D;
        case 'image-float32': return gl.TEXTURE_2D;
        case 'image-float16': return gl.TEXTURE_2D;
        case 'image-depth': return gl.TEXTURE_2D;
    }
    if (isWebGL2(gl)) {
        switch (kind) {
            case 'image-int32': return gl.TEXTURE_2D;
            case 'volume-uint8': return gl.TEXTURE_3D;
            case 'volume-float32': return gl.TEXTURE_3D;
            case 'volume-float16': return gl.TEXTURE_3D;
        }
    }
    throw new Error(`unknown texture kind '${kind}'`);
}

export function getFormat(gl: GLRenderingContext, format: TextureFormat, type: TextureType): number {
    switch (format) {
        case 'alpha':
            if (isWebGL2(gl) && type === 'float') return gl.RED;
            else if (isWebGL2(gl) && type === 'int') return gl.RED_INTEGER;
            else return gl.ALPHA;
        case 'rgb':
            if (isWebGL2(gl) && type === 'int') return gl.RGB_INTEGER;
            return gl.RGB;
        case 'rgba':
            if (isWebGL2(gl) && type === 'int') return gl.RGBA_INTEGER;
            return gl.RGBA;
        case 'depth': return gl.DEPTH_COMPONENT;
    }
}

export function getInternalFormat(gl: GLRenderingContext, format: TextureFormat, type: TextureType): number {
    if (isWebGL2(gl)) {
        switch (format) {
            case 'alpha':
                switch (type) {
                    case 'ubyte': return gl.ALPHA;
                    case 'float': return gl.R32F;
                    case 'fp16': return gl.R16F;
                    case 'int': return gl.R32I;
                }
            case 'rgb':
                switch (type) {
                    case 'ubyte': return gl.RGB;
                    case 'float': return gl.RGB32F;
                    case 'fp16': return gl.RGB16F;
                    case 'int': return gl.RGB32I;
                }
            case 'rgba':
                switch (type) {
                    case 'ubyte': return gl.RGBA;
                    case 'float': return gl.RGBA32F;
                    case 'fp16': return gl.RGBA16F;
                    case 'int': return gl.RGBA32I;
                }
            case 'depth':
                switch (type) {
                    case 'ushort': return gl.DEPTH_COMPONENT16;
                    case 'float': return gl.DEPTH_COMPONENT32F;
                }
        }
    }
    return getFormat(gl, format, type);
}

function getByteCount(format: TextureFormat, type: TextureType, width: number, height: number, depth: number): number {
    const bpe = getFormatSize(format) * getTypeSize(type);
    return bpe * width * height * (depth || 1);
}

function getFormatSize(format: TextureFormat) {
    switch (format) {
        case 'alpha': return 1;
        case 'rgb': return 3;
        case 'rgba': return 4;
        case 'depth': return 4;
    }
}

function getTypeSize(type: TextureType): number {
    switch (type) {
        case 'ubyte': return 1;
        case 'ushort': return 2;
        case 'float': return 4;
        case 'fp16': return 2;
        case 'int': return 4;
    }
}

export function getType(gl: GLRenderingContext, extensions: WebGLExtensions, type: TextureType): number {
    switch (type) {
        case 'ubyte': return gl.UNSIGNED_BYTE;
        case 'ushort': return gl.UNSIGNED_SHORT;
        case 'float': return gl.FLOAT;
        case 'fp16':
            if (extensions.textureHalfFloat) return extensions.textureHalfFloat.HALF_FLOAT;
            else throw new Error('extension "texture_half_float" unavailable');
        case 'int':
            if (isWebGL2(gl)) return gl.INT;
            else throw new Error('texture type "int" requires webgl2');
    }
}

export function getFilter(gl: GLRenderingContext, type: TextureFilter): number {
    switch (type) {
        case 'nearest': return gl.NEAREST;
        case 'linear': return gl.LINEAR;
    }
}

export function getAttachment(gl: GLRenderingContext, extensions: WebGLExtensions, attachment: TextureAttachment): number {
    switch (attachment) {
        case 'depth': return gl.DEPTH_ATTACHMENT;
        case 'stencil': return gl.STENCIL_ATTACHMENT;
        case 'color0': case 0: return gl.COLOR_ATTACHMENT0;
    }
    if (extensions.drawBuffers) {
        switch (attachment) {
            case 'color1': case 1: return extensions.drawBuffers.COLOR_ATTACHMENT1;
            case 'color2': case 2: return extensions.drawBuffers.COLOR_ATTACHMENT2;
            case 'color3': case 3: return extensions.drawBuffers.COLOR_ATTACHMENT3;
            case 'color4': case 4: return extensions.drawBuffers.COLOR_ATTACHMENT4;
            case 'color5': case 5: return extensions.drawBuffers.COLOR_ATTACHMENT5;
            case 'color6': case 6: return extensions.drawBuffers.COLOR_ATTACHMENT6;
            case 'color7': case 7: return extensions.drawBuffers.COLOR_ATTACHMENT7;
        }
    }
    throw new Error('unknown texture attachment');
}

function isImage(x: TextureImage<any> | TextureVolume<any> | HTMLImageElement): x is HTMLImageElement {
    return typeof HTMLImageElement !== 'undefined' && (x instanceof HTMLImageElement);
}

function isTexture2d(x: TextureImage<any> | TextureVolume<any>, target: number, gl: GLRenderingContext): x is TextureImage<any> {
    return target === gl.TEXTURE_2D;
}

function isTexture3d(x: TextureImage<any> | TextureVolume<any>, target: number, gl: WebGL2RenderingContext): x is TextureImage<any> {
    return target === gl.TEXTURE_3D;
}

export interface Texture {
    readonly id: number
    readonly target: number
    readonly format: number
    readonly internalFormat: number
    readonly type: number
    readonly filter: number

    getWidth: () => number
    getHeight: () => number
    getDepth: () => number

    getByteCount: () => number

    define: (width: number, height: number, depth?: number) => void
    /**
     * The `sub` option requires an existing allocation on the GPU, that is, either
     * `define` or `load` without `sub` must have been called before.
     */
    load: (image: TextureImage<any> | TextureVolume<any> | HTMLImageElement, sub?: boolean) => void
    bind: (id: TextureId) => void
    unbind: (id: TextureId) => void
    /** Use `layer` to attach a z-slice of a 3D texture */
    attachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment, layer?: number) => void
    detachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment) => void

    reset: () => void
    destroy: () => void
}

export type TextureId = 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15

export type TextureValues = { [k: string]: ValueCell<TextureValueType> }
export type Textures = [string, Texture][]

function getTexture(gl: GLRenderingContext) {
    const texture = gl.createTexture();
    if (texture === null) {
        throw new Error('Could not create WebGL texture');
    }
    return texture;
}

export function createTexture(gl: GLRenderingContext, extensions: WebGLExtensions, kind: TextureKind, _format: TextureFormat, _type: TextureType, _filter: TextureFilter): Texture {
    const id = getNextTextureId();
    let texture = getTexture(gl);

    // check texture kind and type compatability
    if (
        (kind.endsWith('float32') && _type !== 'float') ||
        (kind.endsWith('float16') && _type !== 'fp16') ||
        (kind.endsWith('uint8') && _type !== 'ubyte') ||
        (kind.endsWith('int32') && _type !== 'int') ||
        (kind.endsWith('depth') && _type !== 'ushort' && _type !== 'float')
    ) {
        throw new Error(`texture kind '${kind}' and type '${_type}' are incompatible`);
    }

    const target = getTarget(gl, kind);
    const filter = getFilter(gl, _filter);
    const format = getFormat(gl, _format, _type);
    const internalFormat = getInternalFormat(gl, _format, _type);
    const type = getType(gl, extensions, _type);

    function init() {
        gl.bindTexture(target, texture);
        gl.texParameteri(target, gl.TEXTURE_MAG_FILTER, filter);
        gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, filter);
        // clamp-to-edge needed for non-power-of-two textures in webgl
        gl.texParameteri(target, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(target, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.bindTexture(target, null);
    }
    init();

    let width = 0, height = 0, depth = 0;
    let loadedData: undefined | TextureImage<any> | TextureVolume<any> | HTMLImageElement;
    let destroyed = false;

    function define(_width: number, _height: number, _depth?: number) {
        if (width === _width && height === _height && depth === (_depth || 0)) return;

        width = _width, height = _height, depth = _depth || 0;
        gl.bindTexture(target, texture);
        if (target === gl.TEXTURE_2D) {
            gl.texImage2D(target, 0, internalFormat, width, height, 0, format, type, null);
        } else if (isWebGL2(gl) && target === gl.TEXTURE_3D && depth !== undefined) {
            gl.texImage3D(target, 0, internalFormat, width, height, depth, 0, format, type, null);
        } else {
            throw new Error('unknown texture target');
        }
    }

    function load(data: TextureImage<any> | TextureVolume<any> | HTMLImageElement, sub = false) {
        gl.bindTexture(target, texture);
        // unpack alignment of 1 since we use textures only for data
        gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
        gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);
        gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, 0);
        if (isImage(data)) {
            gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
            gl.bindTexture(gl.TEXTURE_2D, texture);
            gl.texImage2D(gl.TEXTURE_2D, 0, internalFormat, format, type, data);
        } else if (isTexture2d(data, target, gl)) {
            const _filter = data.filter ? getFilter(gl, data.filter) : filter;
            gl.texParameteri(target, gl.TEXTURE_MAG_FILTER, _filter);
            gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, _filter);
            gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, !!data.flipY);
            if (sub) {
                gl.texSubImage2D(target, 0, 0, 0, data.width, data.height, format, type, data.array);
            } else {
                width = data.width, height = data.height;
                gl.texImage2D(target, 0, internalFormat, width, height, 0, format, type, data.array);
            }
        } else if (isWebGL2(gl) && isTexture3d(data, target, gl)) {
            gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
            if (sub) {
                gl.texSubImage3D(target, 0, 0, 0, 0, data.width, data.height, data.depth, format, type, data.array);
            } else {
                width = data.width, height = data.height, depth = data.depth;
                gl.texImage3D(target, 0, internalFormat, width, height, depth, 0, format, type, data.array);
            }
        } else {
            throw new Error('unknown texture target');
        }
        gl.bindTexture(target, null);
        loadedData = data;
    }

    function attachFramebuffer(framebuffer: Framebuffer, attachment: TextureAttachment, layer?: number) {
        framebuffer.bind();
        if (target === gl.TEXTURE_2D) {
            gl.framebufferTexture2D(gl.FRAMEBUFFER, getAttachment(gl, extensions, attachment), gl.TEXTURE_2D, texture, 0);
        } else if (isWebGL2(gl) && target === gl.TEXTURE_3D) {
            if (layer === undefined) throw new Error('need `layer` to attach 3D texture');
            gl.framebufferTextureLayer(gl.FRAMEBUFFER, getAttachment(gl, extensions, attachment), texture, 0, layer);
        } else {
            throw new Error('unknown/unsupported texture target');
        }
    }

    return {
        id,
        target,
        format,
        internalFormat,
        type,
        filter,

        getWidth: () => width,
        getHeight: () => height,
        getDepth: () => depth,

        getByteCount: () => getByteCount(_format, _type, width, height, depth),

        define,
        load,
        bind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id);
            gl.bindTexture(target, texture);
        },
        unbind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id);
            gl.bindTexture(target, null);
        },
        attachFramebuffer,
        detachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment) => {
            framebuffer.bind();
            if (target === gl.TEXTURE_2D) {
                gl.framebufferTexture2D(gl.FRAMEBUFFER, getAttachment(gl, extensions, attachment), gl.TEXTURE_2D, null, 0);
            } else if (isWebGL2(gl) && target === gl.TEXTURE_3D) {
                gl.framebufferTextureLayer(gl.FRAMEBUFFER, getAttachment(gl, extensions, attachment), null, 0, 0);
            } else {
                throw new Error('unknown texture target');
            }
        },
        reset: () => {
            texture = getTexture(gl);
            init();

            const [_width, _height, _depth] = [width, height, depth];
            width = 0, height = 0, depth = 0; // set to zero to trigger resize
            define(_width, _height, _depth);
            if (loadedData) load(loadedData);
        },
        destroy: () => {
            if (destroyed) return;
            gl.deleteTexture(texture);
            destroyed = true;
        }
    };
}

export function createTextures(ctx: WebGLContext, schema: RenderableSchema, values: TextureValues) {
    const { resources } = ctx;
    const textures: Textures = [];
    Object.keys(schema).forEach(k => {
        const spec = schema[k];
        if (spec.type === 'texture') {
            const value = values[k];
            if (value) {
                if (spec.kind === 'texture') {
                    textures[textures.length] = [k, value.ref.value as Texture];
                } else {
                    const texture = resources.texture(spec.kind, spec.format, spec.dataType, spec.filter);
                    texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
                    textures[textures.length] = [k, texture];
                }
            }
        }
    });
    return textures;
}

/**
 * Loads an image from a url to a textures and triggers update asynchronously.
 * This will not work on node.js without a polyfill for `HTMLImageElement`.
 */
export function loadImageTexture(src: string, cell: ValueCell<Texture>, texture: Texture) {
    const img = new Image();
    img.onload = function () {
        texture.load(img);
        ValueCell.update(cell, texture);
    };
    img.src = src;
}

//

export function createNullTexture(gl?: GLRenderingContext): Texture {
    const target = gl?.TEXTURE_2D ?? 3553;
    return {
        id: getNextTextureId(),
        target,
        format: 0,
        internalFormat: 0,
        type: 0,
        filter: 0,

        getWidth: () => 0,
        getHeight: () => 0,
        getDepth: () => 0,
        getByteCount: () => 0,

        define: () => {},
        load: () => {},
        bind: (id: TextureId) => {
            if (gl) {
                gl.activeTexture(gl.TEXTURE0 + id);
                gl.bindTexture(target, null);
            }
        },
        unbind: (id: TextureId) => {
            if (gl) {
                gl.activeTexture(gl.TEXTURE0 + id);
                gl.bindTexture(target, null);
            }
        },
        attachFramebuffer: () => {},
        detachFramebuffer: () => {},

        reset: () => {},
        destroy: () => {},
    };
}
