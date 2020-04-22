/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    'image-depth': TextureImage<Uint8Array> // TODO should be Uint32Array
    'volume-uint8': TextureVolume<Uint8Array>
    'volume-float32': TextureVolume<Float32Array>
    'texture': Texture
}
export type TextureValueType = ValueOf<TextureKindValue>
export type TextureKind = keyof TextureKindValue
export type TextureType = 'ubyte' | 'ushort' | 'float'
export type TextureFormat = 'alpha' | 'rgb' | 'rgba' | 'depth'
/** Numbers are shortcuts for color attachment */
export type TextureAttachment = 'depth' | 'stencil' | 'color0' | 'color1' | 'color2' | 'color3' | 'color4' | 'color5' | 'color6' | 'color7' | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7
export type TextureFilter = 'nearest' | 'linear'

export function getTarget(gl: GLRenderingContext, kind: TextureKind): number {
    switch (kind) {
        case 'image-uint8': return gl.TEXTURE_2D;
        case 'image-float32': return gl.TEXTURE_2D;
        case 'image-depth': return gl.TEXTURE_2D;
    }
    if (isWebGL2(gl)) {
        switch (kind) {
            case 'volume-uint8': return gl.TEXTURE_3D;
            case 'volume-float32': return gl.TEXTURE_3D;
        }
    }
    throw new Error(`unknown texture kind '${kind}'`);
}

export function getFormat(gl: GLRenderingContext, format: TextureFormat, type: TextureType): number {
    switch (format) {
        case 'alpha':
            if (isWebGL2(gl) && type === 'float') return gl.RED;
            else return gl.ALPHA;
        case 'rgb': return gl.RGB;
        case 'rgba': return gl.RGBA;
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
                }
            case 'rgb':
                switch (type) {
                    case 'ubyte': return gl.RGB;
                    case 'float': return gl.RGB32F;
                }
            case 'rgba':
                switch (type) {
                    case 'ubyte': return gl.RGBA;
                    case 'float': return gl.RGBA32F;
                }
            case 'depth':
                return gl.DEPTH_COMPONENT16;
        }
    }
    return getFormat(gl, format, type);
}

export function getType(gl: GLRenderingContext, type: TextureType): number {
    switch (type) {
        case 'ubyte': return gl.UNSIGNED_BYTE;
        case 'ushort': return gl.UNSIGNED_SHORT;
        case 'float': return gl.FLOAT;
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

    getWidth: () => number
    getHeight: () => number
    getDepth: () => number

    define: (width: number, height: number, depth?: number) => void
    load: (image: TextureImage<any> | TextureVolume<any>) => void
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

type FramebufferAttachment = {
    framebuffer: Framebuffer
    attachment: TextureAttachment
    layer?: number
}

function getTexture(gl: GLRenderingContext) {
    const texture = gl.createTexture();
    if (texture === null) {
        throw new Error('Could not create WebGL texture');
    }
    return texture;
}
// export type TextureProps = { kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter }
export function createTexture(gl: GLRenderingContext, extensions: WebGLExtensions, kind: TextureKind, _format: TextureFormat, _type: TextureType, _filter: TextureFilter): Texture {
    const id = getNextTextureId();
    let texture = getTexture(gl);

    // check texture kind and type compatability
    if (
        (kind.endsWith('float32') && _type !== 'float') ||
        (kind.endsWith('uint8') && _type !== 'ubyte') ||
        (kind.endsWith('depth') && _type !== 'ushort')
    ) {
        throw new Error(`texture kind '${kind}' and type '${_type}' are incompatible`);
    }

    const target = getTarget(gl, kind);
    const filter = getFilter(gl, _filter);
    const format = getFormat(gl, _format, _type);
    const internalFormat = getInternalFormat(gl, _format, _type);
    const type = getType(gl, _type);

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

    let fba: undefined | FramebufferAttachment = undefined;

    let width = 0, height = 0, depth = 0;
    let loadedData: undefined | TextureImage<any> | TextureVolume<any>;
    let destroyed = false;

    function define(_width: number, _height: number, _depth?: number) {
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

    function load(data: TextureImage<any> | TextureVolume<any>) {
        gl.bindTexture(target, texture);
        // unpack alignment of 1 since we use textures only for data
        gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
        gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);
        gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, 0);
        if (isTexture2d(data, target, gl)) {
            gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, !!data.flipY);
            width = data.width, height = data.height;
            gl.texImage2D(target, 0, internalFormat, width, height, 0, format, type, data.array);
        } else if (isWebGL2(gl) && isTexture3d(data, target, gl)) {
            width = data.width, height = data.height, depth = data.depth;
            gl.texImage3D(target, 0, internalFormat, width, height, depth, 0, format, type, data.array);
        } else {
            throw new Error('unknown texture target');
        }
        gl.bindTexture(target, null);
        loadedData = data;
    }

    function attachFramebuffer(framebuffer: Framebuffer, attachment: TextureAttachment, layer?: number) {
        if (fba && fba.framebuffer === framebuffer && fba.attachment === attachment && fba.layer === layer) {
            return;
        }
        framebuffer.bind();
        if (target === gl.TEXTURE_2D) {
            gl.framebufferTexture2D(gl.FRAMEBUFFER, getAttachment(gl, extensions, attachment), gl.TEXTURE_2D, texture, 0);
        } else if (isWebGL2(gl) && target === gl.TEXTURE_3D) {
            if (layer === undefined) throw new Error('need `layer` to attach 3D texture');
            gl.framebufferTextureLayer(gl.FRAMEBUFFER, getAttachment(gl, extensions, attachment), texture, 0, layer);
        } else {
            throw new Error('unknown texture target');
        }
        fba = { framebuffer, attachment, layer };
    }

    return {
        id,
        target,
        format,
        internalFormat,
        type,

        getWidth: () => width,
        getHeight: () => height,
        getDepth: () => depth,

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
            fba = undefined;
        },
        reset: () => {
            texture = getTexture(gl);
            init();

            if (loadedData) {
                load(loadedData);
            } else {
                define(width, height, depth);
            }

            if (fba) {
                // TODO unclear why calling `attachFramebuffer` here does not work reliably after context loss
                // e.g. it still needs to be called in `DrawPass` to work
                fba = undefined;
            }
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
            if (spec.kind === 'texture') {
                textures[textures.length] = [k, values[k].ref.value as Texture];
            } else {
                const texture = resources.texture(spec.kind, spec.format, spec.dataType, spec.filter);
                texture.load(values[k].ref.value as TextureImage<any> | TextureVolume<any>);
                textures[textures.length] = [k, texture];
            }
        }
    });
    return textures;
}