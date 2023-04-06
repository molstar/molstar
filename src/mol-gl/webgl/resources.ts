/**
 * Copyright (c) 2020-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ProgramProps, createProgram, Program } from './program';
import { ShaderType, createShader, Shader, ShaderProps } from './shader';
import { GLRenderingContext } from './compat';
import { Framebuffer, createFramebuffer } from './framebuffer';
import { WebGLExtensions } from './extensions';
import { WebGLState } from './state';
import { AttributeBuffer, UsageHint, ArrayType, AttributeItemSize, createAttributeBuffer, ElementsBuffer, createElementsBuffer, ElementsType, AttributeBuffers } from './buffer';
import { createReferenceCache, ReferenceItem } from '../../mol-util/reference-cache';
import { WebGLStats } from './context';
import { hashString, hashFnv32a } from '../../mol-data/util';
import { DefineValues, ShaderCode } from '../shader-code';
import { RenderableSchema } from '../renderable/schema';
import { createRenderbuffer, Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from './renderbuffer';
import { Texture, TextureKind, TextureFormat, TextureType, TextureFilter, createTexture, CubeFaces, createCubeTexture } from './texture';
import { VertexArray, createVertexArray } from './vertex-array';

function defineValueHash(v: boolean | number | string): number {
    return typeof v === 'boolean' ? (v ? 1 : 0) :
        typeof v === 'number' ? (v * 10000) : hashString(v);
}

function wrapCached<T extends Resource>(resourceItem: ReferenceItem<T>) {
    const wrapped = {
        ...resourceItem.value,
        destroy: () => {
            resourceItem.free();
        }
    };

    return wrapped;
}

//

interface Resource {
    reset: () => void
    destroy: () => void
}

type ResourceName = keyof WebGLStats['resourceCounts']

type ByteCounts = {
    texture: number
    attribute: number
    elements: number
}

export interface WebGLResources {
    attribute: (array: ArrayType, itemSize: AttributeItemSize, divisor: number, usageHint?: UsageHint) => AttributeBuffer
    elements: (array: ElementsType, usageHint?: UsageHint) => ElementsBuffer
    framebuffer: () => Framebuffer
    program: (defineValues: DefineValues, shaderCode: ShaderCode, schema: RenderableSchema) => Program
    renderbuffer: (format: RenderbufferFormat, attachment: RenderbufferAttachment, width: number, height: number) => Renderbuffer
    shader: (type: ShaderType, source: string) => Shader
    texture: (kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter) => Texture,
    cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: () => void) => Texture,
    vertexArray: (program: Program, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer) => VertexArray,

    getByteCounts: () => ByteCounts

    reset: () => void
    destroy: () => void
}

export function createResources(gl: GLRenderingContext, state: WebGLState, stats: WebGLStats, extensions: WebGLExtensions): WebGLResources {
    const sets: { [k in ResourceName]: Set<Resource> } = {
        attribute: new Set<Resource>(),
        elements: new Set<Resource>(),
        framebuffer: new Set<Resource>(),
        program: new Set<Resource>(),
        renderbuffer: new Set<Resource>(),
        shader: new Set<Resource>(),
        texture: new Set<Resource>(),
        cubeTexture: new Set<Resource>(),
        vertexArray: new Set<Resource>(),
    };

    function wrap<T extends Resource>(name: ResourceName, resource: T) {
        sets[name].add(resource);
        stats.resourceCounts[name] += 1;
        return {
            ...resource,
            destroy: () => {
                resource.destroy();
                sets[name].delete(resource);
                stats.resourceCounts[name] -= 1;
            }
        };
    }

    const shaderCache = createReferenceCache(
        (props: ShaderProps) => JSON.stringify(props),
        (props: ShaderProps) => wrap('shader', createShader(gl, props)),
        (shader: Shader) => { shader.destroy(); }
    );

    function getShader(type: ShaderType, source: string) {
        return wrapCached(shaderCache.get({ type, source }));
    }

    const programCache = createReferenceCache(
        (props: ProgramProps) => {
            const array = [props.shaderCode.id];
            const variant = (props.defineValues.dRenderVariant?.ref.value || '') as string;
            Object.keys(props.defineValues).forEach(k => {
                if (!props.shaderCode.ignoreDefine?.(k, variant, props.defineValues)) {
                    array.push(hashString(k), defineValueHash(props.defineValues[k].ref.value));
                }
            });
            return hashFnv32a(array).toString();
        },
        (props: ProgramProps) => wrap('program', createProgram(gl, state, extensions, getShader, props)),
        (program: Program) => { program.destroy(); }
    );

    return {
        attribute: (array: ArrayType, itemSize: AttributeItemSize, divisor: number, usageHint?: UsageHint) => {
            return wrap('attribute', createAttributeBuffer(gl, state, extensions, array, itemSize, divisor, usageHint));
        },
        elements: (array: ElementsType, usageHint?: UsageHint) => {
            return wrap('elements', createElementsBuffer(gl, array, usageHint));
        },
        framebuffer: () => {
            return wrap('framebuffer', createFramebuffer(gl));
        },
        program: (defineValues: DefineValues, shaderCode: ShaderCode, schema: RenderableSchema) => {
            return wrapCached(programCache.get({ defineValues, shaderCode, schema }));
        },
        renderbuffer: (format: RenderbufferFormat, attachment: RenderbufferAttachment, width: number, height: number) => {
            return wrap('renderbuffer', createRenderbuffer(gl, format, attachment, width, height));
        },
        shader: getShader,
        texture: (kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter) => {
            return wrap('texture', createTexture(gl, extensions, kind, format, type, filter));
        },
        cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: () => void) => {
            return wrap('cubeTexture', createCubeTexture(gl, faces, mipmaps, onload));
        },
        vertexArray: (program: Program, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer) => {
            return wrap('vertexArray', createVertexArray(gl, extensions, program, attributeBuffers, elementsBuffer));
        },

        getByteCounts: () => {
            let texture = 0;
            sets.texture.forEach(r => {
                texture += (r as Texture).getByteCount();
            });
            sets.cubeTexture.forEach(r => {
                texture += (r as Texture).getByteCount();
            });

            let attribute = 0;
            sets.attribute.forEach(r => {
                attribute += (r as AttributeBuffer).length * 4;
            });

            let elements = 0;
            sets.elements.forEach(r => {
                elements += (r as ElementsBuffer).length * 4;
            });

            return { texture, attribute, elements };
        },

        reset: () => {
            sets.attribute.forEach(r => r.reset());
            sets.elements.forEach(r => r.reset());
            sets.framebuffer.forEach(r => r.reset());
            sets.renderbuffer.forEach(r => r.reset());
            sets.shader.forEach(r => r.reset());
            sets.program.forEach(r => r.reset());
            sets.vertexArray.forEach(r => r.reset());
            sets.texture.forEach(r => r.reset());
        },
        destroy: () => {
            sets.attribute.forEach(r => r.destroy());
            sets.elements.forEach(r => r.destroy());
            sets.framebuffer.forEach(r => r.destroy());
            sets.renderbuffer.forEach(r => r.destroy());
            sets.shader.forEach(r => r.destroy());
            sets.program.forEach(r => r.destroy());
            sets.vertexArray.forEach(r => r.destroy());
            sets.texture.forEach(r => r.destroy());

            shaderCache.clear();
            programCache.clear();
        }
    };
}