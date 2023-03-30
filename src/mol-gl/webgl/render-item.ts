/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { createAttributeBuffers, ElementsBuffer, AttributeKind } from './buffer';
import { createTextures, Texture } from './texture';
import { WebGLContext, checkError } from './context';
import { ShaderCode, DefineValues } from '../shader-code';
import { Program, Programs } from './program';
import { RenderableSchema, RenderableValues, AttributeSpec, getValueVersions, splitValues, DefineSpec } from '../renderable/schema';
import { idFactory } from '../../mol-util/id-factory';
import { ValueCell } from '../../mol-util';
import { TextureImage, TextureVolume } from '../../mol-gl/renderable/util';
import { checkFramebufferStatus } from './framebuffer';
import { isDebugMode, isTimingMode } from '../../mol-util/debug';
import { VertexArray } from './vertex-array';
import { fillSerial } from '../../mol-util/array';
import { deepClone } from '../../mol-util/object';
import { cloneUniformValues } from './uniform';

const getNextRenderItemId = idFactory();

export type DrawMode = 'points' | 'lines' | 'line-strip' | 'line-loop' | 'triangles' | 'triangle-strip' | 'triangle-fan'

export function getDrawMode(ctx: WebGLContext, drawMode: DrawMode) {
    const { gl } = ctx;
    switch (drawMode) {
        case 'points': return gl.POINTS;
        case 'lines': return gl.LINES;
        case 'line-strip': return gl.LINE_STRIP;
        case 'line-loop': return gl.LINE_LOOP;
        case 'triangles': return gl.TRIANGLES;
        case 'triangle-strip': return gl.TRIANGLE_STRIP;
        case 'triangle-fan': return gl.TRIANGLE_FAN;
    }
}

export interface RenderItem<T extends string> {
    readonly id: number
    readonly materialId: number
    getProgram: (variant: T) => Program

    render: (variant: T, sharedTexturesCount: number) => void
    update: () => void
    destroy: () => void
}

//

const GraphicsRenderVariant = { colorBlended: '', colorWboit: '', colorDpoit: '', pick: '', depth: '', marking: '' };
export type GraphicsRenderVariant = keyof typeof GraphicsRenderVariant
export const GraphicsRenderVariants = Object.keys(GraphicsRenderVariant) as GraphicsRenderVariant[];
export const GraphicsRenderVariantsBlended = GraphicsRenderVariants.filter(v => !['colorWboit', 'colorDpoit'].includes(v));
export const GraphicsRenderVariantsWboit = GraphicsRenderVariants.filter(v => !['colorBlended', 'colorDpoit'].includes(v));
export const GraphicsRenderVariantsDpoit = GraphicsRenderVariants.filter(v => !['colorWboit', 'colorBlended'].includes(v));

const ComputeRenderVariant = { compute: '' };
export type ComputeRenderVariant = keyof typeof ComputeRenderVariant
export const ComputeRenderVariants = Object.keys(ComputeRenderVariant) as ComputeRenderVariant[];

function createProgramVariant(ctx: WebGLContext, variant: string, defineValues: DefineValues, shaderCode: ShaderCode, schema: RenderableSchema) {
    defineValues = { ...defineValues, dRenderVariant: ValueCell.create(variant) };
    if (schema.dRenderVariant === undefined) {
        Object.defineProperty(schema, 'dRenderVariant', { value: DefineSpec('string') });
    }
    return ctx.resources.program(defineValues, shaderCode, schema);
}

//

function createValueChanges() {
    return {
        attributes: false,
        defines: false,
        elements: false,
        textures: false,
    };
}
type ValueChanges = ReturnType<typeof createValueChanges>

function resetValueChanges(valueChanges: ValueChanges) {
    valueChanges.attributes = false;
    valueChanges.defines = false;
    valueChanges.elements = false;
    valueChanges.textures = false;
}

//

export type GraphicsRenderItem = RenderItem<GraphicsRenderVariant>
export function createGraphicsRenderItem(ctx: WebGLContext, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues, materialId: number, variants: GraphicsRenderVariant[]) {
    return createRenderItem(ctx, drawMode, shaderCode, schema, values, materialId, variants);
}

export type ComputeRenderItem = RenderItem<ComputeRenderVariant>
export function createComputeRenderItem(ctx: WebGLContext, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues, materialId = -1) {
    return createRenderItem(ctx, drawMode, shaderCode, schema, values, materialId, ComputeRenderVariants);
}

/**
 * Creates a render item
 *
 * - assumes that `values.drawCount` and `values.instanceCount` exist
 */
export function createRenderItem<T extends string>(ctx: WebGLContext, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues, materialId: number, renderVariants: T[]): RenderItem<T> {
    const id = getNextRenderItemId();
    const { stats, state, resources } = ctx;
    const { instancedArrays, vertexArrayObject } = ctx.extensions;

    // emulate gl_VertexID when needed
    if (values.uVertexCount && !ctx.extensions.noNonInstancedActiveAttribs) {
        const vertexCount = values.uVertexCount.ref.value;
        (values as any).aVertex = ValueCell.create(fillSerial(new Float32Array(vertexCount)));
        (schema as any).aVertex = AttributeSpec('float32', 1, 0);
    }

    const { attributeValues, defineValues, textureValues, materialTextureValues, uniformValues, materialUniformValues, bufferedUniformValues } = splitValues(schema, values);

    const uniformValueEntries = Object.entries(uniformValues);
    const materialUniformValueEntries = Object.entries(materialUniformValues);
    const backBufferUniformValueEntries = Object.entries(bufferedUniformValues);
    const frontBufferUniformValueEntries = Object.entries(cloneUniformValues(bufferedUniformValues));
    const defineValueEntries = Object.entries(defineValues);

    const versions = getValueVersions(values);

    const glDrawMode = getDrawMode(ctx, drawMode);

    const programs: Programs = {};
    for (const k of renderVariants) {
        programs[k] = createProgramVariant(ctx, k, defineValues, shaderCode, schema);
    }

    const textures = createTextures(ctx, schema, textureValues);
    const materialTextures = createTextures(ctx, schema, materialTextureValues);
    const attributeBuffers = createAttributeBuffers(ctx, schema, attributeValues);

    let elementsBuffer: ElementsBuffer | undefined;
    const elements = values.elements;
    if (elements && elements.ref.value) {
        elementsBuffer = resources.elements(elements.ref.value);
    }

    let vertexArray: VertexArray | null = vertexArrayObject ? resources.vertexArray(programs, attributeBuffers, elementsBuffer) : null;

    let drawCount: number = values.drawCount.ref.value;
    let instanceCount: number = values.instanceCount.ref.value;

    stats.drawCount += drawCount;
    stats.instanceCount += instanceCount;
    stats.instancedDrawCount += instanceCount * drawCount;

    const valueChanges = createValueChanges();

    let destroyed = false;
    let currentProgramId = -1;

    return {
        id,
        materialId,
        getProgram: (variant: T) => programs[variant],

        render: (variant: T, sharedTexturesCount: number) => {
            if (drawCount === 0 || instanceCount === 0) return;
            const program = programs[variant];
            if (program.id === currentProgramId && state.currentRenderItemId === id) {
                program.setUniforms(uniformValueEntries);
                program.bindTextures(textures, sharedTexturesCount);
            } else {
                if (program.id !== state.currentProgramId || program.id !== currentProgramId ||
                    materialId === -1 || materialId !== state.currentMaterialId
                ) {
                    // console.log('program.id changed or materialId changed/-1', materialId)
                    if (program.id !== state.currentProgramId) program.use();
                    program.setUniforms(materialUniformValueEntries);
                    program.bindTextures(materialTextures, sharedTexturesCount + textures.length);
                    state.currentMaterialId = materialId;
                    currentProgramId = program.id;
                }
                program.setUniforms(uniformValueEntries);
                program.setUniforms(frontBufferUniformValueEntries);
                program.bindTextures(textures, sharedTexturesCount);
                if (vertexArray) {
                    vertexArray.bind();
                    // need to bind elements buffer explicitly since it is not always recorded in the VAO
                    if (elementsBuffer) elementsBuffer.bind();
                } else {
                    if (elementsBuffer) elementsBuffer.bind();
                    program.bindAttributes(attributeBuffers);
                }
                state.currentRenderItemId = id;
            }
            if (isDebugMode) {
                try {
                    checkFramebufferStatus(ctx.gl);
                } catch (e) {
                    throw new Error(`Framebuffer error rendering item id ${id}: '${e}'`);
                }
            }
            if (elementsBuffer) {
                instancedArrays.drawElementsInstanced(glDrawMode, drawCount, elementsBuffer._dataType, 0, instanceCount);
            } else {
                instancedArrays.drawArraysInstanced(glDrawMode, 0, drawCount, instanceCount);
            }
            if (isTimingMode) {
                stats.calls.drawInstanced += 1;
                stats.calls.counts += instanceCount;
            }
            if (isDebugMode) {
                try {
                    checkError(ctx.gl);
                } catch (e) {
                    throw new Error(`Draw error rendering item id ${id}: '${e}'`);
                }
            }
        },
        update: () => {
            resetValueChanges(valueChanges);

            if (values.aVertex) {
                const vertexCount = values.uVertexCount.ref.value;
                if (values.aVertex.ref.value.length < vertexCount) {
                    ValueCell.update(values.aVertex, fillSerial(new Float32Array(vertexCount)));
                }
            }

            for (let i = 0, il = defineValueEntries.length; i < il; ++i) {
                const [k, value] = defineValueEntries[i];
                if (value.ref.version !== versions[k]) {
                    // console.log('define version changed', k);
                    valueChanges.defines = true;
                    versions[k] = value.ref.version;
                }
            }

            if (valueChanges.defines) {
                // console.log('some defines changed, need to rebuild programs');
                for (const k of renderVariants) {
                    programs[k].destroy();
                    programs[k] = createProgramVariant(ctx, k, defineValues, shaderCode, schema);
                }
            }

            if (values.drawCount.ref.version !== versions.drawCount) {
                // console.log('drawCount version changed');
                stats.drawCount += values.drawCount.ref.value - drawCount;
                stats.instancedDrawCount += instanceCount * values.drawCount.ref.value - instanceCount * drawCount;
                drawCount = values.drawCount.ref.value;
                versions.drawCount = values.drawCount.ref.version;
            }
            if (values.instanceCount.ref.version !== versions.instanceCount) {
                // console.log('instanceCount version changed');
                stats.instanceCount += values.instanceCount.ref.value - instanceCount;
                stats.instancedDrawCount += values.instanceCount.ref.value * drawCount - instanceCount * drawCount;
                instanceCount = values.instanceCount.ref.value;
                versions.instanceCount = values.instanceCount.ref.version;
            }

            for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
                const [k, buffer] = attributeBuffers[i];
                const value = attributeValues[k];
                if (value.ref.version !== versions[k]) {
                    if (buffer.length >= value.ref.value.length) {
                        // console.log('attribute array large enough to update', buffer.id, k, value.ref.id, value.ref.version);
                        buffer.updateSubData(value.ref.value, 0, buffer.length);
                    } else {
                        // console.log('attribute array too small, need to create new attribute', buffer.id, k, value.ref.id, value.ref.version);
                        buffer.destroy();
                        const { itemSize, divisor } = schema[k] as AttributeSpec<AttributeKind>;
                        attributeBuffers[i][1] = resources.attribute(value.ref.value, itemSize, divisor);
                        valueChanges.attributes = true;
                    }
                    versions[k] = value.ref.version;
                }
            }

            if (elementsBuffer && values.elements.ref.version !== versions.elements) {
                if (elementsBuffer.length >= values.elements.ref.value.length) {
                    // console.log('elements array large enough to update', values.elements.ref.id, values.elements.ref.version);
                    elementsBuffer.updateSubData(values.elements.ref.value, 0, elementsBuffer.length);
                } else {
                    // console.log('elements array to small, need to create new elements', values.elements.ref.id, values.elements.ref.version);
                    elementsBuffer.destroy();
                    elementsBuffer = resources.elements(values.elements.ref.value);
                    valueChanges.elements = true;
                }
                versions.elements = values.elements.ref.version;
            }

            if (valueChanges.attributes || valueChanges.defines || valueChanges.elements) {
                // console.log('program/defines or buffers changed, update vao');
                if (vertexArray) vertexArray.destroy();
                vertexArray = vertexArrayObject ? resources.vertexArray(programs, attributeBuffers, elementsBuffer) : null;
            }

            for (let i = 0, il = textures.length; i < il; ++i) {
                const [k, texture] = textures[i];
                const value = textureValues[k];
                if (value.ref.version !== versions[k]) {
                    // update of textures with kind 'texture' is done externally
                    if (schema[k].kind !== 'texture') {
                        // console.log('texture version changed, uploading image', k);
                        texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
                        valueChanges.textures = true;
                    } else {
                        textures[i][1] = value.ref.value as Texture;
                    }
                    versions[k] = value.ref.version;
                }
            }

            for (let i = 0, il = materialTextures.length; i < il; ++i) {
                const [k, texture] = materialTextures[i];
                const value = materialTextureValues[k];
                if (value.ref.version !== versions[k]) {
                    // update of textures with kind 'texture' is done externally
                    if (schema[k].kind !== 'texture') {
                        // console.log('materialTexture version changed, uploading image', k);
                        texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
                        valueChanges.textures = true;
                    } else {
                        materialTextures[i][1] = value.ref.value as Texture;
                    }
                    versions[k] = value.ref.version;
                }
            }

            for (let i = 0, il = backBufferUniformValueEntries.length; i < il; ++i) {
                const [k, uniform] = backBufferUniformValueEntries[i];
                if (uniform.ref.version !== versions[k]) {
                    // console.log('back-buffer uniform version changed, updating front-buffer', k);
                    ValueCell.update(frontBufferUniformValueEntries[i][1], deepClone(uniform.ref.value));
                    versions[k] = uniform.ref.version;
                }
            }
        },
        destroy: () => {
            if (!destroyed) {
                for (const k of renderVariants) {
                    programs[k].destroy();
                }
                if (vertexArray) vertexArray.destroy();
                textures.forEach(([k, texture]) => {
                    // lifetime of textures with kind 'texture' is defined externally
                    if (schema[k].kind !== 'texture') texture.destroy();
                });
                materialTextures.forEach(([k, texture]) => {
                    // lifetime of textures with kind 'texture' is defined externally
                    if (schema[k].kind !== 'texture') texture.destroy();
                });
                attributeBuffers.forEach(([_, buffer]) => buffer.destroy());
                if (elementsBuffer) elementsBuffer.destroy();

                stats.drawCount -= drawCount;
                stats.instanceCount -= instanceCount;
                stats.instancedDrawCount -= instanceCount * drawCount;

                destroyed = true;
            }
        }
    };
}
