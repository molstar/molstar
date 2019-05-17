/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createAttributeBuffers, createElementsBuffer, ElementsBuffer, createAttributeBuffer, AttributeKind } from './buffer';
import { createTextures, Texture } from './texture';
import { WebGLContext, checkError } from './context';
import { ShaderCode } from '../shader-code';
import { Program } from './program';
import { RenderableSchema, RenderableValues, AttributeSpec, getValueVersions, splitValues, Values } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { deleteVertexArray, createVertexArray } from './vertex-array';
import { ValueCell } from 'mol-util';
import { ReferenceItem } from 'mol-util/reference-cache';
import { TextureImage, TextureVolume } from 'mol-gl/renderable/util';
import { checkFramebufferStatus } from './framebuffer';
import { isDebugMode } from 'mol-util/debug';

const getNextRenderItemId = idFactory()

export type DrawMode = 'points' | 'lines' | 'line-strip' | 'line-loop' | 'triangles' | 'triangle-strip' | 'triangle-fan'

export function getDrawMode(ctx: WebGLContext, drawMode: DrawMode) {
    const { gl } = ctx
    switch (drawMode) {
        case 'points': return gl.POINTS
        case 'lines': return gl.LINES
        case 'line-strip': return gl.LINE_STRIP
        case 'line-loop': return gl.LINE_LOOP
        case 'triangles': return gl.TRIANGLES
        case 'triangle-strip': return gl.TRIANGLE_STRIP
        case 'triangle-fan': return gl.TRIANGLE_FAN
    }
}

export interface RenderItem<T extends string> {
    readonly id: number
    readonly materialId: number
    getProgram: (variant: T) => Program

    render: (variant: T) => void
    update: () => Readonly<ValueChanges>
    destroy: () => void
}

//

const GraphicsRenderVariantDefines = {
    'color': {},
    'pickObject': { dColorType: ValueCell.create('objectPicking') },
    'pickInstance': { dColorType: ValueCell.create('instancePicking') },
    'pickGroup': { dColorType: ValueCell.create('groupPicking') },
    'depth': { dColorType: ValueCell.create('depth') }
}
export type GraphicsRenderVariant = keyof typeof GraphicsRenderVariantDefines

const ComputeRenderVariantDefines = {
    'compute': {},
}
export type ComputeRenderVariant = keyof typeof ComputeRenderVariantDefines

type RenderVariantDefines = typeof GraphicsRenderVariantDefines | typeof ComputeRenderVariantDefines

//

type ProgramVariants = { [k: string]: ReferenceItem<Program> }
type VertexArrayVariants = { [k: string]: WebGLVertexArrayObjectOES | null }

interface ValueChanges {
    attributes: boolean
    defines: boolean
    elements: boolean
    textures: boolean
}
function createValueChanges() {
    return {
        attributes: false,
        defines: false,
        elements: false,
        textures: false,
    }
}
function resetValueChanges(valueChanges: ValueChanges) {
    valueChanges.attributes = false
    valueChanges.defines = false
    valueChanges.elements = false
    valueChanges.textures = false
}

//

export type GraphicsRenderItem = RenderItem<keyof typeof GraphicsRenderVariantDefines & string>
export function createGraphicsRenderItem(ctx: WebGLContext, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues, materialId: number) {
    return createRenderItem(ctx, drawMode, shaderCode, schema, values, materialId, GraphicsRenderVariantDefines)
}

export type ComputeRenderItem = RenderItem<keyof typeof ComputeRenderVariantDefines & string>
export function createComputeRenderItem(ctx: WebGLContext, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues, materialId = -1) {
    return createRenderItem(ctx, drawMode, shaderCode, schema, values, materialId, ComputeRenderVariantDefines)
}

/**
 * Creates a render item
 *
 * - assumes that `values.drawCount` and `values.instanceCount` exist
 */
export function createRenderItem<T extends RenderVariantDefines, S extends keyof T & string>(ctx: WebGLContext, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues, materialId: number, renderVariantDefines: T): RenderItem<S> {
    const id = getNextRenderItemId()
    const { stats, state, programCache } = ctx
    const { instancedArrays, vertexArrayObject } = ctx.extensions

    const { attributeValues, defineValues, textureValues, uniformValues, materialUniformValues } = splitValues(schema, values)

    const uniformValueEntries = Object.entries(uniformValues)
    const materialUniformValueEntries = Object.entries(materialUniformValues)
    const defineValueEntries = Object.entries(defineValues)

    const versions = getValueVersions(values)

    const glDrawMode = getDrawMode(ctx, drawMode)

    const programs: ProgramVariants = {}
    Object.keys(renderVariantDefines).forEach(k => {
        const variantDefineValues: Values<RenderableSchema> = (renderVariantDefines as any)[k]
        programs[k] = programCache.get({
            defineValues: { ...defineValues, ...variantDefineValues },
            shaderCode,
            schema
        })
    })

    const textures = createTextures(ctx, schema, textureValues)
    const attributeBuffers = createAttributeBuffers(ctx, schema, attributeValues)

    let elementsBuffer: ElementsBuffer | undefined
    const elements = values.elements
    if (elements && elements.ref.value) {
        elementsBuffer = createElementsBuffer(ctx, elements.ref.value)
    }

    const vertexArrays: VertexArrayVariants = {}
    Object.keys(renderVariantDefines).forEach(k => {
        vertexArrays[k] = createVertexArray(ctx, programs[k].value, attributeBuffers, elementsBuffer)
    })

    let drawCount = values.drawCount.ref.value
    let instanceCount = values.instanceCount.ref.value

    stats.drawCount += drawCount
    stats.instanceCount += instanceCount
    stats.instancedDrawCount += instanceCount * drawCount

    const valueChanges = createValueChanges()

    let destroyed = false
    let currentProgramId = -1

    return {
        id,
        materialId,
        getProgram: (variant: S) => programs[variant].value,

        render: (variant: S) => {
            if (drawCount === 0 || instanceCount === 0) return
            const program = programs[variant].value
            if (program.id === currentProgramId && state.currentRenderItemId === id) {
                program.setUniforms(uniformValueEntries)
                program.bindTextures(textures)
            } else {
                const vertexArray = vertexArrays[variant]
                if (program.id !== state.currentProgramId || program.id !== currentProgramId ||
                    materialId === -1 || materialId !== state.currentMaterialId
                ) {
                    // console.log('program.id changed or materialId changed/-1', materialId)
                    if (program.id !== state.currentProgramId) program.use()
                    program.setUniforms(materialUniformValueEntries)
                    state.currentMaterialId = materialId
                    currentProgramId = program.id
                }
                program.setUniforms(uniformValueEntries)
                program.bindTextures(textures)
                if (vertexArrayObject && vertexArray) {
                    vertexArrayObject.bindVertexArray(vertexArray)
                    // need to bind elements buffer explicitly since it is not always recorded in the VAO
                    if (elementsBuffer) elementsBuffer.bind()
                } else {
                    if (elementsBuffer) elementsBuffer.bind()
                    program.bindAttributes(attributeBuffers)
                }
                state.currentRenderItemId = id
            }
            if (isDebugMode) {
                checkFramebufferStatus(ctx.gl)
            }
            if (elementsBuffer) {
                instancedArrays.drawElementsInstanced(glDrawMode, drawCount, elementsBuffer._dataType, 0, instanceCount);
            } else {
                instancedArrays.drawArraysInstanced(glDrawMode, 0, drawCount, instanceCount)
            }
            if (isDebugMode) {
                try {
                    checkError(ctx.gl)
                } catch (e) {
                    // console.log('shaderCode', shaderCode)
                    // console.log('schema', schema)
                    // console.log('attributeBuffers', attributeBuffers)
                    throw new Error(`Error rendering item id ${id}: '${e}'`)
                }
            }
        },
        update: () => {
            resetValueChanges(valueChanges)

            for (let i = 0, il = defineValueEntries.length; i < il; ++i) {
                const [k, value] = defineValueEntries[i]
                if (value.ref.version !== versions[k]) {
                    // console.log('define version changed', k)
                    valueChanges.defines = true
                    versions[k] = value.ref.version
                }
            }

            if (valueChanges.defines) {
                // console.log('some defines changed, need to rebuild programs')
                Object.keys(renderVariantDefines).forEach(k => {
                    const variantDefineValues: Values<RenderableSchema> = (renderVariantDefines as any)[k]
                    programs[k].free()
                    programs[k] = programCache.get({
                        defineValues: { ...defineValues, ...variantDefineValues },
                        shaderCode,
                        schema
                    })
                })
            }

            if (values.drawCount.ref.version !== versions.drawCount) {
                // console.log('drawCount version changed')
                stats.drawCount += values.drawCount.ref.value - drawCount
                stats.instancedDrawCount += instanceCount * values.drawCount.ref.value - instanceCount * drawCount
                drawCount = values.drawCount.ref.value
                versions.drawCount = values.drawCount.ref.version
            }
            if (values.instanceCount.ref.version !== versions.instanceCount) {
                // console.log('instanceCount version changed')
                stats.instanceCount += values.instanceCount.ref.value - instanceCount
                stats.instancedDrawCount += values.instanceCount.ref.value * drawCount - instanceCount * drawCount
                instanceCount = values.instanceCount.ref.value
                versions.instanceCount = values.instanceCount.ref.version
            }

            for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
                const [k, buffer] = attributeBuffers[i]
                const value = attributeValues[k]
                if (value.ref.version !== versions[k]) {
                    if (buffer.length >= value.ref.value.length) {
                        // console.log('attribute array large enough to update', k, value.ref.id, value.ref.version)
                        buffer.updateData(value.ref.value)
                    } else {
                        // console.log('attribute array to small, need to create new attribute', k, value.ref.id, value.ref.version)
                        buffer.destroy()
                        const { itemSize, divisor } = schema[k] as AttributeSpec<AttributeKind>
                        attributeBuffers[i][1] = createAttributeBuffer(ctx, value.ref.value, itemSize, divisor)
                        valueChanges.attributes = true
                    }
                    versions[k] = value.ref.version
                }
            }

            if (elementsBuffer && values.elements.ref.version !== versions.elements) {
                if (elementsBuffer.length >= values.elements.ref.value.length) {
                    // console.log('elements array large enough to update', values.elements.ref.id, values.elements.ref.version)
                    elementsBuffer.updateData(values.elements.ref.value)
                } else {
                    // console.log('elements array to small, need to create new elements', values.elements.ref.id, values.elements.ref.version)
                    elementsBuffer.destroy()
                    elementsBuffer = createElementsBuffer(ctx, values.elements.ref.value)
                    valueChanges.elements = true
                }
                versions.elements = values.elements.ref.version
            }

            if (valueChanges.attributes || valueChanges.defines || valueChanges.elements) {
                // console.log('program/defines or buffers changed, update vaos')
                const { vertexArrayObject } = ctx.extensions
                if (vertexArrayObject) {
                    Object.keys(renderVariantDefines).forEach(k => {
                        vertexArrayObject.bindVertexArray(vertexArrays[k])
                        if (elementsBuffer && (valueChanges.defines || valueChanges.elements)) {
                            elementsBuffer.bind()
                        }
                        if (valueChanges.attributes || valueChanges.defines) {
                            programs[k].value.bindAttributes(attributeBuffers)
                        }
                        vertexArrayObject.bindVertexArray(null)
                    })
                }
            }

            for (let i = 0, il = textures.length; i < il; ++i) {
                const [k, texture] = textures[i]
                const value = textureValues[k]
                if (value.ref.version !== versions[k]) {
                    // update of textures with kind 'texture' is done externally
                    if (schema[k].kind !== 'texture') {
                        // console.log('texture version changed, uploading image', k)
                        texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>)
                        versions[k] = value.ref.version
                        valueChanges.textures = true
                    } else {
                        textures[i][1] = value.ref.value as Texture
                    }
                }
            }

            return valueChanges
        },
        destroy: () => {
            if (!destroyed) {
                Object.keys(renderVariantDefines).forEach(k => {
                    programs[k].free()
                    deleteVertexArray(ctx, vertexArrays[k])
                })
                textures.forEach(([k, texture]) => {
                    // lifetime of textures with kind 'texture' is defined externally
                    if (schema[k].kind !== 'texture') {
                        texture.destroy()
                    }
                })
                attributeBuffers.forEach(([_, buffer]) => buffer.destroy())
                if (elementsBuffer) elementsBuffer.destroy()
                destroyed = true
            }
        }
    }
}