/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UniformValues } from './uniform';
import { AttributeValues, createAttributeBuffers, createElementsBuffer, ElementsBuffer, createAttributeBuffer, ArrayKind, AttributeBuffers } from './buffer';
import { TextureValues, createTextures } from './texture';
import { Context } from './context';
import { ShaderCode, addShaderDefines, DefineValues } from '../shader-code';
import { Program } from './program';
import { RenderableSchema, RenderableValues, AttributeSpec } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';

const getNextRenderItemId = idFactory()

export type DrawMode = 'points' | 'lines' | 'line-strip' | 'line-loop' | 'triangles' | 'triangle-strip' | 'triangle-fan'

export function getDrawMode(ctx: Context, drawMode: DrawMode) {
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

function splitValues(schema: RenderableSchema, values: RenderableValues) {
    const attributeValues: AttributeValues = {}
    const defineValues: DefineValues = {}
    const textureValues: TextureValues = {}
    const uniformValues: UniformValues = {}
    Object.keys(values).forEach(k => {
        if (schema[k].type === 'attribute') attributeValues[k] = values[k]
        if (schema[k].type === 'define') defineValues[k] = values[k]
        if (schema[k].type === 'texture') textureValues[k] = values[k]
        if (schema[k].type === 'uniform') uniformValues[k] = values[k]
    })
    return { attributeValues, defineValues, textureValues, uniformValues }
}

type Versions<T extends RenderableValues> = { [k in keyof T]: number }
function getValueVersions<T extends RenderableValues>(values: T) {
    const versions: Versions<any> = {}
    Object.keys(values).forEach(k => {
        versions[k] = values[k].ref.version
    })
    return versions as Versions<T>
}

function createVertexArray(ctx: Context, program: Program, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer) {
    const { oesVertexArrayObject } = ctx.extensions
    let vertexArray: WebGLVertexArrayObjectOES | undefined = undefined
    if (oesVertexArrayObject) {
        vertexArray = oesVertexArrayObject.createVertexArrayOES()
        oesVertexArrayObject.bindVertexArrayOES(vertexArray)
        program.bindAttributes(attributeBuffers)
        if (elementsBuffer) elementsBuffer.bind()
        ctx.vaoCount += 1
        oesVertexArrayObject.bindVertexArrayOES(null!)
    }
    return vertexArray
}

function deleteVertexArray(ctx: Context, vertexArray?: WebGLVertexArrayObjectOES) {
    const { oesVertexArrayObject } = ctx.extensions
    if (oesVertexArrayObject && vertexArray) {
        oesVertexArrayObject.deleteVertexArrayOES(vertexArray)
        ctx.vaoCount -= 1
    }
}

export interface RenderItem {
    readonly id: number
    readonly programId: number
    readonly program: Program

    update: () => void
    draw: () => void
    destroy: () => void
}

export function createRenderItem(ctx: Context, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues): RenderItem {
    const id = getNextRenderItemId()
    const { programCache } = ctx
    const { angleInstancedArrays, oesVertexArrayObject } = ctx.extensions

    const { attributeValues, defineValues, textureValues, uniformValues } = splitValues(schema, values)
    const versions = getValueVersions(values)

    const glDrawMode = getDrawMode(ctx, drawMode)
    let programRef = programCache.get(ctx, {
        shaderCode: addShaderDefines(defineValues, shaderCode),
        schema
    })
    let program = programRef.value

    const textures = createTextures(ctx, schema, textureValues)
    const attributeBuffers = createAttributeBuffers(ctx, schema, attributeValues)

    let elementsBuffer: ElementsBuffer | undefined
    const elements = values.elements
    if (elements && elements.ref.value) {
        elementsBuffer = createElementsBuffer(ctx, elements.ref.value)
    }

    let vertexArray: WebGLVertexArrayObjectOES | undefined = createVertexArray(ctx, program, attributeBuffers, elementsBuffer)

    let drawCount = values.drawCount.ref
    let instanceCount = values.instanceCount.ref

    let destroyed = false

    return {
        id,
        get programId () { return program.id },
        get program () { return program },

        draw: () => {
            program.setUniforms(uniformValues)
            if (oesVertexArrayObject && vertexArray) {
                oesVertexArrayObject.bindVertexArrayOES(vertexArray)
            } else {
                program.bindAttributes(attributeBuffers)
                if (elementsBuffer) elementsBuffer.bind()
            }
            program.bindTextures(textures)
            if (elementsBuffer) {
                angleInstancedArrays.drawElementsInstancedANGLE(glDrawMode, drawCount.value, elementsBuffer._dataType, 0, instanceCount.value);
            } else {
                angleInstancedArrays.drawArraysInstancedANGLE(glDrawMode, 0, drawCount.value, instanceCount.value)
            }
        },
        update: () => {
            let defineChange = false
            Object.keys(defineValues).forEach(k => {
                const value = defineValues[k]
                if (value.ref.version === versions[k]) {
                    // console.log('define version unchanged', k)
                } else {
                    console.log('define version changed', k)
                    defineChange = true
                    versions[k] = value.ref.version
                }
            })

            if (defineChange) {
                console.log('some defines changed, need to rebuild program')
                programRef.free()
                // programCache.clear()
                // console.log('programCache.count', programCache.count)
                programRef = programCache.get(ctx, {
                    shaderCode: addShaderDefines(defineValues, shaderCode),
                    schema
                })
                program = programRef.value
            }

            console.log('RenderItem.update', id, values)
            if (values.drawCount.ref.version !== drawCount.version) {
                console.log('drawCount version changed')
                drawCount = values.drawCount.ref
            }
            if (values.instanceCount.ref.version !== instanceCount.version) {
                console.log('instanceCount version changed')
                instanceCount = values.instanceCount.ref
            }

            let bufferChange = false

            Object.keys(attributeValues).forEach(k => {
                const value = attributeValues[k]
                if (value.ref.version === versions[k]) {
                    // console.log('attribute version unchanged', k)
                    return
                }
                const buffer = attributeBuffers[k]
                if (buffer.length >= value.ref.value.length) {
                    console.log('attribute array large enough to update', k)
                    attributeBuffers[k].updateData(value.ref.value)
                } else {
                    console.log('attribute array to small, need to create new attribute', k)
                    attributeBuffers[k].destroy()
                    const spec = schema[k] as AttributeSpec<ArrayKind>
                    attributeBuffers[k] = createAttributeBuffer(ctx, value.ref.value, spec.itemSize, spec.divisor)
                    bufferChange = true
                }
                versions[k] = value.ref.version
            })

            if (elementsBuffer && values.elements.ref.version !== versions.elements) {
                if (elementsBuffer.length >= values.elements.ref.value.length) {
                    console.log('elements array large enough to update')
                    elementsBuffer.updateData(values.elements.ref.value)
                } else {
                    console.log('elements array to small, need to create new elements')
                    elementsBuffer.destroy()
                    elementsBuffer = createElementsBuffer(ctx, values.elements.ref.value)
                    bufferChange = true
                }
                versions.elements = values.elements.ref.version
            }

            if (defineChange || bufferChange) {
                console.log('program/defines or buffers changed, rebuild vao')
                deleteVertexArray(ctx, vertexArray)
                vertexArray = createVertexArray(ctx, program, attributeBuffers, elementsBuffer)
            }

            Object.keys(textureValues).forEach(k => {
                const value = textureValues[k]
                if (value.ref.version === versions[k]) {
                    // console.log('texture version unchanged', k)
                    return
                }
                console.log('texture version changed, uploading image', k)
                textures[k].load(value.ref.value)
            })
        },
        destroy: () => {
            if (destroyed) return
            programRef.free()
            Object.keys(textures).forEach(k => textures[k].destroy())
            Object.keys(attributeBuffers).forEach(k => attributeBuffers[k].destroy())
            if (elementsBuffer) elementsBuffer.destroy()
            deleteVertexArray(ctx, vertexArray)
            destroyed = true
        }
    }
}