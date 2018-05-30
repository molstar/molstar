/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createAttributeBuffers, createElementsBuffer, ElementsBuffer, createAttributeBuffer, ArrayKind } from './buffer';
import { createTextures } from './texture';
import { Context } from './context';
import { ShaderCode, addShaderDefines } from '../shader-code';
import { Program } from './program';
import { RenderableSchema, RenderableValues, AttributeSpec, getValueVersions, splitValues } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { deleteVertexArray, createVertexArray } from './vertex-array';
import { ValueCell } from 'mol-util';

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

export interface RenderItem {
    readonly id: number
    readonly drawProgram: Program
    readonly pickProgram: Program

    update: () => void
    draw: () => void
    pick: () => void
    destroy: () => void
}

export function createRenderItem(ctx: Context, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues): RenderItem {
    const id = getNextRenderItemId()
    const { programCache } = ctx
    const { angleInstancedArrays, oesVertexArrayObject } = ctx.extensions

    const { attributeValues, defineValues, textureValues, uniformValues } = splitValues(schema, values)
    const versions = getValueVersions(values)

    const glDrawMode = getDrawMode(ctx, drawMode)
    let drawProgram = programCache.get(ctx, {
        shaderCode: addShaderDefines(defineValues, shaderCode),
        schema
    })
    let pickProgram = programCache.get(ctx, {
        shaderCode: addShaderDefines({ ...defineValues, dColorType: ValueCell.create('elementPicking') }, shaderCode),
        schema
    })

    const textures = createTextures(ctx, schema, textureValues)
    const attributeBuffers = createAttributeBuffers(ctx, schema, attributeValues)

    let elementsBuffer: ElementsBuffer | undefined
    const elements = values.elements
    if (elements && elements.ref.value) {
        elementsBuffer = createElementsBuffer(ctx, elements.ref.value)
    }

    let drawVertexArray: WebGLVertexArrayObjectOES | undefined = createVertexArray(ctx, drawProgram.value, attributeBuffers, elementsBuffer)
    let pickVertexArray: WebGLVertexArrayObjectOES | undefined = createVertexArray(ctx, pickProgram.value, attributeBuffers, elementsBuffer)

    let drawCount = values.drawCount.ref.value
    let instanceCount = values.instanceCount.ref.value

    let destroyed = false

    function render(program: Program, vertexArray: WebGLVertexArrayObjectOES | undefined) {
        program.setUniforms(uniformValues)
        if (oesVertexArrayObject && vertexArray) {
            oesVertexArrayObject.bindVertexArrayOES(vertexArray)
        } else {
            program.bindAttributes(attributeBuffers)
            if (elementsBuffer) elementsBuffer.bind()
        }
        program.bindTextures(textures)
        if (elementsBuffer) {
            angleInstancedArrays.drawElementsInstancedANGLE(glDrawMode, drawCount, elementsBuffer._dataType, 0, instanceCount);
        } else {
            angleInstancedArrays.drawArraysInstancedANGLE(glDrawMode, 0, drawCount, instanceCount)
        }
    }

    return {
        id,
        get drawProgram () { return drawProgram.value },
        get pickProgram () { return pickProgram.value },

        draw: () => {
            render(drawProgram.value, drawVertexArray)
        },
        pick: () => {
            render(pickProgram.value, pickVertexArray)
        },
        update: () => {
            let defineChange = false
            Object.keys(defineValues).forEach(k => {
                const value = defineValues[k]
                if (value.ref.version !== versions[k]) {
                    console.log('define version changed', k)
                    defineChange = true
                    versions[k] = value.ref.version
                }
            })

            if (defineChange) {
                console.log('some defines changed, need to rebuild program')
                drawProgram.free()
                drawProgram = programCache.get(ctx, {
                    shaderCode: addShaderDefines(defineValues, shaderCode),
                    schema
                })

                pickProgram.free()
                pickProgram = programCache.get(ctx, {
                    shaderCode: addShaderDefines({ ...defineValues, dColorType: ValueCell.create('elementPicking') }, shaderCode),
                    schema
                })
            }

            if (values.drawCount.ref.version !== versions.drawCount) {
                console.log('drawCount version changed')
                drawCount = values.drawCount.ref.value
                versions.drawCount = values.drawCount.ref.version
            }
            if (values.instanceCount.ref.version !== versions.instanceCount) {
                console.log('instanceCount version changed')
                instanceCount = values.instanceCount.ref.value
                versions.instanceCount = values.instanceCount.ref.version
            }

            let bufferChange = false

            Object.keys(attributeValues).forEach(k => {
                const value = attributeValues[k]
                if (value.ref.version !== versions[k]) {
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
                }
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
                deleteVertexArray(ctx, drawVertexArray)
                drawVertexArray = createVertexArray(ctx, drawProgram.value, attributeBuffers, elementsBuffer)
                deleteVertexArray(ctx, pickVertexArray)
                pickVertexArray = createVertexArray(ctx, drawProgram.value, attributeBuffers, elementsBuffer)
            }

            Object.keys(textureValues).forEach(k => {
                const value = textureValues[k]
                if (value.ref.version !== versions[k]) {
                    console.log('texture version changed, uploading image', k)
                    textures[k].load(value.ref.value)
                    versions[k] = value.ref.version
                }
            })
        },
        destroy: () => {
            if (!destroyed) {
                drawProgram.free()
                pickProgram.free()
                Object.keys(textures).forEach(k => textures[k].destroy())
                Object.keys(attributeBuffers).forEach(k => attributeBuffers[k].destroy())
                if (elementsBuffer) elementsBuffer.destroy()
                deleteVertexArray(ctx, drawVertexArray)
                deleteVertexArray(ctx, pickVertexArray)
                destroyed = true
            }
        }
    }
}