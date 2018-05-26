/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UniformValues } from './uniform';
import { AttributeValues, createAttributeBuffers, createElementsBuffer, ElementsBuffer } from './buffer';
import { TextureValues, createTextures } from './texture';
import { Context } from './context';
import { ShaderCode, addShaderDefines, DefineValues } from '../shader-code';
import { Program } from './program';
import { RenderableSchema, RenderableValues } from '../renderable/schema';

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

export interface RenderItem {
    readonly programId: number
    readonly program: Program

    update: () => void
    draw: () => void
    destroy: () => void
}

export function createRenderItem(ctx: Context, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues): RenderItem {
    const { programCache } = ctx
    const { angleInstancedArrays, oesVertexArrayObject } = ctx.extensions

    const { attributeValues, defineValues, textureValues, uniformValues } = splitValues(schema, values)

    const glDrawMode = getDrawMode(ctx, drawMode)
    const programRef = programCache.get(ctx, {
        shaderCode: addShaderDefines(defineValues, shaderCode),
        schema
    })
    const program = programRef.value

    const textures = createTextures(ctx, schema, textureValues)
    const attributeBuffers = createAttributeBuffers(ctx, schema, attributeValues)

    const elements = values.elements

    let vertexArray: WebGLVertexArrayObjectOES
    if (oesVertexArrayObject) {
        vertexArray = oesVertexArrayObject.createVertexArrayOES()
        oesVertexArrayObject.bindVertexArrayOES(vertexArray)
        program.bindAttributes(attributeBuffers)
        ctx.vaoCount += 1
    }

    let elementsBuffer: ElementsBuffer
    if (elements && elements.ref.value) {
        elementsBuffer = createElementsBuffer(ctx, elements.ref.value)
    }

    // needs to come after elements buffer creation to include it in the vao
    if (oesVertexArrayObject) {
        oesVertexArrayObject.bindVertexArrayOES(null!)
    }

    let drawCount = values.drawCount.ref
    let instanceCount = values.instanceCount.ref

    let destroyed = false

    return {
        programId: program.id,
        program,

        draw: () => {
            program.setUniforms(uniformValues)
            if (oesVertexArrayObject) {
                oesVertexArrayObject.bindVertexArrayOES(vertexArray)
            } else {
                program.bindAttributes(attributeBuffers)
                elementsBuffer.bind()
            }
            program.bindTextures(textures)
            if (elementsBuffer) {
                angleInstancedArrays.drawElementsInstancedANGLE(glDrawMode, drawCount.value, elementsBuffer._dataType, 0, instanceCount.value);
            } else {
                angleInstancedArrays.drawArraysInstancedANGLE(glDrawMode, 0, drawCount.value, instanceCount.value)
            }
        },
        update: () => {
            if (values.drawCount.ref.version !== drawCount.version) {
                console.log('drawCount version changed')
                drawCount = values.drawCount.ref
            }
            if (values.instanceCount.ref.version !== instanceCount.version) {
                console.log('instanceCount version changed')
                instanceCount = values.instanceCount.ref
            }

            // Object.keys(attributeValues).forEach(k => {
            //     const value = attributeValues[k]
            //     if (value === undefined) return
            //     const buffer = attributeBuffers[k]
            //     if (buffer.length >= value.length) {
            //         attributeBuffers[k].updateData(value)
            //     } else {

            //     }
            // })
        },
        destroy: () => {
            if (destroyed) return
            programRef.free()
            Object.keys(textures).forEach(k => textures[k].destroy())
            Object.keys(attributeBuffers).forEach(k => attributeBuffers[k].destroy())
            if (elements) {
                elementsBuffer.destroy()
            }
            if (oesVertexArrayObject) {
                oesVertexArrayObject.deleteVertexArrayOES(vertexArray)
                ctx.vaoCount -= 1
            }
            destroyed = true
        }
    }
}