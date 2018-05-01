/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UniformDefs, UniformValues } from './uniform';
import { AttributeDefs, AttributeValues, createAttributeBuffers, createElementsBuffer, ElementsKind, ElementsBuffer } from './buffer';
import { TextureDefs, TextureValues, createTextures } from './texture';
import { Context } from './context';
import { ShaderCode } from '../shader-code';
import { Program } from './program';

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

export type RenderItemProps = {
    shaderCode: ShaderCode

    uniformDefs: UniformDefs
    attributeDefs: AttributeDefs
    textureDefs: TextureDefs

    elementsKind?: ElementsKind
    drawMode: DrawMode
}

export type RenderItemState = {
    uniformValues: UniformValues
    attributeValues: AttributeValues
    textureValues: TextureValues

    elements?: Uint32Array
    drawCount: number
    instanceCount: number
}

export interface RenderItem {
    readonly hash: string
    readonly programId: number
    readonly program: Program

    update: (state: RenderItemState) => void

    draw: () => void
    destroy: () => void
}

export function createRenderItem(ctx: Context, props: RenderItemProps, state: RenderItemState): RenderItem {
    const { programCache } = ctx
    const { angleInstancedArrays, oesVertexArrayObject } = ctx.extensions
    const { shaderCode, uniformDefs, attributeDefs, textureDefs, elementsKind } = props
    const { attributeValues, textureValues, uniformValues, elements } = state

    const hash = JSON.stringify(props)
    const drawMode = getDrawMode(ctx, props.drawMode)
    const programRef = programCache.get(ctx, { shaderCode, uniformDefs, attributeDefs, textureDefs })
    const program = programRef.value

    const textures = createTextures(ctx, textureDefs, textureValues)
    const attributeBuffers = createAttributeBuffers(ctx, attributeDefs, attributeValues)

    let vertexArray: WebGLVertexArrayObjectOES
    if (oesVertexArrayObject) {
        vertexArray = oesVertexArrayObject.createVertexArrayOES()
        oesVertexArrayObject.bindVertexArrayOES(vertexArray)
        program.bindAttributes(attributeBuffers)
        ctx.vaoCount += 1
    }

    let elementsBuffer: ElementsBuffer
    if (elements && elementsKind) {
        elementsBuffer = createElementsBuffer(ctx, elements)
    }

    // needs to come after elements buffer creation to include it in the vao
    if (oesVertexArrayObject) {
        oesVertexArrayObject.bindVertexArrayOES(null!)
    }

    let { drawCount, instanceCount } = state
    let destroyed = false

    return {
        hash,
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
                angleInstancedArrays.drawElementsInstancedANGLE(drawMode, drawCount, elementsBuffer._dataType, 0, instanceCount);
            } else {
                angleInstancedArrays.drawArraysInstancedANGLE(drawMode, 0, drawCount, instanceCount)
            }
        },
        update: (state: RenderItemState) => {
            // TODO
            const { attributeValues } = state
            Object.keys(attributeValues).forEach(k => {
                const value = attributeValues[k]
                if (value !== undefined) attributeBuffers[k].updateData(value)
            })
        },
        destroy: () => {
            if (destroyed) return
            programRef.free()
            Object.keys(textures).forEach(k => textures[k].destroy())
            Object.keys(attributeBuffers).forEach(k => attributeBuffers[k].destroy())
            if (elements && elementsKind) {
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