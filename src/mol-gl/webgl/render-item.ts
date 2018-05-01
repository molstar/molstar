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
    dispose: () => void
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
    }

    let elementsBuffer: ElementsBuffer
    if (elements && elementsKind) {
        elementsBuffer = createElementsBuffer(ctx, elements)
    }

    let { drawCount, instanceCount } = state

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
        dispose: () => {
            // TODO
            programRef.free()
        }
    }
}