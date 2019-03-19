/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from './context';
import { Program } from './program';
import { ElementsBuffer, AttributeBuffers } from './buffer';

export function createVertexArray(ctx: WebGLContext, program: Program, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer) {
    const { vertexArrayObject } = ctx.extensions
    let vertexArray: WebGLVertexArrayObject | null = null
    if (vertexArrayObject) {
        vertexArray = vertexArrayObject.createVertexArray()
        vertexArrayObject.bindVertexArray(vertexArray)
        if (elementsBuffer) elementsBuffer.bind()
        program.bindAttributes(attributeBuffers)
        ctx.vaoCount += 1
        vertexArrayObject.bindVertexArray(null)
    }
    return vertexArray
}

export function updateVertexArray(ctx: WebGLContext, vertexArray: WebGLVertexArrayObject | null, program: Program, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer) {
    const { vertexArrayObject } = ctx.extensions
    if (vertexArrayObject && vertexArray) {
        vertexArrayObject.bindVertexArray(vertexArray)
        if (elementsBuffer) elementsBuffer.bind()
        program.bindAttributes(attributeBuffers)
        vertexArrayObject.bindVertexArray(null)
    }
}

export function deleteVertexArray(ctx: WebGLContext, vertexArray: WebGLVertexArrayObject | null) {
    const { vertexArrayObject } = ctx.extensions
    if (vertexArrayObject && vertexArray) {
        vertexArrayObject.deleteVertexArray(vertexArray)
        ctx.vaoCount -= 1
    }
}