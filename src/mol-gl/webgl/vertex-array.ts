/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context';
import { Program } from './program';
import { AttributeBuffers, ElementsBuffer } from './buffer';

export function createVertexArray(ctx: Context, program: Program, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer) {
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

export function deleteVertexArray(ctx: Context, vertexArray?: WebGLVertexArrayObjectOES) {
    const { oesVertexArrayObject } = ctx.extensions
    if (oesVertexArrayObject && vertexArray) {
        oesVertexArrayObject.deleteVertexArrayOES(vertexArray)
        ctx.vaoCount -= 1
    }
}