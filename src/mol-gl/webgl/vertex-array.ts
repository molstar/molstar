/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './program';
import { ElementsBuffer, AttributeBuffers } from './buffer';
import { WebGLExtensions } from './extensions';
import { idFactory } from '../../mol-util/id-factory';
import { GLRenderingContext } from './compat';

const getNextVertexArrayId = idFactory();

function getVertexArray(extensions: WebGLExtensions): WebGLVertexArrayObject {
    const { vertexArrayObject } = extensions;
    if (!vertexArrayObject) {
        throw new Error('VertexArrayObject not supported');
    }
    const vertexArray = vertexArrayObject.createVertexArray();
    if (!vertexArray) {
        throw new Error('Could not create WebGL vertex array');
    }
    return vertexArray;
}

function getVertexArrayObject(extensions: WebGLExtensions) {
    const { vertexArrayObject } = extensions;
    if (vertexArrayObject === null) {
        throw new Error('VertexArrayObject not supported');
    }
    return vertexArrayObject;
}

export interface VertexArray {
    readonly id: number

    bind: () => void
    update: () => void
    reset: () => void
    destroy: () => void
}

export function createVertexArray(gl: GLRenderingContext, extensions: WebGLExtensions, program: Program, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer): VertexArray {
    const id = getNextVertexArrayId();
    let vertexArray = getVertexArray(extensions);
    let vertexArrayObject = getVertexArrayObject(extensions);

    function update() {
        vertexArrayObject.bindVertexArray(vertexArray);
        if (elementsBuffer) elementsBuffer.bind();
        program.bindAttributes(attributeBuffers);
        vertexArrayObject.bindVertexArray(null);
    }

    update();
    let destroyed = false;

    return {
        id,
        bind: () => {
            vertexArrayObject.bindVertexArray(vertexArray);
        },
        update,
        reset: () => {
            vertexArray = getVertexArray(extensions);
            vertexArrayObject = getVertexArrayObject(extensions);
            update();
        },
        destroy: () => {
            if (destroyed) return;
            if (elementsBuffer) {
                // workaround for ANGLE/Chromium bug
                // - https://bugs.chromium.org/p/angleproject/issues/detail?id=6599
                // - https://bugs.chromium.org/p/chromium/issues/detail?id=1272238
                vertexArrayObject.bindVertexArray(vertexArray);
                gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
            }
            vertexArrayObject.deleteVertexArray(vertexArray);
            destroyed = true;
        }
    };
}