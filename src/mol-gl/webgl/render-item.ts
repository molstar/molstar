/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createAttributeBuffers, createElementsBuffer, ElementsBuffer, createAttributeBuffer, ArrayKind } from './buffer';
import { createTextures } from './texture';
import { WebGLContext } from './context';
import { ShaderCode } from '../shader-code';
import { Program } from './program';
import { RenderableSchema, RenderableValues, AttributeSpec, getValueVersions, splitValues, Values, splitKeys } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { deleteVertexArray, createVertexArray } from './vertex-array';
import { ValueCell } from 'mol-util';
import { ReferenceItem } from 'mol-util/reference-cache';
import { TextureImage, TextureVolume } from 'mol-gl/renderable/util';

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

export interface RenderItem {
    readonly id: number
    getProgram: (variant: RenderVariant) => Program

    render: (variant: RenderVariant) => void
    update: () => Readonly<ValueChanges>
    destroy: () => void
}

const RenderVariantDefines = {
    'draw': {},
    'pickObject': { dColorType: ValueCell.create('objectPicking') },
    'pickInstance': { dColorType: ValueCell.create('instancePicking') },
    'pickGroup': { dColorType: ValueCell.create('groupPicking') }
}
export type RenderVariant = keyof typeof RenderVariantDefines

type ProgramVariants = { [k: string]: ReferenceItem<Program> }
type VertexArrayVariants = { [k: string]: WebGLVertexArrayObjectOES | null }

interface ValueChanges {
    attributes: boolean
    defines: boolean
    elements: boolean
    textures: boolean
    uniforms: boolean
}
function createValueChanges() {
    return {
        attributes: false,
        defines: false,
        elements: false,
        textures: false,
        uniforms: false,
    }
}
function resetValueChanges(valueChanges: ValueChanges) {
    valueChanges.attributes = false
    valueChanges.defines = false
    valueChanges.elements = false
    valueChanges.textures = false
    valueChanges.uniforms = false
}

// TODO make `RenderVariantDefines` a parameter for `createRenderItem`

/**
 * Creates a render item
 *
 * - assumes that `values.drawCount` and `values.instanceCount` exist
 */
export function createRenderItem(ctx: WebGLContext, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues): RenderItem {
    const id = getNextRenderItemId()
    const { programCache } = ctx
    const { instancedArrays, vertexArrayObject } = ctx.extensions

    const { attributeValues, defineValues, textureValues, uniformValues } = splitValues(schema, values)
    const { attributeKeys, defineKeys, textureKeys } = splitKeys(schema)
    const versions = getValueVersions(values)

    const glDrawMode = getDrawMode(ctx, drawMode)

    const programs: ProgramVariants = {}
    Object.keys(RenderVariantDefines).forEach(k => {
        const variantDefineValues: Values<RenderableSchema> = (RenderVariantDefines as any)[k]
        programs[k] = programCache.get(ctx, {
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
    Object.keys(RenderVariantDefines).forEach(k => {
        vertexArrays[k] = createVertexArray(ctx, programs[k].value, attributeBuffers, elementsBuffer)
    })

    let drawCount = values.drawCount.ref.value
    let instanceCount = values.instanceCount.ref.value

    ctx.drawCount += drawCount
    ctx.instanceCount += instanceCount
    ctx.instancedDrawCount += instanceCount * drawCount

    const valueChanges = createValueChanges()

    let destroyed = false

    return {
        id,
        getProgram: (variant: RenderVariant) => programs[variant].value,

        render: (variant: RenderVariant) => {
            if (drawCount === 0 || instanceCount === 0) return
            const program = programs[variant].value
            const vertexArray = vertexArrays[variant]
            program.setUniforms(uniformValues)
            if (vertexArrayObject && vertexArray) {
                vertexArrayObject.bindVertexArray(vertexArray)
                // need to bind elements buffer explicitly since it is not always recorded in the VAO
                if (elementsBuffer) elementsBuffer.bind()
            } else {
                if (elementsBuffer) elementsBuffer.bind()
                program.bindAttributes(attributeBuffers)
            }
            program.bindTextures(textures)
            if (elementsBuffer) {
                instancedArrays.drawElementsInstanced(glDrawMode, drawCount, elementsBuffer._dataType, 0, instanceCount);
            } else {
                instancedArrays.drawArraysInstanced(glDrawMode, 0, drawCount, instanceCount)
            }
        },
        update: () => {
            resetValueChanges(valueChanges)

            for (let i = 0, il = defineKeys.length; i < il; ++i) {
                const k = defineKeys[i]
                const value = defineValues[k]
                if (value.ref.version !== versions[k]) {
                    // console.log('define version changed', k)
                    valueChanges.defines = true
                    versions[k] = value.ref.version
                }
            }

            if (valueChanges.defines) {
                // console.log('some defines changed, need to rebuild programs')
                Object.keys(RenderVariantDefines).forEach(k => {
                    const variantDefineValues: Values<RenderableSchema> = (RenderVariantDefines as any)[k]
                    programs[k].free()
                    programs[k] = programCache.get(ctx, {
                        defineValues: { ...defineValues, ...variantDefineValues },
                        shaderCode,
                        schema
                    })
                })
            }

            if (values.drawCount.ref.version !== versions.drawCount) {
                // console.log('drawCount version changed')
                ctx.drawCount += values.drawCount.ref.value - drawCount
                ctx.instancedDrawCount += instanceCount * values.drawCount.ref.value - instanceCount * drawCount
                drawCount = values.drawCount.ref.value
                versions.drawCount = values.drawCount.ref.version
            }
            if (values.instanceCount.ref.version !== versions.instanceCount) {
                // console.log('instanceCount version changed')
                ctx.instanceCount += values.instanceCount.ref.value - instanceCount
                ctx.instancedDrawCount += values.instanceCount.ref.value * drawCount - instanceCount * drawCount
                instanceCount = values.instanceCount.ref.value
                versions.instanceCount = values.instanceCount.ref.version
            }

            for (let i = 0, il = attributeKeys.length; i < il; ++i) {
                const k = attributeKeys[i]
                const value = attributeValues[k]
                if (value.ref.version !== versions[k]) {
                    const buffer = attributeBuffers[k]
                    if (buffer.length >= value.ref.value.length) {
                        // console.log('attribute array large enough to update', k, value.ref.id, value.ref.version)
                        buffer.updateData(value.ref.value)
                    } else {
                        // console.log('attribute array to small, need to create new attribute', k, value.ref.id, value.ref.version)
                        buffer.destroy()
                        const { itemSize, divisor } = schema[k] as AttributeSpec<ArrayKind>
                        attributeBuffers[k] = createAttributeBuffer(ctx, value.ref.value, itemSize, divisor)
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
                    Object.keys(RenderVariantDefines).forEach(k => {
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

            for (let i = 0, il = textureKeys.length; i < il; ++i) {
                const k = textureKeys[i]
                const value = textureValues[k]
                if (value.ref.version !== versions[k]) {
                    // update of textures with kind 'texture' is done externally
                    if (schema[k].kind !== 'texture') {
                        // console.log('texture version changed, uploading image', k)
                        textures[k].load(value.ref.value as TextureImage<any> | TextureVolume<any>)
                        versions[k] = value.ref.version
                        valueChanges.textures = true
                    }
                }
            }

            return valueChanges
        },
        destroy: () => {
            if (!destroyed) {
                Object.keys(RenderVariantDefines).forEach(k => {
                    programs[k].free()
                    deleteVertexArray(ctx, vertexArrays[k])
                })
                Object.keys(textures).forEach(k => {
                    // lifetime of textures with kind 'texture' is defined externally
                    if (schema[k].kind !== 'texture') {
                        textures[k].destroy()
                    }
                })
                Object.keys(attributeBuffers).forEach(k => attributeBuffers[k].destroy())
                if (elementsBuffer) elementsBuffer.destroy()
                destroyed = true
            }
        }
    }
}