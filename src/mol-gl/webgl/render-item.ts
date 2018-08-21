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
import { RenderableSchema, RenderableValues, AttributeSpec, getValueVersions, splitValues, Values } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { deleteVertexArray, createVertexArray } from './vertex-array';
import { ValueCell } from 'mol-util';
import { ReferenceItem } from 'mol-util/reference-cache';

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
type VertexArrayVariants = { [k: string]: WebGLVertexArrayObjectOES | undefined }

interface ValueChanges {
    attributes: boolean
    defines: boolean
    elements: boolean
    textures: boolean
    uniforms: boolean
}

export function createRenderItem(ctx: Context, drawMode: DrawMode, shaderCode: ShaderCode, schema: RenderableSchema, values: RenderableValues): RenderItem {
    const id = getNextRenderItemId()
    const { programCache } = ctx
    const { angleInstancedArrays, oesVertexArrayObject } = ctx.extensions

    const { attributeValues, defineValues, textureValues, uniformValues } = splitValues(schema, values)
    const versions = getValueVersions(values)

    const glDrawMode = getDrawMode(ctx, drawMode)

    const programs: ProgramVariants = {}
    Object.keys(RenderVariantDefines).forEach(k => {
        const variantDefineValues: Values<RenderableSchema> = (RenderVariantDefines as any)[k]
        programs[k] = programCache.get(ctx, {
            shaderCode: addShaderDefines({ ...defineValues, ...variantDefineValues }, shaderCode),
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

    const valueChanges: ValueChanges = {
        attributes: false,
        defines: false,
        elements: false,
        textures: false,
        uniforms: false
    }

    let destroyed = false

    return {
        id,
        getProgram: (variant: RenderVariant) => programs[variant].value,

        render: (variant: RenderVariant) => {
            if (drawCount === 0 || instanceCount === 0) return
            const program = programs[variant].value
            const vertexArray = vertexArrays[variant]
            program.setUniforms(uniformValues)
            if (oesVertexArrayObject && vertexArray) {
                oesVertexArrayObject.bindVertexArrayOES(vertexArray)
                // TODO need to bind elements buffer explicitely since it is not always recorded in the VAO
                if (elementsBuffer) elementsBuffer.bind()
            } else {
                if (elementsBuffer) elementsBuffer.bind()
                program.bindAttributes(attributeBuffers)
            }
            program.bindTextures(textures)
            if (elementsBuffer) {
                angleInstancedArrays.drawElementsInstancedANGLE(glDrawMode, drawCount, elementsBuffer._dataType, 0, instanceCount);
            } else {
                angleInstancedArrays.drawArraysInstancedANGLE(glDrawMode, 0, drawCount, instanceCount)
            }
        },
        update: () => {
            valueChanges.defines = false
            Object.keys(defineValues).forEach(k => {
                const value = defineValues[k]
                if (value.ref.version !== versions[k]) {
                    // console.log('define version changed', k)
                    valueChanges.defines = true
                    versions[k] = value.ref.version
                }
            })

            if (valueChanges.defines) {
                // console.log('some defines changed, need to rebuild programs')
                Object.keys(RenderVariantDefines).forEach(k => {
                    const variantDefineValues: Values<RenderableSchema> = (RenderVariantDefines as any)[k]
                    programs[k].free()
                    programs[k] = programCache.get(ctx, {
                        shaderCode: addShaderDefines({ ...defineValues, ...variantDefineValues }, shaderCode),
                        schema
                    })
                })
            }

            if (values.drawCount.ref.version !== versions.drawCount) {
                // console.log('drawCount version changed')
                drawCount = values.drawCount.ref.value
                versions.drawCount = values.drawCount.ref.version
            }
            if (values.instanceCount.ref.version !== versions.instanceCount) {
                // console.log('instanceCount version changed')
                instanceCount = values.instanceCount.ref.value
                versions.instanceCount = values.instanceCount.ref.version
            }

            valueChanges.attributes = false
            Object.keys(attributeValues).forEach(k => {
                const value = attributeValues[k]
                if (value.ref.version !== versions[k]) {
                    const buffer = attributeBuffers[k]
                    if (buffer.length >= value.ref.value.length) {
                        // console.log('attribute array large enough to update', k)
                        attributeBuffers[k].updateData(value.ref.value)
                    } else {
                        // console.log('attribute array to small, need to create new attribute', k)
                        attributeBuffers[k].destroy()
                        const spec = schema[k] as AttributeSpec<ArrayKind>
                        attributeBuffers[k] = createAttributeBuffer(ctx, value.ref.value, spec.itemSize, spec.divisor)
                        valueChanges.attributes = true
                    }
                    versions[k] = value.ref.version
                }
            })

            valueChanges.elements = false
            if (elementsBuffer && values.elements.ref.version !== versions.elements) {
                if (elementsBuffer.length >= values.elements.ref.value.length) {
                    // console.log('elements array large enough to update')
                    elementsBuffer.updateData(values.elements.ref.value)
                } else {
                    // console.log('elements array to small, need to create new elements')
                    elementsBuffer.destroy()
                    elementsBuffer = createElementsBuffer(ctx, values.elements.ref.value)
                    valueChanges.elements = true
                }
                versions.elements = values.elements.ref.version
            }

            if (valueChanges.attributes || valueChanges.defines || valueChanges.elements) {
                // console.log('program/defines or buffers changed, rebuild vaos')
                Object.keys(RenderVariantDefines).forEach(k => {
                    deleteVertexArray(ctx, vertexArrays[k])
                    vertexArrays[k] = createVertexArray(ctx, programs[k].value, attributeBuffers, elementsBuffer)
                })
            }

            valueChanges.textures = false
            Object.keys(textureValues).forEach(k => {
                const value = textureValues[k]
                if (value.ref.version !== versions[k]) {
                    // console.log('texture version changed, uploading image', k)
                    textures[k].load(value.ref.value)
                    versions[k] = value.ref.version
                    valueChanges.textures = true
                }
            })

            return valueChanges
        },
        destroy: () => {
            if (!destroyed) {
                Object.keys(RenderVariantDefines).forEach(k => {
                    programs[k].free()
                    deleteVertexArray(ctx, vertexArrays[k])
                })
                Object.keys(textures).forEach(k => textures[k].destroy())
                Object.keys(attributeBuffers).forEach(k => attributeBuffers[k].destroy())
                if (elementsBuffer) elementsBuffer.destroy()
                destroyed = true
            }
        }
    }
}