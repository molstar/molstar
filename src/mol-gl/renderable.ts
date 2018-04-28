/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import PointRenderable from './renderable/point'
import MeshRenderable from './renderable/mesh'
import { Shaders } from './shaders';
import { UniformDefs, UniformValues } from './webgl/uniform';
import { AttributeDefs, AttributeValues, createAttributeBuffers } from './webgl/buffer';
import { TextureDefs, TextureValues, createTextures } from './webgl/texture';
import { Context } from './webgl/context';
import { createProgram } from './webgl/program';

export type RenderableProps = {
    shaders: Shaders
    uniform: UniformDefs
    attribute: AttributeDefs
    texture: TextureDefs
}

export type RenderableState<T extends RenderableProps> = {
    uniform: UniformValues<T['uniform']>
    attribute: AttributeValues<T['attribute']>
    texture: TextureValues<T['texture']>

    drawCount: number
}

export interface Renderable<T extends RenderableProps> {
    readonly hash: string
    readonly programId: number

    loadAttributes: (state: Partial<AttributeValues<T['attribute']>>) => void

    draw: () => void
    dispose: () => void
}

export function createRenderable<T extends RenderableProps>(ctx: Context, props: T, state: RenderableState<T>): Renderable<T> {
    const { gl } = ctx
    const hash = JSON.stringify(props)
    const program = createProgram(ctx, props.shaders, props.uniform, props.attribute, props.texture)
    const attributeBuffers = createAttributeBuffers(ctx, props.attribute, state.attribute)
    const textures = createTextures(gl, props.texture, state.texture)

    function loadAttributes(state: Partial<AttributeValues<T['attribute']>>) {
        Object.keys(state).forEach(k => {
            const value = state[k]
            if (value !== undefined) attributeBuffers[k].updateData(value)
        })
    }

    return {
        hash,
        programId: program.id,

        loadAttributes,

        draw: () => {
            program.setUniforms(state.uniform)
            program.bindAttributes(attributeBuffers)
            program.bindTextures(textures)
            gl.drawArrays(gl.TRIANGLES, 0, state.drawCount);
        },
        dispose: () => {
            // TODO
        }
    }
}

export { PointRenderable, MeshRenderable }