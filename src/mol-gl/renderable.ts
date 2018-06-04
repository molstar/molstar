/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema } from './renderable/schema';
import { RenderVariant, RenderItem } from './webgl/render-item';

export type RenderableState = {
    visible: boolean
    depthMask: boolean
}

export interface Renderable<T extends RenderableValues> {
    readonly values: T
    readonly state: RenderableState

    render: (variant: RenderVariant) => void
    getProgram: (variant: RenderVariant) => Program
    update: () => void
    dispose: () => void
}

export function createRenderable<T extends Values<RenderableSchema>>(renderItem: RenderItem, values: T, state: RenderableState): Renderable<T> {
    return {
        get values () { return values },
        get state () { return state },

        render: (variant: RenderVariant) => renderItem.render(variant),
        getProgram: (variant: RenderVariant) => renderItem.getProgram(variant),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    }
}

export { PointRenderable, PointSchema, PointValues } from './renderable/point'
export { MeshRenderable, MeshSchema, MeshValues } from './renderable/mesh'