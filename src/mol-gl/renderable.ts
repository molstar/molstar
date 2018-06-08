/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema, BaseValues } from './renderable/schema';
import { RenderVariant, RenderItem } from './webgl/render-item';
import { Sphere3D } from 'mol-math/geometry';
import { calculateBoundingSphere } from './renderable/util';

export type RenderableState = {
    visible: boolean
    depthMask: boolean
}

export interface Renderable<T extends RenderableValues & BaseValues> {
    readonly values: T
    readonly state: RenderableState
    readonly boundingSphere: Sphere3D

    render: (variant: RenderVariant) => void
    getProgram: (variant: RenderVariant) => Program
    update: () => void
    dispose: () => void
}

export function createRenderable<T extends Values<RenderableSchema> & BaseValues>(renderItem: RenderItem, values: T, state: RenderableState): Renderable<T> {
    const position = values.aPosition.ref.value
    const transform = values.aTransform.ref.value
    const positionCount = values.drawCount.ref.value / 3 / 3
    const transformCount = values.instanceCount.ref.value

    let boundingSphere: Sphere3D | undefined
    
    return {
        get values () { return values },
        get state () { return state },
        get boundingSphere () {
            if (boundingSphere) return boundingSphere
            boundingSphere = calculateBoundingSphere({ position, positionCount, transform, transformCount })
            return boundingSphere
        },

        render: (variant: RenderVariant) => renderItem.render(variant),
        getProgram: (variant: RenderVariant) => renderItem.getProgram(variant),
        update: () => {
            renderItem.update()
            boundingSphere = undefined
        },
        dispose: () => renderItem.destroy()
    }
}

export { PointRenderable, PointSchema, PointValues } from './renderable/point'
export { MeshRenderable, MeshSchema, MeshValues } from './renderable/mesh'