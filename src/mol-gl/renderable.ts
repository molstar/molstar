/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema } from './renderable/schema';
import { RenderVariant, RenderItem } from './webgl/render-item';
import { Sphere3D } from 'mol-math/geometry';
// import { calculateBoundingSphereFromValues } from './renderable/util';
// import { Sphere } from 'mol-geo/primitive/sphere';
import { Vec3 } from 'mol-math/linear-algebra';

export type RenderableState = {
    visible: boolean
    depthMask: boolean
}

export interface Renderable<T extends RenderableValues> {
    readonly values: T
    readonly state: RenderableState
    readonly boundingSphere: Sphere3D
    readonly opaque: boolean

    render: (variant: RenderVariant) => void
    getProgram: (variant: RenderVariant) => Program
    update: () => void
    dispose: () => void
}

export function createRenderable<T extends Values<RenderableSchema>>(renderItem: RenderItem, values: T, state: RenderableState): Renderable<T> {
    // TODO
    let boundingSphere: Sphere3D = Sphere3D.create(Vec3.zero(), 50)

    return {
        get values () { return values },
        get state () { return state },
        get boundingSphere () {
            return boundingSphere
            // TODO
            // if (boundingSphere) return boundingSphere
            // boundingSphere = calculateBoundingSphereFromValues(values)
            // return boundingSphere
        },
        get opaque () { return values.uAlpha && values.uAlpha.ref.value === 1 },

        render: (variant: RenderVariant) => renderItem.render(variant),
        getProgram: (variant: RenderVariant) => renderItem.getProgram(variant),
        update: () => {
            renderItem.update()
            // TODO
            // const valueChanges = renderItem.update()
            // if (valueChanges.attributes) boundingSphere = undefined
        },
        dispose: () => renderItem.destroy()
    }
}

export { MeshRenderable, MeshSchema, MeshValues } from './renderable/mesh'
export { PointsRenderable, PointsSchema, PointsValues } from './renderable/points'
export { LinesRenderable, LinesSchema, LinesValues } from './renderable/lines'