/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema } from './renderable/schema';
import { RenderVariant, RenderItem } from './webgl/render-item';
import { Sphere3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';

export type RenderableState = {
    visible: boolean
    pickable: boolean
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
    let boundingSphere: Sphere3D = Sphere3D.create(Vec3.zero(), 50)

    return {
        get values () { return values },
        get state () { return state },
        get boundingSphere () {
            if (values.boundingSphere) {
                Sphere3D.copy(boundingSphere, values.boundingSphere.ref.value)
            }
            return boundingSphere
        },
        get opaque () { return values.uAlpha && values.uAlpha.ref.value === 1 },

        render: (variant: RenderVariant) => renderItem.render(variant),
        getProgram: (variant: RenderVariant) => renderItem.getProgram(variant),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    }
}