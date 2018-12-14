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
import { ValueCell } from 'mol-util';
import { idFactory } from 'mol-util/id-factory';

const getNextRenderableId = idFactory()

export type RenderableState = {
    visible: boolean
    pickable: boolean
    opaque: boolean
}

export interface Renderable<T extends RenderableValues> {
    readonly id: number
    readonly values: T
    readonly state: RenderableState
    readonly boundingSphere: Sphere3D
    readonly invariantBoundingSphere: Sphere3D
    readonly z: number

    render: (variant: RenderVariant) => void
    getProgram: (variant: RenderVariant) => Program
    update: () => void
    dispose: () => void
}

export function createRenderable<T extends Values<RenderableSchema>>(renderItem: RenderItem, values: T, state: RenderableState): Renderable<T> {
    const boundingSphere: Sphere3D = Sphere3D.create(Vec3.zero(), 0)
    const invariantBoundingSphere: Sphere3D = Sphere3D.create(Vec3.zero(), 0)

    return {
        id: getNextRenderableId(),
        values,
        state,
        get boundingSphere () {
            if (values.boundingSphere) {
                Sphere3D.copy(boundingSphere, values.boundingSphere.ref.value)
            }
            return boundingSphere
        },
        get invariantBoundingSphere () {
            if (values.invariantBoundingSphere) {
                Sphere3D.copy(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)
            }
            return invariantBoundingSphere
        },
        get z () {
            return boundingSphere.center[2]
        },

        render: (variant: RenderVariant) => {
            if (values.uPickable) {
                ValueCell.updateIfChanged(values.uPickable, state.pickable ? 1 : 0)
            }
            renderItem.render(variant)
        },
        getProgram: (variant: RenderVariant) => renderItem.getProgram(variant),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    }
}