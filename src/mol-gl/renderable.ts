/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema } from './renderable/schema';
import { RenderVariant, RenderItem } from './webgl/render-item';
import { ValueCell } from 'mol-util';
import { idFactory } from 'mol-util/id-factory';
import { clamp } from 'mol-math/interpolate';

const getNextRenderableId = idFactory()

export type RenderableState = {
    visible: boolean
    alphaFactor: number
    pickable: boolean
    opaque: boolean
}

export interface Renderable<T extends RenderableValues> {
    readonly id: number
    readonly values: T
    readonly state: RenderableState

    render: (variant: RenderVariant) => void
    getProgram: (variant: RenderVariant) => Program
    update: () => void
    dispose: () => void
}

export function createRenderable<T extends Values<RenderableSchema>>(renderItem: RenderItem, values: T, state: RenderableState): Renderable<T> {
    return {
        id: getNextRenderableId(),
        values,
        state,

        render: (variant: RenderVariant) => {
            if (values.uAlpha && values.alpha) {
                ValueCell.updateIfChanged(values.uAlpha, clamp(values.alpha.ref.value * state.alphaFactor, 0, 1))
            }
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

//

export interface ComputeRenderable<T extends RenderableValues> {
    readonly id: number
    readonly values: T

    render: () => void
    getProgram: () => Program
    update: () => void
    dispose: () => void
}

export function createComputeRenderable<T extends Values<RenderableSchema>>(renderItem: RenderItem, values: T): ComputeRenderable<T> {
    return {
        id: getNextRenderableId(),
        values,

        render: () => renderItem.render('draw'),
        getProgram: () => renderItem.getProgram('draw'),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    }
}