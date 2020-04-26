/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues, Values, RenderableSchema } from './renderable/schema';
import { GraphicsRenderItem, ComputeRenderItem, GraphicsRenderVariant } from './webgl/render-item';
import { ValueCell } from '../mol-util';
import { idFactory } from '../mol-util/id-factory';
import { clamp } from '../mol-math/interpolate';

const getNextRenderableId = idFactory();

export type RenderableState = {
    visible: boolean
    alphaFactor: number
    pickable: boolean
    opaque: boolean
    writeDepth: boolean,
}

export interface Renderable<T extends RenderableValues> {
    readonly id: number
    readonly materialId: number
    readonly values: T
    readonly state: RenderableState

    render: (variant: GraphicsRenderVariant) => void
    getProgram: (variant: GraphicsRenderVariant) => Program
    update: () => void
    dispose: () => void
}

export function createRenderable<T extends Values<RenderableSchema>>(renderItem: GraphicsRenderItem, values: T, state: RenderableState): Renderable<T> {
    return {
        id: getNextRenderableId(),
        materialId: renderItem.materialId,
        values,
        state,

        render: (variant: GraphicsRenderVariant) => {
            if (values.uAlpha && values.alpha) {
                ValueCell.updateIfChanged(values.uAlpha, clamp(values.alpha.ref.value * state.alphaFactor, 0, 1));
            }
            renderItem.render(variant);
        },
        getProgram: (variant: GraphicsRenderVariant) => renderItem.getProgram(variant),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    };
}

//

export interface ComputeRenderable<T extends RenderableValues> {
    readonly id: number
    readonly values: T

    render: () => void
    update: () => void
    dispose: () => void
}

export function createComputeRenderable<T extends Values<RenderableSchema>>(renderItem: ComputeRenderItem, values: T): ComputeRenderable<T> {
    return {
        id: getNextRenderableId(),
        values,

        render: () => renderItem.render('compute'),
        update: () => renderItem.update(),
        dispose: () => renderItem.destroy()
    };
}