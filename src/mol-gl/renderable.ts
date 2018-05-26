/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Program } from './webgl/program';
import { RenderableValues } from './renderable/schema';

export type RenderableState = {
    visible: boolean
    depthMask: boolean
}

export interface Renderable<T extends RenderableValues> {
    draw: () => void
    values: T
    state: RenderableState
    name: string
    program: Program
    update: () => void
    dispose: () => void
}

export { PointRenderable, PointSchema, PointValues } from './renderable/point'
export { MeshRenderable, MeshSchema, MeshValues } from './renderable/mesh'