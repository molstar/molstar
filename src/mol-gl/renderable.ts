/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import PointRenderable from './renderable/point'
import MeshRenderable from './renderable/mesh'
import { Program } from './webgl/program';

export type BaseProps = {
    objectId: number
    alpha: number
    visible: boolean
    depthMask: boolean

    flatShaded?: boolean
    doubleSided?: boolean
    flipSided?: boolean
}

export interface Renderable<T> {
    draw: () => void
    name: string
    program: Program
    update: (newProps: T) => void
    dispose: () => void
}

export { PointRenderable, MeshRenderable }