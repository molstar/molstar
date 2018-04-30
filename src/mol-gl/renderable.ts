/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import PointRenderable from './renderable/point'
import MeshRenderable from './renderable/mesh'

export interface Renderable<T> {
    draw: () => void
    name: string
    update: (newProps: T) => void
    dispose: () => void
}

export { PointRenderable, MeshRenderable }