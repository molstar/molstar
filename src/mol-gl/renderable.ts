/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import Attribute from './attribute'
import PointRenderable from './renderable/point'
import MeshRenderable from './renderable/mesh'

export type AttributesMutator<T extends AttributesData> = (data: T) => (boolean | void)
export type AttributesData = { [k: string]: Helpers.TypedArray }
export type Attributes<T extends AttributesData> = { [K in keyof T]: Attribute<T[K]> }
export type AttributesBuffers<T extends AttributesData> = { [K in keyof T]: REGL.AttributeConfig }

export interface Renderable {
    draw(): void
    dispose(): void
    stats: REGL.CommandStats
    name: string
    // isPicking: () => boolean
    // isVisible: () => boolean
    // isTransparent: () => boolean
}

export { PointRenderable, MeshRenderable }