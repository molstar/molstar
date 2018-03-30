/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import Attribute from './attribute'
import Point from './renderable/point'

export type AttributesMutator<T extends AttributesData> = (data: T) => (boolean | void)
export type AttributesData = { [k: string]: Helpers.TypedArray }
export type Attributes<T extends AttributesData> = { [K in keyof T]: Attribute<T[K]> }
export type AttributesBuffers<T extends AttributesData> = { [K in keyof T]: REGL.AttributeConfig }

export interface Renderable<T extends AttributesData> {
    draw(): void
}

export { Point }