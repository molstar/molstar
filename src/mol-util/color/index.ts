/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Color from './color'
import { ColorScale } from './scale';

type ColorTable<T extends { [k: string]: number[] }> = { [k in keyof T]: Color[] }
export function ColorTable<T extends { [k: string]: number[] }>(o: T) { return o as ColorTable<T> }

type ColorMap<T extends { [k: string]: number }> = { [k in keyof T]: Color }
export function ColorMap<T extends { [k: string]: number }>(o: T) { return o as ColorMap<T> }

export { Color, ColorScale }
export { ColorBrewer, ColorNames } from './tables'