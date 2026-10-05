/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../param-definition.js';
import { objectForEach } from '../object.js';
import type { ColorMap } from './color.js';

export function getColorMapParams<T extends { [k: string]: number }>(map: ColorMap<T>) {
    const colors: Record<string, PD.Color> = {};
    objectForEach(map, (_, k) => {
        colors[k] = PD.Color(map[k]);
    });
    return colors as { [k in keyof T]: PD.Color };
}