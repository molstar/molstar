/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale } from 'mol-util/color';
import { StructureColorDataProps } from '.';
import { createElementInstanceColor, ColorData } from '../../../util/color-data';

export function elementIndexColorData(props: StructureColorDataProps, colorData?: ColorData) {
    const { group: { units }, elementCount } = props
    const instanceCount = units.length

    const domain = [ 0, instanceCount * elementCount - 1 ]
    const scale = ColorScale.create({ domain })
    return createElementInstanceColor({
        colorFn: (instanceIdx, elementIdx) => scale.color(instanceIdx * elementCount + elementIdx),
        instanceCount,
        elementCount
    }, colorData)
}