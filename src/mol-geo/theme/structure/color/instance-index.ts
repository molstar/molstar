/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale } from 'mol-util/color';
import { StructureColorDataProps } from '.';
import { createInstanceColor } from '../../../util/color-data';

export function instanceIndexColorData(props: StructureColorDataProps) {
    const { group: { units } } = props
    const instanceCount = units.length

    const domain = [ 0, instanceCount - 1 ]
    const scale = ColorScale.create({ domain })
    return createInstanceColor({
        colorFn: scale.color,
        instanceCount
    })
}