/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementGroup, Model } from 'mol-model/structure';

import { StructureColorDataProps } from '.';
import { createAttributeOrElementColor } from '../../../util/color-data';
import { ColorScale } from 'mol-util/color';

function createChainIdMap(model: Model) {
    const { chains } = model.hierarchy
    const { label_asym_id } = chains

    const map = new Map<string, number>()
    let index = 0

    for (let i = 0, il = chains._rowCount; i < il; ++i) {
        const chainId = label_asym_id.value(i)
        if (map.get(chainId) === undefined) {
            map.set(chainId, index)
            index += 1
        }
    }
    return { map, count: index }
}

export function chainIdColorData(props: StructureColorDataProps) {
    const { units, elementGroup, vertexMap } = props
    const unit = units[0]

    const { chains, chainSegments } = unit.model.hierarchy
    const { label_asym_id } = chains
    const { map, count } = createChainIdMap(unit.model)

    const domain = [ 0, count - 1 ]
    const scale = ColorScale.create({ domain })

    return createAttributeOrElementColor(vertexMap, {
        colorFn: (elementIdx: number) => {
            const aI = ElementGroup.getAt(elementGroup, elementIdx);
            const cI = chainSegments.segmentMap[aI]
            const chainId = label_asym_id.value(cI)

            return scale.color(map.get(chainId) || 0)
        },
        vertexMap
    })
}