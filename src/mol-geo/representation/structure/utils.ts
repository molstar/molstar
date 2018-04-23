/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementGroup } from 'mol-model/structure';
import { Mat4 } from 'mol-math/linear-algebra'

import { ColorScale } from 'mol-util/color';
import { createUniformColor, createInstanceColor, createElementInstanceColor } from '../../util/color-data';
import { createUniformSize } from '../../util/size-data';
import { vdwSizeData } from '../../theme/structure/size/vdw';
import VertexMap from '../../shape/vertex-map';
import { ColorTheme, SizeTheme } from '../../theme';
import { elementSymbolColorData } from '../../theme/structure/color/element';
import { OrderedSet } from 'mol-data/int';

export function createTransforms(units: ReadonlyArray<Unit>) {
    const unitCount = units.length
    const transforms = new Float32Array(unitCount * 16)
    for (let i = 0; i < unitCount; i++) {
        Mat4.toArray(units[i].operator.matrix, transforms, i * 16)
    }
    return transforms
}

export function createColors(units: ReadonlyArray<Unit>, elementGroup: ElementGroup, vertexMap: VertexMap, props: ColorTheme) {
    const instanceCount = units.length
    const elementCount = OrderedSet.size(elementGroup.elements)
    switch (props.name) {
        case 'uniform':
            return createUniformColor(props)
        case 'instance-id':
            const instanceDomain = props.domain ? props.domain : [ 0, instanceCount - 1 ]
            const instanceScale = ColorScale.create({ domain: instanceDomain })
            return createInstanceColor({ colorFn: instanceScale.color, instanceCount })
        case 'atom-id':
            const atomDomain = props.domain ? props.domain : [ 0, instanceCount * elementCount - 1 ]
            const atomScale = ColorScale.create({ domain: atomDomain })
            return createElementInstanceColor({
                colorFn: (unitIdx, elementIdx) => atomScale.color(unitIdx * elementCount + elementIdx),
                instanceCount,
                vertexMap
            })
        case 'element-symbol':
            return elementSymbolColorData({ units, elementGroup, vertexMap })
    }
}

export function createSizes(units: ReadonlyArray<Unit>, elementGroup: ElementGroup, vertexMap: VertexMap, props: SizeTheme) {
    switch (props.name) {
        case 'uniform':
            return createUniformSize(props)
        case 'vdw':
            return vdwSizeData({ units, elementGroup, vertexMap })
    }
}