/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementGroup } from 'mol-model/structure';
import { Mat4 } from 'mol-math/linear-algebra'

import { createUniformColor } from '../../util/color-data';
import { createUniformSize } from '../../util/size-data';
import { elementSizeData } from '../../theme/structure/size/element';
import VertexMap from '../../shape/vertex-map';
import { ColorTheme, SizeTheme } from '../../theme';
import { elementIndexColorData, elementSymbolColorData, instanceIndexColorData, chainIdColorData } from '../../theme/structure/color';

export function createTransforms(units: ReadonlyArray<Unit>) {
    const unitCount = units.length
    const transforms = new Float32Array(unitCount * 16)
    for (let i = 0; i < unitCount; i++) {
        Mat4.toArray(units[i].operator.matrix, transforms, i * 16)
    }
    return transforms
}

export function createColors(units: ReadonlyArray<Unit>, elementGroup: ElementGroup, vertexMap: VertexMap, props: ColorTheme) {
    switch (props.name) {
        case 'atom-index':
            return elementIndexColorData({ units, elementGroup, vertexMap })
        case 'chain-id':
            return chainIdColorData({ units, elementGroup, vertexMap })
        case 'element-symbol':
            return elementSymbolColorData({ units, elementGroup, vertexMap })
        case 'instance-index':
            return instanceIndexColorData({ units, elementGroup, vertexMap })
        case 'uniform':
            return createUniformColor(props)
    }
}

export function createSizes(units: ReadonlyArray<Unit>, elementGroup: ElementGroup, vertexMap: VertexMap, props: SizeTheme) {
    switch (props.name) {
        case 'uniform':
            return createUniformSize(props)
        case 'vdw':
            return elementSizeData({ units, elementGroup, vertexMap })
    }
}