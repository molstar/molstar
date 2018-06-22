/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { BaseValues } from 'mol-gl/renderable/schema';
import { BaseProps, MeshProps } from '.';
import { MeshValues, RenderableState } from 'mol-gl/renderable';
import { fillSerial } from 'mol-gl/renderable/util';

type Counts = { drawCount: number, elementCount: number, instanceCount: number }

export function createBaseValues(props: Required<BaseProps>, counts: Counts) {
    return {
        uAlpha: ValueCell.create(props.alpha),
        uInstanceCount: ValueCell.create(counts.instanceCount),
        uElementCount: ValueCell.create(counts.elementCount),
        aInstanceId: ValueCell.create(fillSerial(new Float32Array(counts.instanceCount))),
        drawCount: ValueCell.create(counts.drawCount),
        instanceCount: ValueCell.create(counts.instanceCount),
    }
}

export function createMeshValues(props: Required<MeshProps>, counts: Counts) {
    return {
        ...createBaseValues(props, counts),
        dDoubleSided: ValueCell.create(props.doubleSided),
        dFlatShaded: ValueCell.create(props.flatShaded),
        dFlipSided: ValueCell.create(props.flipSided),
        dUseFog: ValueCell.create(props.useFog),
    }
}

export function createRenderableState(props: Required<BaseProps>): RenderableState {
    return {
        visible: props.visible,
        depthMask: props.depthMask
    }
}

export function updateBaseValues(values: BaseValues, props: Required<BaseProps>) {
    ValueCell.updateIfChanged(values.uAlpha, props.alpha)
}

export function updateMeshValues(values: MeshValues, props: Required<MeshProps>) {
    updateBaseValues(values, props)
    ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided)
    ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded)
    ValueCell.updateIfChanged(values.dFlipSided, props.flipSided)
    ValueCell.updateIfChanged(values.dUseFog, props.useFog)
}

export function updateRenderableState(state: RenderableState, props: Required<BaseProps>) {
    state.visible = props.visible
    state.depthMask = props.depthMask
}