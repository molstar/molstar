/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual, VisualUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { UnitsLinesVisual, UnitsLinesParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubesLines } from '../../../util/marching-cubes/algorithm';
import { Lines } from '../../../geometry/lines/lines';
import { GaussianDensityProps, GaussianDensityParams } from 'mol-model/structure/structure/unit/gaussian-density';
import { paramDefaultValues, SelectParam, NumberParam, BooleanParam } from 'mol-util/parameter';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';

async function createGaussianWireframe(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianDensityProps, lines?: Lines): Promise<Lines> {
    const { smoothness } = props
    const { transform, field, idField } = await unit.computeGaussianDensity(props, ctx)

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    }
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx)

    Lines.transformImmediate(wireframe, transform)

    return wireframe
}

export const GaussianWireframeParams = {
    ...UnitsLinesParams,
    ...GaussianDensityParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 2, 0, 10, 0.1),
    lineSizeAttenuation: BooleanParam('Line Size Attenuation', '', false),
}
export const DefaultGaussianWireframeProps = paramDefaultValues(GaussianWireframeParams)
export type GaussianWireframeProps = typeof DefaultGaussianWireframeProps

export function GaussianWireframeVisual(): UnitsVisual<GaussianWireframeProps> {
    return UnitsLinesVisual<GaussianWireframeProps>({
        defaultProps: DefaultGaussianWireframeProps,
        createGeometry: createGaussianWireframe,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianWireframeProps, currentProps: GaussianWireframeProps) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true
            if (newProps.ignoreCache !== currentProps.ignoreCache) state.createGeometry = true
        }
    })
}