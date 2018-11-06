/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { UnitsLinesVisual, UnitsLinesParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { GaussianDensityProps, GaussianDensityParams } from 'mol-model/structure/structure/unit/gaussian-density';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { Lines } from 'mol-geo/geometry/lines/lines';
import { computeMarchingCubesLines } from 'mol-geo/util/marching-cubes/algorithm';
import { VisualContext } from 'mol-repr';

async function createGaussianWireframe(ctx: VisualContext, unit: Unit, structure: Structure, props: GaussianDensityProps, lines?: Lines): Promise<Lines> {
    const { smoothness } = props
    const { transform, field, idField } = await unit.computeGaussianDensity(props, ctx.runtime)

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    }
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime)

    Lines.transformImmediate(wireframe, transform)

    return wireframe
}

export const GaussianWireframeParams = {
    ...UnitsLinesParams,
    ...GaussianDensityParams,
    sizeTheme: PD.SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: PD.NumberParam('Size Value', '', 2, 0, 10, 0.1),
    lineSizeAttenuation: PD.BooleanParam('Line Size Attenuation', '', false),
}
export const DefaultGaussianWireframeProps = PD.paramDefaultValues(GaussianWireframeParams)
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