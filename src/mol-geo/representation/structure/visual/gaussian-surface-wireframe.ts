/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual, VisualUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { UnitsLinesVisual, DefaultUnitsLinesProps } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubesLines } from '../../../util/marching-cubes/algorithm';
import { Lines } from '../../../geometry/lines/lines';
import { GaussianDensityProps, DefaultGaussianDensityProps } from 'mol-model/structure/structure/unit/gaussian-density';

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

export const DefaultGaussianWireframeProps = {
    ...DefaultUnitsLinesProps,
    ...DefaultGaussianDensityProps,
}
export type GaussianWireframeProps = typeof DefaultGaussianWireframeProps

export function GaussianWireframeVisual(): UnitsVisual<GaussianWireframeProps> {
    return UnitsLinesVisual<GaussianWireframeProps>({
        defaultProps: DefaultGaussianWireframeProps,
        createLines: createGaussianWireframe,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianWireframeProps, currentProps: GaussianWireframeProps) => {
            if (newProps.resolutionFactor !== currentProps.resolutionFactor) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
        }
    })
}