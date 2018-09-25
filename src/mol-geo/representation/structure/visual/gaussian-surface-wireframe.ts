/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual, VisualUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../geometry/mesh/mesh';
import { UnitsLinesVisual, DefaultUnitsLinesProps } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubes } from '../../../util/marching-cubes/algorithm';
import { Lines } from '../../../geometry/lines/lines';
import { LinesBuilder } from '../../../geometry/lines/lines-builder';
import { GaussianDensityProps, DefaultGaussianDensityProps } from 'mol-model/structure/structure/unit/gaussian-density';

async function createGaussianWireframe(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianDensityProps, lines?: Lines): Promise<Lines> {
    const { smoothness } = props
    const { transform, field, idField } = await unit.computeGaussianDensity(props, ctx)

    const surface = await computeMarchingCubes({
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    }).runAsChild(ctx)

    Mesh.transformImmediate(surface, transform)

    const vb = surface.vertexBuffer.ref.value
    const ib = surface.indexBuffer.ref.value
    const gb = surface.groupBuffer.ref.value

    const builder = LinesBuilder.create(surface.triangleCount * 3, surface.triangleCount / 10, lines)

    // TODO avoid duplicate lines and move to '../../../geometry/lines/lines' as Lines.fromMesh()
    for (let i = 0, il = surface.triangleCount * 3; i < il; i += 3) {
        const i0 = ib[i], i1 = ib[i + 1], i2 = ib[i + 2];
        const x0 = vb[i0 * 3], y0 = vb[i0 * 3 + 1], z0 = vb[i0 * 3 + 2];
        const x1 = vb[i1 * 3], y1 = vb[i1 * 3 + 1], z1 = vb[i1 * 3 + 2];
        const x2 = vb[i2 * 3], y2 = vb[i2 * 3 + 1], z2 = vb[i2 * 3 + 2];
        builder.add(x0, y0, z0, x1, y1, z1, gb[i0])
        builder.add(x0, y0, z0, x2, y2, z2, gb[i0])
        builder.add(x1, y1, z1, x2, y2, z2, gb[i1])
    }

    const l = builder.getLines();
    console.log(l)
    return l
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