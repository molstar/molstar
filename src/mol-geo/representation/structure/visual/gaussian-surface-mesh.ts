/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual, VisualUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../geometry/mesh/mesh';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { computeMarchingCubes } from '../../../util/marching-cubes/algorithm';
import { SizeThemeProps } from 'mol-view/theme/size';
import { Color } from 'mol-util/color';
import { computeGaussianDensity } from './util/gaussian';
import { ColorThemeProps } from 'mol-view/theme/color';

export interface GaussianSurfaceMeshProps {
    sizeTheme: SizeThemeProps

    resolutionFactor: number
    probeRadius: number
    isoValue: number
}

async function createGaussianSurfaceMesh(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const { isoValue } = props

    console.time('surface density')
    const { transform, field } = await computeGaussianDensity(unit, structure, props).runAsChild(ctx)
    console.timeEnd('surface density')

    console.time('surface mc')
    const surface = await computeMarchingCubes({
        isoLevel: Math.exp(-isoValue),
        scalarField: field,
        oldSurface: mesh
    }).runAsChild(ctx)
    console.timeEnd('surface mc')

    Mesh.transformImmediate(surface, transform)
    Mesh.computeNormalsImmediate(surface)

    return surface;
}

export const DefaultGaussianSurfaceProps = {
    ...DefaultUnitsMeshProps,
    linearSegments: 8,
    radialSegments: 12,
    aspectRatio: 5,
    arrowFactor: 1.5,

    flipSided: true,
    // flatShaded: true,
    alpha: 0.7,
    colorTheme: { name: 'uniform' as 'uniform', value: Color(0xDDDDDD) } as ColorThemeProps,

    resolutionFactor: 7,
    probeRadius: 0,
    isoValue: 1.5,
}
export type GaussianSurfaceProps = typeof DefaultGaussianSurfaceProps

export function GaussianSurfaceVisual(): UnitsVisual<GaussianSurfaceProps> {
    return UnitsMeshVisual<GaussianSurfaceProps>({
        defaultProps: DefaultGaussianSurfaceProps,
        createMesh: createGaussianSurfaceMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianSurfaceProps, currentProps: GaussianSurfaceProps) => {
            if (newProps.resolutionFactor !== currentProps.resolutionFactor) state.createGeometry = true
            if (newProps.probeRadius !== currentProps.probeRadius) state.createGeometry = true
            if (newProps.isoValue !== currentProps.isoValue) state.createGeometry = true
        }
    })
}