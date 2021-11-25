/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { MolecularSurfaceCalculationParams } from '../../../mol-math/geometry/molecular-surface';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { computeUnitMolecularSurface } from './util/molecular-surface';
import { computeMarchingCubesMesh } from '../../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { CommonSurfaceParams, getUnitExtraRadius } from './util/common';
import { Sphere3D } from '../../../mol-math/geometry';
import { MeshValues } from '../../../mol-gl/renderable/mesh';
import { Texture } from '../../../mol-gl/webgl/texture';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { applyMeshColorSmoothing } from '../../../mol-geo/geometry/mesh/color-smoothing';
import { ColorSmoothingParams, getColorSmoothingProps } from '../../../mol-geo/geometry/base';

export const MolecularSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...MolecularSurfaceCalculationParams,
    ...CommonSurfaceParams,
    ...ColorSmoothingParams,
};
export type MolecularSurfaceMeshParams = typeof MolecularSurfaceMeshParams
export type MolecularSurfaceMeshProps = PD.Values<MolecularSurfaceMeshParams>

type MolecularSurfaceMeta = {
    resolution?: number
    colorTexture?: Texture
}

//

async function createMolecularSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: MolecularSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const { transform, field, idField, resolution } = await computeUnitMolecularSurface(structure, unit, props).runInContext(ctx.runtime);

    const params = {
        isoLevel: props.probeRadius,
        scalarField: field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

    if (props.includeParent) {
        const iterations = Math.ceil(2 / props.resolution);
        Mesh.smoothEdges(surface, { iterations, maxNewEdgeLength: Math.sqrt(2) });
    }

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) Mesh.uniformTriangleGroup(surface);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.probeRadius + getUnitExtraRadius(unit));
    surface.setBoundingSphere(sphere);
    (surface.meta as MolecularSurfaceMeta).resolution = resolution;

    return surface;
}

export function MolecularSurfaceMeshVisual(materialId: number): UnitsVisual<MolecularSurfaceMeshParams> {
    return UnitsMeshVisual<MolecularSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(MolecularSurfaceMeshParams),
        createGeometry: createMolecularSurfaceMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<MolecularSurfaceMeshParams>, currentProps: PD.Values<MolecularSurfaceMeshParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.probeRadius !== currentProps.probeRadius) state.createGeometry = true;
            if (newProps.probePositions !== currentProps.probePositions) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        processValues: (values: MeshValues, geometry: Mesh, props: PD.Values<MolecularSurfaceMeshParams>, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as MolecularSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp.resolution, csp.stride, webgl, colorTexture);
                (geometry.meta as MolecularSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as MolecularSurfaceMeta).colorTexture?.destroy();
        }
    }, materialId);
}