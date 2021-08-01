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
import { applyMeshColorSmoothing, ColorSmoothingParams, getColorSmoothingProps } from './util/color';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';

export const MolecularSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...MolecularSurfaceCalculationParams,
    ...CommonSurfaceParams,
    ...ColorSmoothingParams,
    clipSphere: PD.MappedStatic('off', {
        on: PD.Group({
            center: PD.Vec3(Vec3()),
            radius: PD.Numeric(1, { min: 0, max: 100, step: 0.1 })
        }),
        off: PD.Group({})
    }),
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

    if (props.clipSphere.name === 'on') {
        const { center, radius } = props.clipSphere.params;
        const c = Vec3.transformMat4(Vec3(), center, Mat4.invert(Mat4(), transform));
        const r = radius / Mat4.getMaxScaleOnAxis(transform);
        Mesh.trimByPositionTest(surface, p => Vec3.distance(p, c) < r);

        const iterations = Math.ceil(2 / props.resolution);
        Mesh.smoothEdges(surface, iterations);
    }

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) Mesh.uniformTriangleGroup(surface);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.probeRadius + getUnitExtraRadius(unit));
    surface.setBoundingSphere(sphere);
    (surface.meta.resolution as MolecularSurfaceMeta['resolution']) = resolution;

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

            if (newProps.clipSphere.name !== currentProps.clipSphere.name) {
                state.createGeometry = true;
            } else if (newProps.clipSphere.name === 'on' && currentProps.clipSphere.name === 'on') {
                if (!Vec3.exactEquals(newProps.clipSphere.params.center, currentProps.clipSphere.params.center)) state.createGeometry = true;
                if (newProps.clipSphere.params.radius !== currentProps.clipSphere.params.radius) state.createGeometry = true;
            }
        },
        processValues: (values: MeshValues, geometry: Mesh, props: PD.Values<MolecularSurfaceMeshParams>, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as MolecularSurfaceMeta;
            const csp = getColorSmoothingProps(props, theme, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp.resolution, csp.stride, webgl, colorTexture);
                (geometry.meta.colorTexture as MolecularSurfaceMeta['colorTexture']) = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as MolecularSurfaceMeta).colorTexture?.destroy();
        }
    }, materialId);
}