/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsTextureMeshParams, UnitsVisual, UnitsMeshVisual, UnitsTextureMeshVisual } from '../units-visual';
import { GaussianDensityParams, computeUnitGaussianDensity, GaussianDensityTextureProps, computeUnitGaussianDensityTexture2d, GaussianDensityProps, computeStructureGaussianDensity } from './util/gaussian';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '../../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { TextureMesh } from '../../../mol-geo/geometry/texture-mesh/texture-mesh';
import { calcActiveVoxels } from '../../../mol-gl/compute/marching-cubes/active-voxels';
import { createHistogramPyramid } from '../../../mol-gl/compute/histogram-pyramid/reduction';
import { createIsosurfaceBuffers } from '../../../mol-gl/compute/marching-cubes/isosurface';
import { Sphere3D } from '../../../mol-math/geometry';
import { ComplexVisual, ComplexMeshParams, ComplexMeshVisual } from '../complex-visual';

export const GaussianSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...UnitsTextureMeshParams,
    ...GaussianDensityParams,
};
export type GaussianSurfaceMeshParams = typeof GaussianSurfaceMeshParams

export function getGaussianSurfaceVisual(webgl?: WebGLContext) {
    return webgl && webgl.extensions.drawBuffers ? GaussianSurfaceTextureMeshVisual : GaussianSurfaceMeshVisual;
}

//

async function createGaussianSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props;
    const { transform, field, idField } = await computeUnitGaussianDensity(structure, unit, props, ctx.webgl).runInContext(ctx.runtime);

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) Mesh.uniformTriangleGroup(surface);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.radiusOffset);
    surface.setBoundingSphere(sphere);

    return surface;
}

export function GaussianSurfaceMeshVisual(materialId: number): UnitsVisual<GaussianSurfaceMeshParams> {
    return UnitsMeshVisual<GaussianSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(GaussianSurfaceMeshParams),
        createGeometry: createGaussianSurfaceMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianSurfaceMeshParams>, currentProps: PD.Values<GaussianSurfaceMeshParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true;
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true;
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
        }
    }, materialId);
}

//

export const StructureGaussianSurfaceMeshParams = {
    ...ComplexMeshParams,
    ...GaussianDensityParams,
    ignoreHydrogens: PD.Boolean(false),
};
export type StructureGaussianSurfaceMeshParams = typeof StructureGaussianSurfaceMeshParams

async function createStructureGaussianSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props;
    const { transform, field, idField } = await computeStructureGaussianDensity(structure, props, ctx.webgl).runInContext(ctx.runtime);

    const params = {
        isoLevel: Math.exp(-smoothness),
        scalarField: field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) Mesh.uniformTriangleGroup(surface);

    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, props.radiusOffset);
    surface.setBoundingSphere(sphere);

    return surface;
}

export function StructureGaussianSurfaceMeshVisual(materialId: number): ComplexVisual<StructureGaussianSurfaceMeshParams> {
    return ComplexMeshVisual<StructureGaussianSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(StructureGaussianSurfaceMeshParams),
        createGeometry: createStructureGaussianSurfaceMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianSurfaceMeshParams>, currentProps: PD.Values<GaussianSurfaceMeshParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true;
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true;
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
        }
    }, materialId);
}

//

async function createGaussianSurfaceTextureMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityTextureProps, textureMesh?: TextureMesh): Promise<TextureMesh> {
    if (!ctx.webgl) throw new Error('webgl context required to create gaussian surface texture-mesh');
    const isoLevel = Math.exp(-props.smoothness);

    const densityTextureData = await computeUnitGaussianDensityTexture2d(structure, unit, props, ctx.webgl).runInContext(ctx.runtime);
    // console.log(densityTextureData)
    // console.log('vertexGroupTexture', readTexture(ctx.webgl, densityTextureData.texture))
    // ctx.webgl.waitForGpuCommandsCompleteSync()

    const activeVoxelsTex = calcActiveVoxels(ctx.webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, isoLevel, densityTextureData.gridTexScale);
    // ctx.webgl.waitForGpuCommandsCompleteSync()

    const compacted = createHistogramPyramid(ctx.webgl, activeVoxelsTex, densityTextureData.gridTexScale);
    // ctx.webgl.waitForGpuCommandsCompleteSync()

    const gv = createIsosurfaceBuffers(ctx.webgl, activeVoxelsTex, densityTextureData.texture, compacted, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.transform, isoLevel, textureMesh ? textureMesh.vertexGroupTexture.ref.value : undefined, textureMesh ? textureMesh.normalTexture.ref.value : undefined);
    // ctx.webgl.waitForGpuCommandsCompleteSync()

    // const boundingSphere = Sphere3D.zero()
    // Sphere3D.addVec3(boundingSphere, boundingSphere, densityTextureData.gridDimension)
    const boundingSphere = Sphere3D.fromBox3D(Sphere3D(), densityTextureData.bbox);
    const surface = TextureMesh.create(gv.vertexCount, 1, gv.vertexGroupTexture, gv.normalTexture, boundingSphere, textureMesh);

    // ctx.webgl.waitForGpuCommandsCompleteSync()
    return surface;
}

export function GaussianSurfaceTextureMeshVisual(materialId: number): UnitsVisual<GaussianSurfaceMeshParams> {
    return UnitsTextureMeshVisual<GaussianSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(GaussianSurfaceMeshParams),
        createGeometry: createGaussianSurfaceTextureMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianSurfaceMeshParams>, currentProps: PD.Values<GaussianSurfaceMeshParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true;
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
        }
    }, materialId);
}