/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsTextureMeshParams, UnitsVisual, UnitsMeshVisual, UnitsTextureMeshVisual } from '../units-visual';
import { GaussianDensityParams, computeUnitGaussianDensity, computeUnitGaussianDensityTexture2d, GaussianDensityProps, computeStructureGaussianDensity, computeStructureGaussianDensityTexture2d } from './util/gaussian';
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
import { ComplexVisual, ComplexMeshParams, ComplexMeshVisual, ComplexTextureMeshVisual, ComplexTextureMeshParams } from '../complex-visual';
import { getUnitExtraRadius, getStructureExtraRadius } from './util/common';
import { WebGLContext } from '../../../mol-gl/webgl/context';

const SharedParams = {
    ...GaussianDensityParams,
    ignoreHydrogens: PD.Boolean(false),
    useGpu: PD.Boolean(false),
};

export const GaussianSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...UnitsTextureMeshParams,
    ...SharedParams,
};
export type GaussianSurfaceMeshParams = typeof GaussianSurfaceMeshParams

export const StructureGaussianSurfaceMeshParams = {
    ...ComplexMeshParams,
    ...ComplexTextureMeshParams,
    ...SharedParams,
};
export type StructureGaussianSurfaceMeshParams = typeof StructureGaussianSurfaceMeshParams

function gpuSupport(webgl: WebGLContext) {
    return webgl.extensions.colorBufferFloat && webgl.extensions.textureFloat && webgl.extensions.blendMinMax && webgl.extensions.drawBuffers;
}

export function GaussianSurfaceVisual(materialId: number, props?: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) {
    return props?.useGpu && webgl && gpuSupport(webgl)
        ? GaussianSurfaceTextureMeshVisual(materialId)
        : GaussianSurfaceMeshVisual(materialId);
}

export function StructureGaussianSurfaceVisual(materialId: number, props?: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) {
    return props?.useGpu && webgl && gpuSupport(webgl)
        ? StructureGaussianSurfaceTextureMeshVisual(materialId)
        : StructureGaussianSurfaceMeshVisual(materialId);
}

//

async function createGaussianSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props;
    const { transform, field, idField, radiusFactor } = await computeUnitGaussianDensity(structure, unit, props, ctx.webgl).runInContext(ctx.runtime);

    const params = {
        isoLevel: Math.exp(-smoothness) / radiusFactor,
        scalarField: field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) Mesh.uniformTriangleGroup(surface);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.radiusOffset + getUnitExtraRadius(unit));
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
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
        },
        mustRecreate: (props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return props.useGpu && !!webgl;
        }
    }, materialId);
}

//

async function createStructureGaussianSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props;
    const { transform, field, idField, radiusFactor } = await computeStructureGaussianDensity(structure, props, ctx.webgl).runInContext(ctx.runtime);

    const params = {
        isoLevel: Math.exp(-smoothness) / radiusFactor,
        scalarField: field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) Mesh.uniformTriangleGroup(surface);

    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, props.radiusOffset + getStructureExtraRadius(structure));
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
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
        },
        mustRecreate: (props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return props.useGpu && !!webgl;
        }
    }, materialId);
}

//

const GaussianSurfaceName = 'gaussian-surface';

async function createGaussianSurfaceTextureMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, textureMesh?: TextureMesh): Promise<TextureMesh> {
    if (!ctx.webgl) throw new Error('webgl context required to create gaussian surface texture-mesh');

    const { namedTextures, resources, extensions: { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } } = ctx.webgl;
    if (!namedTextures[GaussianSurfaceName]) {
        namedTextures[GaussianSurfaceName] = colorBufferHalfFloat && textureHalfFloat
            ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
            : colorBufferFloat && textureFloat
                ? resources.texture('image-float32', 'rgba', 'float', 'linear')
                : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    }

    // console.time('computeUnitGaussianDensityTexture2d');
    const densityTextureData = await computeUnitGaussianDensityTexture2d(structure, unit, true, props, ctx.webgl, namedTextures[GaussianSurfaceName]).runInContext(ctx.runtime);
    // console.log(densityTextureData);
    // console.log('vertexGroupTexture', readTexture(ctx.webgl, densityTextureData.texture));
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('computeUnitGaussianDensityTexture2d');

    const isoLevel = Math.exp(-props.smoothness) / densityTextureData.radiusFactor;

    // console.time('calcActiveVoxels');
    const activeVoxelsTex = calcActiveVoxels(ctx.webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, isoLevel, densityTextureData.gridTexScale);
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('calcActiveVoxels');

    // console.time('createHistogramPyramid');
    const compacted = createHistogramPyramid(ctx.webgl, activeVoxelsTex, densityTextureData.gridTexScale, densityTextureData.gridTexDim);
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('createHistogramPyramid');

    // console.time('createIsosurfaceBuffers');
    const gv = createIsosurfaceBuffers(ctx.webgl, activeVoxelsTex, densityTextureData.texture, compacted, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.transform, isoLevel, textureMesh ? textureMesh.vertexGroupTexture.ref.value : undefined, textureMesh ? textureMesh.normalTexture.ref.value : undefined);
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('createIsosurfaceBuffers');

    const boundingSphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.radiusOffset + getStructureExtraRadius(structure));
    const surface = TextureMesh.create(gv.vertexCount, 1, gv.vertexGroupTexture, gv.normalTexture, boundingSphere, textureMesh);
    // console.log({
    //     renderables: ctx.webgl.namedComputeRenderables,
    //     framebuffers: ctx.webgl.namedFramebuffers,
    //     textures: ctx.webgl.namedTextures,
    // });
    // ctx.webgl.waitForGpuCommandsCompleteSync();
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
        },
        mustRecreate: (props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return !props.useGpu || !webgl;
        }
    }, materialId);
}

//

async function createStructureGaussianSurfaceTextureMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, textureMesh?: TextureMesh): Promise<TextureMesh> {
    if (!ctx.webgl) throw new Error('webgl context required to create structure gaussian surface texture-mesh');

    const { namedTextures, resources, extensions: { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } } = ctx.webgl;
    if (!namedTextures[GaussianSurfaceName]) {
        namedTextures[GaussianSurfaceName] = colorBufferHalfFloat && textureHalfFloat
            ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
            : colorBufferFloat && textureFloat
                ? resources.texture('image-float32', 'rgba', 'float', 'linear')
                : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    }

    // console.time('computeUnitGaussianDensityTexture2d');
    const densityTextureData = await computeStructureGaussianDensityTexture2d(structure, true, props, ctx.webgl, namedTextures[GaussianSurfaceName]).runInContext(ctx.runtime);
    // console.log(densityTextureData);
    // console.log('vertexGroupTexture', readTexture(ctx.webgl, densityTextureData.texture));
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('computeUnitGaussianDensityTexture2d');

    const isoLevel = Math.exp(-props.smoothness) / densityTextureData.radiusFactor;

    // console.time('calcActiveVoxels');
    const activeVoxelsTex = calcActiveVoxels(ctx.webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, isoLevel, densityTextureData.gridTexScale);
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('calcActiveVoxels');

    // console.time('createHistogramPyramid');
    const compacted = createHistogramPyramid(ctx.webgl, activeVoxelsTex, densityTextureData.gridTexScale, densityTextureData.gridTexDim);
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('createHistogramPyramid');

    // console.time('createIsosurfaceBuffers');
    const gv = createIsosurfaceBuffers(ctx.webgl, activeVoxelsTex, densityTextureData.texture, compacted, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.transform, isoLevel, textureMesh ? textureMesh.vertexGroupTexture.ref.value : undefined, textureMesh ? textureMesh.normalTexture.ref.value : undefined);
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    // console.timeEnd('createIsosurfaceBuffers');

    const boundingSphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, props.radiusOffset + getStructureExtraRadius(structure));
    const surface = TextureMesh.create(gv.vertexCount, 1, gv.vertexGroupTexture, gv.normalTexture, boundingSphere, textureMesh);
    // console.log({
    //     renderables: ctx.webgl.namedComputeRenderables,
    //     framebuffers: ctx.webgl.namedFramebuffers,
    //     textures: ctx.webgl.namedTextures,
    // });
    // ctx.webgl.waitForGpuCommandsCompleteSync();
    return surface;
}

export function StructureGaussianSurfaceTextureMeshVisual(materialId: number): ComplexVisual<StructureGaussianSurfaceMeshParams> {
    return ComplexTextureMeshVisual<StructureGaussianSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(StructureGaussianSurfaceMeshParams),
        createGeometry: createStructureGaussianSurfaceTextureMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureGaussianSurfaceMeshParams>, currentProps: PD.Values<StructureGaussianSurfaceMeshParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true;
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true;
            if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
            if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
        },
        mustRecreate: (props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return !props.useGpu || !webgl;
        }
    }, materialId);
}