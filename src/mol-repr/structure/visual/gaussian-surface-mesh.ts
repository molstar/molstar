/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsTextureMeshParams, UnitsVisual, UnitsMeshVisual, UnitsTextureMeshVisual, StructureGroup } from '../units-visual';
import { GaussianDensityParams, computeUnitGaussianDensity, computeUnitGaussianDensityTexture2d, GaussianDensityProps, computeStructureGaussianDensity, computeStructureGaussianDensityTexture2d } from './util/gaussian';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '../../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { TextureMesh } from '../../../mol-geo/geometry/texture-mesh/texture-mesh';
import { extractIsosurface } from '../../../mol-gl/compute/marching-cubes/isosurface';
import { Sphere3D } from '../../../mol-math/geometry';
import { ComplexVisual, ComplexMeshParams, ComplexMeshVisual, ComplexTextureMeshVisual, ComplexTextureMeshParams } from '../complex-visual';
import { getUnitExtraRadius, getStructureExtraRadius, getVolumeSliceInfo } from './util/common';
import { WebGLContext } from '../../../mol-gl/webgl/context';

const SharedParams = {
    ...GaussianDensityParams,
    ignoreHydrogens: PD.Boolean(false),
    tryUseGpu: PD.Boolean(true),
};
type SharedParams = typeof SharedParams

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

function suitableForGpu(structure: Structure, props: PD.Values<SharedParams>, webgl: WebGLContext) {
    // lower resolutions are about as fast on CPU vs integrated GPU,
    // very low resolutions have artifacts when calculated on GPU
    if (props.resolution > 1) return false;
    // the GPU is much more memory contraint, especially true for integrated GPUs,
    // being conservative here still allows for small and medium sized assemblies
    const d = webgl.maxTextureSize / 3;
    const { areaCells, maxAreaCells } = getVolumeSliceInfo(structure.boundary.box, props.resolution, d * d);
    return areaCells < maxAreaCells;
}

export function GaussianSurfaceVisual(materialId: number, structure: Structure, props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) {
    if (props.tryUseGpu && webgl && gpuSupport(webgl) && suitableForGpu(structure, props, webgl)) {
        return GaussianSurfaceTextureMeshVisual(materialId);
    }
    return GaussianSurfaceMeshVisual(materialId);
}

export function StructureGaussianSurfaceVisual(materialId: number, structure: Structure, props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) {
    if (props.tryUseGpu && webgl && gpuSupport(webgl) && suitableForGpu(structure, props, webgl)) {
        return StructureGaussianSurfaceTextureMeshVisual(materialId);
    }
    return StructureGaussianSurfaceMeshVisual(materialId);
}

//

async function createGaussianSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props;
    const { transform, field, idField, radiusFactor } = await computeUnitGaussianDensity(structure, unit, props).runInContext(ctx.runtime);

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
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return props.tryUseGpu && !!webgl && suitableForGpu(structureGroup.structure, props, webgl);
        }
    }, materialId);
}

//

async function createStructureGaussianSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness } = props;
    const { transform, field, idField, radiusFactor } = await computeStructureGaussianDensity(structure, props).runInContext(ctx.runtime);

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
        mustRecreate: (structure: Structure, props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return props.tryUseGpu && !!webgl && suitableForGpu(structure, props, webgl);
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

    const buffer = textureMesh?.doubleBuffer.get();
    const gv = extractIsosurface(ctx.webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.gridTexScale, densityTextureData.transform, isoLevel, true, buffer?.vertex, buffer?.group, buffer?.normal);

    const boundingSphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, props.radiusOffset + getStructureExtraRadius(structure));
    const surface = TextureMesh.create(gv.vertexCount, 1, gv.vertexTexture, gv.groupTexture, gv.normalTexture, boundingSphere, textureMesh);

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
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return !props.tryUseGpu || !webgl || !suitableForGpu(structureGroup.structure, props, webgl);
        },
        dispose: (geometry: TextureMesh) => {
            geometry.vertexTexture.ref.value.destroy();
            geometry.groupTexture.ref.value.destroy();
            geometry.normalTexture.ref.value.destroy();
            geometry.doubleBuffer.destroy();
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

    const buffer = textureMesh?.doubleBuffer.get();
    const gv = extractIsosurface(ctx.webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.gridTexScale, densityTextureData.transform, isoLevel, true, buffer?.vertex, buffer?.group, buffer?.normal);

    const boundingSphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, props.radiusOffset + getStructureExtraRadius(structure));
    const surface = TextureMesh.create(gv.vertexCount, 1, gv.vertexTexture, gv.groupTexture, gv.normalTexture, boundingSphere, textureMesh);

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
        mustRecreate: (structure: Structure, props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return !props.tryUseGpu || !webgl || !suitableForGpu(structure, props, webgl);
        },
        dispose: (geometry: TextureMesh) => {
            geometry.vertexTexture.ref.value.destroy();
            geometry.groupTexture.ref.value.destroy();
            geometry.normalTexture.ref.value.destroy();
            geometry.doubleBuffer.destroy();
        }
    }, materialId);
}