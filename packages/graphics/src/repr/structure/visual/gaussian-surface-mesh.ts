/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { UnitsMeshParams, UnitsTextureMeshParams, type UnitsVisual, UnitsMeshVisual, UnitsTextureMeshVisual } from '../units-visual.js';
import { GaussianDensityParams, computeUnitGaussianDensity, computeUnitGaussianDensityTexture2d, type GaussianDensityProps, computeStructureGaussianDensity, computeStructureGaussianDensityTexture2d } from './util/gaussian.js';
import type { VisualContext } from '../../visual.js';
import type { Unit, Structure } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element.js';
import type { VisualUpdateState } from '../../util.js';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import { extractIsosurface } from '@molstar/graphics/gl/compute/marching-cubes/isosurface';
import { Sphere3D } from '@molstar/core/math/geometry';
import { type ComplexVisual, ComplexMeshParams, ComplexMeshVisual, ComplexTextureMeshVisual, ComplexTextureMeshParams } from '../complex-visual.js';
import { getVolumeSliceInfo, type StructureGroup } from './util/common.js';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import type { TextureMeshValues } from '@molstar/graphics/gl/renderable/texture-mesh';
import type { Texture } from '@molstar/graphics/gl/webgl/texture';
import { applyMeshColorSmoothing } from '@molstar/graphics/geo/geometry/mesh/color-smoothing';
import { applyTextureMeshColorSmoothing } from '@molstar/graphics/geo/geometry/texture-mesh/color-smoothing';
import { ColorSmoothingParams, getColorSmoothingProps } from '@molstar/graphics/geo/geometry/base';
import { Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { isTimingMode } from '@molstar/core/util/debug';
import { ValueCell } from '@molstar/core/util/value-cell';

const SharedParams = {
    ...GaussianDensityParams,
    ...ColorSmoothingParams,
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

function gpuSupport(webgl: WebGLContext): boolean {
    return !!(webgl.extensions.colorBufferFloat && webgl.extensions.textureFloat && webgl.extensions.textureFloatLinear && webgl.extensions.blendMinMax && webgl.extensions.drawBuffers);
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

function useGpu(structure: Structure, props: PD.Values<SharedParams>, webgl?: WebGLContext): boolean {
    return !props.includeParent && props.floodfill === 'off' && props.tryUseGpu && !!webgl && gpuSupport(webgl) && suitableForGpu(structure, props, webgl);
}

export function GaussianSurfaceVisual(materialId: number, structure: Structure, props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) {
    if (useGpu(structure, props, webgl)) {
        return GaussianSurfaceTextureMeshVisual(materialId);
    }
    return GaussianSurfaceMeshVisual(materialId);
}

export function StructureGaussianSurfaceVisual(materialId: number, structure: Structure, props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) {
    if (useGpu(structure, props, webgl)) {
        return StructureGaussianSurfaceTextureMeshVisual(materialId);
    }
    return StructureGaussianSurfaceMeshVisual(materialId);
}

type GaussianSurfaceMeta = {
    resolution?: number
    colorTexture?: Texture
}

//

async function createGaussianSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness, floodfill, includeParent, radiusOffset } = props;
    const { transform, field, idField, radiusFactor, resolution, maxRadius } = await computeUnitGaussianDensity(structure, unit, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const params = {
        isoLevel,
        scalarField: floodfill !== 'off' ? Tensor.createFloodfilled(field, isoLevel, floodfill) : field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);
    (surface.meta.resolution as GaussianSurfaceMeta['resolution']) = resolution;

    if (includeParent) {
        const iterations = Math.ceil(2 / resolution);
        Mesh.smoothEdges(surface, { iterations, maxNewEdgeLength: Math.sqrt(2) });
    }

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius + extraRadius);
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
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.radiusOffset !== currentProps.radiusOffset ||
                newProps.smoothness !== currentProps.smoothness ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return useGpu(structureGroup.structure, props, webgl);
        },
        processValues: (values: MeshValues, geometry: Mesh, props: PD.Values<GaussianSurfaceMeshParams>, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta.colorTexture as GaussianSurfaceMeta['colorTexture']) = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
        }
    }, materialId);
}

//

async function createStructureGaussianSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, mesh?: Mesh): Promise<Mesh> {
    const { smoothness, floodfill, includeParent, radiusOffset } = props;
    const { transform, field, idField, radiusFactor, resolution, maxRadius } = await computeStructureGaussianDensity(structure, theme.size, props).runInContext(ctx.runtime);

    const isoLevel = Math.exp(-smoothness) / radiusFactor;
    const params = {
        isoLevel,
        scalarField: floodfill !== 'off' ? Tensor.createFloodfilled(field, isoLevel, floodfill) : field,
        idField
    };
    const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);
    (surface.meta.resolution as GaussianSurfaceMeta['resolution']) = resolution;

    if (includeParent) {
        const iterations = Math.ceil(2 / resolution);
        Mesh.smoothEdges(surface, { iterations, maxNewEdgeLength: Math.sqrt(2) });
    }

    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius + extraRadius);
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
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.radiusOffset !== currentProps.radiusOffset ||
                newProps.smoothness !== currentProps.smoothness ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        mustRecreate: (structure: Structure, props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return useGpu(structure, props, webgl);
        },
        processValues: (values: MeshValues, geometry: Mesh, props: PD.Values<GaussianSurfaceMeshParams>, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta.colorTexture as GaussianSurfaceMeta['colorTexture']) = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
        }
    }, materialId);
}

//

const GaussianSurfaceName = 'gaussian-surface';

function createGaussianSurfaceTextureMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, textureMesh?: TextureMesh): TextureMesh {
    const { webgl } = ctx;
    if (!webgl) throw new Error('webgl context required to create gaussian surface texture-mesh');

    const { namedTextures, resources, extensions: { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } } = webgl;
    if (!namedTextures[GaussianSurfaceName]) {
        namedTextures[GaussianSurfaceName] = colorBufferHalfFloat && textureHalfFloat
            ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
            : colorBufferFloat && textureFloat
                ? resources.texture('image-float32', 'rgba', 'float', 'linear')
                : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    }

    const axisOrder = Vec3.create(0, 1, 2);
    const groupCount = unit.elements.length;

    const create = (textureMesh?: TextureMesh) => {
        if (isTimingMode) webgl.timer.mark('createGaussianSurfaceTextureMesh');
        const densityTextureData = computeUnitGaussianDensityTexture2d(structure, unit, theme.size, true, props, webgl, namedTextures[GaussianSurfaceName]);
        const isoLevel = Math.exp(-props.smoothness) / densityTextureData.radiusFactor;

        const buffer = textureMesh?.doubleBuffer.get();
        const gv = extractIsosurface(webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.gridDataDim, densityTextureData.gridTexScale, densityTextureData.transform, isoLevel, false, true, axisOrder, true, buffer?.vertex, buffer?.group, buffer?.normal);
        if (isTimingMode) webgl.timer.markEnd('createGaussianSurfaceTextureMesh');

        const extraRadius = props.radiusOffset * (1 + Math.exp(-props.smoothness));
        const boundingSphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, densityTextureData.maxRadius + extraRadius);
        const surface = TextureMesh.create(gv.vertexCount, groupCount, gv.vertexTexture, gv.groupTexture, gv.normalTexture, boundingSphere, textureMesh);
        (surface.meta as GaussianSurfaceMeta).resolution = densityTextureData.resolution;
        return surface;
    };

    const surface = create(textureMesh);
    surface.meta.webgl = webgl;
    surface.meta.reset = () => {
        create(surface);
    };

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
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.radiusOffset !== currentProps.radiusOffset ||
                newProps.smoothness !== currentProps.smoothness ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<GaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return !useGpu(structureGroup.structure, props, webgl);
        },
        processValues: (values: TextureMeshValues, geometry: TextureMesh, props: PD.Values<GaussianSurfaceMeshParams>, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp && webgl) {
                applyTextureMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta as GaussianSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: TextureMesh) => {
            geometry.vertexTexture.ref.value.destroy();
            geometry.groupTexture.ref.value.destroy();
            geometry.normalTexture.ref.value.destroy();
            geometry.doubleBuffer.destroy();

            (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
        }
    }, materialId);
}

//

function createStructureGaussianSurfaceTextureMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, textureMesh?: TextureMesh): TextureMesh {
    const { webgl } = ctx;
    if (!webgl) throw new Error('webgl context required to create structure gaussian surface texture-mesh');

    const { namedTextures, resources, extensions: { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } } = webgl;
    if (!namedTextures[GaussianSurfaceName]) {
        namedTextures[GaussianSurfaceName] = colorBufferHalfFloat && textureHalfFloat
            ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
            : colorBufferFloat && textureFloat
                ? resources.texture('image-float32', 'rgba', 'float', 'linear')
                : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    }

    const axisOrder = Vec3.create(0, 1, 2);
    const groupCount = structure.elementCount;

    const create = (textureMesh?: TextureMesh) => {
        if (isTimingMode) webgl.timer.mark('createStructureGaussianSurfaceTextureMesh');
        const densityTextureData = computeStructureGaussianDensityTexture2d(structure, theme.size, true, props, webgl, namedTextures[GaussianSurfaceName]);
        const isoLevel = Math.exp(-props.smoothness) / densityTextureData.radiusFactor;

        const buffer = textureMesh?.doubleBuffer.get();
        const gv = extractIsosurface(webgl, densityTextureData.texture, densityTextureData.gridDim, densityTextureData.gridTexDim, densityTextureData.gridDataDim, densityTextureData.gridTexScale, densityTextureData.transform, isoLevel, false, true, axisOrder, true, buffer?.vertex, buffer?.group, buffer?.normal);
        if (isTimingMode) webgl.timer.markEnd('createStructureGaussianSurfaceTextureMesh');

        const extraRadius = props.radiusOffset * (1 + Math.exp(-props.smoothness));
        const boundingSphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, densityTextureData.maxRadius + extraRadius);
        const surface = TextureMesh.create(gv.vertexCount, groupCount, gv.vertexTexture, gv.groupTexture, gv.normalTexture, boundingSphere, textureMesh);
        (surface.meta as GaussianSurfaceMeta).resolution = densityTextureData.resolution;
        return surface;
    };

    const surface = create(textureMesh);
    surface.meta.webgl = webgl;
    surface.meta.reset = () => {
        create(surface);
    };

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
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.radiusOffset !== currentProps.radiusOffset ||
                newProps.smoothness !== currentProps.smoothness ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        mustRecreate: (structure: Structure, props: PD.Values<StructureGaussianSurfaceMeshParams>, webgl?: WebGLContext) => {
            return !useGpu(structure, props, webgl);
        },
        processValues: (values: TextureMeshValues, geometry: TextureMesh, props: PD.Values<GaussianSurfaceMeshParams>, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp && webgl) {
                applyTextureMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta as GaussianSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: TextureMesh) => {
            geometry.vertexTexture.ref.value.destroy();
            geometry.groupTexture.ref.value.destroy();
            geometry.normalTexture.ref.value.destroy();
            geometry.doubleBuffer.destroy();

            (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
        }
    }, materialId);
}