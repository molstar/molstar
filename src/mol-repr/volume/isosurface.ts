/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Grid, Volume } from '../../mol-model/volume';
import { VisualContext } from '../visual';
import { Theme, ThemeRegistryContext } from '../../mol-theme/theme';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh, computeMarchingCubesLines } from '../../mol-geo/util/marching-cubes/algorithm';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider, VolumeKey } from './representation';
import { VisualUpdateState } from '../util';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../representation';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { Tensor, Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { fillSerial } from '../../mol-util/array';
import { createVolumeCellLocationIterator, createVolumeTexture2d, eachVolumeLoci, getVolumeTexture2dLayout } from './util';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { extractIsosurface } from '../../mol-gl/compute/marching-cubes/isosurface';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { Texture } from '../../mol-gl/webgl/texture';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { ValueCell } from '../../mol-util/value-cell';

export const VolumeIsosurfaceParams = {
    isoValue: Volume.IsoValueParam,
};
export type VolumeIsosurfaceParams = typeof VolumeIsosurfaceParams
export type VolumeIsosurfaceProps = PD.Values<VolumeIsosurfaceParams>

export const VolumeIsosurfaceTextureParams = {
    isoValue: Volume.IsoValueParam,
    tryUseGpu: PD.Boolean(true),
    gpuDataType: PD.Select('byte', PD.arrayToOptions(['byte', 'float', 'halfFloat'] as const), { hideIf: p => !p.tryUseGpu }),
};
export type VolumeIsosurfaceGpuParams = typeof VolumeIsosurfaceTextureParams
export type VolumeIsosurfaceGpuProps = PD.Values<VolumeIsosurfaceGpuParams>

function gpuSupport(webgl: WebGLContext) {
    return webgl.extensions.colorBufferFloat && webgl.extensions.textureFloat && webgl.extensions.drawBuffers;
}

const Padding = 1;

function suitableForGpu(volume: Volume, webgl: WebGLContext) {
    // small volumes are about as fast or faster on CPU vs integrated GPU
    if (volume.grid.cells.data.length < Math.pow(10, 3)) return false;
    // the GPU is much more memory contraint, especially true for integrated GPUs,
    // fallback to CPU for large volumes
    const gridDim = volume.grid.cells.space.dimensions as Vec3;
    const { powerOfTwoSize } = getVolumeTexture2dLayout(gridDim, Padding);
    return powerOfTwoSize <= webgl.maxTextureSize / 2;
}

export function IsosurfaceVisual(materialId: number, volume: Volume, key: number, props: PD.Values<IsosurfaceMeshParams>, webgl?: WebGLContext) {
    if (props.tryUseGpu && webgl && gpuSupport(webgl) && suitableForGpu(volume, webgl)) {
        return IsosurfaceTextureMeshVisual(materialId);
    }
    return IsosurfaceMeshVisual(materialId);
}

function getLoci(volume: Volume, props: VolumeIsosurfaceProps) {
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return Volume.Isosurface.Loci(volume, props.isoValue, instances);
}

function getIsosurfaceLoci(pickingId: PickingId, volume: Volume, key: number, props: VolumeIsosurfaceProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;
    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        if (granularity === 'volume') {
            return Volume.Loci(volume, instances);
        } else if (granularity === 'object') {
            return Volume.Isosurface.Loci(volume, props.isoValue, instances);
        } else {
            const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
            return Volume.Cell.Loci(volume, [{ indices, instances }]);
        }
    }
    return EmptyLoci;
}

export function eachIsosurface(loci: Loci, volume: Volume, key: number, props: VolumeIsosurfaceProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, { isoValue: props.isoValue }, apply);
}

//

export async function createVolumeIsosurfaceMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeIsosurfaceProps, mesh?: Mesh) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const ids = fillSerial(new Int32Array(volume.grid.cells.data.length));

    const surface = await computeMarchingCubesMesh({
        isoLevel: Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue,
        scalarField: volume.grid.cells,
        idField: Tensor.create(volume.grid.cells.space, Tensor.Data1(ids))
    }, mesh).runAsChild(ctx.runtime);

    const transform = Grid.getGridToCartesianTransform(volume.grid);
    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        // 2nd arg means not to split triangles based on group id. Splitting triangles
        // is too expensive if each cell has its own group id as is the case here.
        Mesh.uniformTriangleGroup(surface, false);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    surface.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));

    return surface;
}

export const IsosurfaceMeshParams = {
    ...Mesh.Params,
    ...TextureMesh.Params,
    ...VolumeIsosurfaceParams,
    ...VolumeIsosurfaceTextureParams,
    quality: { ...Mesh.Params.quality, isEssential: false },
};
export type IsosurfaceMeshParams = typeof IsosurfaceMeshParams

export function IsosurfaceMeshVisual(materialId: number): VolumeVisual<IsosurfaceMeshParams> {
    return VolumeVisual<Mesh, IsosurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceMeshParams),
        createGeometry: createVolumeIsosurfaceMesh,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getIsosurfaceLoci,
        eachLocation: eachIsosurface,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<IsosurfaceMeshParams>, currentProps: PD.Values<IsosurfaceMeshParams>) => {
            if (!Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats)) state.createGeometry = true;
        },
        geometryUtils: Mesh.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<IsosurfaceMeshParams>, webgl?: WebGLContext) => {
            return props.tryUseGpu && !!webgl && suitableForGpu(volumekey.volume, webgl);
        }
    }, materialId);
}

//

namespace VolumeIsosurfaceTexture {
    const name = 'volume-isosurface-texture';
    export const descriptor = CustomPropertyDescriptor({ name });
    export function clear(volume: Volume) {
        delete volume._propertyData[name];
    }
    export function get(volume: Volume, webgl: WebGLContext, props: VolumeIsosurfaceGpuProps) {
        const transform = Grid.getGridToCartesianTransform(volume.grid);
        const gridDimension = Vec3.clone(volume.grid.cells.space.dimensions as Vec3);
        const { width, height, powerOfTwoSize: texDim } = getVolumeTexture2dLayout(gridDimension, Padding);
        const gridTexDim = Vec3.create(width, height, 0);
        const gridTexScale = Vec2.create(width / texDim, height / texDim);
        // console.log({ texDim, width, height, gridDimension });

        if (texDim > webgl.maxTextureSize / 2) {
            throw new Error('volume too large for gpu isosurface extraction');
        }

        const dataType = props.gpuDataType === 'halfFloat' && !webgl.extensions.textureHalfFloat ? 'float' : props.gpuDataType;

        if (volume._propertyData[name]?.dataType !== dataType) {
            const texture = dataType === 'byte'
                ? webgl.resources.texture('image-uint8', 'alpha', 'ubyte', 'linear')
                : dataType === 'halfFloat'
                    ? webgl.resources.texture('image-float16', 'alpha', 'fp16', 'linear')
                    : webgl.resources.texture('image-float32', 'alpha', 'float', 'linear');
            volume._propertyData[name] = { texture, dataType };
            texture.define(texDim, texDim);
            // load volume into sub-section of texture
            texture.load(createVolumeTexture2d(volume, 'data', Padding, dataType), true);
            volume.customProperties.add(descriptor);
            volume.customProperties.assets(descriptor, [{ dispose: () => texture.destroy() }]);
        }

        gridDimension[0] += Padding;
        gridDimension[1] += Padding;

        return {
            texture: volume._propertyData[name].texture as Texture,
            transform,
            gridDimension,
            gridTexDim,
            gridTexScale
        };
    }
}

function createVolumeIsosurfaceTextureMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeIsosurfaceGpuProps, textureMesh?: TextureMesh) {
    const { webgl } = ctx;
    if (!webgl) throw new Error('webgl context required to create volume isosurface texture-mesh');

    if (volume.grid.cells.data.length <= 1) {
        return TextureMesh.createEmpty(textureMesh);
    }

    const { max, min } = volume.grid.stats;
    const diff = max - min;
    const value = Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue;
    const isoLevel = ((value - min) / diff);

    const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
    const groupCount = volume.grid.cells.data.length;
    const boundingSphere = Volume.getBoundingSphere(volume); // getting isosurface bounding-sphere is too expensive here

    const create = (textureMesh?: TextureMesh) => {
        const { texture, gridDimension, gridTexDim, gridTexScale, transform } = VolumeIsosurfaceTexture.get(volume, webgl, props);

        const buffer = textureMesh?.doubleBuffer.get();
        const gv = extractIsosurface(webgl, texture, gridDimension, gridTexDim, gridTexScale, transform, isoLevel, value < 0, false, axisOrder, true, buffer?.vertex, buffer?.group, buffer?.normal);

        return TextureMesh.create(gv.vertexCount, groupCount, gv.vertexTexture, gv.groupTexture, gv.normalTexture, boundingSphere, textureMesh);
    };

    const surface = create(textureMesh);
    surface.meta.webgl = webgl;
    surface.meta.reset = () => {
        VolumeIsosurfaceTexture.clear(volume);
        create(surface);
    };

    return surface;
}

export function IsosurfaceTextureMeshVisual(materialId: number): VolumeVisual<IsosurfaceMeshParams> {
    return VolumeVisual<TextureMesh, IsosurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceMeshParams),
        createGeometry: createVolumeIsosurfaceTextureMesh,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getIsosurfaceLoci,
        eachLocation: eachIsosurface,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<IsosurfaceMeshParams>, currentProps: PD.Values<IsosurfaceMeshParams>) => {
            if (!Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats)) state.createGeometry = true;
            if (newProps.gpuDataType !== currentProps.gpuDataType) state.createGeometry = true;
        },
        geometryUtils: TextureMesh.Utils,
        mustRecreate: (volumeKey: VolumeKey, props: PD.Values<IsosurfaceMeshParams>, webgl?: WebGLContext) => {
            return !props.tryUseGpu || !webgl || !suitableForGpu(volumeKey.volume, webgl);
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

export async function createVolumeIsosurfaceWireframe(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeIsosurfaceProps, lines?: Lines) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const ids = fillSerial(new Int32Array(volume.grid.cells.data.length));

    const wireframe = await computeMarchingCubesLines({
        isoLevel: Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue,
        scalarField: volume.grid.cells,
        idField: Tensor.create(volume.grid.cells.space, Tensor.Data1(ids))
    }, lines).runAsChild(ctx.runtime);

    const transform = Grid.getGridToCartesianTransform(volume.grid);
    Lines.transform(wireframe, transform);

    wireframe.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));

    return wireframe;
}

export const IsosurfaceWireframeParams = {
    ...Lines.Params,
    ...VolumeIsosurfaceParams,
    quality: { ...Lines.Params.quality, isEssential: false },
    sizeFactor: PD.Numeric(3, { min: 0, max: 10, step: 0.1 }),
};
export type IsosurfaceWireframeParams = typeof IsosurfaceWireframeParams

export function IsosurfaceWireframeVisual(materialId: number): VolumeVisual<IsosurfaceWireframeParams> {
    return VolumeVisual<Lines, IsosurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceWireframeParams),
        createGeometry: createVolumeIsosurfaceWireframe,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getIsosurfaceLoci,
        eachLocation: eachIsosurface,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<IsosurfaceWireframeParams>, currentProps: PD.Values<IsosurfaceWireframeParams>) => {
            if (!Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats)) state.createGeometry = true;
        },
        geometryUtils: Lines.Utils
    }, materialId);
}

//

const IsosurfaceVisuals = {
    'solid': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceMeshParams>) => VolumeRepresentation('Isosurface mesh', ctx, getParams, IsosurfaceVisual, getLoci),
    'wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceWireframeParams>) => VolumeRepresentation('Isosurface wireframe', ctx, getParams, IsosurfaceWireframeVisual, getLoci),
};

export const IsosurfaceParams = {
    ...IsosurfaceMeshParams,
    ...IsosurfaceWireframeParams,
    visuals: PD.MultiSelect(['solid'], PD.objectToOptions(IsosurfaceVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type IsosurfaceParams = typeof IsosurfaceParams
export function getIsosurfaceParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(IsosurfaceParams);
    p.isoValue = Volume.createIsoValueParam(Volume.IsoValue.relative(2), volume.grid.stats);
    return p;
}

export type IsosurfaceRepresentation = VolumeRepresentation<IsosurfaceParams>
export function IsosurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceParams>): IsosurfaceRepresentation {
    return Representation.createMulti('Isosurface', ctx, getParams, Representation.StateBuilder, IsosurfaceVisuals as unknown as Representation.Def<Volume, IsosurfaceParams>);
}

export const IsosurfaceRepresentationProvider = VolumeRepresentationProvider({
    name: 'isosurface',
    label: 'Isosurface',
    description: 'Displays a triangulated isosurface of volumetric data.',
    factory: IsosurfaceRepresentation,
    getParams: getIsosurfaceParams,
    defaultValues: PD.getDefaultValues(IsosurfaceParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume)
});