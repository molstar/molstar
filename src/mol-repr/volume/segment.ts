/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { computeMarchingCubesMeshWebGPU } from '../../mol-gl/webgpu/marching-cubes';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Grid, Volume } from '../../mol-model/volume';
import { VisualContext } from '../visual';
import { Theme, ThemeRegistryContext } from '../../mol-theme/theme';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '../../mol-geo/util/marching-cubes/algorithm';
import { VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { VisualUpdateState } from '../util';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../representation';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Mat4, Tensor, Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { createSegmentTexture2d, createSegmentSampler, eachVolumeLoci, getVolumeTexture2dLayout } from './util';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { ValueCell } from '../../mol-util/value-cell';
import { extractIsosurface } from '../../mol-gl/compute/marching-cubes/isosurface';
import { Box3D } from '../../mol-math/geometry/primitives/box3d';
import { SortedArray } from '../../mol-data/int/sorted-array';
import { Interval } from '../../mol-data/int/interval';
import { OrderedSet } from '../../mol-data/int/ordered-set';
import { VolumeKey, VolumeVisual } from './visual';

export const VolumeSegmentParams = {
    segments: PD.Converted(
        (v: number[]) => v.map(x => `${x}`),
        (v: string[]) => v.map(x => parseInt(x)),
        PD.MultiSelect<string>(['0'], PD.arrayToOptions(['0']), {
            isEssential: true
        })
    )
};
export type VolumeSegmentParams = typeof VolumeSegmentParams
export type VolumeSegmentProps = PD.Values<VolumeSegmentParams>

function gpuSupport(webgl: WebGLContext): boolean {
    return !!(webgl.extensions.colorBufferFloat && webgl.extensions.textureFloat && webgl.extensions.drawBuffers);
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

const _translate = Mat4();
function getSegmentTransform(grid: Grid, segmentBox: Box3D) {
    const transform = Grid.getGridToCartesianTransform(grid);
    const translate = Mat4.fromTranslation(_translate, segmentBox.min);
    return Mat4.mul(Mat4(), transform, translate);
}

function useGpu(volume: Volume, props: PD.Values<SegmentMeshParams>, webgl?: WebGLContext): boolean {
    return props.tryUseGpu && !!webgl && gpuSupport(webgl) && suitableForGpu(volume, webgl);
}

export function SegmentVisual(materialId: number, volume: Volume, key: number, props: PD.Values<SegmentMeshParams>, webgl?: WebGLContext) {
    if (useGpu(volume, props, webgl)) {
        return SegmentTextureMeshVisual(materialId);
    }
    return SegmentMeshVisual(materialId);
}

function getLoci(volume: Volume, props: VolumeSegmentProps) {
    const segments = SortedArray.ofUnsortedArray<Volume.SegmentIndex>(props.segments);
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return Volume.Segment.Loci(volume, [{ segments, instances }]);
}

function getSegmentLoci(pickingId: PickingId, volume: Volume, key: number, props: VolumeSegmentProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;

    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        if (granularity === 'volume') {
            return Volume.Loci(volume, instances);
        } else if (granularity === 'object' || groupId === PickingId.Null) {
            const segments = OrderedSet.ofSingleton(key as Volume.SegmentIndex);
            return Volume.Segment.Loci(volume, [{ segments, instances }]);
        } else {
            const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
            return Volume.Cell.Loci(volume, [{ indices, instances }]);
        }
    }
    return EmptyLoci;
}

export function eachSegment(loci: Loci, volume: Volume, key: number, props: VolumeSegmentProps, apply: (interval: Interval) => boolean) {
    const segments = SortedArray.ofSingleton(key);
    return eachVolumeLoci(loci, volume, { segments }, apply);
}

//

function getSegmentCells(set: number[], bbox: Box3D, cells: Tensor): Tensor {
    const sample = createSegmentSampler(cells, set);

    const dim = Box3D.size(Vec3(), bbox);
    const [xn, yn, zn] = dim;

    const [minx, miny, minz] = bbox.min;

    const axisOrder = [...cells.space.axisOrderSlowToFast];
    const segmentSpace = Tensor.Space(dim, axisOrder, Uint8Array);
    const segmentCells = Tensor.create(segmentSpace, segmentSpace.create());

    const segData = segmentCells.data;
    const segSet = segmentSpace.set;

    for (let z = 0; z < zn; ++z) {
        for (let y = 0; y < yn; ++y) {
            for (let x = 0; x < xn; ++x) {
                segSet(segData, x, y, z, sample(x + minx, y + miny, z + minz));
            }
        }
    }

    return segmentCells;
}

export async function createVolumeSegmentMesh(ctx: VisualContext, volume: Volume, key: Volume.SegmentIndex, theme: Theme, props: VolumeSegmentProps & { tryUseGpu?: boolean }, mesh?: Mesh) {
    const segmentation = Volume.Segmentation.get(volume);
    if (!segmentation) throw new Error('missing volume segmentation');

    ctx.runtime.update({ message: 'Marching cubes...' });

    const bbox = Box3D.clone(segmentation.bounds[key]);
    Box3D.expand(bbox, bbox, Vec3.create(2, 2, 2));

    const set = Array.from(segmentation.segments.get(key)!.values());
    const cells = getSegmentCells(set, bbox, volume.grid.cells);
    // Segment geometry is cropped, but voxel loci refer to the source grid.
    const ids = new Int32Array(cells.data.length), coordinate = [0, 0, 0];
    for (let i = 0; i < ids.length; i++) {
        cells.space.getCoords(i, coordinate);
        for (let a = 0; a < 3; a++) coordinate[a] = Math.max(0, Math.min(volume.grid.cells.space.dimensions[a] - 1, coordinate[a] + bbox.min[a]));
        ids[i] = volume.grid.cells.space.dataOffset(...coordinate);
    }
    const params = { isoLevel: 128, scalarField: cells, idField: Tensor.create(cells.space, Tensor.Data1(ids)) };
    const useNative = ctx.webgpu && props.tryUseGpu !== false && cells.space.dimensions.every(d => d >= 2) && cells.data.length * 4 <= Math.min(ctx.webgpu.device.limits.maxStorageBufferBindingSize, ctx.webgpu.device.limits.maxBufferSize);
    const surface = useNative
        ? await computeMarchingCubesMeshWebGPU(ctx.runtime, ctx.webgpu!, params, mesh)
        : await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

    const transform = getSegmentTransform(volume.grid, bbox);
    Mesh.transform(surface, transform);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        // 2nd arg means not to split triangles based on group id. Splitting triangles
        // is too expensive if each cell has its own group id as is the case here.
        Mesh.uniformTriangleGroup(surface, false);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    const segments = OrderedSet.ofSingleton(key as Volume.SegmentIndex);
    surface.setBoundingSphere(Volume.Segment.getBoundingSphere(volume, [{ segments, instances }]));

    return surface;
}

export const SegmentMeshParams = {
    ...Mesh.Params,
    ...TextureMesh.Params,
    ...VolumeSegmentParams,
    quality: { ...Mesh.Params.quality, isEssential: false },
    tryUseGpu: PD.Boolean(true),
};
export type SegmentMeshParams = typeof SegmentMeshParams

export function SegmentMeshVisual(materialId: number): VolumeVisual<SegmentMeshParams> {
    return VolumeVisual<Mesh, SegmentMeshParams>({
        defaultProps: PD.getDefaultValues(SegmentMeshParams),
        createGeometry: createVolumeSegmentMesh,
        createLocationIterator: (volume: Volume, key: number) => {
            const l = Volume.Segment.Location(volume, key);
            return LocationIterator(volume.grid.cells.data.length, volume.instances ? volume.instances.length : 1, 1, () => l);
        },
        getLoci: getSegmentLoci,
        eachLocation: eachSegment,
        setUpdateState: (state: VisualUpdateState, newVolume: Volume, currentVolume: Volume, newProps: PD.Values<SegmentMeshParams>, currentProps: PD.Values<SegmentMeshParams>) => {
            state.createGeometry = newProps.tryUseGpu !== currentProps.tryUseGpu;
        },
        geometryUtils: Mesh.Utils,
        mustRecreate: (volumeKey: VolumeKey, props: PD.Values<SegmentMeshParams>, webgl?: WebGLContext) => {
            return useGpu(volumeKey.volume, props, webgl);
        }
    }, materialId);
}

//

const SegmentTextureName = 'segment-texture';

function getSegmentTexture(volume: Volume, segment: Volume.SegmentIndex, webgl: WebGLContext) {
    const segmentation = Volume.Segmentation.get(volume);
    if (!segmentation) throw new Error('missing volume segmentation');

    const { resources } = webgl;

    const bbox = Box3D.clone(segmentation.bounds[segment]);
    Box3D.expand(bbox, bbox, Vec3.create(2, 2, 2));

    const transform = getSegmentTransform(volume.grid, bbox);
    const gridDimension = Box3D.size(Vec3(), bbox);
    const { width, height, powerOfTwoSize: texDim } = getVolumeTexture2dLayout(gridDimension, Padding);
    const gridTexDim = Vec3.create(width, height, 0);
    const gridDataDim = Vec3.clone(gridDimension);
    const gridTexScale = Vec2.create(width / texDim, height / texDim);
    // console.log({ texDim, width, height, gridDimension });

    if (texDim > webgl.maxTextureSize / 2) {
        throw new Error('volume too large for gpu segment extraction');
    }

    if (!webgl.namedTextures[SegmentTextureName]) {
        webgl.namedTextures[SegmentTextureName] = resources.texture('image-uint8', 'alpha', 'ubyte', 'linear');
    }
    const texture = webgl.namedTextures[SegmentTextureName];

    texture.define(texDim, texDim);
    // load volume into sub-section of texture
    const set = Array.from(segmentation.segments.get(segment)!.values());
    texture.load(createSegmentTexture2d(volume, set, bbox, Padding), true);

    gridDimension[0] += Padding;
    gridDimension[1] += Padding;

    return {
        texture,
        transform,
        gridDimension,
        gridTexDim,
        gridDataDim,
        gridTexScale
    };
}

async function createVolumeSegmentTextureMesh(ctx: VisualContext, volume: Volume, segment: Volume.SegmentIndex, theme: Theme, props: VolumeSegmentProps, textureMesh?: TextureMesh) {
    if (!ctx.webgl) throw new Error('webgl context required to create volume segment texture-mesh');

    if (volume.grid.cells.data.length <= 1) {
        return TextureMesh.createEmpty(textureMesh);
    }

    const { texture, gridDimension, gridTexDim, gridDataDim, gridTexScale, transform } = getSegmentTexture(volume, segment, ctx.webgl);

    const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
    const buffer = textureMesh?.doubleBuffer.get();
    const gv = extractIsosurface(ctx.webgl, texture, gridDimension, gridTexDim, gridDataDim, gridTexScale, transform, 0.5, false, false, axisOrder, true, buffer?.vertex, buffer?.group, buffer?.normal);

    const groupCount = volume.grid.cells.data.length;
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    const segments = OrderedSet.ofSingleton(segment as Volume.SegmentIndex);
    const boundingSphere = Volume.Segment.getBoundingSphere(volume, [{ segments, instances }]);
    const surface = TextureMesh.create(gv.vertexCount, groupCount, gv.vertexTexture, gv.groupTexture, gv.normalTexture, boundingSphere, textureMesh);

    return surface;
}

export function SegmentTextureMeshVisual(materialId: number): VolumeVisual<SegmentMeshParams> {
    return VolumeVisual<TextureMesh, SegmentMeshParams>({
        defaultProps: PD.getDefaultValues(SegmentMeshParams),
        createGeometry: createVolumeSegmentTextureMesh,
        createLocationIterator: (volume: Volume, segment: number) => {
            const l = Volume.Segment.Location(volume, segment);
            return LocationIterator(volume.grid.cells.data.length, volume.instances ? volume.instances.length : 1, 1, () => l);
        },
        getLoci: getSegmentLoci,
        eachLocation: eachSegment,
        setUpdateState: (state: VisualUpdateState, newVolume: Volume, currentVolume: Volume, newProps: PD.Values<SegmentMeshParams>, currentProps: PD.Values<SegmentMeshParams>) => {
            state.createGeometry = newProps.tryUseGpu !== currentProps.tryUseGpu;
        },
        geometryUtils: TextureMesh.Utils,
        mustRecreate: (volumeKey: VolumeKey, props: PD.Values<SegmentMeshParams>, webgl?: WebGLContext) => {
            return !useGpu(volumeKey.volume, props, webgl);
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

function getSegments(props: VolumeSegmentProps, _volume: Volume): SortedArray {
    return SortedArray.ofUnsortedArray(props.segments);
}

const SegmentVisuals = {
    'segment': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, SegmentMeshParams>) => VolumeRepresentation('Segment mesh', ctx, getParams, SegmentVisual, getLoci, getSegments),
};

export const SegmentParams = {
    ...SegmentMeshParams,
    visuals: PD.MultiSelect(['segment'], PD.objectToOptions(SegmentVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type SegmentParams = typeof SegmentParams
export function getSegmentParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(SegmentParams);

    const segmentation = Volume.Segmentation.get(volume);
    if (segmentation) {
        const segments = Array.from(segmentation.segments.keys());
        p.segments = PD.Converted(
            (v: number[]) => v.map(x => `${x}`),
            (v: string[]) => v.map(x => parseInt(x)),
            PD.MultiSelect(segments.map(x => `${x}`), PD.arrayToOptions(segments.map(x => `${x}`)), {
                isEssential: true
            })
        );
    }
    return p;
}

export type SegmentRepresentation = VolumeRepresentation<SegmentParams>
export function SegmentRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, SegmentParams>): SegmentRepresentation {
    return Representation.createMulti('Segment', ctx, getParams, Representation.StateBuilder, SegmentVisuals as unknown as Representation.Def<Volume, SegmentParams>);
}

export const SegmentRepresentationProvider = VolumeRepresentationProvider({
    name: 'segment',
    label: 'Segment',
    description: 'Displays a triangulated segment of volumetric data.',
    factory: SegmentRepresentation,
    getParams: getSegmentParams,
    defaultValues: PD.getDefaultValues(SegmentParams),
    defaultColorTheme: { name: 'volume-segment' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !!Volume.Segmentation.get(volume)
});