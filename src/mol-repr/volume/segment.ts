/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Grid, Volume } from '../../mol-model/volume';
import { VisualContext } from '../visual';
import { Theme, ThemeRegistryContext } from '../../mol-theme/theme';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '../../mol-geo/util/marching-cubes/algorithm';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider, VolumeKey } from './representation';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { VisualUpdateState } from '../util';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../representation';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, SortedArray } from '../../mol-data/int';
import { Mat4, Tensor, Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { fillSerial } from '../../mol-util/array';
import { createSegmentTexture2d, eachVolumeLoci, getVolumeTexture2dLayout } from './util';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { ValueCell } from '../../mol-util/value-cell';
import { extractIsosurface } from '../../mol-gl/compute/marching-cubes/isosurface';
import { Box3D } from '../../mol-math/geometry/primitives/box3d';

export const VolumeSegmentParams = {
    segments: PD.Converted(
        (v: number[]) => v.map(x => `${x}`),
        (v: string[]) => v.map(x => parseInt(x)),
        PD.MultiSelect(['0'], PD.arrayToOptions(['0']), {
            isEssential: true
        })
    )
};
export type VolumeSegmentParams = typeof VolumeSegmentParams
export type VolumeSegmentProps = PD.Values<VolumeSegmentParams>

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

const _translate = Mat4();
function getSegmentTransform(grid: Grid, segmentBox: Box3D) {
    const transform = Grid.getGridToCartesianTransform(grid);
    const translate = Mat4.fromTranslation(_translate, segmentBox.min);
    return Mat4.mul(Mat4(), transform, translate);
}

export function SegmentVisual(materialId: number, volume: Volume, key: number, props: PD.Values<SegmentMeshParams>, webgl?: WebGLContext) {
    if (props.tryUseGpu && webgl && gpuSupport(webgl) && suitableForGpu(volume, webgl)) {
        return SegmentTextureMeshVisual(materialId);
    }
    return SegmentMeshVisual(materialId);
}

function getLoci(volume: Volume, props: VolumeSegmentProps) {
    return Volume.Segment.Loci(volume, props.segments);
}

function getSegmentLoci(pickingId: PickingId, volume: Volume, key: number, props: VolumeSegmentProps, id: number) {
    const { objectId, groupId } = pickingId;

    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        if (granularity === 'volume') {
            return Volume.Loci(volume);
        } else if (granularity === 'object') {
            return Volume.Segment.Loci(volume, [key]);
        } else {
            return Volume.Cell.Loci(volume, Interval.ofSingleton(groupId as Volume.CellIndex));
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
    const data = cells.data;
    const o = cells.space.dataOffset;

    const dim = Box3D.size(Vec3(), bbox);
    const [xn, yn, zn] = dim;
    const xn1 = xn - 1;
    const yn1 = yn - 1;
    const zn1 = zn - 1;

    const [minx, miny, minz] = bbox.min;
    const [maxx, maxy, maxz] = bbox.max;

    const axisOrder = [...cells.space.axisOrderSlowToFast];
    const segmentSpace = Tensor.Space(dim, axisOrder, Uint8Array);
    const segmentCells = Tensor.create(segmentSpace, segmentSpace.create());

    const segData = segmentCells.data;
    const segSet = segmentSpace.set;

    for (let z = 0; z < zn; ++z) {
        for (let y = 0; y < yn; ++y) {
            for (let x = 0; x < xn; ++x) {
                const v0 = set.includes(data[o(x + minx, y + miny, z + minz)]) ? 255 : 0;
                const xp = set.includes(data[o(Math.min(xn1 + maxx, x + 1 + minx), y + miny, z + minz)]) ? 255 : 0;
                const xn = set.includes(data[o(Math.max(0, x - 1 + minx), y + miny, z + minz)]) ? 255 : 0;
                const yp = set.includes(data[o(x + minx, Math.min(yn1 + maxy, y + 1 + miny), z + minz)]) ? 255 : 0;
                const yn = set.includes(data[o(x + minx, Math.max(0, y - 1 + miny), z + minz)]) ? 255 : 0;
                const zp = set.includes(data[o(x + minx, y + miny, Math.min(zn1 + maxz, z + 1 + minz))]) ? 255 : 0;
                const zn = set.includes(data[o(x + minx, y + miny, Math.max(0, z - 1 + minz))]) ? 255 : 0;

                segSet(segData, x, y, z, Math.round((v0 + v0 + xp + xn + yp + yn + zp + zn) / 8));
            }
        }
    }

    return segmentCells;
}

export async function createVolumeSegmentMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeSegmentProps, mesh?: Mesh) {
    const segmentation = Volume.Segmentation.get(volume);
    if (!segmentation) throw new Error('missing volume segmentation');

    ctx.runtime.update({ message: 'Marching cubes...' });

    const bbox = Box3D.clone(segmentation.bounds[key]);
    Box3D.expand(bbox, bbox, Vec3.create(2, 2, 2));

    const set = Array.from(segmentation.segments.get(key)!.values());
    const cells = getSegmentCells(set, bbox, volume.grid.cells);
    const ids = fillSerial(new Int32Array(cells.data.length));

    const surface = await computeMarchingCubesMesh({
        isoLevel: 128,
        scalarField: cells,
        idField: Tensor.create(cells.space, Tensor.Data1(ids))
    }, mesh).runAsChild(ctx.runtime);

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

    surface.setBoundingSphere(Volume.Segment.getBoundingSphere(volume, [key]));

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
            return LocationIterator(volume.grid.cells.data.length, 1, 1, () => l);
        },
        getLoci: getSegmentLoci,
        eachLocation: eachSegment,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<SegmentMeshParams>, currentProps: PD.Values<SegmentMeshParams>) => {
        },
        geometryUtils: Mesh.Utils,
        mustRecreate: (volumeKey: VolumeKey, props: PD.Values<SegmentMeshParams>, webgl?: WebGLContext) => {
            return props.tryUseGpu && !!webgl && suitableForGpu(volumeKey.volume, webgl);
        }
    }, materialId);
}

//

const SegmentTextureName = 'segment-texture';

function getSegmentTexture(volume: Volume, segment: number, webgl: WebGLContext) {
    const segmentation = Volume.Segmentation.get(volume);
    if (!segmentation) throw new Error('missing volume segmentation');

    const { resources } = webgl;

    const bbox = Box3D.clone(segmentation.bounds[segment]);
    Box3D.expand(bbox, bbox, Vec3.create(2, 2, 2));

    const transform = getSegmentTransform(volume.grid, bbox);
    const gridDimension = Box3D.size(Vec3(), bbox);
    const { width, height, powerOfTwoSize: texDim } = getVolumeTexture2dLayout(gridDimension, Padding);
    const gridTexDim = Vec3.create(width, height, 0);
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
        gridTexScale
    };
}

async function createVolumeSegmentTextureMesh(ctx: VisualContext, volume: Volume, segment: number, theme: Theme, props: VolumeSegmentProps, textureMesh?: TextureMesh) {
    if (!ctx.webgl) throw new Error('webgl context required to create volume segment texture-mesh');

    if (volume.grid.cells.data.length <= 1) {
        return TextureMesh.createEmpty(textureMesh);
    }

    const { texture, gridDimension, gridTexDim, gridTexScale, transform } = getSegmentTexture(volume, segment, ctx.webgl);

    const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
    const buffer = textureMesh?.doubleBuffer.get();
    const gv = extractIsosurface(ctx.webgl, texture, gridDimension, gridTexDim, gridTexScale, transform, 0.5, false, false, axisOrder, true, buffer?.vertex, buffer?.group, buffer?.normal);

    const groupCount = volume.grid.cells.data.length;
    const surface = TextureMesh.create(gv.vertexCount, groupCount, gv.vertexTexture, gv.groupTexture, gv.normalTexture, Volume.Segment.getBoundingSphere(volume, [segment]), textureMesh);

    return surface;
}

export function SegmentTextureMeshVisual(materialId: number): VolumeVisual<SegmentMeshParams> {
    return VolumeVisual<TextureMesh, SegmentMeshParams>({
        defaultProps: PD.getDefaultValues(SegmentMeshParams),
        createGeometry: createVolumeSegmentTextureMesh,
        createLocationIterator: (volume: Volume, segment: number) => {
            const l = Volume.Segment.Location(volume, segment);
            return LocationIterator(volume.grid.cells.data.length, 1, 1, () => l);
        },
        getLoci: getSegmentLoci,
        eachLocation: eachSegment,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<SegmentMeshParams>, currentProps: PD.Values<SegmentMeshParams>) => {
        },
        geometryUtils: TextureMesh.Utils,
        mustRecreate: (volumeKey: VolumeKey, props: PD.Values<SegmentMeshParams>, webgl?: WebGLContext) => {
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

function getSegments(props: VolumeSegmentProps): SortedArray {
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