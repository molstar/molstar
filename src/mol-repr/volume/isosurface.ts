/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { NullLocation } from '../../mol-model/location';
import { VisualUpdateState } from '../util';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../representation';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval } from '../../mol-data/int';
import { Tensor } from '../../mol-math/linear-algebra';
import { fillSerial } from '../../mol-util/array';
import { eachVolumeLoci } from './util';

export const VolumeIsosurfaceParams = {
    isoValue: Volume.IsoValueParam
};
export type VolumeIsosurfaceParams = typeof VolumeIsosurfaceParams
export type VolumeIsosurfaceProps = PD.Values<VolumeIsosurfaceParams>

function getLoci(volume: Volume, props: VolumeIsosurfaceProps) {
    return Volume.Isosurface.Loci(volume, props.isoValue);
}

function getIsosurfaceLoci(pickingId: PickingId, volume: Volume, props: VolumeIsosurfaceProps, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        return Volume.Cell.Loci(volume, Interval.ofSingleton(groupId as Volume.CellIndex));
    }
    return EmptyLoci;
}

export function eachIsosurface(loci: Loci, volume: Volume, props: VolumeIsosurfaceProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, props.isoValue, apply);
}

//

export async function createVolumeIsosurfaceMesh(ctx: VisualContext, volume: Volume, theme: Theme, props: VolumeIsosurfaceProps, mesh?: Mesh) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const ids = fillSerial(new Int32Array(volume.grid.cells.data.length));

    const surface = await computeMarchingCubesMesh({
        isoLevel: Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue,
        scalarField: volume.grid.cells,
        idField: Tensor.create(volume.grid.cells.space, Tensor.Data1(ids))
    }, mesh).runAsChild(ctx.runtime);

    const transform = Grid.getGridToCartesianTransform(volume.grid);
    Mesh.transform(surface, transform);

    surface.setBoundingSphere(Volume.getBoundingSphere(volume));

    return surface;
}

export const IsosurfaceMeshParams = {
    ...Mesh.Params,
    quality: { ...Mesh.Params.quality, isEssential: false },
    ...VolumeIsosurfaceParams
};
export type IsosurfaceMeshParams = typeof IsosurfaceMeshParams

export function IsosurfaceMeshVisual(materialId: number): VolumeVisual<IsosurfaceMeshParams> {
    return VolumeVisual<Mesh, IsosurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceMeshParams),
        createGeometry: createVolumeIsosurfaceMesh,
        createLocationIterator: (volume: Volume) => LocationIterator(volume.grid.cells.data.length, 1, () => NullLocation),
        getLoci: getIsosurfaceLoci,
        eachLocation: eachIsosurface,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<IsosurfaceMeshParams>, currentProps: PD.Values<IsosurfaceMeshParams>) => {
            if (!Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats)) state.createGeometry = true;
        },
        geometryUtils: Mesh.Utils
    }, materialId);
}

//

export async function createVolumeIsosurfaceWireframe(ctx: VisualContext, volume: Volume, theme: Theme, props: VolumeIsosurfaceProps, lines?: Lines) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const ids = fillSerial(new Int32Array(volume.grid.cells.data.length));

    const wireframe = await computeMarchingCubesLines({
        isoLevel: Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue,
        scalarField: volume.grid.cells,
        idField: Tensor.create(volume.grid.cells.space, Tensor.Data1(ids))
    }, lines).runAsChild(ctx.runtime);

    const transform = Grid.getGridToCartesianTransform(volume.grid);
    Lines.transform(wireframe, transform);

    wireframe.setBoundingSphere(Volume.getBoundingSphere(volume));

    return wireframe;
}

export const IsosurfaceWireframeParams = {
    ...Lines.Params,
    quality: { ...Lines.Params.quality, isEssential: false },
    sizeFactor: PD.Numeric(1.5, { min: 0, max: 10, step: 0.1 }),
    ...VolumeIsosurfaceParams
};
export type IsosurfaceWireframeParams = typeof IsosurfaceWireframeParams

export function IsosurfaceWireframeVisual(materialId: number): VolumeVisual<IsosurfaceWireframeParams> {
    return VolumeVisual<Lines, IsosurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceWireframeParams),
        createGeometry: createVolumeIsosurfaceWireframe,
        createLocationIterator: (volume: Volume) => LocationIterator(volume.grid.cells.data.length, 1, () => NullLocation),
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
    'solid': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceMeshParams>) => VolumeRepresentation('Isosurface mesh', ctx, getParams, IsosurfaceMeshVisual, getLoci),
    'wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceWireframeParams>) => VolumeRepresentation('Isosurface wireframe', ctx, getParams, IsosurfaceWireframeVisual, getLoci),
};

export const IsosurfaceParams = {
    ...IsosurfaceMeshParams,
    ...IsosurfaceWireframeParams,
    visuals: PD.MultiSelect(['solid'], PD.objectToOptions(IsosurfaceVisuals)),
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
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume)
});