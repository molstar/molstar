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
import { toPrecision } from '../../mol-util/number';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { Tensor } from '../../mol-math/linear-algebra';
import { fillSerial } from '../../mol-util/array';

const defaultStats: Grid['stats'] = { min: -1, max: 1, mean: 0, sigma: 0.1  };
export function createIsoValueParam(defaultValue: Volume.IsoValue, stats?: Grid['stats']) {
    const sts = stats || defaultStats;
    const { min, max, mean, sigma } = sts;

    // using ceil/floor could lead to "ouf of bounds" when converting
    const relMin = (min - mean) / sigma;
    const relMax = (max - mean) / sigma;

    let def = defaultValue;
    if (defaultValue.kind === 'absolute') {
        if (defaultValue.absoluteValue < min) def = Volume.IsoValue.absolute(min);
        else if (defaultValue.absoluteValue > max) def = Volume.IsoValue.absolute(max);
    } else {
        if (defaultValue.relativeValue < relMin) def = Volume.IsoValue.relative(relMin);
        else if (defaultValue.relativeValue > relMax) def = Volume.IsoValue.relative(relMax);
    }

    return PD.Conditioned(
        def,
        {
            'absolute': PD.Converted(
                (v: Volume.IsoValue) => Volume.IsoValue.toAbsolute(v, Grid.One.stats).absoluteValue,
                (v: number) => Volume.IsoValue.absolute(v),
                PD.Numeric(mean, { min, max, step: toPrecision(sigma / 100, 2) })
            ),
            'relative': PD.Converted(
                (v: Volume.IsoValue) => Volume.IsoValue.toRelative(v, Grid.One.stats).relativeValue,
                (v: number) => Volume.IsoValue.relative(v),
                PD.Numeric(Math.min(1, relMax), { min: relMin, max: relMax, step: toPrecision(Math.round(((max - min) / sigma)) / 100, 2) })
            )
        },
        (v: Volume.IsoValue) => v.kind === 'absolute' ? 'absolute' : 'relative',
        (v: Volume.IsoValue, c: 'absolute' | 'relative') => c === 'absolute' ? Volume.IsoValue.toAbsolute(v, sts) : Volume.IsoValue.toRelative(v, sts),
        { isEssential: true }
    );
}

export const IsoValueParam = createIsoValueParam(Volume.IsoValue.relative(2));
type IsoValueParam = typeof IsoValueParam

export const VolumeIsosurfaceParams = {
    isoValue: IsoValueParam
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

function eachIsosurface(loci: Loci, volume: Volume, props: VolumeIsosurfaceProps, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (Volume.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (apply(Interval.ofLength(volume.grid.cells.data.length))) changed = true;
    } else if (Volume.Isosurface.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (!Volume.IsoValue.areSame(loci.isoValue, props.isoValue, volume.grid.stats)) return false;
        if (apply(Interval.ofLength(volume.grid.cells.data.length))) changed = true;
    } else if (Volume.Cell.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (Interval.is(loci.indices)) {
            if (apply(loci.indices)) changed = true;
        } else {
            OrderedSet.forEach(loci.indices, v => {
                if (apply(Interval.ofSingleton(v))) changed = true;
            });
        }
    }
    return changed;
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
    ctx.runtime.update({ message: 'Transforming mesh...' });
    Mesh.transform(surface, transform);
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
    p.isoValue = createIsoValueParam(Volume.IsoValue.relative(2), volume.grid.stats);
    return p;
}

export type IsosurfaceRepresentation = VolumeRepresentation<IsosurfaceParams>
export function IsosurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceParams>): IsosurfaceRepresentation {
    return Representation.createMulti('Isosurface', ctx, getParams, Representation.StateBuilder, IsosurfaceVisuals as unknown as Representation.Def<Volume, IsosurfaceParams>);
}

export const IsosurfaceRepresentationProvider = VolumeRepresentationProvider({
    name: 'isosurface',
    label: 'Isosurface',
    description: 'Displays an isosurface of volumetric data.',
    factory: IsosurfaceRepresentation,
    getParams: getIsosurfaceParams,
    defaultValues: PD.getDefaultValues(IsosurfaceParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume)
});