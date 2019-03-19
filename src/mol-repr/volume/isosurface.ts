/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData, VolumeIsoValue } from 'mol-model/volume'
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { EmptyLoci } from 'mol-model/loci';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh, computeMarchingCubesLines } from 'mol-geo/util/marching-cubes/algorithm';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { VisualUpdateState } from 'mol-repr/util';
import { RepresentationContext, RepresentationParamsGetter, Representation } from 'mol-repr/representation';
import { Theme, ThemeRegistryContext } from 'mol-theme/theme';
import { VisualContext } from 'mol-repr/visual';
import { NullLocation } from 'mol-model/location';
import { Lines } from 'mol-geo/geometry/lines/lines';

const defaultStats: VolumeData['dataStats'] = { min: -1, max: 1, mean: 0, sigma: 0.1, histogram: { bins: [-1, 1], counts: [1] }  };
export function createIsoValueParam(defaultValue: VolumeIsoValue, stats?: VolumeData['dataStats']) {
    const sts = stats || defaultStats;
    const { min, max, mean, sigma, histogram } = sts;

    // using ceil/floor could lead to "ouf of bounds" when converting
    const relMin = (min - mean) / sigma;
    const relMax = (max - mean) / sigma;

    let def = defaultValue;
    if (defaultValue.kind === 'absolute') {
        if (defaultValue.absoluteValue < min) def = VolumeIsoValue.absolute(min);
        else if (defaultValue.absoluteValue > max) def = VolumeIsoValue.absolute(max);
    } else {
        if (defaultValue.relativeValue < relMin) def = VolumeIsoValue.relative(relMin);
        else if (defaultValue.relativeValue > relMax) def = VolumeIsoValue.relative(relMax);
    }

    return PD.Conditioned(
        def,
        {
            'absolute': PD.Converted(
                (v: VolumeIsoValue) => VolumeIsoValue.toAbsolute(v, VolumeData.One.dataStats).absoluteValue,
                (v: number) => VolumeIsoValue.absolute(v),
                histogram ? PD.Histogram(mean, histogram) : PD.Numeric(mean, { min, max, step: sigma / 100 })
            ),
            'relative': PD.Converted(
                (v: VolumeIsoValue) => VolumeIsoValue.toRelative(v, VolumeData.One.dataStats).relativeValue,
                (v: number) => VolumeIsoValue.relative(v),
                PD.Numeric(Math.min(1, relMax), { min: relMin, max: relMax, step: Math.round(((max - min) / sigma)) / 100 })
            )
        },
        (v: VolumeIsoValue) => v.kind === 'absolute' ? 'absolute' : 'relative',
        (v: VolumeIsoValue, c: 'absolute' | 'relative') => c === 'absolute' ? VolumeIsoValue.toAbsolute(v, sts) : VolumeIsoValue.toRelative(v, sts)
    );
}

export const IsoValueParam = createIsoValueParam(VolumeIsoValue.relative(2));
type IsoValueParam = typeof IsoValueParam

export const VolumeIsosurfaceParams = {
    isoValue: IsoValueParam
}
export type VolumeIsosurfaceParams = typeof VolumeIsosurfaceParams
export type VolumeIsosurfaceProps = PD.Values<VolumeIsosurfaceParams>

//

export async function createVolumeIsosurfaceMesh(ctx: VisualContext, volume: VolumeData, theme: Theme, props: VolumeIsosurfaceProps, mesh?: Mesh) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const surface = await computeMarchingCubesMesh({
        isoLevel: VolumeIsoValue.toAbsolute(props.isoValue, volume.dataStats).absoluteValue,
        scalarField: volume.data
    }, mesh).runAsChild(ctx.runtime);

    const transform = VolumeData.getGridToCartesianTransform(volume);
    ctx.runtime.update({ message: 'Transforming mesh...' });
    Mesh.transformImmediate(surface, transform);
    Mesh.computeNormalsImmediate(surface)

    return surface;
}

export const IsosurfaceMeshParams = {
    ...Mesh.Params,
    ...VolumeIsosurfaceParams
}
export type IsosurfaceMeshParams = typeof IsosurfaceMeshParams

export function IsosurfaceMeshVisual(materialId: number): VolumeVisual<IsosurfaceMeshParams> {
    return VolumeVisual<Mesh, IsosurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceMeshParams),
        createGeometry: createVolumeIsosurfaceMesh,
        createLocationIterator: (volume: VolumeData) => LocationIterator(1, 1, () => NullLocation),
        getLoci: () => EmptyLoci,
        eachLocation: () => false,
        setUpdateState: (state: VisualUpdateState, volume: VolumeData, newProps: PD.Values<IsosurfaceMeshParams>, currentProps: PD.Values<IsosurfaceMeshParams>) => {
            if (!VolumeIsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.dataStats)) state.createGeometry = true
        },
        geometryUtils: Mesh.Utils
    }, materialId)
}

//

export async function createVolumeIsosurfaceWireframe(ctx: VisualContext, volume: VolumeData, theme: Theme, props: VolumeIsosurfaceProps, lines?: Lines) {
    ctx.runtime.update({ message: 'Marching cubes...' });

    const wireframe = await computeMarchingCubesLines({
        isoLevel: VolumeIsoValue.toAbsolute(props.isoValue, volume.dataStats).absoluteValue,
        scalarField: volume.data
    }, lines).runAsChild(ctx.runtime)

    const transform = VolumeData.getGridToCartesianTransform(volume);
    Lines.transformImmediate(wireframe, transform)

    return wireframe;
}

export const IsosurfaceWireframeParams = {
    ...Lines.Params,
    ...VolumeIsosurfaceParams
}
export type IsosurfaceWireframeParams = typeof IsosurfaceWireframeParams

export function IsosurfaceWireframeVisual(materialId: number): VolumeVisual<IsosurfaceWireframeParams> {
    return VolumeVisual<Lines, IsosurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(IsosurfaceWireframeParams),
        createGeometry: createVolumeIsosurfaceWireframe,
        createLocationIterator: (volume: VolumeData) => LocationIterator(1, 1, () => NullLocation),
        getLoci: () => EmptyLoci,
        eachLocation: () => false,
        setUpdateState: (state: VisualUpdateState, volume: VolumeData, newProps: PD.Values<IsosurfaceWireframeParams>, currentProps: PD.Values<IsosurfaceWireframeParams>) => {
            if (!VolumeIsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.dataStats)) state.createGeometry = true
        },
        geometryUtils: Lines.Utils
    }, materialId)
}

//

const IsosurfaceVisuals = {
    'solid': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceMeshParams>) => VolumeRepresentation('Isosurface mesh', ctx, getParams, IsosurfaceMeshVisual),
    'wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceWireframeParams>) => VolumeRepresentation('Isosurface wireframe', ctx, getParams, IsosurfaceWireframeVisual),
}
type IsosurfaceVisualName = keyof typeof IsosurfaceVisuals
const IsosurfaceVisualOptions = Object.keys(IsosurfaceVisuals).map(name => [name, name] as [IsosurfaceVisualName, string])

export const IsosurfaceParams = {
    ...IsosurfaceMeshParams,
    ...IsosurfaceWireframeParams,
    visuals: PD.MultiSelect<IsosurfaceVisualName>(['solid'], IsosurfaceVisualOptions),
}
export type IsosurfaceParams = typeof IsosurfaceParams
export function getIsosurfaceParams(ctx: ThemeRegistryContext, volume: VolumeData) {
    const p = PD.clone(IsosurfaceParams);
    p.isoValue = createIsoValueParam(VolumeIsoValue.relative(2), volume.dataStats);
    return p
}

export type IsosurfaceRepresentation = VolumeRepresentation<IsosurfaceParams>
export function IsosurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<VolumeData, IsosurfaceParams>): IsosurfaceRepresentation {
    return Representation.createMulti('Isosurface', ctx, getParams, Representation.StateBuilder, IsosurfaceVisuals as unknown as Representation.Def<VolumeData, IsosurfaceParams>)
}

export const IsosurfaceRepresentationProvider: VolumeRepresentationProvider<IsosurfaceParams> = {
    label: 'Isosurface',
    description: 'Displays an isosurface of volumetric data.',
    factory: IsosurfaceRepresentation,
    getParams: getIsosurfaceParams,
    defaultValues: PD.getDefaultValues(IsosurfaceParams),
    defaultColorTheme: 'uniform',
    defaultSizeTheme: 'uniform'
}