/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { ColorTheme } from '../../../mol-theme/color';
import { UniformSizeTheme } from '../../../mol-theme/size/uniform';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Particle, ParticleList } from '../../../mol-model/particles/particle-list';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { VisualUpdateState } from '../../util';
import { VisualContext } from '../../visual';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../representation';
import { ParticleRepresentation, ParticleRepresentationProvider } from '../representation';
import { ParticleVisual } from '../visual';
import { Location } from '../../../mol-model/location';
import { WebGLContext } from '../../../mol-gl/webgl/context';

// ---- Axes visual (3 orientation lines, instanced per particle) ---------------

const AxisLengthOptions = { min: 0.1, max: 1000, step: 0.1 } as const;

export const AxesParticlesParams = {
    ...Lines.Params,
    axisLength: PD.Numeric(10, AxisLengthOptions, { description: 'Length of the particle orientation axes.' }),
    xColor: PD.Color(ColorNames.red),
    yColor: PD.Color(ColorNames.green),
    zColor: PD.Color(ColorNames.blue),
};
export type AxesParticlesParams = typeof AxesParticlesParams;
export type AxesParticlesProps = PD.Values<AxesParticlesParams>;

type AxesLocation = Particle.Location & { groupIndex: number }

function axisColorTheme(props: AxesParticlesProps): ColorTheme<{}> {
    const colors = [props.xColor, props.yColor, props.zColor] as const;
    const theme: ColorTheme<{}> = {
        factory: (() => theme) as any,
        granularity: 'group',
        color: (location: Location) => {
            const g = (location as any).groupIndex as number | undefined;
            if (g !== undefined && g >= 0 && g < 3) return colors[g];
            return Color(0xCCCCCC);
        },
        props: {},
    };
    return theme;
}

function createAxesGeometry(_ctx: VisualContext, _particles: ParticleList, _theme: Theme, props: AxesParticlesProps, lines?: Lines): Lines {
    const { axisLength } = props;
    const builder = LinesBuilder.create(3, 3, lines);
    builder.add(0, 0, 0, axisLength, 0, 0, 0); // X axis, group 0
    builder.add(0, 0, 0, 0, axisLength, 0, 1); // Y axis, group 1
    builder.add(0, 0, 0, 0, 0, axisLength, 2); // Z axis, group 2
    return builder.getLines();
}

function createAxesLocationIterator(particles: ParticleList, _geometry: Lines): LocationIterator {
    const { count } = particles;
    const location = Object.assign(Particle.Location(particles, 0), { groupIndex: 0 }) as AxesLocation;
    return LocationIterator(3, count, 1, (groupIndex, instanceIndex) => {
        location.index = instanceIndex;
        location.groupIndex = groupIndex;
        return location;
    });
}

function getAxesLoci(pickingId: PickingId, particles: ParticleList, _props: AxesParticlesProps, id: number, _geometry: Lines): Loci {
    const { objectId, instanceId } = pickingId;
    if (id !== objectId) return EmptyLoci;
    if (instanceId < 0 || instanceId >= particles.count) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofSingleton(instanceId as any));
}

function eachAxesLoci(loci: Loci, particles: ParticleList, _props: AxesParticlesProps, apply: (interval: Interval) => boolean, _geometry: Lines): boolean {
    if (!Particle.isLoci(loci)) return false;
    if (Particle.isLociEmpty(loci)) return false;
    if (loci.particles !== particles) return false;
    let changed = false;
    OrderedSet.forEach(loci.indices, idx => {
        // Each particle maps to 3 consecutive group slots (one per axis).
        const start = idx * 3;
        if (apply(Interval.ofBounds(start, start + 3))) changed = true;
    });
    return changed;
}

export function AxesParticlesVisual(materialId: number, _particles?: ParticleList, _props?: PD.Values<AxesParticlesParams>, _webgl?: WebGLContext): ParticleVisual<AxesParticlesParams> {
    return ParticleVisual<Lines, AxesParticlesParams>({
        defaultProps: PD.getDefaultValues(AxesParticlesParams),
        createGeometry: createAxesGeometry,
        createLocationIterator: createAxesLocationIterator,
        getLoci: getAxesLoci,
        eachLocation: eachAxesLoci,
        setUpdateState: (state: VisualUpdateState, _np: ParticleList, _cp: ParticleList, newProps: AxesParticlesProps, currentProps: AxesParticlesProps) => {
            state.createGeometry = newProps.axisLength !== currentProps.axisLength;
            if (newProps.xColor !== currentProps.xColor ||
                newProps.yColor !== currentProps.yColor ||
                newProps.zColor !== currentProps.zColor) {
                state.updateColor = true;
            }
        },
        overrideTheme: (_theme: Theme, props: AxesParticlesProps): Theme => ({
            color: axisColorTheme(props),
            size: UniformSizeTheme({} as any, { value: 1 }),
        }),
        geometryUtils: Lines.Utils,
    }, materialId);
}

// ---- OrientationParticles representation ------------------------------------

function getAllLoci(particles: ParticleList): Loci {
    const { count } = particles;
    if (count === 0) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofBounds(0, count) as any);
}

const OrientationParticlesVisualKinds = ['axes'] as const;
type OrientationParticlesVisualKind = typeof OrientationParticlesVisualKinds[number];

const OrientationVisuals = {
    'axes': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, AxesParticlesParams>) =>
        ParticleRepresentation('Axes', ctx, getParams, AxesParticlesVisual, getAllLoci),
};

export const OrientationParticlesParams = {
    ...AxesParticlesParams,
    visuals: PD.MultiSelect(['axes'] as OrientationParticlesVisualKind[], PD.arrayToOptions(OrientationParticlesVisualKinds as unknown as OrientationParticlesVisualKind[])),
};
export type OrientationParticlesParams = typeof OrientationParticlesParams;
export type OrientationParticlesProps = PD.Values<OrientationParticlesParams>;

export function getOrientationParticlesParams(_ctx: ThemeRegistryContext, data: ParticleList) {
    const hasRotations = !!data.rotations;
    return {
        ...OrientationParticlesParams,
        visuals: PD.MultiSelect(
            hasRotations ? (['axes'] as OrientationParticlesVisualKind[]) : ([] as OrientationParticlesVisualKind[]),
            PD.arrayToOptions(OrientationParticlesVisualKinds as unknown as OrientationParticlesVisualKind[]),
        ),
    };
}

export type OrientationParticlesRepresentation = ParticleRepresentation<OrientationParticlesParams>

export function OrientationParticlesRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, OrientationParticlesParams>): OrientationParticlesRepresentation {
    return Representation.createMulti('Orientation Particles', ctx, getParams, Representation.StateBuilder, OrientationVisuals as unknown as Representation.Def<ParticleList, OrientationParticlesParams>);
}

export const OrientationParticlesRepresentationProvider = ParticleRepresentationProvider({
    name: 'orientation',
    label: 'Orientation',
    description: 'Displays particle orientation axes.',
    factory: OrientationParticlesRepresentation,
    getParams: getOrientationParticlesParams,
    defaultValues: PD.getDefaultValues(OrientationParticlesParams),
    defaultColorTheme: { name: 'particle-entity' },
    defaultSizeTheme: { name: 'particle-size' },
    isApplicable: (data: ParticleList) => data.count > 0 && !!data.rotations,
});
