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
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { Mat4 } from '../../../mol-math/linear-algebra/3d/mat4';
import { Quat } from '../../../mol-math/linear-algebra/3d/quat';

// ---- Axes visual (3 orientation lines per particle, baked/non-instanced) -----

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

const _m = Mat4();
const _q = Quat();
const _pos = Vec3();
const _end = Vec3();

function createAxesGeometry(_ctx: VisualContext, particles: ParticleList, _theme: Theme, props: AxesParticlesProps, lines?: Lines): Lines {
    const { count, coordinates, rotations } = particles;
    const { axisLength } = props;
    const segmentCount = Math.max(1, count * 3);
    const builder = LinesBuilder.create(segmentCount, Math.max(1, Math.ceil(segmentCount / 10)), lines);

    for (let i = 0; i < count; ++i) {
        Vec3.fromArray(_pos, coordinates, i * 3);

        if (rotations) {
            const o = i * 4;
            Quat.set(_q, rotations[o], rotations[o + 1], rotations[o + 2], rotations[o + 3]);
            Mat4.fromQuat(_m, _q);
        } else {
            Mat4.setIdentity(_m);
        }
        const basis = Mat4.extractBasis(_m);

        Vec3.scaleAndAdd(_end, _pos, basis.x, axisLength);
        builder.add(_pos[0], _pos[1], _pos[2], _end[0], _end[1], _end[2], i * 3 + 0);
        Vec3.scaleAndAdd(_end, _pos, basis.y, axisLength);
        builder.add(_pos[0], _pos[1], _pos[2], _end[0], _end[1], _end[2], i * 3 + 1);
        Vec3.scaleAndAdd(_end, _pos, basis.z, axisLength);
        builder.add(_pos[0], _pos[1], _pos[2], _end[0], _end[1], _end[2], i * 3 + 2);
    }

    const result = builder.getLines();
    result.setBoundingSphere(Particle.getBoundary(particles).sphere);
    return result;
}

function createAxesLocationIterator(particles: ParticleList, _geometry: Lines): LocationIterator {
    const { count } = particles;
    const location = Object.assign(Particle.Location(particles, 0), { groupIndex: 0 }) as AxesLocation;
    return LocationIterator(count * 3, 1, 1, (groupIndex: number) => {
        location.index = Math.floor(groupIndex / 3);
        location.groupIndex = groupIndex % 3;
        return location;
    }, true /* nonInstanceable */);
}

function getAxesLoci(pickingId: PickingId, particles: ParticleList, _props: AxesParticlesProps, id: number, _geometry: Lines): Loci {
    const { objectId, groupId } = pickingId;
    if (id !== objectId) return EmptyLoci;
    const particleIdx = Math.floor(groupId / 3);
    if (groupId < 0 || particleIdx >= particles.count) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofSingleton(particleIdx as any));
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
        setUpdateState: (state: VisualUpdateState, newParticles: ParticleList, currentParticles: ParticleList, newProps: AxesParticlesProps, currentProps: AxesParticlesProps) => {
            // axis lines are baked from particle coordinates/rotations, so any new particle data requires rebuilding the geometry
            if (newParticles !== currentParticles || newProps.axisLength !== currentProps.axisLength) {
                state.createGeometry = true;
            }
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
