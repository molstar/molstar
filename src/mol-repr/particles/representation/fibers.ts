/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder, StripLinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { addTube } from '../../../mol-geo/geometry/mesh/builder/tube';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Particle, ParticleList } from '../../../mol-model/particles/particle-list';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { VisualUpdateState } from '../../util';
import { VisualContext } from '../../visual';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../representation';
import { ParticleRepresentation, ParticleRepresentationProvider } from '../representation';
import { ParticleVisual } from '../visual';
import { computeFrenetFrames } from '../../../mol-math/linear-algebra/3d/frenet-frames';
import { WebGLContext } from '../../../mol-gl/webgl/context';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3spline = Vec3.spline;
const v3toArray = Vec3.toArray;

const Tension = 0.5; // standard Catmull-Rom tension (matches StandardTension in polymer)

// ---- Shared helpers ---------------------------------------------------------

function fiberInterpolatedPointCount(particles: ParticleList, linearSegments: number): number {
    const { fibers } = particles;
    if (!fibers) return 0;
    let count = 0;
    for (let f = 0; f < fibers.count; ++f) {
        const n = fibers.offsets[f + 1] - fibers.offsets[f];
        if (n >= 2) count += (n - 1) * linearSegments + 1;
    }
    return count;
}

function createFibersLocationIterator(particles: ParticleList, _geometry: Lines | Mesh): LocationIterator {
    const { fibers } = particles;
    const groupCount = fibers ? fibers.count : 0;
    const location = Particle.Location(particles, 0);
    return LocationIterator(groupCount, 1, 1, (groupIndex: number) => {
        if (fibers && groupIndex < fibers.count) {
            location.index = fibers.indices[fibers.offsets[groupIndex]];
        }
        return location;
    }, true /* nonInstanceable */);
}

function getFibersLoci(pickingId: PickingId, particles: ParticleList, _props: any, id: number, _geometry: Lines | Mesh): Loci {
    const { objectId, groupId } = pickingId;
    if (id !== objectId) return EmptyLoci;
    const { fibers } = particles;
    if (!fibers || groupId < 0 || groupId >= fibers.count) return EmptyLoci;
    const start = fibers.offsets[groupId];
    const end = fibers.offsets[groupId + 1];
    const indices: number[] = [];
    for (let i = start; i < end; ++i) {
        indices.push(fibers.indices[i]);
    }
    indices.sort((a, b) => a - b);
    return Particle.Loci(particles, OrderedSet.ofSortedArray(indices));
}

function eachFibers(loci: Loci, particles: ParticleList, _props: any, apply: (interval: Interval) => boolean, _geometry: Lines | Mesh): boolean {
    if (!Particle.isLoci(loci)) return false;
    if (Particle.isLociEmpty(loci)) return false;
    if (loci.particles !== particles) return false;
    const { fibers } = particles;
    if (!fibers) return false;
    let changed = false;
    for (let f = 0; f < fibers.count; ++f) {
        const start = fibers.offsets[f];
        const end = fibers.offsets[f + 1];
        for (let i = start; i < end; ++i) {
            if (OrderedSet.has(loci.indices, fibers.indices[i])) {
                if (apply(Interval.ofSingleton(f))) changed = true;
                break;
            }
        }
    }
    return changed;
}

function getAllFibersLoci(particles: ParticleList): Loci {
    const { count, fibers } = particles;
    if (count === 0 || !fibers || fibers.count === 0) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofBounds(0, count) as any);
}

// ---- Curve interpolation ----------------------------------------------------

const _p0 = Vec3(), _p1 = Vec3(), _p2 = Vec3(), _p3 = Vec3(), _pt = Vec3();

/**
 * Build Catmull-Rom interpolated curve points for a single fiber (fStart..fEnd exclusive).
 * Output: (n-1)*linearSegments + 1 points written into `curvePoints` starting at `offset`.
 * Returns the number of points written.
 */
function buildFiberCurvePoints(
    particles: ParticleList,
    fStart: number, fEnd: number,
    linearSegments: number,
    curvePoints: Float32Array,
    offset: number
): number {
    const n = fEnd - fStart;
    const { coordinates, fibers } = particles;

    const positions: Vec3[] = new Array(n);
    for (let k = 0; k < n; ++k) {
        const pIdx = fibers!.indices[fStart + k];
        positions[k] = Vec3.fromArray(Vec3(), coordinates, pIdx * 3);
    }

    let outIdx = offset;
    for (let k = 0; k < n - 1; ++k) {
        Vec3.copy(_p1, positions[k]);
        Vec3.copy(_p2, positions[k + 1]);
        Vec3.copy(_p0, k > 0 ? positions[k - 1] : _p1);
        Vec3.copy(_p3, k < n - 2 ? positions[k + 2] : _p2);

        for (let s = 0; s < linearSegments; ++s) {
            v3spline(_pt, _p0, _p1, _p2, _p3, s / linearSegments, Tension);
            v3toArray(_pt, curvePoints, outIdx * 3);
            ++outIdx;
        }
    }
    // Final particle position
    v3toArray(positions[n - 1], curvePoints, outIdx * 3);
    ++outIdx;

    return outIdx - offset;
}

/**
 * Build linearly interpolated size values (width/height) for a single fiber.
 * Written into `widthValues` and `heightValues` starting at `offset`.
 */
function buildFiberSizeValues(
    particles: ParticleList,
    fStart: number, fEnd: number,
    linearSegments: number,
    tubeSizeFactor: number,
    theme: Theme,
    widthValues: Float32Array,
    heightValues: Float32Array,
    offset: number
): void {
    const n = fEnd - fStart;
    const { fibers } = particles;
    const location = Particle.Location(particles);

    const sizes = new Float32Array(n);
    for (let k = 0; k < n; ++k) {
        location.index = fibers!.indices[fStart + k];
        sizes[k] = theme.size.size(location) * tubeSizeFactor;
    }

    let outIdx = offset;
    for (let k = 0; k < n - 1; ++k) {
        const s1 = sizes[k], s2 = sizes[k + 1];
        for (let s = 0; s < linearSegments; ++s) {
            const size = s1 + (s2 - s1) * (s / linearSegments);
            widthValues[outIdx] = size;
            heightValues[outIdx] = size;
            ++outIdx;
        }
    }
    widthValues[outIdx] = sizes[n - 1];
    heightValues[outIdx] = sizes[n - 1];
}

// ---- Lines visual -----------------------------------------------------------

export const FibersLinesParams = {
    ...Lines.Params,
    linearSegments: PD.Numeric(8, { min: 1, max: 48, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    useLineStrips: PD.Boolean(true),
};
export type FibersLinesParams = typeof FibersLinesParams;
export type FibersLinesProps = PD.Values<FibersLinesParams>;

function createParticleFibersLines(_ctx: VisualContext, particles: ParticleList, _theme: Theme, props: FibersLinesProps, lines?: Lines): Lines {
    const { fibers } = particles;
    const { linearSegments, useLineStrips } = props;
    const totalPoints = fiberInterpolatedPointCount(particles, linearSegments);

    if (useLineStrips) {
        const builder = StripLinesBuilder.create(Math.max(1, totalPoints), Math.max(1, Math.ceil(totalPoints / 10)), lines);

        if (fibers) {
            const fiberCurvePoints = new Float32Array(totalPoints * 3);
            let globalOffset = 0;
            for (let f = 0; f < fibers.count; ++f) {
                const fStart = fibers.offsets[f];
                const fEnd = fibers.offsets[f + 1];
                if (fEnd - fStart < 2) continue;

                const nPoints = buildFiberCurvePoints(particles, fStart, fEnd, linearSegments, fiberCurvePoints, globalOffset);
                builder.start(f);
                for (let i = 0; i < nPoints; ++i) {
                    const oi = (globalOffset + i) * 3;
                    builder.add(fiberCurvePoints[oi], fiberCurvePoints[oi + 1], fiberCurvePoints[oi + 2]);
                }
                builder.end();
                globalOffset += nPoints;
            }
        }

        const result = builder.getLines();
        result.setBoundingSphere(Particle.getBoundary(particles).sphere);
        return result;
    } else {
        const builder = LinesBuilder.create(Math.max(1, totalPoints), Math.max(1, Math.ceil(totalPoints / 10)), lines);
        const a = Vec3(), b = Vec3();

        if (fibers) {
            const fiberCurvePoints = new Float32Array(totalPoints * 3);
            let globalOffset = 0;
            for (let f = 0; f < fibers.count; ++f) {
                const fStart = fibers.offsets[f];
                const fEnd = fibers.offsets[f + 1];
                if (fEnd - fStart < 2) continue;

                const nPoints = buildFiberCurvePoints(particles, fStart, fEnd, linearSegments, fiberCurvePoints, globalOffset);
                Vec3.fromArray(a, fiberCurvePoints, globalOffset * 3);
                for (let i = 1; i < nPoints; ++i) {
                    Vec3.fromArray(b, fiberCurvePoints, (globalOffset + i) * 3);
                    builder.addVec(a, b, f);
                    Vec3.copy(a, b);
                }
                globalOffset += nPoints;
            }
        }

        const result = builder.getLines();
        result.setBoundingSphere(Particle.getBoundary(particles).sphere);
        return result;
    }
}

export function FibersLinesVisual(materialId: number, _particles?: ParticleList, _props?: PD.Values<FibersLinesParams>, _webgl?: WebGLContext): ParticleVisual<FibersLinesParams> {
    return ParticleVisual<Lines, FibersLinesParams>({
        defaultProps: PD.getDefaultValues(FibersLinesParams),
        createGeometry: createParticleFibersLines,
        createLocationIterator: createFibersLocationIterator,
        getLoci: getFibersLoci,
        eachLocation: eachFibers,
        setUpdateState: (state: VisualUpdateState, newParticles: ParticleList, currentParticles: ParticleList, newProps: FibersLinesProps, currentProps: FibersLinesProps) => {
            // fiber curves are baked from particle coordinates, so any new particle data requires rebuilding the geometry
            if (newParticles !== currentParticles ||
                newProps.linearSegments !== currentProps.linearSegments ||
                newProps.useLineStrips !== currentProps.useLineStrips) {
                state.createGeometry = true;
            }
        },
        geometryUtils: Lines.Utils,
    }, materialId);
}

// ---- Tube mesh visual -------------------------------------------------------

export const FibersTubeMeshParams = {
    ...Mesh.Params,
    tubeSizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    linearSegments: PD.Numeric(8, { min: 1, max: 48, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    radialSegments: PD.Numeric(8, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
};
export type FibersTubeMeshParams = typeof FibersTubeMeshParams;
export type FibersTubeMeshProps = PD.Values<FibersTubeMeshParams>;

function addFiberTube(
    builderState: MeshBuilder.State,
    particles: ParticleList,
    fStart: number, fEnd: number,
    fiberIndex: number,
    linearSegments: number,
    radialSegments: number,
    tubeSizeFactor: number,
    theme: Theme
): void {
    const n = fEnd - fStart;
    if (n < 2) return;

    const totalPoints = (n - 1) * linearSegments + 1;
    const curvePoints = new Float32Array(totalPoints * 3);
    const normalVectors = new Float32Array(totalPoints * 3);
    const binormalVectors = new Float32Array(totalPoints * 3);
    const widthValues = new Float32Array(totalPoints);
    const heightValues = new Float32Array(totalPoints);

    buildFiberCurvePoints(particles, fStart, fEnd, linearSegments, curvePoints, 0);
    buildFiberSizeValues(particles, fStart, fEnd, linearSegments, tubeSizeFactor, theme, widthValues, heightValues, 0);

    builderState.currentGroup = fiberIndex;
    computeFrenetFrames(curvePoints, normalVectors, binormalVectors, totalPoints);
    addTube(builderState, curvePoints, normalVectors, binormalVectors, totalPoints - 1, radialSegments, widthValues, heightValues, true, true, 'elliptical');
}

function createParticleFibersTubeMesh(_ctx: VisualContext, particles: ParticleList, theme: Theme, props: FibersTubeMeshProps, mesh?: Mesh): Mesh {
    const { fibers } = particles;
    const { tubeSizeFactor, linearSegments, radialSegments } = props;
    const totalPoints = fiberInterpolatedPointCount(particles, linearSegments);
    const vertexCount = Math.max(1, totalPoints * radialSegments * 2);
    const builderState = MeshBuilder.createState(vertexCount, Math.ceil(vertexCount / 10), mesh);

    if (fibers) {
        for (let f = 0; f < fibers.count; ++f) {
            const fStart = fibers.offsets[f];
            const fEnd = fibers.offsets[f + 1];
            if (fEnd - fStart < 2) continue;
            addFiberTube(builderState, particles, fStart, fEnd, f, linearSegments, radialSegments, tubeSizeFactor, theme);
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    m.setBoundingSphere(Particle.getBoundary(particles).sphere);
    return m;
}

export function FibersTubeMeshVisual(materialId: number, _particles?: ParticleList, _props?: PD.Values<FibersTubeMeshParams>, _webgl?: WebGLContext): ParticleVisual<FibersTubeMeshParams> {
    return ParticleVisual<Mesh, FibersTubeMeshParams>({
        defaultProps: PD.getDefaultValues(FibersTubeMeshParams),
        createGeometry: createParticleFibersTubeMesh,
        createLocationIterator: createFibersLocationIterator,
        getLoci: getFibersLoci,
        eachLocation: eachFibers,
        setUpdateState: (state: VisualUpdateState, newParticles: ParticleList, currentParticles: ParticleList, newProps: FibersTubeMeshProps, currentProps: FibersTubeMeshProps) => {
            // fiber tubes are baked from particle coordinates, so any new particle data requires rebuilding the geometry
            if (newParticles !== currentParticles ||
                newProps.tubeSizeFactor !== currentProps.tubeSizeFactor ||
                newProps.linearSegments !== currentProps.linearSegments ||
                newProps.radialSegments !== currentProps.radialSegments) {
                state.createGeometry = true;
            }
        },
        geometryUtils: Mesh.Utils,
    }, materialId);
}

// ---- Combined Fibers representation -----------------------------------------

const FibersVisuals = {
    'lines': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, FibersLinesParams>) =>
        ParticleRepresentation('Fibers lines', ctx, getParams, FibersLinesVisual, getAllFibersLoci),
    'tube-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, FibersTubeMeshParams>) =>
        ParticleRepresentation('Fibers tube-mesh', ctx, getParams, FibersTubeMeshVisual, getAllFibersLoci),
};

export const FibersParams = {
    ...FibersLinesParams,
    ...FibersTubeMeshParams,
    visuals: PD.MultiSelect(['lines'], PD.objectToOptions(FibersVisuals)),
};
export type FibersParams = typeof FibersParams;

export function getFibersParams(_ctx: ThemeRegistryContext, _particles: ParticleList) {
    return PD.clone(FibersParams);
}

export type FibersRepresentation = ParticleRepresentation<FibersParams>;
export function FibersRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, FibersParams>): FibersRepresentation {
    return Representation.createMulti('Fibers', ctx, getParams, Representation.StateBuilder, FibersVisuals as unknown as Representation.Def<ParticleList, FibersParams>);
}

export const FibersRepresentationProvider = ParticleRepresentationProvider({
    name: 'fibers',
    label: 'Fibers',
    description: 'Displays particle fibers as lines or tubes.',
    factory: FibersRepresentation,
    getParams: getFibersParams,
    defaultValues: PD.getDefaultValues(FibersParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (data: ParticleList) => data.count > 0 && !!data.fibers && data.fibers.count > 0,
});
