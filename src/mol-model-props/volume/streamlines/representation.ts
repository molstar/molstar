/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder, StripLinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { addTube } from '../../../mol-geo/geometry/mesh/builder/tube';
import { Volume, Grid } from '../../../mol-model/volume';
import { VolumeRepresentation, VolumeRepresentationProvider } from '../../../mol-repr/volume/representation';
import { VisualUpdateState } from '../../../mol-repr/util';
import { VisualContext } from '../../../mol-repr/visual';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../../../mol-repr/representation';
import { CommonStreamlinesParams, createStreamlinesLocationIterator, eachStreamlines, getStreamlinesLoci, getStreamlinesVisualLoci, Streamline, StreamlinesLocation, StreamlinesIndex, Streamlines, streamlinePassesFilter } from './shared';
import { Sphere3D } from '../../../mol-math/geometry';
import { StreamlinesProvider } from '../streamlines';
import { CustomProperty } from '../../common/custom-property';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { computeFrenetFrames } from '../../../mol-math/linear-algebra/3d/frenet-frames';
import { VolumeVisual } from '../../../mol-repr/volume/visual';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4 = Vec3.transformMat4;
const v3transformMat4Offset = Vec3.transformMat4Offset;
const v3toArray = Vec3.toArray;

export const StreamlinesLinesParams = {
    ...CommonStreamlinesParams,
    ...Lines.Params,
    useLineStrips: PD.Boolean(true),
};
export type StreamlinesLinesParams = typeof StreamlinesLinesParams
export type StreamlinesLinesProps = PD.Values<StreamlinesLinesParams>

export function VolumeStreamlinesLinesVisual(materialId: number): VolumeVisual<StreamlinesLinesParams> {
    return VolumeVisual<Lines, StreamlinesLinesParams>({
        defaultProps: PD.getDefaultValues(StreamlinesLinesParams),
        createGeometry: createVolumeStreamlinesLines,
        createLocationIterator: createStreamlinesLocationIterator,
        getLoci: getStreamlinesLoci,
        eachLocation: eachStreamlines,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<StreamlinesLinesParams>, currentProps: PD.Values<StreamlinesLinesParams>) => {
            const streamlinesHash = StreamlinesProvider.get(volume).version;
            if ((state.info.streamlinesHash as number) !== streamlinesHash) {
                if (state.info.streamlinesHash !== undefined) {
                    state.createGeometry = true;
                    state.updateLocation = true;
                }
                state.info.streamlinesHash = streamlinesHash;
            }
            if (newProps.anchorEnabled !== currentProps.anchorEnabled ||
                !Vec3.equals(newProps.anchorCenter, currentProps.anchorCenter) ||
                newProps.anchorRadius !== currentProps.anchorRadius) {
                state.createGeometry = true;
            }
            if (newProps.dashEnabled !== currentProps.dashEnabled ||
                newProps.dashPoints !== currentProps.dashPoints ||
                newProps.dashShift !== currentProps.dashShift) {
                state.createGeometry = true;
            }
            if (newProps.useLineStrips !== currentProps.useLineStrips) {
                state.createGeometry = true;
            }
        },
        initUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<StreamlinesLinesParams>, newTheme: Theme) => {
            const streamlinesHash = StreamlinesProvider.get(volume).version;
            state.info.streamlinesHash = streamlinesHash;
        },
        geometryUtils: Lines.Utils,
    }, materialId);
}

function streamlinePointCount(streamlines: Streamlines): number {
    let count = 0;
    for (const streamline of streamlines) {
        count += streamline.length;
    }
    return count;
}

function createVolumeStreamlinesLines(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: StreamlinesLinesProps, lines?: Lines): Lines {
    const { cells: { space } } = volume.grid;
    const gridDimension = space.dimensions as Vec3;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const streamlines = StreamlinesProvider.get(volume).value!;
    const pointCount = streamlinePointCount(streamlines);

    const { dashEnabled, dashPoints, dashShift } = props;
    const cycleLength = dashPoints * 2;

    let builder: LinesBuilder | StripLinesBuilder;

    if (props.useLineStrips) {
        const _builder = StripLinesBuilder.create(pointCount, Math.ceil(pointCount / 10), lines);
        builder = _builder;

        const b = Vec3();
        for (let s = 0, sl = streamlines.length; s < sl; ++s) {
            const l = streamlines[s];
            if (!streamlinePassesFilter(l, gridToCartn, props)) continue;

            if (dashEnabled) {
                let inDash = false;
                for (let i = 0, il = l.length; i < il; ++i) {
                    const inCycle = i % cycleLength;
                    const shouldDraw = dashShift ? (inCycle >= dashPoints) : (inCycle < dashPoints);

                    if (shouldDraw) {
                        if (!inDash) {
                            _builder.start(s);
                            inDash = true;
                        }
                        v3transformMat4(b, l[i], gridToCartn);
                        _builder.addVec(b);
                    } else if (inDash) {
                        v3transformMat4(b, l[i], gridToCartn);
                        _builder.addVec(b);
                        _builder.end();
                        inDash = false;
                    }
                }
                if (inDash) _builder.end();
            } else {
                _builder.start(s);
                for (let i = 0, il = l.length; i < il; ++i) {
                    v3transformMat4(b, l[i], gridToCartn);
                    _builder.addVec(b);
                }
                _builder.end();
            }
        }
    } else {
        const _builder = LinesBuilder.create(pointCount, Math.ceil(pointCount / 10), lines);
        builder = _builder;

        const a = Vec3(), b = Vec3();
        for (let s = 0, sl = streamlines.length; s < sl; ++s) {
            const l = streamlines[s];
            if (!streamlinePassesFilter(l, gridToCartn, props)) continue;

            if (dashEnabled) {
                Vec3.transformMat4(a, l[0], gridToCartn);
                for (let i = 1, il = l.length; i < il; ++i) {
                    Vec3.transformMat4(b, l[i], gridToCartn);
                    const inCycle = (i - 1) % (dashPoints * 2);
                    if (dashShift ? (inCycle >= dashPoints) : (inCycle < dashPoints)) {
                        _builder.addVec(a, b, s);
                    }
                    Vec3.copy(a, b);
                }
            } else {
                Vec3.transformMat4(a, l[0], gridToCartn);
                for (let i = 1, il = l.length; i < il; ++i) {
                    Vec3.transformMat4(b, l[i], gridToCartn);
                    _builder.addVec(a, b, s);
                    Vec3.copy(a, b);
                }
            }
        }
    }

    const result = builder.getLines();
    result.setBoundingSphere(Sphere3D.fromDimensionsAndTransform(Sphere3D(), gridDimension, gridToCartn));
    return result;
}

//

export const StreamlinesTubeMeshParams = {
    ...CommonStreamlinesParams,
    ...Mesh.Params,
    tubeSizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    radialSegments: PD.Numeric(8, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
};
export type StreamlinesTubeMeshParams = typeof StreamlinesTubeMeshParams
export type StreamlinesTubeMeshProps = PD.Values<StreamlinesTubeMeshParams>

export function VolumeStreamlinesTubeMeshVisual(materialId: number): VolumeVisual<StreamlinesTubeMeshParams> {
    return VolumeVisual<Mesh, StreamlinesTubeMeshParams>({
        defaultProps: PD.getDefaultValues(StreamlinesTubeMeshParams),
        createGeometry: createVolumeStreamlinesTubeMesh,
        createLocationIterator: createStreamlinesLocationIterator,
        getLoci: getStreamlinesLoci,
        eachLocation: eachStreamlines,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<StreamlinesTubeMeshParams>, currentProps: PD.Values<StreamlinesTubeMeshParams>) => {
            const streamlinesHash = StreamlinesProvider.get(volume).version;
            if ((state.info.streamlinesHash as number) !== streamlinesHash) {
                if (state.info.streamlinesHash !== undefined) {
                    state.createGeometry = true;
                    state.updateLocation = true;
                }
                state.info.streamlinesHash = streamlinesHash;
            }
            if (newProps.tubeSizeFactor !== currentProps.tubeSizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments) {
                state.createGeometry = true;
            }
            if (newProps.anchorEnabled !== currentProps.anchorEnabled ||
                !Vec3.equals(newProps.anchorCenter, currentProps.anchorCenter) ||
                newProps.anchorRadius !== currentProps.anchorRadius) {
                state.createGeometry = true;
            }
            if (newProps.dashEnabled !== currentProps.dashEnabled ||
                newProps.dashPoints !== currentProps.dashPoints ||
                newProps.dashShift !== currentProps.dashShift) {
                state.createGeometry = true;
            }
        },
        initUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<StreamlinesTubeMeshParams>, newTheme: Theme) => {
            const streamlinesHash = StreamlinesProvider.get(volume).version;
            state.info.streamlinesHash = streamlinesHash;
        },
        geometryUtils: Mesh.Utils,
    }, materialId);
}

function createVolumeStreamlinesTubeMesh(ctx: VisualContext, volume: Volume, _key: number, theme: Theme, props: StreamlinesTubeMeshProps, mesh?: Mesh): Mesh {
    const { cells: { space } } = volume.grid;
    const gridDimension = space.dimensions as Vec3;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const streamlines = StreamlinesProvider.get(volume).value!;

    const { tubeSizeFactor, radialSegments, dashEnabled, dashPoints, dashShift } = props;

    // Estimate vertex count
    const pointCount = streamlinePointCount(streamlines);
    const vertexCount = pointCount * radialSegments * 2;
    const builderState = MeshBuilder.createState(vertexCount, Math.ceil(vertexCount / 10), mesh);

    const location = StreamlinesLocation(streamlines, volume);
    for (let s = 0, sl = streamlines.length; s < sl; ++s) {
        const l = streamlines[s];
        if (!streamlinePassesFilter(l, gridToCartn, props)) continue;

        builderState.currentGroup = s;
        location.element.index = s as StreamlinesIndex;
        const tubeSize = theme.size.size(location) * tubeSizeFactor;

        if (dashEnabled) {
            addDashedStreamlineTube(builderState, l, gridToCartn, radialSegments, tubeSize, dashPoints, dashShift);
        } else {
            addStreamlineTube(builderState, l, gridToCartn, radialSegments, tubeSize);
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    m.setBoundingSphere(Sphere3D.fromDimensionsAndTransform(Sphere3D(), gridDimension, gridToCartn));
    return m;
}

/**
 * Add a tube along a streamline path
 */
function addStreamlineTube(state: MeshBuilder.State, streamline: Streamline, gridToCartn: Mat4, radialSegments: number, tubeSize: number) {
    const n = streamline.length;
    if (n < 2) return;

    const linearSegments = n - 1;
    const curvePoints = new Float32Array(n * 3);
    const normalVectors = new Float32Array(n * 3);
    const binormalVectors = new Float32Array(n * 3);
    const widthValues = new Float32Array(n);
    const heightValues = new Float32Array(n);

    for (let i = 0; i < n; ++i) {
        const p = streamline[i];
        v3transformMat4Offset(curvePoints, p, gridToCartn, i * 3, 0, 0);
        widthValues[i] = tubeSize;
        heightValues[i] = tubeSize;
    }

    computeFrenetFrames(curvePoints, normalVectors, binormalVectors, n);
    addTube(state, curvePoints, normalVectors, binormalVectors, linearSegments, radialSegments, widthValues, heightValues, true, true, 'elliptical');
}

/**
 * Add dashed tube segments along a streamline path.
 * Uses point count to determine dash/gap boundaries.
 */
function addDashedStreamlineTube(state: MeshBuilder.State, streamline: Streamline, gridToCartn: Mat4, radialSegments: number, tubeSize: number, dashPoints: number, dashShift: boolean) {
    const n = streamline.length;
    if (n < 2) return;

    const allPoints: Vec3[] = [];
    for (let i = 0; i < n; ++i) {
        allPoints.push(v3transformMat4(Vec3(), streamline[i], gridToCartn));
    }

    const cycleLength = dashPoints * 2;
    let i = dashShift ? dashPoints : 0;
    while (i < n - 1) {
        const dashStart = i;
        const dashEnd = Math.min(i + dashPoints, n - 1);
        if (dashEnd > dashStart) {
            emitTubeSegment(state, allPoints, dashStart, dashEnd, radialSegments, tubeSize);
        }
        i += cycleLength;
    }
}

/**
 * Emit a tube segment from startIdx to endIdx (inclusive) with caps on both ends.
 */
function emitTubeSegment(state: MeshBuilder.State, points: Vec3[], startIdx: number, endIdx: number, radialSegments: number, tubeSize: number) {
    const segmentLength = endIdx - startIdx + 1;
    if (segmentLength < 2) return;

    const curvePoints = new Float32Array(segmentLength * 3);
    const normalVectors = new Float32Array(segmentLength * 3);
    const binormalVectors = new Float32Array(segmentLength * 3);
    const widthValues = new Float32Array(segmentLength);
    const heightValues = new Float32Array(segmentLength);

    for (let i = 0; i < segmentLength; ++i) {
        v3toArray(points[startIdx + i], curvePoints, i * 3);
        widthValues[i] = tubeSize;
        heightValues[i] = tubeSize;
    }

    computeFrenetFrames(curvePoints, normalVectors, binormalVectors, segmentLength);
    addTube(state, curvePoints, normalVectors, binormalVectors, segmentLength - 1, radialSegments, widthValues, heightValues, true, true, 'elliptical');
}

//

const StreamlinesVisuals = {
    'lines': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, StreamlinesLinesParams>) => VolumeRepresentation('Streamlines lines', ctx, getParams, VolumeStreamlinesLinesVisual, getStreamlinesVisualLoci),
    'tube-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, StreamlinesTubeMeshParams>) => VolumeRepresentation('Streamlines tube-mesh', ctx, getParams, VolumeStreamlinesTubeMeshVisual, getStreamlinesVisualLoci),
};

export const StreamlinesParams = {
    ...StreamlinesLinesParams,
    ...StreamlinesTubeMeshParams,
    visuals: PD.MultiSelect(['lines'], PD.objectToOptions(StreamlinesVisuals)),
    density: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
};
export type StreamlinesParams = typeof StreamlinesParams;

export function getStreamlinesParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(StreamlinesParams);
    return p;
}

export type StreamlinesRepresentation = VolumeRepresentation<StreamlinesParams>
export function StreamlinesRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, StreamlinesParams>): StreamlinesRepresentation {
    return Representation.createMulti('Streamlines', ctx, getParams, Representation.StateBuilder, StreamlinesVisuals as unknown as Representation.Def<Volume, StreamlinesParams>);
}

export const StreamlinesRepresentationProvider = VolumeRepresentationProvider({
    name: 'streamlines',
    label: 'Streamlines',
    description: 'Displays streamlines.',
    factory: StreamlinesRepresentation,
    getParams: getStreamlinesParams,
    defaultValues: PD.getDefaultValues(StreamlinesParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, volume: Volume) => StreamlinesProvider.attach(ctx, volume, void 0, true),
        detach: (data) => StreamlinesProvider.ref(data, false)
    },
});
