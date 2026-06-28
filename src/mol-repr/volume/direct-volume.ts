/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Vec2, Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { Box3D } from '../../mol-math/geometry';
import { Grid, Volume } from '../../mol-model/volume';
import { RuntimeContext } from '../../mol-task';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { DirectVolume } from '../../mol-geo/geometry/direct-volume/direct-volume';
import { VisualContext } from '../visual';
import { Theme, ThemeRegistryContext } from '../../mol-theme/theme';
import { VolumeRepresentation, VolumeRepresentationProvider } from './representation';
import { VisualUpdateState } from '../util';
import { RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Loci, EmptyLoci } from '../../mol-model/loci';
import { PickingId } from '../../mol-geo/geometry/picking';
import { createVolumeCellLocationIterator, createVolumeTexture2d, createVolumeTexture3d, eachVolumeLoci, getVolumeTexture2dLayout } from './util';
import { Texture } from '../../mol-gl/webgl/texture';
import { Interval } from '../../mol-data/int/interval';
import { OrderedSet } from '../../mol-data/int/ordered-set';
import { VolumeVisual } from './visual';
import { clamp } from '../../mol-math/interpolate';

function getBoundingBox(gridDimension: Vec3, transform: Mat4) {
    const bbox = Box3D();
    Box3D.add(bbox, gridDimension);
    Box3D.transform(bbox, bbox, transform);
    return bbox;
}

// 2d volume texture

export function createDirectVolume2d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, props: PD.Values<DirectVolumeParams>, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
    const { width, height } = getVolumeTexture2dLayout(gridDimension);
    if (Math.max(width, height) > webgl.maxTextureSize / 2) {
        throw new Error('volume too large for direct-volume rendering');
    }

    const dataType = props.dataType === 'halfFloat' && !webgl.extensions.textureHalfFloat ? 'float' : props.dataType;

    const textureImage = createVolumeTexture2d(volume, 'normals', 0, dataType);
    // debugTexture(createImageData(textureImage.array, textureImage.width, textureImage.height), 1/3)
    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const bbox = getBoundingBox(gridDimension, transform);

    let texture: Texture;
    if (directVolume && directVolume.dataType.ref.value === dataType) {
        texture = directVolume.gridTexture.ref.value;
    } else {
        texture = dataType === 'byte'
            ? webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear')
            : dataType === 'halfFloat'
                ? webgl.resources.texture('image-float16', 'rgba', 'fp16', 'linear')
                : webgl.resources.texture('image-float32', 'rgba', 'float', 'linear');
    }
    texture.load(textureImage);

    const { unitToCartn, cellDim } = getUnitToCartn(volume.grid);
    const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
    return DirectVolume.create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, volume.grid.stats, false, axisOrder, dataType, directVolume);
}

// 3d volume texture

function getUnitToCartn(grid: Grid) {
    if (grid.transform.kind === 'matrix') {
        return {
            unitToCartn: Mat4.mul(Mat4(),
                grid.transform.matrix,
                Mat4.fromScaling(Mat4(), grid.cells.space.dimensions as Vec3)
            ),
            cellDim: Mat4.getScaling(Vec3(), grid.transform.matrix)
        };
    }
    const box = grid.transform.fractionalBox;
    const size = Box3D.size(Vec3(), box);
    return {
        unitToCartn: Mat4.mul3(Mat4(),
            grid.transform.cell.fromFractional,
            Mat4.fromTranslation(Mat4(), box.min),
            Mat4.fromScaling(Mat4(), size)
        ),
        cellDim: Vec3.div(Vec3(), grid.transform.cell.size, grid.cells.space.dimensions as Vec3)
    };
}

export function createDirectVolume3d(ctx: RuntimeContext, webgl: WebGLContext, volume: Volume, props: PD.Values<DirectVolumeParams>, directVolume?: DirectVolume) {
    const gridDimension = volume.grid.cells.space.dimensions as Vec3;
    if (Math.max(...gridDimension) > webgl.max3dTextureSize / 2) {
        throw new Error('volume too large for direct-volume rendering');
    }

    const dataType = props.dataType === 'halfFloat' && !webgl.extensions.textureHalfFloat ? 'float' : props.dataType;

    const textureVolume = createVolumeTexture3d(volume, dataType);
    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const bbox = getBoundingBox(gridDimension, transform);

    let texture: Texture;
    if (directVolume && directVolume.dataType.ref.value === dataType) {
        texture = directVolume.gridTexture.ref.value;
    } else {
        texture = dataType === 'byte'
            ? webgl.resources.texture('volume-uint8', 'rgba', 'ubyte', 'linear')
            : dataType === 'halfFloat'
                ? webgl.resources.texture('volume-float16', 'rgba', 'fp16', 'linear')
                : webgl.resources.texture('volume-float32', 'rgba', 'float', 'linear');
    }
    texture.load(textureVolume);

    const { unitToCartn, cellDim } = getUnitToCartn(volume.grid);
    const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
    return DirectVolume.create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, volume.grid.stats, false, axisOrder, dataType, directVolume);
}

//

export async function createDirectVolume(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: PD.Values<DirectVolumeParams>, directVolume?: DirectVolume) {
    const { runtime, webgl } = ctx;
    if (webgl === undefined) throw new Error('DirectVolumeVisual requires `webgl` in VisualContext');

    return webgl.isWebGL2 ?
        createDirectVolume3d(runtime, webgl, volume, props, directVolume) :
        createDirectVolume2d(runtime, webgl, volume, props, directVolume);
}

function getLoci(volume: Volume, props: PD.Values<DirectVolumeParams>) {
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return Volume.Loci(volume, instances);
}

export function getDirectVolumeLoci(pickingId: PickingId, volume: Volume, key: number, props: DirectVolumeProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;
    if (id === objectId) {
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
        return Volume.Cell.Loci(volume, [{ indices, instances }]);
    }
    return EmptyLoci;
}

export function eachDirectVolume(loci: Loci, volume: Volume, key: number, props: DirectVolumeProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, undefined, apply);
}

//

export const DirectVolumeParams = {
    ...DirectVolume.Params,
    quality: { ...DirectVolume.Params.quality, isEssential: false },
    dataType: PD.Select('byte', PD.arrayToOptions(['byte', 'float', 'halfFloat'] as const)),
};
export type DirectVolumeParams = typeof DirectVolumeParams
export function getDirectVolumeParams(ctx: ThemeRegistryContext, volume: Volume) {
    const params = PD.clone(DirectVolumeParams);
    params.controlPoints.getVolume = () => volume;
    params.controlPoints.defaultValue = computeRampControlPoints(volume.grid);
    params.controlPoints.getPresets = () => buildPresets(volume.grid);
    return params;
}
export type DirectVolumeProps = PD.Values<DirectVolumeParams>

export function DirectVolumeVisual(materialId: number): VolumeVisual<DirectVolumeParams> {
    return VolumeVisual<DirectVolume, DirectVolumeParams>({
        defaultProps: PD.getDefaultValues(DirectVolumeParams),
        createGeometry: createDirectVolume,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getDirectVolumeLoci,
        eachLocation: eachDirectVolume,
        setUpdateState: (state: VisualUpdateState, newVolume: Volume, currentVolume: Volume, newProps: PD.Values<DirectVolumeParams>, currentProps: PD.Values<DirectVolumeParams>) => {
            state.createGeometry = newProps.dataType !== currentProps.dataType;
        },
        geometryUtils: DirectVolume.Utils,
        dispose: (geometry: DirectVolume) => {
            geometry.gridTexture.ref.value.destroy();
        },
    }, materialId);
}

export function DirectVolumeRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, DirectVolumeParams>): VolumeRepresentation<DirectVolumeParams> {
    return VolumeRepresentation('Direct Volume', ctx, getParams, DirectVolumeVisual, getLoci);
}

export const DirectVolumeRepresentationProvider = VolumeRepresentationProvider({
    name: 'direct-volume',
    label: 'Direct Volume',
    description: 'Direct rendering of volumetric data.',
    factory: DirectVolumeRepresentation,
    getParams: getDirectVolumeParams,
    defaultValues: PD.getDefaultValues(DirectVolumeParams),
    defaultColorTheme: { name: 'volume-value' },
    defaultSizeTheme: { name: 'uniform' },
    locationKinds: ['position-location', 'direct-location'],
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume)
});

//

const PeakAlpha = 0.05;
const PeakRamp = 0.005;
const MinPeakX = 0.05;
const MaxPeakX = 0.95;

/**
 * Compute data-aware default control points from the volume's histogram.
 *
 * Strategy:
 *   - signed volumes (min < 0): symmetric ramps on both sides of 0
 *     (−3σ → −σ on the negative side, +σ → +3σ on the positive side).
 *   - unsigned volumes: a smooth ramp from μ → μ+σ → min(p99, μ+3σ)
 *
 * Robust stats are computed from the cached histogram with the bin
 * containing 0 excluded (typical for masked/padded density maps).
 */
function computeRampControlPoints(grid: Grid): Vec2[] {
    const { min, max } = grid.stats;
    const range = max - min;
    if (range <= 0 || !isFinite(range)) {
        return fallbackControlPoints();
    }

    const robust = Grid.getRobustStats(grid, { ignoreZero: true });
    const { mean, sigma, p99 } = robust;
    if (!isFinite(sigma) || sigma <= 0 || robust.count === 0) {
        return fallbackControlPoints();
    }

    const toX = (v: number) => clamp((v - min) / range, 0.001, 0.999);

    if (min < 0) {
        // Signed volume: mirrored ramps around 0.
        const points: Vec2[] = [];
        const negOuter = toX(Math.min(mean - 3 * sigma, -sigma));
        const negInner = toX(-sigma);
        const zero = toX(0);
        const posInner = toX(sigma);
        const posOuter = toX(Math.max(mean + 3 * sigma, sigma));
        const xs = [negOuter, negInner, zero, posInner, posOuter].sort((a, b) => a - b);
        // Negative ramp: alpha PeakAlpha at outer, falls to 0 at zero.
        points.push(Vec2.create(xs[0], PeakAlpha));
        points.push(Vec2.create(xs[1], PeakAlpha * 0.5));
        points.push(Vec2.create(xs[2], 0));
        // Positive ramp: 0 at zero, rises to PeakAlpha at outer.
        points.push(Vec2.create(xs[3], PeakAlpha * 0.5));
        points.push(Vec2.create(xs[4], PeakAlpha));
        return points;
    }

    // Unsigned: ramp from μ up to ~p99 (capped at μ+3σ).
    const xStart = toX(mean);
    const xMid = toX(mean + sigma);
    const xEnd = toX(Math.min(mean + 3 * sigma, p99));
    const ramp = buildRampPoints(xStart, xMid, xEnd, PeakAlpha);

    return ramp;
}

/**
 * Build a ramp described by three monotonic x positions:
 *   xStart -> alpha 0
 *   xMid   -> alpha maxAlpha/2
 *   xEnd   -> alpha maxAlpha (then plateau briefly and close to 0)
 */
function buildRampPoints(xStart: number, xMid: number, xEnd: number, maxAlpha: number): Vec2[] {
    const eps = 0.002;
    const a = clamp(xStart, 0.001, 0.997);
    let b = clamp(xMid, a + eps, 0.998);
    let c = clamp(xEnd, b + eps, 0.999);
    // Avoid collapsed ramps when sigma is tiny relative to range.
    if (c - a < 3 * eps) {
        b = a + eps;
        c = a + 2 * eps;
    }
    const points: Vec2[] = [];
    points.push(Vec2.create(a, 0));
    points.push(Vec2.create(b, maxAlpha * 0.5));
    points.push(Vec2.create(c, maxAlpha));
    // Plateau plus close so the upper end of the volume range stays opaque
    // until very near the top of the data range.
    const dEnd = clamp(c + eps, c + 0.0001, 0.9995);
    points.push(Vec2.create(dEnd, maxAlpha));
    return points;
}

function fallbackControlPoints(): Vec2[] {
    return [
        Vec2.create(0.19, 0.0), Vec2.create(0.2, PeakAlpha), Vec2.create(0.25, PeakAlpha), Vec2.create(0.26, 0.0),
        Vec2.create(0.79, 0.0), Vec2.create(0.8, PeakAlpha), Vec2.create(0.85, PeakAlpha), Vec2.create(0.86, 0.0),
    ];
}

/**
 * Build a single narrow peak (iso-surface-like) at `mean + sigmaOffset*sigma`.
 * `halfWidth` is in normalized [0,1] x-space; `alpha` is the peak height.
 */
function buildSinglePeak(grid: Grid, sigmaOffset: number, halfWidth: number, alpha: number): Vec2[] {
    const { min, max } = grid.stats;
    const range = max - min;
    if (range <= 0 || !isFinite(range)) {
        return fallbackControlPoints();
    }
    const robust = Grid.getRobustStats(grid, { ignoreZero: true });
    const { mean, sigma } = robust;
    if (!isFinite(sigma) || sigma <= 0 || robust.count === 0) {
        return fallbackControlPoints();
    }
    const cx = clamp((mean + sigmaOffset * sigma - min) / range, MinPeakX, MaxPeakX);
    const dx = clamp(halfWidth, PeakRamp + 0.001, 0.2);
    const x0 = clamp(cx - dx, 0.001, 0.999);
    const x1 = clamp(cx - dx + PeakRamp, 0.001, 0.999);
    const x2 = clamp(cx + dx - PeakRamp, 0.001, 0.999);
    const x3 = clamp(cx + dx, 0.001, 0.999);
    if (!(x0 < x1 && x1 < x2 && x2 < x3)) return fallbackControlPoints();
    return [
        Vec2.create(x0, 0),
        Vec2.create(x1, alpha),
        Vec2.create(x2, alpha),
        Vec2.create(x3, 0),
    ];
}

/** Sharp iso-surface-like peak around mean+2σ. */
function computeSharpSurfaceControlPoints(grid: Grid): Vec2[] {
    return buildSinglePeak(grid, 2, 0.012, PeakAlpha);
}

/** Narrow high-contour peak at mean+3σ with stronger alpha. */
function computeHighContourControlPoints(grid: Grid): Vec2[] {
    return buildSinglePeak(grid, 3, 0.015, 0.08);
}

/**
 * Wide soft ramp starting at the mean, plateauing around mean+σ; reveals
 * weak density (low-resolution maps, partial occupancy).
 */
function computeLowContourControlPoints(grid: Grid): Vec2[] {
    const { min, max } = grid.stats;
    const range = max - min;
    if (range <= 0 || !isFinite(range)) {
        return fallbackControlPoints();
    }
    const robust = Grid.getRobustStats(grid, { ignoreZero: true });
    const { mean, sigma } = robust;
    if (!isFinite(sigma) || sigma <= 0 || robust.count === 0) {
        return fallbackControlPoints();
    }
    const toX = (v: number) => clamp((v - min) / range, 0.001, 0.999);
    const xStart = toX(mean - 0.5 * sigma);
    const xMid = toX(mean + 0.5 * sigma);
    const xEnd = toX(mean + 1.5 * sigma);
    return buildRampPoints(xStart, xMid, xEnd, PeakAlpha);
}

/**
 * Tomogram-style preset: emphasizes low values (typical of cryo-ET where
 * matter is darker than background prior to inversion). Builds a ramp
 * descending from p1 → mean−σ.
 */
function computeTomogramControlPoints(grid: Grid): Vec2[] {
    const { min, max } = grid.stats;
    const range = max - min;
    if (range <= 0 || !isFinite(range)) return fallbackControlPoints();

    const robust = Grid.getRobustStats(grid, { ignoreZero: true });
    const { mean, sigma, p1 } = robust;
    if (!isFinite(sigma) || sigma <= 0 || robust.count === 0) return fallbackControlPoints();

    const toX = (v: number) => clamp((v - min) / range, 0.001, 0.999);
    // Mirror of the unsigned ramp, descending from low values.
    const xStart = toX(Math.max(mean - 3 * sigma, p1));
    const xMid = toX(mean - sigma);
    const xEnd = toX(mean);
    const a = clamp(xStart, 0.001, 0.997);
    const b = clamp(xMid, a + 0.002, 0.998);
    const c = clamp(xEnd, b + 0.002, 0.999);
    return [
        Vec2.create(a, PeakAlpha),
        Vec2.create(b, PeakAlpha * 0.5),
        Vec2.create(c, 0),
    ];
}

/**
 * Symmetric narrow peaks at ±3σ for signed difference maps (Fo−Fc, mFo−DFc).
 * Falls back to the default control points when the volume is unsigned.
 */
function computeDifferenceMapControlPoints(grid: Grid): Vec2[] {
    const { min, max } = grid.stats;
    const range = max - min;
    if (range <= 0 || !isFinite(range) || min >= 0) {
        return fallbackControlPoints();
    }
    const robust = Grid.getRobustStats(grid, { ignoreZero: true });
    const { mean, sigma } = robust;
    if (!isFinite(sigma) || sigma <= 0 || robust.count === 0) {
        return fallbackControlPoints();
    }
    const toX = (v: number) => clamp((v - min) / range, MinPeakX, MaxPeakX);
    const xs = [toX(mean - 3 * sigma), toX(mean + 3 * sigma)];
    xs.sort((a, b) => a - b);

    const points: Vec2[] = [];
    let dx = 0.012;
    if (xs[1] - xs[0] < 2 * (dx + PeakRamp)) {
        dx = Math.max(PeakRamp + 0.001, (xs[1] - xs[0]) / 2 - PeakRamp);
    }
    for (const cx of xs) {
        const x = clamp(cx, MinPeakX, MaxPeakX);
        const x0 = clamp(x - dx, 0.001, 0.999);
        const x1 = clamp(x - dx + PeakRamp, 0.001, 0.999);
        const x2 = clamp(x + dx - PeakRamp, 0.001, 0.999);
        const x3 = clamp(x + dx, 0.001, 0.999);
        points.push(Vec2.create(x0, 0));
        points.push(Vec2.create(x1, PeakAlpha));
        points.push(Vec2.create(x2, PeakAlpha));
        points.push(Vec2.create(x3, 0));
    }
    return points;
}

/**
 * Linear ramp from p1 → p99 with monotonically increasing alpha. Suited for
 * data normalized to [0,1] like occupancy or probability/mask volumes.
 */
function computeOccupancyControlPoints(grid: Grid): Vec2[] {
    const { min, max } = grid.stats;
    const range = max - min;
    if (range <= 0 || !isFinite(range)) return fallbackControlPoints();

    const robust = Grid.getRobustStats(grid, { ignoreZero: true });
    const { p1, p99 } = robust;
    if (!isFinite(p99 - p1) || p99 <= p1) return fallbackControlPoints();

    const toX = (v: number) => clamp((v - min) / range, 0.001, 0.999);
    const a = toX(p1);
    const c = toX(p99);
    const b = clamp((a + c) * 0.5, a + 0.002, c - 0.002);
    const points: Vec2[] = [
        Vec2.create(a, 0),
        Vec2.create(b, PeakAlpha * 0.5),
        Vec2.create(c, PeakAlpha),
    ];
    // Hold alpha near the top of the range so the most-occupied voxels stay
    // visible instead of fading out above p99.
    if (c < 0.997) points.push(Vec2.create(0.999, PeakAlpha));
    return points;
}

/**
 * Soft, low-alpha ramp covering [μ, p99] for cloudy/accumulation-style
 * rendering (MO densities, mesoscale fields).
 */
function computeWideVolumetricControlPoints(grid: Grid): Vec2[] {
    const { min, max } = grid.stats;
    const range = max - min;
    if (range <= 0 || !isFinite(range)) return fallbackControlPoints();

    const robust = Grid.getRobustStats(grid, { ignoreZero: true });
    const { mean, sigma, p99 } = robust;
    if (!isFinite(sigma) || sigma <= 0) return fallbackControlPoints();

    const toX = (v: number) => clamp((v - min) / range, 0.001, 0.999);
    const xStart = toX(Math.max(mean - sigma, min));
    const xMid = toX(mean + sigma);
    const xEnd = toX(p99);
    return buildRampPoints(xStart, xMid, xEnd, 0.02);
}

/**
 * Build the list of presets exposed via the LineGraph control. Presets that
 * don't apply to the current volume (e.g. difference-map preset on an
 * unsigned volume) are filtered out rather than producing degenerate ramps.
 */
function buildPresets(grid: Grid): [Vec2[], string][] {
    const cached = (grid as any)._directVolumeControlPointsPresets as [Vec2[], string][] | undefined;
    if (cached) return cached;
    const { min, max } = grid.stats;
    const isSigned = min < 0;
    const presets: [Vec2[], string][] = [
        [computeRampControlPoints(grid), 'Ramp (auto)'],
        [computeSharpSurfaceControlPoints(grid), 'Sharp surface (~+2σ)'],
        [computeHighContourControlPoints(grid), 'High contour (~+3σ)'],
        [computeLowContourControlPoints(grid), 'Low contour (~+1σ)'],
    ];
    // Tomogram preset: include when the data has meaningful negative spread
    // (signed maps, or unsigned but strongly left-skewed).
    if (isSigned || (isFinite(max - min) && grid.stats.mean - min > max - grid.stats.mean)) {
        presets.push([computeTomogramControlPoints(grid), 'Tomogram (inverted)']);
    }
    if (isSigned) {
        presets.push([computeDifferenceMapControlPoints(grid), 'Difference map (±3σ)']);
    }
    presets.push([computeOccupancyControlPoints(grid), 'Occupancy / probability']);
    presets.push([computeWideVolumetricControlPoints(grid), 'Wide volumetric']);
    (grid as any)._directVolumeControlPointsPresets = presets;
    return presets;
}
