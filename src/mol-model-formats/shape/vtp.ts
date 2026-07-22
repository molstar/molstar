/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { RuntimeContext, Task } from '../../mol-task';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { VtpFile } from '../../mol-io/reader/vtp/schema';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Shape } from '../../mol-model/shape';
import { ChunkedArray } from '../../mol-data/util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ColorNames } from '../../mol-util/color/names';
import { ColorScale } from '../../mol-util/color/scale';
import { ColorListOptionsScale, ColorListName } from '../../mol-util/color/lists';
import { ValueCell } from '../../mol-util/value-cell';
import { deepClone } from '../../mol-util/object';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Grid, Volume } from '../../mol-model/volume';

export interface VtpData {
    source: VtpFile,
    transforms?: Mat4[],
}

/**
 * Enumerates the volumes selectable by the `volume` color option, e.g. from the plugin state.
 * Injected by the caller so this format module stays free of any plugin-layer dependency.
 */
export type VolumeRefOptions = PD.ValueRef<Volume>['getOptions'];

const _identityTransforms = [Mat4.identity()];

// --- Params ---

function fmtStat(x: number): string {
    if (!isFinite(x)) return x.toString();
    const abs = Math.abs(x);
    if (abs === 0) return '0';
    if (abs < 0.001 || abs >= 1e5) return x.toExponential(2);
    return x.toPrecision(3);
}

function computeAttrStatsText(vtpFile: VtpFile, attrKey: string): string {
    let rawVals: ArrayLike<number> | undefined;
    let nComp = 1;
    let name = '';
    if (attrKey.startsWith('point:')) {
        name = attrKey.slice(6);
        const arr = vtpFile.pointData.get(name);
        if (!arr) return '';
        rawVals = arr.values.toArray(); nComp = arr.numberOfComponents;
    } else if (attrKey.startsWith('cell:')) {
        name = attrKey.slice(5);
        const arr = vtpFile.cellData.get(name);
        if (!arr) return '';
        rawVals = arr.values.toArray(); nComp = arr.numberOfComponents;
    }
    if (!rawVals) return '';
    const n = rawVals.length / nComp;
    let min = Infinity, max = -Infinity, sum = 0, sumSq = 0;
    for (let i = 0; i < n; i++) {
        let v: number;
        if (nComp === 1) {
            v = rawVals[i];
        } else {
            let mag2 = 0;
            for (let c = 0; c < nComp; c++) mag2 += rawVals[i * nComp + c] ** 2;
            v = Math.sqrt(mag2);
        }
        if (v < min) min = v;
        if (v > max) max = v;
        sum += v; sumSq += v * v;
    }
    const mean = sum / n;
    const std = Math.sqrt(Math.max(0, sumSq / n - mean * mean));
    const typeStr = nComp > 1 ? `${nComp}-comp magnitude` : 'scalar';
    return `${name}  (${typeStr}, n=${n})\nmin: ${fmtStat(min)}   max: ${fmtStat(max)}\nmean: ${fmtStat(mean)}   σ: ${fmtStat(std)}`;
}

function buildAttrOptions(vtpFile?: VtpFile): PD.SelectOption<string>[] {
    if (!vtpFile) return [['', 'No attributes']];
    const opts: PD.SelectOption<string>[] = [];
    for (const [name, arr] of vtpFile.pointData) {
        const nComp = arr.numberOfComponents;
        opts.push([`point:${name}`, nComp > 1 ? `Point: ${name} (mag)` : `Point: ${name}`]);
    }
    for (const [name, arr] of vtpFile.cellData) {
        const nComp = arr.numberOfComponents;
        opts.push([`cell:${name}`, nComp > 1 ? `Cell: ${name} (mag)` : `Cell: ${name}`]);
    }
    return opts.length > 0 ? opts : [['', 'No attributes']];
}

function defaultAttrKey(vtpFile?: VtpFile): string {
    if (!vtpFile) return '';
    // Prefer 1-component arrays first (direct scalar mapping)
    for (const [name, arr] of vtpFile.cellData) {
        if (arr.numberOfComponents === 1) return `cell:${name}`;
    }
    for (const [name, arr] of vtpFile.pointData) {
        if (arr.numberOfComponents === 1) return `point:${name}`;
    }
    if (vtpFile.cellData.size > 0) return `cell:${vtpFile.cellData.keys().next().value}`;
    if (vtpFile.pointData.size > 0) return `point:${vtpFile.pointData.keys().next().value}`;
    return '';
}

function scalarRange(values: ArrayLike<number>): [number, number] {
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < values.length; i++) {
        if (values[i] < min) min = values[i];
        if (values[i] > max) max = values[i];
    }
    return [min === Infinity ? 0 : min, max === -Infinity ? 1 : max];
}

/**
 * Robust domain derived from the value distribution: keep the central `coveragePct`% of values and
 * clip the extreme tails, so a handful of outliers don't wash out the color range. This reproduces
 * the tighter min/max that ParaView's "separate color map" produces (max pulled closer to the mean).
 * Uses linearly-interpolated percentiles on a sorted copy of the values.
 */
function distributionRange(values: ArrayLike<number>, coveragePct: number): [number, number] {
    // Drop non-finite values so a single NaN (e.g. masked/no-data vertices) can't corrupt the sort;
    // keeps parity with the NaN-robust `auto`/scalarRange path.
    let m = 0;
    for (let i = 0; i < values.length; i++) if (Number.isFinite(values[i])) m++;
    if (m === 0) return [0, 1];
    const sorted = new Float64Array(m);
    for (let i = 0, k = 0; i < values.length; i++) if (Number.isFinite(values[i])) sorted[k++] = values[i];
    sorted.sort();
    const n = m;
    if (n === 1) return [sorted[0], sorted[0]];
    const tail = (100 - Math.max(0, Math.min(100, coveragePct))) / 200; // fraction dropped from each side
    const quantile = (f: number) => {
        const idx = f * (n - 1);
        const i = Math.floor(idx);
        return i + 1 < n ? sorted[i] + (sorted[i + 1] - sorted[i]) * (idx - i) : sorted[i];
    };
    return [quantile(tail), quantile(1 - tail)];
}

/**
 * IRAF/DS9 "zscale" auto range: subsample the values, fit a line to the sorted samples with
 * iterative k-sigma rejection, then derive [z1, z2] from the fitted slope around the median.
 * `contrast` (0.25 default) sets the stretch — lower spreads the mid-range over more of the color
 * scale. Reproduces the contrast astronomy viewers use for a dominant background plus faint structure.
 * Port of astropy's ZScaleInterval algorithm.
 */
function zscaleRange(values: ArrayLike<number>, contrast: number): [number, number] {
    const n = values.length;
    if (n === 0) return [0, 1];
    const nsamples = 1000;
    const stride = Math.max(1, Math.floor(n / nsamples));
    const sample: number[] = [];
    // skip non-finite values so masked/no-data vertices can't corrupt the sort/fit
    for (let i = 0; i < n && sample.length < nsamples; i += stride) if (Number.isFinite(values[i])) sample.push(values[i]);
    sample.sort((a, b) => a - b);
    const npix = sample.length;
    if (npix === 0) return [0, 1];
    const vmin = sample[0], vmax = sample[npix - 1];
    if (npix < 5) return [vmin, vmax];

    const krej = 2.5, maxIter = 5;
    const minpix = Math.max(5, Math.floor(npix * 0.5));
    const grow = Math.max(1, Math.floor(npix * 0.01) >> 1);
    const bad = new Uint8Array(npix);
    let ngood = npix, lastNgood = npix + 1;
    let slope = 0, intercept = 0;

    for (let iter = 0; iter < maxIter; iter++) {
        if (ngood >= lastNgood || ngood < minpix) break;
        // least-squares line fit y = intercept + slope*x over good pixels (x = sorted index)
        let s0 = 0, sx = 0, sxx = 0, sxy = 0, sy = 0;
        for (let i = 0; i < npix; i++) {
            if (bad[i]) continue;
            s0++; sx += i; sxx += i * i; sxy += i * sample[i]; sy += sample[i];
        }
        const denom = s0 * sxx - sx * sx;
        if (denom === 0) break;
        slope = (s0 * sxy - sx * sy) / denom;
        intercept = (sy - slope * sx) / s0;

        let sumSq = 0, cnt = 0;
        for (let i = 0; i < npix; i++) {
            if (bad[i]) continue;
            const r = sample[i] - (intercept + slope * i);
            sumSq += r * r; cnt++;
        }
        const threshold = krej * Math.sqrt(sumSq / Math.max(1, cnt));

        // reject points far from the fit, then dilate the rejection mask by `grow`
        const grown = new Uint8Array(npix);
        for (let i = 0; i < npix; i++) {
            const r = sample[i] - (intercept + slope * i);
            if (bad[i] || r < -threshold || r > threshold) {
                for (let j = Math.max(0, i - grow), jl = Math.min(npix - 1, i + grow); j <= jl; j++) grown[j] = 1;
            }
        }
        bad.set(grown);
        lastNgood = ngood;
        ngood = 0;
        for (let i = 0; i < npix; i++) if (!bad[i]) ngood++;
    }

    if (ngood < minpix) return [vmin, vmax];
    const s = contrast > 0 ? slope / contrast : slope;
    const center = (npix - 1) >> 1;
    const median = npix % 2 ? sample[(npix - 1) / 2] : (sample[npix / 2 - 1] + sample[npix / 2]) / 2;
    const z1 = Math.max(vmin, median - center * s);
    const z2 = Math.min(vmax, median + (npix - 1 - center) * s);
    return z2 > z1 ? [z1, z2] : [vmin, vmax];
}

type MappingName = 'linear' | 'log';

/**
 * Build the value→scale transform for a color mapping, plus the transformed domain endpoints.
 * The color scale is created over [dlo, dhi] and looked up with `f(value)`, so a monotonic curve
 * (log) reuses the linear ColorScale machinery.
 */
function makeMappingTransform(mapping: MappingName, lo: number, hi: number): { f: (v: number) => number, dlo: number, dhi: number } {
    if (mapping === 'log') {
        const safeLo = lo > 0 ? lo : (hi > 0 ? hi * 1e-6 : 1e-6);
        const safeHi = hi > safeLo ? hi : safeLo * 10;
        return { f: (v: number) => Math.log10(Math.max(v, safeLo)), dlo: Math.log10(safeLo), dhi: Math.log10(safeHi) };
    }
    return { f: (v: number) => v, dlo: lo, dhi: hi };
}


const DomainHelpText = 'auto: full [min, max]. custom: fixed range. best: robust range from the value distribution (clips outlier tails). zscale: IRAF/DS9 contrast auto-range.';

/** Color-scale params shared by the `attribute` and `volume` color options: a color list, a domain (limits) and a mapping (scale function). */
function scaleColoringParams(getUsedDomainText?: () => string) {
    return {
        colors: PD.Select('viridis' as ColorListName, ColorListOptionsScale),
        domain: PD.MappedStatic('auto', {
            custom: PD.Interval([-1, 1], { step: 0.001 }),
            auto: PD.Group({
                symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
            }),
            best: PD.Group({
                coverage: PD.Numeric(95, { min: 50, max: 100, step: 0.5 }, { description: 'Central percentage of values the domain covers; the most extreme values in the remaining tails are clipped so outliers do not wash out the color range.' }),
                symmetric: PD.Boolean(false, { description: 'If true the range is made symmetric as [-|max|, |max|].' })
            }),
            zscale: PD.Group({
                contrast: PD.Numeric(0.25, { min: 0.05, max: 1, step: 0.01 }, { description: 'IRAF/DS9 zscale contrast; lower values stretch the mid-range over more of the color scale.' })
            })
        }, { help: () => {
            const used = getUsedDomainText?.();
            return { description: used ? `${DomainHelpText}\n${used}` : DomainHelpText };
        } }),
        mapping: PD.Select<MappingName>('linear', [
            ['linear', 'Linear'],
            ['log', 'Log']
        ], { description: 'Scale function mapping values onto the color list. log needs a positive domain (non-positive values are clamped to the low bound).' })
    };
}

export function createVtpShapeParams(vtpFile?: VtpFile, getStats?: () => string, getVolumeOptions?: VolumeRefOptions, getUsedDomainText?: () => string) {
    const attrOptions = buildAttrOptions(vtpFile);
    const defKey = defaultAttrKey(vtpFile);
    const defStats = vtpFile && defKey ? computeAttrStatsText(vtpFile, defKey) : '';
    const hasAttrs = (vtpFile?.pointData.size ?? 0) + (vtpFile?.cellData.size ?? 0) > 0;

    return {
        ...Mesh.Params,
        doubleSided: { ...Mesh.Params.doubleSided, defaultValue: true },
        interior: { ...Mesh.Params.interior, defaultValue: { ...Mesh.Params.interior.defaultValue, colorStrength: 0 } },
        colorTheme: PD.MappedStatic(hasAttrs ? 'attribute' : 'uniform', {
            attribute: PD.Group({
                name: PD.Select(defKey, attrOptions, { help: () => {
                    return { description: getStats ? getStats() : defStats };
                } }),
                ...scaleColoringParams(getUsedDomainText)
            }),
            'external-volume': PD.Group({
                volume: PD.ValueRef<Volume>(getVolumeOptions ?? (() => []), (ref, getData) => getData(ref)),
                ...scaleColoringParams(getUsedDomainText)
            }, { label: 'External Volume', description: 'Color by trilinearly sampling a selected volume at each vertex.' }),
            uniform: PD.Group({
                color: PD.Color(ColorNames.grey, { label: 'Uniform Color' })
            })
        }, { isEssential: true }),
        scale: PD.Numeric(1, { min: 0.01, max: 100, step: 0.01 }, { label: 'Scale', isEssential: true, description: 'Uniform scale factor applied to the mesh.' }),
    };
}

export const VtpShapeParams = createVtpShapeParams();
export type VtpShapeParams = typeof VtpShapeParams;

// --- Mesh building ---

async function buildMesh(ctx: RuntimeContext, vtpFile: VtpFile, mesh?: Mesh): Promise<Mesh> {
    const { positions, connectivity, numberOfPoints, numberOfTriangles } = vtpFile;

    const builderState = MeshBuilder.createState(numberOfPoints, numberOfTriangles, mesh);
    const { vertices, indices, groups } = builderState;

    const chunkSize = 50000;

    for (let i = 0, il = numberOfPoints; i < il; i += chunkSize) {
        const end = Math.min(i + chunkSize, il);
        for (let v = i; v < end; v++) {
            ChunkedArray.add3(vertices, positions[3 * v], positions[3 * v + 1], positions[3 * v + 2]);
            ChunkedArray.add(groups, v);
        }
        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'Building VTP mesh...', current: i, max: numberOfPoints + numberOfTriangles });
        }
    }

    for (let i = 0, il = numberOfTriangles; i < il; i += chunkSize) {
        const end = Math.min(i + chunkSize, il);
        for (let t = i; t < end; t++) {
            ChunkedArray.add3(indices, connectivity[3 * t], connectivity[3 * t + 1], connectivity[3 * t + 2]);
        }
        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'Building VTP mesh...', current: numberOfPoints + i, max: numberOfPoints + numberOfTriangles });
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    // Requires consistent CCW winding. Files with alternating CW/CCW faces
    // (e.g. some lat/lon sphere exports) will show pinwheel shading artifacts.
    Mesh.computeNormals(m);
    ValueCell.updateIfChanged(m.varyingGroup, true);
    return m;
}

// --- Color lookup ---

function cellToVertexAverage(vtpFile: VtpFile, cellValues: ArrayLike<number>): Float64Array {
    const { connectivity, triangleCellIndex, numberOfPoints } = vtpFile;
    const nTris = connectivity.length / 3;
    const sum = new Float64Array(numberOfPoints);
    const count = new Uint32Array(numberOfPoints);
    for (let i = 0; i < nTris; i++) {
        const val = cellValues[triangleCellIndex[i]];
        const v0 = connectivity[3 * i];
        const v1 = connectivity[3 * i + 1];
        const v2 = connectivity[3 * i + 2];
        sum[v0] += val; count[v0]++;
        sum[v1] += val; count[v1]++;
        sum[v2] += val; count[v2]++;
    }
    const smoothed = new Float64Array(numberOfPoints);
    for (let v = 0; v < numberOfPoints; v++) {
        smoothed[v] = count[v] > 0 ? sum[v] / count[v] : 0;
    }
    return smoothed;
}

/**
 * Average multi-component cell vectors to vertices, then return the magnitude of each vertex vector.
 * Matching smgui _cell_to_vertex_interpolation_vector: averaging direction vectors, not their magnitudes,
 * so that crease/boundary vertices (where adjacent face normals cancel) get low magnitude values.
 */
function cellToVertexAverageMag(vtpFile: VtpFile, cellVectors: ArrayLike<number>, nComp: number): Float64Array {
    const { connectivity, triangleCellIndex, numberOfPoints } = vtpFile;
    const nTris = connectivity.length / 3;
    const sum = new Float64Array(numberOfPoints * nComp);
    const count = new Uint32Array(numberOfPoints);
    for (let i = 0; i < nTris; i++) {
        const ci = triangleCellIndex[i];
        const v0 = connectivity[3 * i];
        const v1 = connectivity[3 * i + 1];
        const v2 = connectivity[3 * i + 2];
        for (let c = 0; c < nComp; c++) {
            const val = cellVectors[ci * nComp + c];
            sum[v0 * nComp + c] += val;
            sum[v1 * nComp + c] += val;
            sum[v2 * nComp + c] += val;
        }
        count[v0]++; count[v1]++; count[v2]++;
    }
    const result = new Float64Array(numberOfPoints);
    for (let v = 0; v < numberOfPoints; v++) {
        if (count[v] > 0) {
            let mag2 = 0;
            for (let c = 0; c < nComp; c++) {
                const avg = sum[v * nComp + c] / count[v];
                mag2 += avg * avg;
            }
            result[v] = Math.sqrt(mag2);
        }
    }
    return result;
}

function vecMag(values: ArrayLike<number>, baseIdx: number, nComp: number): number {
    let mag2 = 0;
    for (let c = 0; c < nComp; c++) mag2 += values[baseIdx + c] ** 2;
    return Math.sqrt(mag2);
}

interface VertexResult {
    values: ArrayLike<number>;
    isMagnitude: boolean; // true for multi-component (magnitude) attributes
}

/**
 * Compute per-vertex scalar values from a VTP attribute, matching smgui exactly:
 *   - CellData scalar: average cell scalars to vertices
 *   - CellData multi-component: average raw vectors to vertices first, then compute magnitude
 *     (matches smgui _cell_to_vertex_interpolation_vector → norm; crease vertices cancel → low magnitude)
 *   - PointData scalar: use values directly
 *   - PointData multi-component: compute per-vertex magnitude directly
 * Returns null for unknown or missing attributes.
 */
function computeVertexValues(vtpFile: VtpFile, attribute: string): VertexResult | null {
    if (!attribute) return null;

    if (attribute.startsWith('cell:')) {
        const arr = vtpFile.cellData.get(attribute.slice(5));
        if (!arr) return null;
        const nComp = arr.numberOfComponents;
        const vals = arr.values.toArray();
        if (nComp === 1) {
            return { values: cellToVertexAverage(vtpFile, vals), isMagnitude: false };
        }
        // Average vectors to vertices first, then take magnitude — crease vertices get lower values
        return { values: cellToVertexAverageMag(vtpFile, vals, nComp), isMagnitude: true };
    }

    if (attribute.startsWith('point:')) {
        const arr = vtpFile.pointData.get(attribute.slice(6));
        if (!arr) return null;
        const nComp = arr.numberOfComponents;
        const vals = arr.values.toArray();
        if (nComp === 1) return { values: vals, isMagnitude: false };
        // multi-component point data: per-vertex magnitude
        const nPts = vals.length / nComp;
        const mag = new Float64Array(nPts);
        for (let i = 0; i < nPts; i++) mag[i] = vecMag(vals, i * nComp, nComp);
        return { values: mag, isMagnitude: true };
    }

    return null;
}

/** Resolve the volume behind a `ValueRef`; returns undefined for an unset ref or before reconciliation binds it. */
function resolveVolumeRef(ref: PD.ValueRef<Volume>['defaultValue']): Volume | undefined {
    try {
        return ref.getValue();
    } catch {
        return undefined;
    }
}

/**
 * Sample a volume at each vertex position via trilinear interpolation, reusing the same
 * `Grid.makeGetTrilinearlyInterpolated` primitive as `ExternalVolumeColorTheme`. The result flows
 * through the shared color-scale path in `makeColorFn`, just like a per-vertex attribute.
 *
 * Sampling assumes VTP positions share the volume's Cartesian frame — the same contract
 * ExternalVolumeColorTheme imposes on any surface mesh. Vertices outside the grid (NaN) are clamped
 * to 0 so the auto domain and color scale stay well-defined. The `scale` param is applied via instance
 * transforms only, so at scale !== 1 the sampled colors won't follow the displayed geometry.
 */
function computeVolumeVertexValues(vtpFile: VtpFile, volume: Volume): VertexResult {
    const { positions, numberOfPoints } = vtpFile;
    const getInterpolated = Grid.makeGetTrilinearlyInterpolated(volume.grid, 'none');
    const values = new Float64Array(numberOfPoints);
    const p = Vec3();
    for (let i = 0; i < numberOfPoints; i++) {
        Vec3.set(p, positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
        const v = getInterpolated(p);
        values[i] = Number.isNaN(v) ? 0 : v;
    }
    return { values, isMagnitude: false };
}

/**
 * Build a color function from pre-computed per-vertex scalar values (from a VTP attribute or a sampled volume).
 *   - domain `auto`: full [min, max] (matches smgui). Magnitude attributes whose values are nearly identical
 *     (e.g. PointData unit normals, range < 1% of max) fall back to [0, max].
 *   - domain `custom`: fixed [lo, hi].
 *   - domain `best`: robust range from the value distribution — clips outlier tails so a few extreme values
 *     don't wash out the color range (reproduces ParaView's tighter "separate color map" min/max).
 *   - `logScale`: maps values on a log10 scale within the domain, boosting low-end contrast for long-tailed data.
 */
interface ColorFnResult {
    color: (gid: number) => Color;
    /** Data-space [min, max] actually mapped to the color scale (before any log transform); undefined for uniform. */
    usedDomain?: [number, number];
}

function makeColorFn(colorTheme: PD.Values<VtpShapeParams>['colorTheme'], result: VertexResult | null): ColorFnResult {
    if (colorTheme.name === 'uniform') return { color: () => colorTheme.params.color };
    if (!result) return { color: () => ColorNames.grey };

    const { values, isMagnitude } = result;
    const { colors, domain, mapping } = colorTheme.params;

    let lo: number, hi: number;
    if (domain.name === 'custom') {
        [lo, hi] = domain.params;
    } else if (domain.name === 'best') {
        [lo, hi] = distributionRange(values, domain.params.coverage);
        if (domain.params.symmetric) {
            const absMax = Math.max(Math.abs(lo), Math.abs(hi));
            lo = -absMax; hi = absMax;
        }
    } else if (domain.name === 'zscale') {
        [lo, hi] = zscaleRange(values, domain.params.contrast);
    } else {
        [lo, hi] = scalarRange(values);
        if (isMagnitude && hi > 0 && (hi - lo) / hi < 0.01) lo = 0; // degenerate unit-vector case
        if (domain.params.symmetric) {
            const absMax = Math.max(Math.abs(lo), Math.abs(hi));
            lo = -absMax; hi = absMax;
        }
    }
    if (hi - lo < 1e-9) hi = lo + 1e-9; // prevent degenerate domain

    // Apply the scale function (linear/log/symlog/asinh): color scale is built over the transformed
    // domain and looked up with the transformed value, so all mappings reuse the linear ColorScale.
    const { f, dlo, dhi } = makeMappingTransform(mapping, lo, hi);
    const scaleDomain: [number, number] = dhi > dlo ? [dlo, dhi] : [dlo, dlo + 1e-9];
    const colorScale = ColorScale.create({ listOrName: colors, domain: scaleDomain });
    return { color: (gid: number) => colorScale.color(f(values[gid])), usedDomain: [lo, hi] };
}

/** Human-readable summary of the range actually mapped to colors, appended to the domain help text. */
function formatUsedDomain(d?: [number, number]): string {
    if (!d) return '';
    return `current range: [${fmtStat(d[0])}, ${fmtStat(d[1])}]`;
}

function createShape(vtpData: VtpData, mesh: Mesh, colorFn: (gid: number) => Color, props: PD.Values<VtpShapeParams>): Shape<Mesh> {
    const scaleT = Mat4.fromUniformScaling(Mat4(), props.scale);
    const baseTransforms = vtpData.transforms ?? _identityTransforms;
    const transforms = baseTransforms.map(t => Mat4.mul(Mat4(), t, scaleT));
    const { numberOfPoints } = vtpData.source;
    return Shape.create(
        'vtp-mesh', vtpData.source, mesh,
        colorFn,
        () => 1,
        (gid: number) => gid >= 0 ? `Vertex ${gid}` : '',
        transforms,
        numberOfPoints
    );
}

function makeShapeGetter() {
    let _vtpData: VtpData | undefined;
    let _props: PD.Values<VtpShapeParams> | undefined;
    let _shape: Shape<Mesh> | undefined;
    let _mesh: Mesh | undefined;
    let _colorFn: ((gid: number) => Color) | undefined;
    let _vertexResult: VertexResult | null = null;
    let _lastStatsAttr: string | undefined;
    let _lastStatsText: string = '';
    let _lastUsedDomain: [number, number] | undefined;

    const getShape = async (ctx: RuntimeContext, vtpData: VtpData, props: PD.Values<VtpShapeParams>, shape?: Shape<Mesh>) => {
        const params = VtpShapeParams;
        const needsNewMesh = !_vtpData || _vtpData.source !== vtpData.source;

        const newAttrKey = props.colorTheme.name === 'attribute' ? props.colorTheme.params.name : null;
        const oldAttrKey = _props
            ? (_props.colorTheme.name === 'attribute' ? _props.colorTheme.params.name : null)
            : undefined;
        // Per-vertex values come from a VTP attribute or a sampled volume; recompute when either the
        // selected attribute or the referenced volume changes (a changed ref covers switching modes too).
        const newVolRef = props.colorTheme.name === 'external-volume' ? props.colorTheme.params.volume.ref : undefined;
        const oldVolRef = _props?.colorTheme.name === 'external-volume' ? _props.colorTheme.params.volume.ref : undefined;
        const resultChanged = needsNewMesh || newAttrKey !== oldAttrKey || newVolRef !== oldVolRef;

        const colorThemeChanged = !_props || !PD.isParamEqual(params.colorTheme, _props.colorTheme, props.colorTheme);
        const needsNewColor = needsNewMesh || resultChanged || colorThemeChanged;

        const needsNewShape = needsNewMesh || needsNewColor ||
            _vtpData?.transforms !== vtpData.transforms ||
            !_props || _props.scale !== props.scale;

        if (needsNewMesh) {
            _mesh = await buildMesh(ctx, vtpData.source, shape?.geometry);
        }

        if (resultChanged) {
            if (newAttrKey) {
                _vertexResult = computeVertexValues(vtpData.source, newAttrKey);
            } else if (props.colorTheme.name === 'external-volume') {
                const volume = resolveVolumeRef(props.colorTheme.params.volume);
                _vertexResult = volume ? computeVolumeVertexValues(vtpData.source, volume) : null;
            } else {
                _vertexResult = null;
            }
            if (newAttrKey && newAttrKey !== _lastStatsAttr) {
                _lastStatsText = computeAttrStatsText(vtpData.source, newAttrKey);
                _lastStatsAttr = newAttrKey;
            } else if (!newAttrKey) {
                _lastStatsText = '';
                _lastStatsAttr = undefined;
            }
        }

        if (needsNewColor) {
            const cf = makeColorFn(props.colorTheme, _vertexResult);
            _colorFn = cf.color;
            _lastUsedDomain = cf.usedDomain;
        }

        if (needsNewShape) {
            _shape = createShape(vtpData, _mesh!, _colorFn!, props);
        }

        _vtpData = vtpData;
        if (needsNewMesh || needsNewColor || needsNewShape) {
            _props = deepClone(props);
        }
        return _shape!;
    };

    const getStats = () => ({ text: _lastStatsText, attr: _lastStatsAttr });
    const getUsedDomain = () => _lastUsedDomain;
    return { getShape, getStats, getUsedDomain };
}

export function shapeFromVtp(source: VtpFile, params?: { transforms?: Mat4[], getVolumeOptions?: VolumeRefOptions }) {
    return Task.create<ShapeProvider<VtpData, Mesh, VtpShapeParams>>('Shape Provider', async () => {
        const getter = makeShapeGetter();
        const provider: ShapeProvider<VtpData, Mesh, VtpShapeParams> = {
            label: 'VTP Mesh',
            data: { source, transforms: params?.transforms },
            params: createVtpShapeParams(source, () => getter.getStats().text, params?.getVolumeOptions, () => formatUsedDomain(getter.getUsedDomain())),
            getShape: getter.getShape,
            geometryUtils: Mesh.Utils,
        };
        return provider;
    });
}
