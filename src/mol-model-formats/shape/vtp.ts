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

export interface VtpData {
    source: VtpFile,
    transforms?: Mat4[],
}

const _identityTransforms = [Mat4.identity()];

// --- Params ---

const UNIFORM_KEY = 'uniform::';

function fmtStat(x: number): string {
    if (!isFinite(x)) return x.toString();
    const abs = Math.abs(x);
    if (abs === 0) return '0';
    if (abs < 0.001 || abs >= 1e5) return x.toExponential(2);
    return x.toPrecision(3);
}

function computeAttrStatsText(vtpFile: VtpFile, attrKey: string): string {
    if (attrKey === UNIFORM_KEY) return '';
    let rawVals: Float64Array | undefined;
    let nComp = 1;
    let name = '';
    if (attrKey.startsWith('point:')) {
        name = attrKey.slice(6);
        const arr = vtpFile.pointData.get(name);
        if (!arr) return '';
        rawVals = arr.values; nComp = arr.desc.numberOfComponents;
    } else if (attrKey.startsWith('cell:')) {
        name = attrKey.slice(5);
        const arr = vtpFile.cellData.get(name);
        if (!arr) return '';
        rawVals = arr.values; nComp = arr.desc.numberOfComponents;
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
    const opts: PD.SelectOption<string>[] = [[UNIFORM_KEY, 'Uniform color']];
    if (!vtpFile) return opts;
    for (const [name, arr] of vtpFile.pointData) {
        const nComp = arr.desc.numberOfComponents;
        opts.push([`point:${name}`, nComp > 1 ? `Point: ${name} (mag)` : `Point: ${name}`]);
    }
    for (const [name, arr] of vtpFile.cellData) {
        const nComp = arr.desc.numberOfComponents;
        opts.push([`cell:${name}`, nComp > 1 ? `Cell: ${name} (mag)` : `Cell: ${name}`]);
    }
    return opts;
}

function defaultAttrKey(vtpFile?: VtpFile): string {
    if (!vtpFile) return UNIFORM_KEY;
    // Prefer 1-component arrays first (direct scalar mapping)
    for (const [name, arr] of vtpFile.cellData) {
        if (arr.desc.numberOfComponents === 1) return `cell:${name}`;
    }
    for (const [name, arr] of vtpFile.pointData) {
        if (arr.desc.numberOfComponents === 1) return `point:${name}`;
    }
    if (vtpFile.cellData.size > 0) return `cell:${vtpFile.cellData.keys().next().value}`;
    if (vtpFile.pointData.size > 0) return `point:${vtpFile.pointData.keys().next().value}`;
    return UNIFORM_KEY;
}

function scalarRange(values: Float64Array): [number, number] {
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < values.length; i++) {
        if (values[i] < min) min = values[i];
        if (values[i] > max) max = values[i];
    }
    return [min === Infinity ? 0 : min, max === -Infinity ? 1 : max];
}

function magnitudeRange(values: Float64Array, nComp: number): [number, number] {
    const n = values.length / nComp;
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < n; i++) {
        let mag2 = 0;
        for (let c = 0; c < nComp; c++) mag2 += values[i * nComp + c] ** 2;
        const mag = Math.sqrt(mag2);
        if (mag < min) min = mag;
        if (mag > max) max = mag;
    }
    return [min === Infinity ? 0 : min, max === -Infinity ? 1 : max];
}

function defaultDomain(vtpFile: VtpFile | undefined, key: string): [number, number] {
    if (!vtpFile || key === UNIFORM_KEY) return [0, 1];
    let arr;
    if (key.startsWith('cell:')) arr = vtpFile.cellData.get(key.slice(5));
    else if (key.startsWith('point:')) arr = vtpFile.pointData.get(key.slice(6));
    if (!arr) return [0, 1];
    const nComp = arr.desc.numberOfComponents;
    if (nComp > 1) return magnitudeRange(arr.values, nComp);
    // Use VTK-provided range when available; otherwise compute from data.
    // Many writers (Python VTK, SurfaceMorphometrics) omit RangeMin/RangeMax.
    if (arr.desc.hasRange) return [arr.desc.rangeMin, arr.desc.rangeMax];
    return scalarRange(arr.values);
}

export function createVtpShapeParams(vtpFile?: VtpFile) {
    const attrOptions = buildAttrOptions(vtpFile);
    const defKey = defaultAttrKey(vtpFile);
    const [defMin, defMax] = defaultDomain(vtpFile, defKey);
    const defStats = vtpFile ? computeAttrStatsText(vtpFile, defKey) : '';

    return {
        ...Mesh.Params,
        doubleSided: { ...Mesh.Params.doubleSided, defaultValue: true },
        interior: { ...Mesh.Params.interior, defaultValue: { ...Mesh.Params.interior.defaultValue, colorStrength: 0 } },
        attribute: PD.Select(defKey, attrOptions, { label: 'Color Attribute' }),
        colormap: PD.Select('viridis' as ColorListName, ColorListOptionsScale, { label: 'Colormap' }),
        domainMin: PD.Numeric(defMin, {}, { label: 'Domain Min', isHidden: true }),
        domainMax: PD.Numeric(defMax, {}, { label: 'Domain Max', isHidden: true }),
        statsText: PD.Text(defStats, { label: 'Statistics', multiline: true }),
        uniformColor: PD.Color(ColorNames.grey, { label: 'Uniform Color' }),
        scale: PD.Numeric(1, { min: 0.01, max: 100, step: 0.01 }, { label: 'Scale', description: 'Uniform scale factor applied to the mesh.' }),
    };
}

export const VtpShapeParams = createVtpShapeParams();
export type VtpShapeParams = typeof VtpShapeParams;

// --- Mesh building ---
// Vertices are intentionally NOT shared (3 unique entries per triangle).
// This allows each vertex to carry an independent group ID, which is required
// for flat per-cell CellData coloring: gid/3 = triangle index, gid%3 = vertex
// within that triangle. Sharing vertices would break CellData flat coloring
// since a vertex at a cell boundary could only hold one of its cells' values.

async function buildMesh(ctx: RuntimeContext, vtpFile: VtpFile, mesh?: Mesh): Promise<Mesh> {
    const { positions, connectivity, numberOfTriangles } = vtpFile;
    const nVerts = numberOfTriangles * 3;

    const builderState = MeshBuilder.createState(nVerts, numberOfTriangles, mesh);
    const { vertices, indices, groups } = builderState;

    const chunkSize = 50000;

    for (let ti = 0, il = numberOfTriangles; ti < il; ti += chunkSize) {
        const end = Math.min(ti + chunkSize, il);
        for (let i = ti; i < end; i++) {
            const v0 = connectivity[3 * i];
            const v1 = connectivity[3 * i + 1];
            const v2 = connectivity[3 * i + 2];

            const gBase = 3 * i;
            ChunkedArray.add3(vertices, positions[3 * v0], positions[3 * v0 + 1], positions[3 * v0 + 2]);
            ChunkedArray.add(groups, gBase);
            ChunkedArray.add3(vertices, positions[3 * v1], positions[3 * v1 + 1], positions[3 * v1 + 2]);
            ChunkedArray.add(groups, gBase + 1);
            ChunkedArray.add3(vertices, positions[3 * v2], positions[3 * v2 + 1], positions[3 * v2 + 2]);
            ChunkedArray.add(groups, gBase + 2);

            ChunkedArray.add3(indices, gBase, gBase + 1, gBase + 2);
        }

        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'Building VTP mesh...', current: ti, max: numberOfTriangles });
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    Mesh.computeNormals(m);
    ValueCell.updateIfChanged(m.varyingGroup, true);
    return m;
}

// --- Color lookup ---

function cellToVertexAverage(vtpFile: VtpFile, cellValues: Float64Array): Float64Array {
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

// Average multi-component cell vectors to vertices, then return the magnitude of each vertex vector.
// Matching smgui _cell_to_vertex_interpolation_vector: averaging direction vectors, not their magnitudes,
// so that crease/boundary vertices (where adjacent face normals cancel) get low magnitude values.
function cellToVertexAverageMag(vtpFile: VtpFile, cellVectors: Float64Array, nComp: number): Float64Array {
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

function vecMag(values: Float64Array, baseIdx: number, nComp: number): number {
    let mag2 = 0;
    for (let c = 0; c < nComp; c++) mag2 += values[baseIdx + c] ** 2;
    return Math.sqrt(mag2);
}

interface VertexResult {
    values: Float64Array;
    isMagnitude: boolean; // true for multi-component (magnitude) attributes
}

// Compute per-vertex scalar values from a VTP attribute, matching smgui exactly:
//   - CellData scalar: average cell scalars to vertices
//   - CellData multi-component: average raw vectors to vertices first, then compute magnitude
//     (matches smgui _cell_to_vertex_interpolation_vector → norm; crease vertices cancel → low magnitude)
//   - PointData scalar: use values directly
//   - PointData multi-component: compute per-vertex magnitude directly
// Returns null for UNIFORM_KEY or missing attributes.
function computeVertexValues(vtpFile: VtpFile, attribute: string): VertexResult | null {
    if (attribute === UNIFORM_KEY) return null;

    if (attribute.startsWith('cell:')) {
        const arr = vtpFile.cellData.get(attribute.slice(5));
        if (!arr) return null;
        const nComp = arr.desc.numberOfComponents;
        if (nComp === 1) {
            return { values: cellToVertexAverage(vtpFile, arr.values), isMagnitude: false };
        }
        // Average vectors to vertices first, then take magnitude — crease vertices get lower values
        return { values: cellToVertexAverageMag(vtpFile, arr.values, nComp), isMagnitude: true };
    }

    if (attribute.startsWith('point:')) {
        const arr = vtpFile.pointData.get(attribute.slice(6));
        if (!arr) return null;
        const nComp = arr.desc.numberOfComponents;
        if (nComp === 1) return { values: arr.values, isMagnitude: false };
        // multi-component point data: per-vertex magnitude
        const nPts = arr.values.length / nComp;
        const mag = new Float64Array(nPts);
        for (let i = 0; i < nPts; i++) mag[i] = vecMag(arr.values, i * nComp, nComp);
        return { values: mag, isMagnitude: true };
    }

    return null;
}

// Build a color function from pre-computed per-vertex scalar values.
// Domain is [min, max] (matching smgui). For magnitude attributes where all values are
// nearly identical (e.g. PointData unit normals, range < 1% of max), fall back to [0, max]
// so the entire colormap is used and values map to the high end (yellow in viridis).
function makeColorFn(vtpFile: VtpFile, colormap: ColorListName, uniformColor: Color, attribute: string, result: VertexResult | null): (gid: number) => Color {
    if (attribute === UNIFORM_KEY || !result) return () => uniformColor;

    const { values, isMagnitude } = result;
    let [lo, hi] = scalarRange(values);
    if (isMagnitude && hi > 0 && (hi - lo) / hi < 0.01) lo = 0; // degenerate unit-vector case
    if (hi - lo < 1e-9) hi = lo + 1e-9; // prevent degenerate domain

    const colorScale = ColorScale.create({ listOrName: colormap, domain: [lo, hi] });
    const { connectivity } = vtpFile;
    return (gid: number) => colorScale.color(values[connectivity[gid]]);
}

function createShape(vtpData: VtpData, mesh: Mesh, colorFn: (gid: number) => Color, scale: number): Shape<Mesh> {
    const scaleT = Mat4.fromUniformScaling(Mat4(), scale);
    const baseTransforms = vtpData.transforms ?? _identityTransforms;
    const transforms = baseTransforms.map(t => Mat4.mul(Mat4(), t, scaleT));
    return Shape.create(
        'vtp-mesh', vtpData.source, mesh,
        colorFn,
        () => 1,
        (gid: number) => {
            if (gid < 0) return '';
            return `Triangle ${Math.floor(gid / 3)}`;
        },
        transforms
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

    const getShape = async (ctx: RuntimeContext, vtpData: VtpData, props: PD.Values<VtpShapeParams>, shape?: Shape<Mesh>) => {
        const needsNewMesh = !_vtpData || _vtpData.source !== vtpData.source;
        const attributeChanged = !_props || _props.attribute !== props.attribute || needsNewMesh;
        const needsNewColor = needsNewMesh || attributeChanged ||
            !_props ||
            _props.colormap !== props.colormap ||
            _props.uniformColor !== props.uniformColor;
        const needsNewShape = needsNewMesh || needsNewColor ||
            _vtpData?.transforms !== vtpData.transforms ||
            !_props || _props.scale !== props.scale;

        if (needsNewMesh) {
            _mesh = await buildMesh(ctx, vtpData.source, shape?.geometry);
        }

        if (attributeChanged) {
            _vertexResult = computeVertexValues(vtpData.source, props.attribute);
            if (props.attribute !== _lastStatsAttr) {
                _lastStatsText = computeAttrStatsText(vtpData.source, props.attribute);
                _lastStatsAttr = props.attribute;
            }
        }

        if (needsNewColor) {
            _colorFn = makeColorFn(vtpData.source, props.colormap as ColorListName, props.uniformColor, props.attribute, _vertexResult);
        }

        if (needsNewShape) {
            _shape = createShape(vtpData, _mesh!, _colorFn!, props.scale);
        }

        _vtpData = vtpData;
        if (needsNewMesh || needsNewColor || needsNewShape) {
            _props = deepClone(props);
        }
        return _shape!;
    };

    const getStats = () => ({ text: _lastStatsText, attr: _lastStatsAttr });
    return { getShape, getStats };
}

export function shapeFromVtp(source: VtpFile, params?: { transforms?: Mat4[] }) {
    return Task.create<ShapeProvider<VtpData, Mesh, VtpShapeParams>>('Shape Provider', async () => {
        const getter = makeShapeGetter();
        const provider: ShapeProvider<VtpData, Mesh, VtpShapeParams> & { onParamsUpdate?: (props: PD.Values<VtpShapeParams>) => Record<string, unknown> | null } = {
            label: 'VTP Mesh',
            data: { source, transforms: params?.transforms },
            params: createVtpShapeParams(source),
            getShape: getter.getShape,
            geometryUtils: Mesh.Utils,
            onParamsUpdate(props: PD.Values<VtpShapeParams>) {
                const { text, attr } = getter.getStats();
                if (!attr || props.attribute !== attr) return null;
                if ((props as any).statsText === text) return null; // already up to date
                return { statsText: text };
            }
        };
        return provider;
    });
}
