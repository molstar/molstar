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


export function createVtpShapeParams(vtpFile?: VtpFile, getStats?: () => string) {
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
                colors: PD.Select('viridis' as ColorListName, ColorListOptionsScale),
                domain: PD.MappedStatic('auto', {
                    custom: PD.Interval([-1, 1], { step: 0.001 }),
                    auto: PD.Group({
                        symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
                    })
                })
            }),
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

/**
 * Build a color function from pre-computed per-vertex scalar values.
 * Auto domain matches smgui [min, max]. For magnitude attributes where all values are
 * nearly identical (e.g. PointData unit normals, range < 1% of max), fall back to [0, max].
 */
function makeColorFn(vtpFile: VtpFile, colorTheme: PD.Values<VtpShapeParams>['colorTheme'], result: VertexResult | null): (gid: number) => Color {
    if (colorTheme.name === 'uniform') return () => colorTheme.params.color;
    if (!result) return () => ColorNames.grey;

    const { values, isMagnitude } = result;
    const { colors, domain } = colorTheme.params;

    let lo: number, hi: number;
    if (domain.name === 'custom') {
        [lo, hi] = domain.params;
    } else {
        [lo, hi] = scalarRange(values);
        if (isMagnitude && hi > 0 && (hi - lo) / hi < 0.01) lo = 0; // degenerate unit-vector case
        if (domain.params.symmetric) {
            const absMax = Math.max(Math.abs(lo), Math.abs(hi));
            lo = -absMax; hi = absMax;
        }
    }
    if (hi - lo < 1e-9) hi = lo + 1e-9; // prevent degenerate domain

    const colorScale = ColorScale.create({ listOrName: colors, domain: [lo, hi] });
    return (gid: number) => colorScale.color(values[gid]);
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

    const getShape = async (ctx: RuntimeContext, vtpData: VtpData, props: PD.Values<VtpShapeParams>, shape?: Shape<Mesh>) => {
        const params = VtpShapeParams;
        const needsNewMesh = !_vtpData || _vtpData.source !== vtpData.source;

        const newAttrKey = props.colorTheme.name === 'attribute' ? props.colorTheme.params.name : null;
        const oldAttrKey = _props
            ? (_props.colorTheme.name === 'attribute' ? _props.colorTheme.params.name : null)
            : undefined;
        const attributeChanged = needsNewMesh || newAttrKey !== oldAttrKey;

        const colorThemeChanged = !_props || !PD.isParamEqual(params.colorTheme, _props.colorTheme, props.colorTheme);
        const needsNewColor = needsNewMesh || attributeChanged || colorThemeChanged;

        const needsNewShape = needsNewMesh || needsNewColor ||
            _vtpData?.transforms !== vtpData.transforms ||
            !_props || _props.scale !== props.scale;

        if (needsNewMesh) {
            _mesh = await buildMesh(ctx, vtpData.source, shape?.geometry);
        }

        if (attributeChanged) {
            _vertexResult = newAttrKey ? computeVertexValues(vtpData.source, newAttrKey) : null;
            if (newAttrKey && newAttrKey !== _lastStatsAttr) {
                _lastStatsText = computeAttrStatsText(vtpData.source, newAttrKey);
                _lastStatsAttr = newAttrKey;
            } else if (!newAttrKey) {
                _lastStatsText = '';
                _lastStatsAttr = undefined;
            }
        }

        if (needsNewColor) {
            _colorFn = makeColorFn(vtpData.source, props.colorTheme, _vertexResult);
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
    return { getShape, getStats };
}

export function shapeFromVtp(source: VtpFile, params?: { transforms?: Mat4[] }) {
    return Task.create<ShapeProvider<VtpData, Mesh, VtpShapeParams>>('Shape Provider', async () => {
        const getter = makeShapeGetter();
        const provider: ShapeProvider<VtpData, Mesh, VtpShapeParams> = {
            label: 'VTP Mesh',
            data: { source, transforms: params?.transforms },
            params: createVtpShapeParams(source, () => getter.getStats().text),
            getShape: getter.getShape,
            geometryUtils: Mesh.Utils,
        };
        return provider;
    });
}
