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

function buildAttrOptions(vtpFile?: VtpFile): PD.SelectOption<string>[] {
    const opts: PD.SelectOption<string>[] = [[UNIFORM_KEY, 'Uniform color']];
    if (!vtpFile) return opts;
    for (const [name, arr] of vtpFile.pointData) {
        const nComp = arr.desc.numberOfComponents;
        const label = nComp > 1 ? `Point: ${name} (mag)` : `Point: ${name}`;
        opts.push([`point:${name}`, label]);
    }
    for (const [name, arr] of vtpFile.cellData) {
        const nComp = arr.desc.numberOfComponents;
        const label = nComp > 1 ? `Cell: ${name} (mag)` : `Cell: ${name}`;
        opts.push([`cell:${name}`, label]);
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

    return {
        ...Mesh.Params,
        doubleSided: { ...Mesh.Params.doubleSided, defaultValue: true },
        interior: { ...Mesh.Params.interior, defaultValue: { ...Mesh.Params.interior.defaultValue, colorStrength: 0 } },
        attribute: PD.Select(defKey, attrOptions, { label: 'Color Attribute' }),
        smoothColor: PD.Boolean(false, { label: 'Smooth Color', description: 'Interpolate per-cell data to vertices for smooth color gradients (no effect for per-vertex attributes).' }),
        colormap: PD.Select('viridis' as ColorListName, ColorListOptionsScale, { label: 'Colormap' }),
        domainMin: PD.Numeric(defMin, {}, { label: 'Domain Min' }),
        domainMax: PD.Numeric(defMax, {}, { label: 'Domain Max' }),
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

function vecMag(values: Float64Array, baseIdx: number, nComp: number): number {
    let mag2 = 0;
    for (let c = 0; c < nComp; c++) mag2 += values[baseIdx + c] ** 2;
    return Math.sqrt(mag2);
}

function makeColorFn(vtpFile: VtpFile, props: PD.Values<VtpShapeParams>): (gid: number) => Color {
    const { attribute, uniformColor, colormap, domainMin, domainMax, smoothColor } = props;

    if (attribute === UNIFORM_KEY) {
        return () => uniformColor;
    }

    const colorScale = ColorScale.create({ listOrName: colormap, domain: [domainMin, domainMax] });

    if (attribute.startsWith('cell:')) {
        const name = attribute.slice(5);
        const arr = vtpFile.cellData.get(name);
        if (!arr) return () => uniformColor;
        const nComp = arr.desc.numberOfComponents;
        const { triangleCellIndex } = vtpFile;
        if (nComp === 1) {
            if (smoothColor) {
                const smoothed = cellToVertexAverage(vtpFile, arr.values);
                // Domain must be recomputed from smoothed values: averaging compresses the
                // raw [min, max] range, so using the raw domain would place all smoothed
                // values in the middle of the colormap.
                const [sMin, sMax] = scalarRange(smoothed);
                const smoothScale = ColorScale.create({ listOrName: colormap, domain: [sMin, sMax] });
                const { connectivity } = vtpFile;
                return (gid: number) => smoothScale.color(smoothed[connectivity[gid]]);
            }
            return (gid: number) => colorScale.color(arr.values[triangleCellIndex[Math.floor(gid / 3)]]);
        }
        // multi-component: magnitude per cell, no smooth (averaging magnitudes not meaningful)
        const { values } = arr;
        return (gid: number) => colorScale.color(vecMag(values, triangleCellIndex[Math.floor(gid / 3)] * nComp, nComp));
    }

    if (attribute.startsWith('point:')) {
        const name = attribute.slice(6);
        const arr = vtpFile.pointData.get(name);
        if (!arr) return () => uniformColor;
        const nComp = arr.desc.numberOfComponents;
        const { connectivity, numberOfTriangles } = vtpFile;

        if (nComp === 1) {
            if (smoothColor) {
                // Per-vertex lookup: shows fine-grained per-vertex variation.
                return (gid: number) => colorScale.color(arr.values[connectivity[gid]]);
            }
            // Default: flat per-triangle average of the 3 corner vertex values.
            // For smoothly-varying per-vertex data (e.g. v_n, curvature), adjacent
            // triangles have nearly identical averages → visually smooth surface.
            // (True Gouraud shading is not available through the group-ID texture
            // mechanism used here, so flat-averaged is the best approximation.)
            const triAvg = new Float64Array(numberOfTriangles);
            for (let i = 0; i < numberOfTriangles; i++) {
                triAvg[i] = (arr.values[connectivity[3 * i]] + arr.values[connectivity[3 * i + 1]] + arr.values[connectivity[3 * i + 2]]) / 3;
            }
            const [avgMin, avgMax] = scalarRange(triAvg);
            const avgScale = ColorScale.create({ listOrName: colormap, domain: [avgMin, avgMax] });
            return (gid: number) => avgScale.color(triAvg[Math.floor(gid / 3)]);
        }

        // multi-component: magnitude
        const { values } = arr;
        if (smoothColor) {
            // Per-vertex magnitude (fine-grained variation).
            return (gid: number) => colorScale.color(vecMag(values, connectivity[gid] * nComp, nComp));
        }
        // Default: flat per-triangle average of vertex magnitudes.
        const triAvg = new Float64Array(numberOfTriangles);
        for (let i = 0; i < numberOfTriangles; i++) {
            const m0 = vecMag(values, connectivity[3 * i] * nComp, nComp);
            const m1 = vecMag(values, connectivity[3 * i + 1] * nComp, nComp);
            const m2 = vecMag(values, connectivity[3 * i + 2] * nComp, nComp);
            triAvg[i] = (m0 + m1 + m2) / 3;
        }
        const [avgMin, avgMax] = scalarRange(triAvg);
        const avgScale = ColorScale.create({ listOrName: colormap, domain: [avgMin, avgMax] });
        return (gid: number) => avgScale.color(triAvg[Math.floor(gid / 3)]);
    }

    return () => uniformColor;
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
    // Effective domain persisted in the closure so that secondary param changes
    // (e.g. smoothColor toggle) reuse the correct domain rather than the stale
    // domainMin/domainMax stored in the state transform's props.
    let _computedDomainMin = 0;
    let _computedDomainMax = 1;

    const getShape = async (ctx: RuntimeContext, vtpData: VtpData, props: PD.Values<VtpShapeParams>, shape?: Shape<Mesh>) => {
        const needsNewMesh = !_vtpData || _vtpData.source !== vtpData.source;
        const needsNewColor = needsNewMesh || !_props ||
            _props.attribute !== props.attribute ||
            _props.smoothColor !== props.smoothColor ||
            _props.colormap !== props.colormap ||
            _props.domainMin !== props.domainMin ||
            _props.domainMax !== props.domainMax ||
            _props.uniformColor !== props.uniformColor;
        const needsNewShape = needsNewMesh || needsNewColor ||
            _vtpData?.transforms !== vtpData.transforms ||
            !_props || _props.scale !== props.scale;

        if (needsNewMesh) {
            _mesh = await buildMesh(ctx, vtpData.source, shape?.geometry);
        }

        if (needsNewColor) {
            // Recompute domain only when the attribute (or underlying data) changes.
            // Store it in the closure so subsequent toggles (smoothColor, colormap, etc.)
            // continue to use the correct domain for the current attribute.
            const attributeChanged = !_props || _props.attribute !== props.attribute || needsNewMesh;
            if (attributeChanged) {
                [_computedDomainMin, _computedDomainMax] = defaultDomain(vtpData.source, props.attribute);
            }
            _colorFn = makeColorFn(vtpData.source, { ...props, domainMin: _computedDomainMin, domainMax: _computedDomainMax });
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

    return getShape;
}

export function shapeFromVtp(source: VtpFile, params?: { transforms?: Mat4[] }) {
    return Task.create<ShapeProvider<VtpData, Mesh, VtpShapeParams>>('Shape Provider', async () => {
        return {
            label: 'VTP Mesh',
            data: { source, transforms: params?.transforms },
            params: createVtpShapeParams(source),
            getShape: makeShapeGetter(),
            geometryUtils: Mesh.Utils
        };
    });
}
