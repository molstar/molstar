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

export type VtpData = {
    source: VtpFile,
    transforms?: Mat4[],
}

const _identityTransforms = [Mat4.identity()];

// --- Params ---

const UNIFORM_KEY = 'uniform::';

function buildAttrOptions(vtpFile?: VtpFile): PD.SelectOption<string>[] {
    const opts: PD.SelectOption<string>[] = [[UNIFORM_KEY, 'Uniform color']];
    if (!vtpFile) return opts;
    for (const [name] of vtpFile.pointData) {
        opts.push([`point:${name}`, `Point: ${name}`]);
    }
    for (const [name] of vtpFile.cellData) {
        opts.push([`cell:${name}`, `Cell: ${name}`]);
    }
    return opts;
}

function defaultAttrKey(vtpFile?: VtpFile): string {
    if (!vtpFile) return UNIFORM_KEY;
    if (vtpFile.cellData.size > 0) return `cell:${vtpFile.cellData.keys().next().value}`;
    if (vtpFile.pointData.size > 0) return `point:${vtpFile.pointData.keys().next().value}`;
    return UNIFORM_KEY;
}

function defaultDomain(vtpFile: VtpFile | undefined, key: string): [number, number] {
    if (!vtpFile || key === UNIFORM_KEY) return [0, 1];
    if (key.startsWith('cell:')) {
        const arr = vtpFile.cellData.get(key.slice(5));
        if (arr) return [arr.desc.rangeMin, arr.desc.rangeMax];
    } else if (key.startsWith('point:')) {
        const arr = vtpFile.pointData.get(key.slice(6));
        if (arr) return [arr.desc.rangeMin, arr.desc.rangeMax];
    }
    return [0, 1];
}

export function createVtpShapeParams(vtpFile?: VtpFile) {
    const attrOptions = buildAttrOptions(vtpFile);
    const defKey = defaultAttrKey(vtpFile);
    const [defMin, defMax] = defaultDomain(vtpFile, defKey);

    return {
        ...Mesh.Params,
        doubleSided: { ...Mesh.Params.doubleSided, defaultValue: true },
        attribute: PD.Select(defKey, attrOptions, { label: 'Color Attribute' }),
        smoothColor: PD.Boolean(true, { label: 'Smooth Color', description: 'Interpolate cell data to vertices for smooth color gradients.' }),
        colormap: PD.Select('red-yellow-blue' as ColorListName, ColorListOptionsScale, { label: 'Colormap' }),
        domainMin: PD.Numeric(defMin, {}, { label: 'Domain Min' }),
        domainMax: PD.Numeric(defMax, {}, { label: 'Domain Max' }),
        uniformColor: PD.Color(ColorNames.grey, { label: 'Uniform Color' }),
        scale: PD.Numeric(1, { min: 0.01, max: 100, step: 0.01 }, { label: 'Scale', description: 'Uniform scale factor applied to the mesh.' }),
    };
}

export const VtpShapeParams = createVtpShapeParams();
export type VtpShapeParams = typeof VtpShapeParams;

// --- Mesh building ---

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
    const { connectivity, numberOfPoints } = vtpFile;
    const nTris = connectivity.length / 3;
    const sum = new Float64Array(numberOfPoints);
    const count = new Uint32Array(numberOfPoints);
    for (let i = 0; i < nTris; i++) {
        const val = cellValues[i];
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
        if (smoothColor) {
            const smoothed = cellToVertexAverage(vtpFile, arr.values);
            const { connectivity } = vtpFile;
            return (gid: number) => colorScale.color(smoothed[connectivity[gid]]);
        }
        return (gid: number) => colorScale.color(arr.values[Math.floor(gid / 3)]);
    }

    if (attribute.startsWith('point:')) {
        const name = attribute.slice(6);
        const arr = vtpFile.pointData.get(name);
        if (!arr) return () => uniformColor;
        const { values } = arr;
        const { connectivity } = vtpFile;
        return (gid: number) => colorScale.color(values[connectivity[gid]]);
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
            _colorFn = makeColorFn(vtpData.source, props);
        }

        if (needsNewShape) {
            _shape = createShape(vtpData, _mesh!, _colorFn!, props.scale);
        }

        _vtpData = vtpData;
        _props = deepClone(props);
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
