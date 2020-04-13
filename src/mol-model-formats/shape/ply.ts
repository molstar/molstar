/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Schäfer, Marco <marco.schaefer@uni-tuebingen.de>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext, Task } from '../../mol-task';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { PlyFile, PlyTable, PlyList } from '../../mol-io/reader/ply/schema';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Shape } from '../../mol-model/shape';
import { ChunkedArray } from '../../mol-data/util';
import { arrayMax, fillSerial } from '../../mol-util/array';
import { Column } from '../../mol-data/db';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ColorNames } from '../../mol-util/color/names';
import { deepClone } from '../../mol-util/object';
import { stringToWords } from '../../mol-util/string';

// TODO support 'edge' element, see https://www.mathworks.com/help/vision/ug/the-ply-format.html
// TODO support missing face element

function createPlyShapeParams(plyFile?: PlyFile) {
    const vertex = plyFile && plyFile.getElement('vertex') as PlyTable;
    const material = plyFile && plyFile.getElement('material') as PlyTable;

    const defaultValues = { group: '', vRed: '', vGreen: '', vBlue: '', mRed: '', mGreen: '', mBlue: '' };

    const groupOptions: [string, string][] = [['', '']];
    const colorOptions: [string, string][] = [['', '']];
    if (vertex) {
        for (let i = 0, il = vertex.propertyNames.length; i < il; ++i) {
            const name = vertex.propertyNames[i];
            const type = vertex.propertyTypes[i];
            if (
                type === 'uchar' || type === 'uint8' ||
                type === 'ushort' || type === 'uint16' ||
                type === 'uint' || type === 'uint32'
            ) groupOptions.push([ name, name ]);
            if (type === 'uchar' || type === 'uint8') colorOptions.push([ name, name ]);
        }

        // TODO hardcoded as convenience for data provided by MegaMol
        if (vertex.propertyNames.includes('atomid')) defaultValues.group = 'atomid';
        else if (vertex.propertyNames.includes('material_index')) defaultValues.group = 'material_index';

        if (vertex.propertyNames.includes('red')) defaultValues.vRed = 'red';
        if (vertex.propertyNames.includes('green')) defaultValues.vGreen = 'green';
        if (vertex.propertyNames.includes('blue')) defaultValues.vBlue = 'blue';
    }

    const materialOptions: [string, string][] = [['', '']];
    if (material) {
        for (let i = 0, il = material.propertyNames.length; i < il; ++i) {
            const name = material.propertyNames[i];
            const type = material.propertyTypes[i];
            if (type === 'uchar' || type === 'uint8') materialOptions.push([ name, name ]);
        }

        if (material.propertyNames.includes('red')) defaultValues.mRed = 'red';
        if (material.propertyNames.includes('green')) defaultValues.mGreen = 'green';
        if (material.propertyNames.includes('blue')) defaultValues.mBlue = 'blue';
    }

    const defaultColoring = defaultValues.vRed && defaultValues.vGreen && defaultValues.vBlue ? 'vertex' :
        defaultValues.mRed && defaultValues.mGreen && defaultValues.mBlue ? 'material' : 'uniform';

    return {
        ...Mesh.Params,

        coloring: PD.MappedStatic(defaultColoring, {
            vertex: PD.Group({
                red: PD.Select(defaultValues.vRed, colorOptions, { label: 'Red Property' }),
                green: PD.Select(defaultValues.vGreen, colorOptions, { label: 'Green Property' }),
                blue: PD.Select(defaultValues.vBlue, colorOptions, { label: 'Blue Property' }),
            }, { isFlat: true }),
            material: PD.Group({
                red: PD.Select(defaultValues.mRed, materialOptions, { label: 'Red Property' }),
                green: PD.Select(defaultValues.mGreen, materialOptions, { label: 'Green Property' }),
                blue: PD.Select(defaultValues.mBlue, materialOptions, { label: 'Blue Property' }),
            }, { isFlat: true }),
            uniform: PD.Group({
                color: PD.Color(ColorNames.grey)
            }, { isFlat: true })
        }),
        grouping: PD.MappedStatic(defaultValues.group ? 'vertex' : 'none', {
            vertex: PD.Group({
                group: PD.Select(defaultValues.group, groupOptions, { label: 'Group Property' }),
            }, { isFlat: true }),
            none: PD.Group({ })
        }),
    };
}

export const PlyShapeParams = createPlyShapeParams();
export type PlyShapeParams = typeof PlyShapeParams

function addVerticesRange(begI: number, endI: number, state: MeshBuilder.State, vertex: PlyTable, groupIds: ArrayLike<number>) {
    const { vertices, normals, groups } = state;

    const x = vertex.getProperty('x');
    const y = vertex.getProperty('y');
    const z = vertex.getProperty('z');
    if (!x || !y || !z) throw new Error('missing coordinate properties');

    const nx = vertex.getProperty('nx');
    const ny = vertex.getProperty('ny');
    const nz = vertex.getProperty('nz');

    const hasNormals = !!nx && !!ny && !!nz;

    for (let i = begI; i < endI; ++i) {
        ChunkedArray.add3(vertices, x.value(i), y.value(i), z.value(i));
        if (hasNormals) ChunkedArray.add3(normals, nx!.value(i), ny!.value(i), nz!.value(i));
        ChunkedArray.add(groups, groupIds[i]);
    }
}

function addFacesRange(begI: number, endI: number, state: MeshBuilder.State, face: PlyList) {
    const { indices } = state;

    for (let i = begI; i < endI; ++i) {
        const { entries, count } = face.value(i);
        if (count === 3) {
            // triangle
            ChunkedArray.add3(indices, entries[0], entries[1], entries[2]);
        } else if (count === 4) {
            // quadrilateral
            ChunkedArray.add3(indices, entries[2], entries[1], entries[0]);
            ChunkedArray.add3(indices, entries[2], entries[0], entries[3]);
        }
    }
}

async function getMesh(ctx: RuntimeContext, vertex: PlyTable, face: PlyList, groupIds: ArrayLike<number>, mesh?: Mesh) {
    const builderState = MeshBuilder.createState(vertex.rowCount, vertex.rowCount / 4, mesh);

    const x = vertex.getProperty('x');
    const y = vertex.getProperty('y');
    const z = vertex.getProperty('z');
    if (!x || !y || !z) throw new Error('missing coordinate properties');

    const nx = vertex.getProperty('nx');
    const ny = vertex.getProperty('ny');
    const nz = vertex.getProperty('nz');

    const hasNormals = !!nx && !!ny && !!nz;
    const updateChunk = 100000;

    for (let i = 0, il = vertex.rowCount; i < il; i += updateChunk) {
        addVerticesRange(i, Math.min(i + updateChunk, il), builderState, vertex, groupIds);

        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'adding ply mesh vertices', current: i, max: il });
        }
    }

    for (let i = 0, il = face.rowCount; i < il; i += updateChunk) {
        addFacesRange(i, Math.min(i + updateChunk, il), builderState, face);

        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'adding ply mesh faces', current: i, max: il });
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    if (!hasNormals) Mesh.computeNormals(m);

    return m;
}

const int = Column.Schema.int;

type Grouping = { ids: ArrayLike<number>, map: ArrayLike<number>, label: string }
function getGrouping(vertex: PlyTable, props: PD.Values<PlyShapeParams>): Grouping {
    const { grouping } = props;
    const { rowCount } = vertex;
    const column = grouping.name === 'vertex' ? vertex.getProperty(grouping.params.group) : undefined;
    const label = grouping.name === 'vertex' ? stringToWords(grouping.params.group) : 'Vertex';

    const ids = column ? column.toArray({ array: Uint32Array }) : fillSerial(new Uint32Array(rowCount));
    const maxId = column ? arrayMax(ids) : rowCount - 1; // assumes uint ids
    const map = new Uint32Array(maxId + 1);
    for (let i = 0, il = ids.length; i < il; ++i) map[ids[i]] = i;
    return { ids, map, label  };
}

type Coloring = { kind: 'vertex' | 'material' | 'uniform', red: Column<number>, green: Column<number>, blue: Column<number> }
function getColoring(vertex: PlyTable, material: PlyTable | undefined, props: PD.Values<PlyShapeParams>): Coloring {
    const { coloring } = props;
    const { rowCount } = vertex;

    let red: Column<number>, green: Column<number>, blue: Column<number>;
    if (coloring.name === 'vertex') {
        red = vertex.getProperty(coloring.params.red) || Column.ofConst(127, rowCount, int);
        green = vertex.getProperty(coloring.params.green) || Column.ofConst(127, rowCount, int);
        blue = vertex.getProperty(coloring.params.blue) || Column.ofConst(127, rowCount, int);
    } else if (coloring.name === 'material') {
        red = (material && material.getProperty(coloring.params.red)) || Column.ofConst(127, rowCount, int);
        green = (material && material.getProperty(coloring.params.green)) || Column.ofConst(127, rowCount, int);
        blue = (material && material.getProperty(coloring.params.blue)) || Column.ofConst(127, rowCount, int);
    } else {
        const [r, g, b] = Color.toRgb(coloring.params.color);
        red = Column.ofConst(r, rowCount, int);
        green = Column.ofConst(g, rowCount, int);
        blue = Column.ofConst(b, rowCount, int);
    }
    return { kind: coloring.name, red, green, blue };
}

function createShape(plyFile: PlyFile, mesh: Mesh, coloring: Coloring, grouping: Grouping) {
    const { kind, red, green, blue } = coloring;
    const { ids, map, label } = grouping;
    return Shape.create(
        'ply-mesh', plyFile, mesh,
        (groupId: number) => {
            const idx = kind === 'material' ? groupId : map[groupId];
            return Color.fromRgb(red.value(idx), green.value(idx), blue.value(idx));
        },
        () => 1, // size: constant
        (groupId: number) => {
            return `${label} ${ids[groupId]}`;
        }
    );
}

function makeShapeGetter() {
    let _plyFile: PlyFile | undefined;
    let _props: PD.Values<PlyShapeParams> | undefined;

    let _shape: Shape<Mesh>;
    let _mesh: Mesh;
    let _coloring: Coloring;
    let _grouping: Grouping;

    const getShape = async (ctx: RuntimeContext, plyFile: PlyFile, props: PD.Values<PlyShapeParams>, shape?: Shape<Mesh>) => {

        const vertex = plyFile.getElement('vertex') as PlyTable;
        if (!vertex) throw new Error('missing vertex element');

        const face = plyFile.getElement('face') as PlyList;
        if (!face) throw new Error('missing face element');

        const material = plyFile.getElement('material') as PlyTable;

        let newMesh = false;
        let newColor = false;

        if (!_plyFile || _plyFile !== plyFile) {
            newMesh = true;
        }

        if (!_props || !PD.isParamEqual(PlyShapeParams.grouping, _props.grouping, props.grouping)) {
            newMesh = true;
        }

        if (!_props || !PD.isParamEqual(PlyShapeParams.coloring, _props.coloring, props.coloring)) {
            newColor = true;
        }

        if (newMesh) {
            _coloring = getColoring(vertex, material, props);
            _grouping = getGrouping(vertex, props);
            _mesh = await getMesh(ctx, vertex, face, _grouping.ids, shape && shape.geometry);
            _shape = createShape(plyFile, _mesh, _coloring, _grouping);
        } else if (newColor) {
            _coloring = getColoring(vertex, material, props);
            _shape = createShape(plyFile, _mesh, _coloring, _grouping);
        }

        _plyFile = plyFile;
        _props = deepClone(props);

        return _shape;
    };
    return getShape;
}

export function shapeFromPly(source: PlyFile, params?: {}) {
    return Task.create<ShapeProvider<PlyFile, Mesh, PlyShapeParams>>('Shape Provider', async ctx => {
        return {
            label: 'Mesh',
            data: source,
            params: createPlyShapeParams(source),
            getShape: makeShapeGetter(),
            geometryUtils: Mesh.Utils
        };
    });
}