/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { AssemblySymmetryValue, AssemblySymmetryProvider, AssemblySymmetry } from './prop';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3, Mat4, Mat3 } from '../../../mol-math/linear-algebra';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { RuntimeContext } from '../../../mol-task';
import { Shape } from '../../../mol-model/shape';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../../../mol-repr/shape/representation';
import { MarkerActions } from '../../../mol-util/marker-action';
import { Prism, PrismCage } from '../../../mol-geo/primitive/prism';
import { Wedge, WedgeCage } from '../../../mol-geo/primitive/wedge';
import { Primitive, transformPrimitive } from '../../../mol-geo/primitive/primitive';
import { memoize1 } from '../../../mol-util/memoize';
import { polygon } from '../../../mol-geo/primitive/polygon';
import { ColorMap, Color } from '../../../mol-util/color';
import { TableLegend } from '../../../mol-util/legend';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { Cage, transformCage, cloneCage } from '../../../mol-geo/primitive/cage';
import { OctahedronCage } from '../../../mol-geo/primitive/octahedron';
import { TetrahedronCage } from '../../../mol-geo/primitive/tetrahedron';
import { IcosahedronCage } from '../../../mol-geo/primitive/icosahedron';
import { degToRad, radToDeg } from '../../../mol-math/misc';
import { Mutable } from '../../../mol-util/type-helpers';
import { equalEps } from '../../../mol-math/linear-algebra/3d/common';
import { Structure } from '../../../mol-model/structure';
import { isInteger } from '../../../mol-util/number';
import { Sphere3D } from '../../../mol-math/geometry';

const OrderColors = ColorMap({
    '2': ColorNames.deepskyblue,
    '3': ColorNames.lime,
    'N': ColorNames.red,
});
const OrderColorsLegend = TableLegend(Object.keys(OrderColors).map(name => {
    return [name, (OrderColors as any)[name] as Color] as [string, Color];
}));

function axesColorHelp(value: { name: string, params: {} }) {
    return value.name === 'byOrder'
        ? { description: 'Color axes by their order', legend: OrderColorsLegend }
        : {};
}

const SharedParams = {
    ...Mesh.Params,
    scale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 }),
};

const AxesParams = {
    ...SharedParams,
    axesColor: PD.MappedStatic('byOrder', {
        byOrder: PD.EmptyGroup(),
        uniform: PD.Group({
            colorValue: PD.Color(ColorNames.orange),
        }, { isFlat: true })
    }, { help: axesColorHelp }),
};
type AxesParams = typeof AxesParams

const CageParams = {
    ...SharedParams,
    cageColor: PD.Color(ColorNames.orange),
};
type CageParams = typeof CageParams

const AssemblySymmetryVisuals = {
    // cage should come before 'axes' so that the representative loci uses the cage shape
    'cage': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CageParams>) => ShapeRepresentation(getCageShape, Mesh.Utils, { modifyState: s => ({ ...s, markerActions: MarkerActions.Highlighting }) }),
    'axes': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, AxesParams>) => ShapeRepresentation(getAxesShape, Mesh.Utils, { modifyState: s => ({ ...s, markerActions: MarkerActions.Highlighting }) }),
};

export const AssemblySymmetryParams = {
    ...AxesParams,
    ...CageParams,
    visuals: PD.MultiSelect(['axes', 'cage'], PD.objectToOptions(AssemblySymmetryVisuals)),
};
export type AssemblySymmetryParams = typeof AssemblySymmetryParams
export type AssemblySymmetryProps = PD.Values<AssemblySymmetryParams>

//

function getAssemblyName(s: Structure) {
    const id = s.units[0].conformation.operator.assembly?.id || '';
    return isInteger(id) ? `Assembly ${id}` : id;
}

const t = Mat4.identity();
const tmpV = Vec3();
const tmpCenter = Vec3();
const tmpScale = Vec3();

const getOrderPrimitive = memoize1((order: number): Primitive | undefined => {
    if (order < 2) {
        return Prism(polygon(48, false));
    } else if (order === 2) {
        const lens = Prism(polygon(48, false));
        const m = Mat4.identity();
        Mat4.scale(m, m, Vec3.create(1, 0.35, 1));
        transformPrimitive(lens, m);
        return lens;
    } else if (order === 3) {
        return Wedge();
    } else {
        return Prism(polygon(order, false));
    }
});

function getAxesMesh(data: AssemblySymmetryValue, props: PD.Values<AxesParams>, mesh?: Mesh) {
    const { scale } = props;

    const { rotation_axes } = data;
    if (!AssemblySymmetry.isRotationAxes(rotation_axes)) return Mesh.createEmpty(mesh);

    const { start, end } = rotation_axes[0];
    const radius = (Vec3.distance(start, end) / 500) * scale;

    Vec3.set(tmpScale, radius * 7, radius * 7, radius * 0.4);

    const cylinderProps = { radiusTop: radius, radiusBottom: radius };
    const builderState = MeshBuilder.createState(256, 128, mesh);

    builderState.currentGroup = 0;
    Vec3.scale(tmpCenter, Vec3.add(tmpCenter, start, end), 0.5);

    for (let i = 0, il = rotation_axes.length; i < il; ++i) {
        const { order, start, end } = rotation_axes[i];
        builderState.currentGroup = i;
        addCylinder(builderState, start, end, 1, cylinderProps);

        const primitive = getOrderPrimitive(order);
        if (primitive) {
            Vec3.scale(tmpCenter, Vec3.add(tmpCenter, start, end), 0.5);
            if (Vec3.dot(Vec3.unitY, Vec3.sub(tmpV, start, tmpCenter)) === 0) {
                Mat4.targetTo(t, start, tmpCenter, Vec3.unitY);
            } else {
                Mat4.targetTo(t, start, tmpCenter, Vec3.unitX);
            }
            Mat4.scale(t, t, tmpScale);

            Mat4.setTranslation(t, start);
            MeshBuilder.addPrimitive(builderState, t, primitive);
            Mat4.setTranslation(t, end);
            MeshBuilder.addPrimitive(builderState, t, primitive);
        }
    }
    return MeshBuilder.getMesh(builderState);
}

function getAxesShape(ctx: RuntimeContext, data: Structure, props: AssemblySymmetryProps, shape?: Shape<Mesh>) {
    const assemblySymmetry = AssemblySymmetryProvider.get(data).value!;
    const geo = getAxesMesh(assemblySymmetry, props, shape && shape.geometry);
    const getColor = (groupId: number) => {
        if (props.axesColor.name === 'byOrder') {
            const { rotation_axes } = assemblySymmetry;
            const order = rotation_axes![groupId]?.order;
            if (order === 2) return OrderColors[2];
            else if (order === 3) return OrderColors[3];
            else return OrderColors.N;
        } else {
            return props.axesColor.params.colorValue;
        }
    };
    const getLabel = (groupId: number) => {
        const { type, symbol, kind, rotation_axes } = assemblySymmetry;
        const order = rotation_axes![groupId]?.order;
        return [
            `<small>${data.model.entryId}</small>`,
            `<small>${getAssemblyName(data)}</small>`,
            `Axis ${groupId + 1} with Order ${order} of ${type} ${kind} (${symbol})`
        ].join(' | ');
    };
    return Shape.create('Axes', data, geo, getColor, () => 1, getLabel);
}

//

const getSymbolCage = memoize1((symbol: string): Cage | undefined => {
    if (symbol.startsWith('D') || symbol.startsWith('C')) {
        // z axis is prism axis, x/y axes cut through edge midpoints
        const fold = parseInt(symbol.substr(1));
        let cage: Cage;
        if (fold === 2) {
            cage = PrismCage(polygon(4, false));
        } else if (fold === 3) {
            cage = WedgeCage();
        } else if (fold > 3) {
            cage = PrismCage(polygon(fold, false));
        } else {
            return;
        }
        if (fold % 2 === 0) {
            return cage;
        } else {
            const m = Mat4.identity();
            Mat4.rotate(m, m, 1 / fold * Math.PI / 2, Vec3.unitZ);
            return transformCage(cloneCage(cage), m);
        }
    } else if (symbol === 'O') {
        // x/y/z axes cut through order 4 vertices
        return OctahedronCage();
    } else if (symbol === 'I') {
        // z axis cut through order 5 vertex
        // x axis cut through edge midpoint
        const cage = IcosahedronCage();
        const m = Mat4.identity();
        Mat4.rotate(m, m, degToRad(31.7), Vec3.unitX);
        return transformCage(cloneCage(cage), m);
    } else if (symbol === 'T') {
        // x/y/z axes cut through edge midpoints
        return TetrahedronCage();
    }
});

function getSymbolScale(symbol: string) {
    if (symbol.startsWith('D') || symbol.startsWith('C')) {
        return 0.75;
    } else if (symbol === 'O') {
        return 1.2;
    } else if (symbol === 'I') {
        return 0.25;
    } else if (symbol === 'T') {
        return 0.8;
    }
    return 1;
}

function setSymbolTransform(t: Mat4, symbol: string, axes: AssemblySymmetry.RotationAxes, size: number, structure: Structure) {
    const eye = Vec3();
    const target = Vec3();
    const dir = Vec3();
    const up = Vec3();
    let pair: Mutable<AssemblySymmetry.RotationAxes> | undefined = undefined;

    if (symbol.startsWith('C')) {
        pair = [axes[0]];
    } else if (symbol.startsWith('D')) {
        const fold = parseInt(symbol.substr(1));
        if (fold === 2) {
            pair = axes.filter(a => a.order === 2);
        } else if (fold >= 3) {
            const aN = axes.filter(a => a.order === fold)[0];
            const a2 = axes.filter(a => a.order === 2)[1];
            pair = [aN, a2];
        }
    } else if (symbol === 'O') {
        pair = axes.filter(a => a.order === 4);
    } else if (symbol === 'I') {
        const a5 = axes.filter(a => a.order === 5)[0];
        const a5dir = Vec3.sub(Vec3(), a5.end, a5.start);
        pair = [a5];
        for (const a of axes.filter(a => a.order === 3)) {
            const d = radToDeg(Vec3.angle(Vec3.sub(up, a.end, a.start), a5dir));
            if (!pair[1] && (equalEps(d, 100.81, 0.1) || equalEps(d, 79.19, 0.1))) {
                pair[1] = a;
                break;
            }
        }
    } else if (symbol === 'T') {
        pair = axes.filter(a => a.order === 2);
    }

    Mat4.setIdentity(t);
    if (pair) {
        const [aA, aB] = pair;
        Vec3.scale(eye, Vec3.add(eye, aA.end, aA.start), 0.5);
        Vec3.copy(target, aA.end);
        if (aB) {
            Vec3.sub(up, aB.end, aB.start);
            Vec3.sub(dir, eye, target);
            if (Vec3.dot(dir, up) < 0) Vec3.negate(up, up);
            Mat4.targetTo(t, eye, target, up);

            if (symbol.startsWith('D')) {
                const { sphere } = structure.lookup3d.boundary;
                let sizeXY = (sphere.radius * 2) * 0.8; // fallback for missing extrema
                if (Sphere3D.hasExtrema(sphere)) {
                    const n = Mat3.directionTransform(Mat3(), t);
                    const dirs = unitCircleDirections.map(d => Vec3.transformMat3(Vec3(), d, n));
                    sizeXY = getMaxProjectedDistance(sphere.extrema, dirs, sphere.center) * 1.6;
                }
                Mat4.scale(t, t, Vec3.create(sizeXY, sizeXY, Vec3.distance(aA.start, aA.end) * 0.9));
            } else {
                Mat4.scaleUniformly(t, t, size * getSymbolScale(symbol));
            }
        } else {
            if (Vec3.dot(Vec3.unitY, Vec3.sub(tmpV, aA.end, aA.start)) === 0) {
                Vec3.copy(up, Vec3.unitY);
            } else {
                Vec3.copy(up, Vec3.unitX);
            }
            Mat4.targetTo(t, eye, target, up);

            const { sphere } = structure.lookup3d.boundary;
            let sizeXY = (sphere.radius * 2) * 0.8; // fallback for missing extrema
            if (Sphere3D.hasExtrema(sphere)) {
                const n = Mat3.directionTransform(Mat3(), t);
                const dirs = unitCircleDirections.map(d => Vec3.transformMat3(Vec3(), d, n));
                sizeXY = getMaxProjectedDistance(sphere.extrema, dirs, sphere.center);
            }
            Mat4.scale(t, t, Vec3.create(sizeXY, sizeXY, size * 0.9));
        }
    }
}

const unitCircleDirections = (function() {
    const dirs: Vec3[] = [];
    const circle = polygon(12, false, 1);
    for (let i = 0, il = circle.length; i < il; i += 3) {
        dirs.push(Vec3.fromArray(Vec3(), circle, i));
    }
    return dirs;
})();
const tmpProj = Vec3();

function getMaxProjectedDistance(points: Vec3[], directions: Vec3[], center: Vec3) {
    let maxDist = 0;
    for (const p of points) {
        for (const d of directions) {
            Vec3.projectPointOnVector(tmpProj, p, d, center);
            const dist = Vec3.distance(tmpProj, center);
            if (dist > maxDist) maxDist = dist;
        }
    }
    return maxDist;
}

function getCageMesh(data: Structure, props: PD.Values<CageParams>, mesh?: Mesh) {
    const assemblySymmetry = AssemblySymmetryProvider.get(data).value!;
    const { scale } = props;

    const { rotation_axes, symbol } = assemblySymmetry;
    if (!AssemblySymmetry.isRotationAxes(rotation_axes)) return Mesh.createEmpty(mesh);

    const structure = AssemblySymmetry.getStructure(data, assemblySymmetry);

    const cage = getSymbolCage(symbol);
    if (!cage) return Mesh.createEmpty(mesh);

    const { start, end } = rotation_axes[0];
    const size = Vec3.distance(start, end);
    const radius = (size / 500) * scale;

    const builderState = MeshBuilder.createState(256, 128, mesh);
    builderState.currentGroup = 0;
    setSymbolTransform(t, symbol, rotation_axes, size, structure);
    Vec3.scale(tmpCenter, Vec3.add(tmpCenter, start, end), 0.5);
    Mat4.setTranslation(t, tmpCenter);
    MeshBuilder.addCage(builderState, t, cage, radius, 1, 8);

    return MeshBuilder.getMesh(builderState);
}

function getCageShape(ctx: RuntimeContext, data: Structure, props: AssemblySymmetryProps, shape?: Shape<Mesh>) {
    const assemblySymmetry = AssemblySymmetryProvider.get(data).value!;
    const geo = getCageMesh(data, props, shape && shape.geometry);
    const getColor = (groupId: number) => {
        return props.cageColor;
    };
    const getLabel = (groupId: number) => {
        const { type, symbol, kind } = assemblySymmetry;
        data.model.entryId;
        return [
            `<small>${data.model.entryId}</small>`,
            `<small>${getAssemblyName(data)}</small>`,
            `Cage of ${type} ${kind} (${symbol})`
        ].join(' | ');
    };
    return Shape.create('Cage', data, geo, getColor, () => 1, getLabel);
}

//

export type AssemblySymmetryRepresentation = Representation<Structure, AssemblySymmetryParams>
export function AssemblySymmetryRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, AssemblySymmetryParams>): AssemblySymmetryRepresentation {
    return Representation.createMulti('Assembly Symmetry', ctx, getParams, Representation.StateBuilder, AssemblySymmetryVisuals as unknown as Representation.Def<Structure, AssemblySymmetryParams>);
}