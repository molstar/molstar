/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, Symmetry } from '../../../mol-model/structure';
import { ShapeRepresentation } from '../representation';
import { Shape } from '../../../mol-model/shape';
import { ColorNames } from '../../../mol-util/color/names';
import { RuntimeContext } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { BoxCage } from '../../../mol-geo/primitive/box';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { transformCage, cloneCage } from '../../../mol-geo/primitive/cage';
import { Sphere3D } from '../../../mol-math/geometry';
import { RepresentationParamsGetter, Representation, RepresentationContext } from '../../representation';

const translate05 = Mat4.fromTranslation(Mat4(), Vec3.create(0.5, 0.5, 0.5));
const unitCage = transformCage(cloneCage(BoxCage()), translate05);

const tmpRef = Vec3();
const tmpTranslate = Mat4();

interface UnitcellData {
    symmetry: Symmetry
    ref: Vec3
}

const CellParams = {
    ...Mesh.Params,
    cellColor: PD.Color(ColorNames.orange),
    cellScale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 })
};
type MeshParams = typeof CellParams

const UnitcellVisuals = {
    'mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<UnitcellData, MeshParams>) => ShapeRepresentation(getUnitcellShape, Mesh.Utils),
};

export const UnitcellParams = {
    ...CellParams
};
export type UnitcellParams = typeof UnitcellParams
export type UnitcellProps = PD.Values<UnitcellParams>

function getUnitcellMesh(data: UnitcellData, props: UnitcellProps, mesh?: Mesh) {
    const state = MeshBuilder.createState(256, 128, mesh);

    const { fromFractional } = data.symmetry.spacegroup.cell;

    Vec3.floor(tmpRef, data.ref);
    Mat4.fromTranslation(tmpTranslate, tmpRef);
    const cellCage = transformCage(cloneCage(unitCage), tmpTranslate);

    const radius = (Math.cbrt(data.symmetry.spacegroup.cell.volume) / 300) * props.cellScale;
    state.currentGroup = 1;
    MeshBuilder.addCage(state, fromFractional, cellCage, radius, 2, 20);

    const cpA = Vec3.create(0, 0, 0);
    Vec3.transformMat4(cpA, Vec3.add(cpA, cpA, tmpRef), fromFractional);
    const cpB = Vec3.create(1, 1, 1);
    Vec3.transformMat4(cpB, Vec3.add(cpB, cpB, tmpRef), fromFractional);
    const cpC = Vec3.create(1, 0, 0);
    Vec3.transformMat4(cpC, Vec3.add(cpC, cpC, tmpRef), fromFractional);
    const cpD = Vec3.create(0, 1, 1);
    Vec3.transformMat4(cpD, Vec3.add(cpD, cpD, tmpRef), fromFractional);

    const cpE = Vec3.create(0, 0, 1);
    Vec3.transformMat4(cpE, Vec3.add(cpE, cpE, tmpRef), fromFractional);
    const cpF = Vec3.create(1, 0, 1);
    Vec3.transformMat4(cpF, Vec3.add(cpF, cpF, tmpRef), fromFractional);
    const cpG = Vec3.create(1, 1, 0);
    Vec3.transformMat4(cpG, Vec3.add(cpG, cpG, tmpRef), fromFractional);
    const cpH = Vec3.create(0, 1, 0);
    Vec3.transformMat4(cpH, Vec3.add(cpH, cpH, tmpRef), fromFractional);

    const center = Vec3();
    Vec3.add(center, cpA, cpB);
    Vec3.scale(center, center, 0.5);
    const d = Math.max(Vec3.distance(cpA, cpB), Vec3.distance(cpC, cpD));
    const sphere = Sphere3D.create(center, d / 2);
    Sphere3D.setExtrema(sphere, [cpA, cpB, cpC, cpD, cpE, cpF, cpG, cpH]);
    Sphere3D.expand(sphere, sphere, radius);

    const m = MeshBuilder.getMesh(state);
    m.setBoundingSphere(sphere);
    return m;
}

function getUnitcellShape(ctx: RuntimeContext, data: UnitcellData, props: UnitcellProps, shape?: Shape<Mesh>) {
    const geo = getUnitcellMesh(data, props, shape && shape.geometry);
    const label = Symmetry.getUnitcellLabel(data.symmetry);
    return Shape.create(label, data, geo, () => props.cellColor, () => 1, () => label);
}

//

export function getUnitcellData(model: Model, symmetry: Symmetry) {
    return {
        symmetry,
        ref: Vec3.transformMat4(Vec3(), Model.getCenter(model), symmetry.spacegroup.cell.toFractional)
    };
}

export type UnitcellRepresentation = Representation<UnitcellData, UnitcellParams>
export function UnitcellRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<UnitcellData, UnitcellParams>): UnitcellRepresentation {
    return Representation.createMulti('Unit Cell', ctx, getParams, Representation.StateBuilder, UnitcellVisuals as unknown as Representation.Def<UnitcellData, UnitcellParams>);
}