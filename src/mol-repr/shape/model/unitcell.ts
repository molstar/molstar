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

const CellRef = {
    origin: 'Origin',
    model: 'Model'
};

const CellAttachment = {
    corner: 'Corner',
    center: 'Center'
};

const CellParams = {
    ...Mesh.Params,
    cellColor: PD.Color(ColorNames.orange),
    cellScale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 }),
    ref: PD.Select('model', PD.objectToOptions(CellRef), { isEssential: true }),
    attachment: PD.Select('corner', PD.objectToOptions(CellAttachment), { isEssential: true }),
};
type CellParams = typeof CellParams

const UnitcellVisuals = {
    'mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<UnitcellData, CellParams>) => ShapeRepresentation(getUnitcellShape, Mesh.Utils),
};

export const UnitcellParams = {
    ...CellParams
};
export type UnitcellParams = typeof UnitcellParams
export type UnitcellProps = PD.Values<UnitcellParams>

function getUnitcellMesh(data: UnitcellData, props: UnitcellProps, mesh?: Mesh) {
    const state = MeshBuilder.createState(256, 128, mesh);

    const { fromFractional } = data.symmetry.spacegroup.cell;

    Vec3.copy(tmpRef, data.ref);
    if (props.attachment === 'center') {
        Vec3.trunc(tmpRef, tmpRef);
        Vec3.subScalar(tmpRef, tmpRef, 0.5);
    } else {
        Vec3.floor(tmpRef, tmpRef);
    }
    Mat4.fromTranslation(tmpTranslate, tmpRef);
    const cellCage = transformCage(cloneCage(unitCage), tmpTranslate);

    const radius = (Math.cbrt(data.symmetry.spacegroup.cell.volume) / 300) * props.cellScale;
    state.currentGroup = 1;
    MeshBuilder.addCage(state, fromFractional, cellCage, radius, 2, 20);

    const sphere = Sphere3D.fromDimensionsAndTransform(Sphere3D(), Vec3.unit, fromFractional);
    Vec3.transformMat4(tmpRef, tmpRef, fromFractional);
    Sphere3D.translate(sphere, sphere, tmpRef);
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

export function getUnitcellData(model: Model, symmetry: Symmetry, props: UnitcellProps) {
    const ref = Vec3();
    if (props.ref === 'model') {
        Vec3.transformMat4(ref, Model.getCenter(model), symmetry.spacegroup.cell.toFractional);
    }
    return { symmetry, ref };
}

export type UnitcellRepresentation = Representation<UnitcellData, UnitcellParams>
export function UnitcellRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<UnitcellData, UnitcellParams>): UnitcellRepresentation {
    return Representation.createMulti('Unit Cell', ctx, getParams, Representation.StateBuilder, UnitcellVisuals as unknown as Representation.Def<UnitcellData, UnitcellParams>);
}