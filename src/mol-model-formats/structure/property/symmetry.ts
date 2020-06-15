/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { Spacegroup, SpacegroupCell, SymmetryOperator } from '../../../mol-math/geometry';
import { Tensor, Vec3, Mat3 } from '../../../mol-math/linear-algebra';
import { Symmetry } from '../../../mol-model/structure/model/properties/symmetry';
import { createAssemblies } from './assembly';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { FormatPropertyProvider } from '../common/property';
import { Table } from '../../../mol-data/db';

export { ModelSymmetry };

namespace ModelSymmetry {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'model_symmetry',
    };

    export const Provider = FormatPropertyProvider.create<Symmetry>(Descriptor);

    type Data = {
        symmetry: Table<mmCIF_Schema['symmetry']>
        cell: Table<mmCIF_Schema['cell']>
        struct_ncs_oper: Table<mmCIF_Schema['struct_ncs_oper']>
        atom_sites: Table<mmCIF_Schema['atom_sites']>
        pdbx_struct_assembly: Table<mmCIF_Schema['pdbx_struct_assembly']>
        pdbx_struct_assembly_gen: Table<mmCIF_Schema['pdbx_struct_assembly_gen']>
        pdbx_struct_oper_list: Table<mmCIF_Schema['pdbx_struct_oper_list']>
    }

    export function fromData(data: Data): Symmetry {
        const assemblies = createAssemblies(data.pdbx_struct_assembly, data.pdbx_struct_assembly_gen, data.pdbx_struct_oper_list);
        const spacegroup = getSpacegroup(data.symmetry, data.cell);
        const isNonStandardCrystalFrame = checkNonStandardCrystalFrame(data.atom_sites, spacegroup);
        return { assemblies, spacegroup, isNonStandardCrystalFrame, ncsOperators: getNcsOperators(data.struct_ncs_oper) };
    }
}

function checkNonStandardCrystalFrame(atom_sites: Table<mmCIF_Schema['atom_sites']>, spacegroup: Spacegroup) {
    if (atom_sites._rowCount === 0) return false;
    // TODO: parse atom_sites transform and check if it corresponds to the toFractional matrix
    return false;
}

function getSpacegroupNameOrNumber(symmetry: Table<mmCIF_Schema['symmetry']>) {
    const groupNumber = symmetry['Int_Tables_number'].value(0);
    const groupName = symmetry['space_group_name_H-M'].value(0);
    if (!symmetry['Int_Tables_number'].isDefined) return groupName;
    if (!symmetry['space_group_name_H-M'].isDefined) return groupNumber;
    return groupName;
}

function getSpacegroup(symmetry: Table<mmCIF_Schema['symmetry']>, cell: Table<mmCIF_Schema['cell']>): Spacegroup {
    if (symmetry._rowCount === 0 || cell._rowCount === 0) return Spacegroup.ZeroP1;
    const nameOrNumber = getSpacegroupNameOrNumber(symmetry);
    const spaceCell = SpacegroupCell.create(nameOrNumber,
        Vec3.create(cell.length_a.value(0), cell.length_b.value(0), cell.length_c.value(0)),
        Vec3.scale(Vec3.zero(), Vec3.create(cell.angle_alpha.value(0), cell.angle_beta.value(0), cell.angle_gamma.value(0)), Math.PI / 180));

    return Spacegroup.create(spaceCell);
}

function getNcsOperators(struct_ncs_oper: Table<mmCIF_Schema['struct_ncs_oper']>) {
    if (struct_ncs_oper._rowCount === 0) return void 0;
    const { id, matrix, vector } = struct_ncs_oper;

    const matrixSpace = mmCIF_Schema.struct_ncs_oper.matrix.space, vectorSpace = mmCIF_Schema.struct_ncs_oper.vector.space;

    const opers: SymmetryOperator[] = [];
    for (let i = 0; i < struct_ncs_oper._rowCount; i++) {
        const m = Tensor.toMat3(Mat3(), matrixSpace, matrix.value(i));
        const v = Tensor.toVec3(Vec3(), vectorSpace, vector.value(i));
        if (!SymmetryOperator.checkIfRotationAndTranslation(m, v)) continue;
        // ignore non-identity 'given' NCS operators
        if (struct_ncs_oper.code.value(i) === 'given' && !Mat3.isIdentity(m) && !Vec3.isZero(v)) continue;
        const ncsId = id.value(i);
        opers[opers.length] = SymmetryOperator.ofRotationAndOffset(`ncs_${ncsId}`, m, v, ncsId);
    }
    return opers;
}