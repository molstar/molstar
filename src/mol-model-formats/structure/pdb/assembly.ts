/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifCategory, CifField } from '../../../mol-io/reader/cif';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { Mat4 } from '../../../mol-math/linear-algebra';
import { Tokens } from '../../../mol-io/reader/common/text/tokenizer';

export function parseCryst1(id: string, record: string): CifCategory[] {
    // COLUMNS       DATA TYPE      CONTENTS
    // --------------------------------------------------------------------------------
    //  1 -  6       Record name    "CRYST1"
    //  7 - 15       Real(9.3)      a (Angstroms)
    // 16 - 24       Real(9.3)      b (Angstroms)
    // 25 - 33       Real(9.3)      c (Angstroms)
    // 34 - 40       Real(7.2)      alpha (degrees)
    // 41 - 47       Real(7.2)      beta (degrees)
    // 48 - 54       Real(7.2)      gamma (degrees)
    // 56 - 66       LString        Space group
    // 67 - 70       Integer        Z value

    const get = (s: number, l: number) => (record.substr(s, l) || '').trim();

    const cell: CifCategory.Fields<mmCIF_Schema['cell']> = {
        entry_id: CifField.ofString(id),
        length_a: CifField.ofString(get(6, 9)),
        length_b: CifField.ofString(get(15, 9)),
        length_c: CifField.ofString(get(24, 9)),
        angle_alpha: CifField.ofString(get(33, 7)),
        angle_beta: CifField.ofString(get(40, 7)),
        angle_gamma: CifField.ofString(get(47, 7)),
        Z_PDB: CifField.ofString(get(66, 4)),
        pdbx_unique_axis: CifField.ofString('?')
    };
    const symmetry: CifCategory.Fields<mmCIF_Schema['symmetry']> = {
        entry_id: CifField.ofString(id),
        'space_group_name_H-M': CifField.ofString(get(55, 11)),
        Int_Tables_number: CifField.ofString('?'),
        cell_setting: CifField.ofString('?'),
        space_group_name_Hall: CifField.ofString('?')
    };
    return [CifCategory.ofFields('cell', cell), CifCategory.ofFields('symmetry', symmetry)];
}

interface PdbAssembly {
    id: string,
    details: string,
    groups: { chains: string[], operators: { id: number, matrix: Mat4 }[] }[]
}

function PdbAssembly(id: string, details: string): PdbAssembly {
    return { id, details, groups: [] };
}

export function parseRemark350(lines: Tokens, lineStart: number, lineEnd: number): CifCategory[] {
    const assemblies: PdbAssembly[] = [];

    // Read the assemblies
    let current: PdbAssembly, group: PdbAssembly['groups'][0], matrix: Mat4, operId = 1, asmId = 1;
    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);
    for (let i = lineStart; i < lineEnd; i++) {
        let line = getLine(i);
        if (line.substr(11, 12) === 'BIOMOLECULE:') {
            const id = line.substr(23).trim();
            let details = `Biomolecule ${id}`;
            line = getLine(i + 1);
            if (line.substr(11, 30) !== 'APPLY THE FOLLOWING TO CHAINS:') {
                i++;
                details = line.substr(11).trim();
            }
            current = PdbAssembly(id, details);
            assemblies.push(current);
        } else if (line.substr(13, 5) === 'BIOMT') {
            const biomt = line.split(/\s+/);
            const row = parseInt(line[18]) - 1;

            if (row === 0) {
                matrix = Mat4.identity();
                group!.operators.push({ id: operId++, matrix });
            }

            Mat4.setValue(matrix!, row, 0, parseFloat(biomt[4]));
            Mat4.setValue(matrix!, row, 1, parseFloat(biomt[5]));
            Mat4.setValue(matrix!, row, 2, parseFloat(biomt[6]));
            Mat4.setValue(matrix!, row, 3, parseFloat(biomt[7]));
        } else if (
            line.substr(11, 30) === 'APPLY THE FOLLOWING TO CHAINS:' ||
            line.substr(11, 30) === '                   AND CHAINS:') {

            if (line.substr(11, 5) === 'APPLY') {
                group = { chains: [], operators: [] };
                current!.groups.push(group);
            }

            const chainList = line.substr(41, 30).split(',');
            for (let j = 0, jl = chainList.length; j < jl; ++j) {
                const c = chainList[j].trim();
                if (c) group!.chains.push(c);
            }
        } else if (line.substr(11, 33) === 'APPLYING THE FOLLOWING TO CHAINS:') {
            // variant in older PDB format version
            current = PdbAssembly(`${asmId}`, `Biomolecule ${asmId}`);
            assemblies.push(current);
            asmId += 1;

            group = { chains: [], operators: [] };
            current!.groups.push(group);

            i++;
            line = getLine(i);

            const chainList = line.substr(11, 69).split(',');
            for (let j = 0, jl = chainList.length; j < jl; ++j) {
                const c = chainList[j].trim();
                if (c) group!.chains.push(c);
            }
        }
    }

    if (assemblies.length === 0) return [];

    // Generate CIF

    // pdbx_struct_assembly
    const pdbx_struct_assembly: CifCategory.SomeFields<mmCIF_Schema['pdbx_struct_assembly']> = {
        id: CifField.ofStrings(assemblies.map(a => a.id)),
        details: CifField.ofStrings(assemblies.map(a => a.details))
    };


    // pdbx_struct_assembly_gen
    const pdbx_struct_assembly_gen_rows: { [P in keyof CifCategory.Fields<mmCIF_Schema['pdbx_struct_assembly_gen']>]: string }[] = [];
    for (const asm of assemblies) {
        for (const group of asm.groups) {
            pdbx_struct_assembly_gen_rows.push({
                assembly_id: asm.id,
                oper_expression: group.operators.map(o => o.id).join(','),
                asym_id_list: group.chains.join(',')
            });
        }
    }
    const pdbx_struct_assembly_gen: CifCategory.Fields<mmCIF_Schema['pdbx_struct_assembly_gen']> = {
        assembly_id: CifField.ofStrings(pdbx_struct_assembly_gen_rows.map(r => r.assembly_id)),
        oper_expression: CifField.ofStrings(pdbx_struct_assembly_gen_rows.map(r => r.oper_expression)),
        asym_id_list: CifField.ofStrings(pdbx_struct_assembly_gen_rows.map(r => r.asym_id_list))
    };

    // pdbx_struct_oper_list
    const pdbx_struct_oper_list_rows: { [P in keyof CifCategory.Fields<mmCIF_Schema['pdbx_struct_oper_list']>]?: string }[] = [];
    for (const asm of assemblies) {
        for (const group of asm.groups) {
            for (const oper of group.operators) {
                const row = {
                    id: '' + oper.id,
                    type: '?',
                    name: '?',
                    symmetry_operation: '?'
                } as (typeof pdbx_struct_oper_list_rows)[0] as any;
                for (let i = 0; i < 3; i++) {
                    for (let j = 0; j < 3; j++) {
                        row[`matrix[${i + 1}][${j + 1}]`] = '' + Mat4.getValue(oper.matrix, i, j);
                    }
                    row[`vector[${i + 1}]`] = '' + Mat4.getValue(oper.matrix, i, 3);
                }
                pdbx_struct_oper_list_rows.push(row);
            }
        }
    }

    const pdbx_struct_oper_list: CifCategory.SomeFields<mmCIF_Schema['pdbx_struct_oper_list']> = {
        id: CifField.ofStrings(pdbx_struct_oper_list_rows.map(r => r.id!)),
        type: CifField.ofStrings(pdbx_struct_oper_list_rows.map(r => r.type!)),
        name: CifField.ofStrings(pdbx_struct_oper_list_rows.map(r => r.name!)),
        symmetry_operation: CifField.ofStrings(pdbx_struct_oper_list_rows.map(r => r.symmetry_operation!))
    };
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            const k = `matrix[${i + 1}][${j + 1}]`;
            (pdbx_struct_oper_list as any)[k] = CifField.ofStrings(pdbx_struct_oper_list_rows.map(r => (r as any)[k]!));
        }
        const k = `vector[${i + 1}]`;
        (pdbx_struct_oper_list as any)[k] = CifField.ofStrings(pdbx_struct_oper_list_rows.map(r => (r as any)[k]!));
    }

    return [
        CifCategory.ofFields('pdbx_struct_assembly', pdbx_struct_assembly),
        CifCategory.ofFields('pdbx_struct_assembly_gen', pdbx_struct_assembly_gen),
        CifCategory.ofFields('pdbx_struct_oper_list', pdbx_struct_oper_list)
    ];
}

export function parseMtrix(lines: Tokens, lineStart: number, lineEnd: number): CifCategory[] {
    const matrices: Mat4[] = [];
    let matrix: Mat4;

    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);
    for (let i = lineStart; i < lineEnd; i++) {
        let line = getLine(i);

        const ncs = line.split(/\s+/);
        const row = parseInt(line[5]) - 1;

        if (row === 0) {
            matrix = Mat4.identity();
            matrices.push(matrix);
        }

        Mat4.setValue(matrix!, row, 0, parseFloat(ncs[2]));
        Mat4.setValue(matrix!, row, 1, parseFloat(ncs[3]));
        Mat4.setValue(matrix!, row, 2, parseFloat(ncs[4]));
        Mat4.setValue(matrix!, row, 3, parseFloat(ncs[5]));
    }

    if (matrices.length === 0) return [];

    const struct_ncs_oper_rows: { [P in keyof CifCategory.Fields<mmCIF_Schema['struct_ncs_oper']>]?: string }[] = [];
    let id = 1;
    for (const oper of matrices) {
        const row = {
            id: 'ncsop' + (id++),
            code: '.',
            details: '.'
        } as (typeof struct_ncs_oper_rows)[0] as any;
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                row[`matrix[${i + 1}][${j + 1}]`] = '' + Mat4.getValue(oper, i, j);
            }
            row[`vector[${i + 1}]`] = '' + Mat4.getValue(oper, i, 3);
        }
        struct_ncs_oper_rows.push(row);
    }

    const struct_ncs_oper: CifCategory.SomeFields<mmCIF_Schema['struct_ncs_oper']> = {
        id: CifField.ofStrings(struct_ncs_oper_rows.map(r => r.id!)),
        code: CifField.ofStrings(struct_ncs_oper_rows.map(r => r.code!)),
        details: CifField.ofStrings(struct_ncs_oper_rows.map(r => r.details!)),
    };
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            const k = `matrix[${i + 1}][${j + 1}]`;
            (struct_ncs_oper as any)[k] = CifField.ofStrings(struct_ncs_oper_rows.map(r => (r as any)[k]!));
        }
        const k = `vector[${i + 1}]`;
        (struct_ncs_oper as any)[k] = CifField.ofStrings(struct_ncs_oper_rows.map(r => (r as any)[k]!));
    }

    return [CifCategory.ofFields('struct_ncs_oper', struct_ncs_oper)];
}