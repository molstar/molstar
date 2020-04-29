/**
 * Copyright (c) 2017-2020 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../../mol-model/structure/model/model';
import { Structure } from '../../../../mol-model/structure';
import { BondType } from '../../../../mol-model/structure/model/types';
import { Column, Table } from '../../../../mol-data/db';
import { CustomPropertyDescriptor } from '../../../../mol-model/custom-property';
import { mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif';
import { SortedArray } from '../../../../mol-data/int';
import { CifWriter } from '../../../../mol-io/writer/cif';
import { ElementIndex, ResidueIndex } from '../../../../mol-model/structure/model/indexing';
import { getInterBondOrderFromTable } from '../../../../mol-model/structure/model/properties/atomic/bonds';
import { FormatPropertyProvider } from '../../common/property';

export interface StructConn {
    readonly data: Table<mmCIF_Schema['struct_conn']>
    readonly byAtomIndex: Map<ElementIndex, ReadonlyArray<StructConn.Entry>>
    readonly entries: ReadonlyArray<StructConn.Entry>
}

export namespace StructConn {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'struct_conn',
        cifExport: {
            prefix: '',
            categories: [{
                name: 'struct_conn',
                instance(ctx) {
                    const p = Provider.get(ctx.firstModel);
                    if (!p || p.entries.length === 0) return CifWriter.Category.Empty;

                    const structure = ctx.structures[0];

                    const indices: number[] = [];
                    for (const e of p.entries) {
                        if (hasAtom(structure, e.partnerA.atomIndex) &&
                                hasAtom(structure, e.partnerB.atomIndex)) {
                            indices[indices.length] = e.rowIndex;
                        }
                    }

                    return CifWriter.Category.ofTable(p.data, indices);
                }
            }]
        }
    };

    export const Provider = FormatPropertyProvider.create<StructConn>(Descriptor);

    function hasAtom({ units }: Structure, element: ElementIndex) {
        for (let i = 0, _i = units.length; i < _i; i++) {
            if (SortedArray.indexOf(units[i].elements, element) >= 0) return true;
        }
        return false;
    }

    export function getAtomIndexFromEntries(entries: StructConn['entries']) {
        const m = new Map();
        for (const e of entries) {
            const { partnerA: { atomIndex: iA }, partnerB: { atomIndex: iB } } = e;
            if (m.has(iA)) m.get(iA)!.push(e);
            else m.set(iA, [e]);

            if (m.has(iB)) m.get(iB)!.push(e);
            else m.set(iB, [e]);
        }
        return m;
    }

    export interface Entry {
        rowIndex: number,
        distance: number,
        order: number,
        flags: number,
        partnerA: { residueIndex: ResidueIndex, atomIndex: ElementIndex, symmetry: string },
        partnerB: { residueIndex: ResidueIndex, atomIndex: ElementIndex, symmetry: string }
    }

    export function getEntriesFromStructConn(struct_conn: Table<mmCIF_Schema['struct_conn']>, model: Model): StructConn['entries'] {
        const { conn_type_id, pdbx_dist_value, pdbx_value_order } = struct_conn;
        const p1 = {
            label_asym_id: struct_conn.ptnr1_label_asym_id,
            label_seq_id: struct_conn.ptnr1_label_seq_id,
            auth_seq_id: struct_conn.ptnr1_auth_seq_id,
            label_atom_id: struct_conn.ptnr1_label_atom_id,
            label_alt_id: struct_conn.pdbx_ptnr1_label_alt_id,
            ins_code: struct_conn.pdbx_ptnr1_PDB_ins_code,
            symmetry: struct_conn.ptnr1_symmetry
        };
        const p2: typeof p1 = {
            label_asym_id: struct_conn.ptnr2_label_asym_id,
            label_seq_id: struct_conn.ptnr2_label_seq_id,
            auth_seq_id: struct_conn.ptnr2_auth_seq_id,
            label_atom_id: struct_conn.ptnr2_label_atom_id,
            label_alt_id: struct_conn.pdbx_ptnr2_label_alt_id,
            ins_code: struct_conn.pdbx_ptnr2_PDB_ins_code,
            symmetry: struct_conn.ptnr2_symmetry
        };

        const _p = (row: number, ps: typeof p1) => {
            if (ps.label_asym_id.valueKind(row) !== Column.ValueKind.Present) return void 0;
            const asymId = ps.label_asym_id.value(row);
            const entityIndex = model.atomicHierarchy.index.findEntity(asymId);
            if (entityIndex < 0) return void 0;
            const residueIndex = model.atomicHierarchy.index.findResidue(
                model.entities.data.id.value(entityIndex),
                asymId,
                ps.auth_seq_id.value(row),
                ps.ins_code.value(row)
            );
            if (residueIndex < 0) return void 0;
            const atomName = ps.label_atom_id.value(row);
            // turns out "mismat" records might not have atom name value
            if (!atomName) return void 0;
            const atomIndex = model.atomicHierarchy.index.findAtomOnResidue(residueIndex, atomName, ps.label_alt_id.value(row));
            if (atomIndex < 0) return void 0;
            return { residueIndex, atomIndex, symmetry: ps.symmetry.value(row) };
        };

        const entries: StructConn.Entry[] = [];
        for (let i = 0; i < struct_conn._rowCount; i++) {
            const partnerA = _p(i, p1);
            const partnerB = _p(i, p2);
            if (partnerA === undefined || partnerB === undefined) continue;

            const type = conn_type_id.value(i);
            const orderType = (pdbx_value_order.value(i) || '').toLowerCase();
            let flags = BondType.Flag.None;
            let order = 1;

            switch (orderType) {
                case 'sing': order = 1; break;
                case 'doub': order = 2; break;
                case 'trip': order = 3; break;
                case 'quad': order = 4; break;
                default:
                    order = getInterBondOrderFromTable(
                        struct_conn.ptnr1_label_comp_id.value(i),
                        struct_conn.ptnr1_label_atom_id.value(i),
                        struct_conn.ptnr2_label_comp_id.value(i),
                        struct_conn.ptnr2_label_atom_id.value(i)
                    );
            }

            switch (type) {
                case 'covale':
                    flags = BondType.Flag.Covalent;
                    break;
                case 'disulf': flags = BondType.Flag.Covalent | BondType.Flag.Disulfide; break;
                case 'hydrog':
                    flags = BondType.Flag.HydrogenBond;
                    break;
                case 'metalc': flags = BondType.Flag.MetallicCoordination; break;
            }

            entries.push({
                rowIndex: i, flags, order, distance: pdbx_dist_value.value(i), partnerA, partnerB
            });
        }

        return entries;
    }
}