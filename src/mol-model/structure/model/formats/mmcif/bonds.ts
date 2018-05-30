/**
 * Copyright (c) 2017-2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Model from '../../model'
import { BondType } from '../../types'
import { findEntityIdByAsymId, findAtomIndexByLabelName } from './util'
import { Column } from 'mol-data/db'
import { IntraUnitBonds } from '../../../structure/unit/bonds';

export class StructConn implements IntraUnitBonds.StructConn {
    private _residuePairIndex: Map<string, StructConn.Entry[]> | undefined = void 0;
    private _atomIndex: Map<number, StructConn.Entry[]> | undefined = void 0;

    private static _resKey(rA: number, rB: number) {
        if (rA < rB) return `${rA}-${rB}`;
        return `${rB}-${rA}`;
    }

    private getResiduePairIndex() {
        if (this._residuePairIndex) return this._residuePairIndex;
        this._residuePairIndex = new Map();
        for (const e of this.entries) {
            const ps = e.partners;
            const l = ps.length;
            for (let i = 0; i < l - 1; i++) {
                for (let j = i + i; j < l; j++) {
                    const key = StructConn._resKey(ps[i].residueIndex, ps[j].residueIndex);
                    if (this._residuePairIndex.has(key)) {
                        this._residuePairIndex.get(key)!.push(e);
                    } else {
                        this._residuePairIndex.set(key, [e]);
                    }
                }
            }
        }
        return this._residuePairIndex;
    }

    private getAtomIndex() {
        if (this._atomIndex) return this._atomIndex;
        this._atomIndex = new Map();
        for (const e of this.entries) {
            for (const p of e.partners) {
                const key = p.atomIndex;
                if (this._atomIndex.has(key)) {
                    this._atomIndex.get(key)!.push(e);
                } else {
                    this._atomIndex.set(key, [e]);
                }
            }
        }
        return this._atomIndex;
    }

    private static _emptyEntry = [];

    getResidueEntries(residueAIndex: number, residueBIndex: number): ReadonlyArray<StructConn.Entry> {
        return this.getResiduePairIndex().get(StructConn._resKey(residueAIndex, residueBIndex)) || StructConn._emptyEntry;
    }

    getAtomEntries(atomIndex: number): ReadonlyArray<StructConn.Entry> {
        return this.getAtomIndex().get(atomIndex) || StructConn._emptyEntry;
    }

    constructor(public entries: StructConn.Entry[]) {
    }
}

export namespace StructConn {
    export interface Entry extends IntraUnitBonds.StructConnEntry {
        distance: number,
        order: number,
        flags: number,
        partners: { residueIndex: number, atomIndex: number, symmetry: string }[]
    }

    type StructConnType =
        | 'covale'
        | 'covale_base'
        | 'covale_phosphate'
        | 'covale_sugar'
        | 'disulf'
        | 'hydrog'
        | 'metalc'
        | 'mismat'
        | 'modres'
        | 'saltbr'

    export function create(model: Model): StructConn | undefined {
        if (model.sourceData.kind !== 'mmCIF') return
        const { struct_conn } = model.sourceData.data;
        if (!struct_conn._rowCount) return void 0;

        const { conn_type_id, pdbx_dist_value, pdbx_value_order } = struct_conn;
        const p1 = {
            label_asym_id: struct_conn.ptnr1_label_asym_id,
            label_comp_id: struct_conn.ptnr1_label_comp_id,
            label_seq_id: struct_conn.ptnr1_label_seq_id,
            label_atom_id: struct_conn.ptnr1_label_atom_id,
            label_alt_id: struct_conn.pdbx_ptnr1_label_alt_id,
            ins_code: struct_conn.pdbx_ptnr1_PDB_ins_code,
            symmetry: struct_conn.ptnr1_symmetry
        };
        const p2: typeof p1 = {
            label_asym_id: struct_conn.ptnr2_label_asym_id,
            label_comp_id: struct_conn.ptnr2_label_comp_id,
            label_seq_id: struct_conn.ptnr2_label_seq_id,
            label_atom_id: struct_conn.ptnr2_label_atom_id,
            label_alt_id: struct_conn.pdbx_ptnr2_label_alt_id,
            ins_code: struct_conn.pdbx_ptnr2_PDB_ins_code,
            symmetry: struct_conn.ptnr2_symmetry
        };

        const _p = (row: number, ps: typeof p1) => {
            if (ps.label_asym_id.valueKind(row) !== Column.ValueKind.Present) return void 0;
            const asymId = ps.label_asym_id.value(row)
            const residueIndex = model.atomicHierarchy.findResidueKey(
                findEntityIdByAsymId(model, asymId),
                ps.label_comp_id.value(row),
                asymId,
                ps.label_seq_id.value(row),
                ps.ins_code.value(row)
            );
            if (residueIndex < 0) return void 0;
            const atomName = ps.label_atom_id.value(row);
            // turns out "mismat" records might not have atom name value
            if (!atomName) return void 0;
            const atomIndex = findAtomIndexByLabelName(model, residueIndex, atomName, ps.label_alt_id.value(row));
            if (atomIndex < 0) return void 0;
            return { residueIndex, atomIndex, symmetry: ps.symmetry.value(row) || '1_555' };
        }

        const _ps = (row: number) => {
            const ret = [];
            let p = _p(row, p1);
            if (p) ret.push(p);
            p = _p(row, p2);
            if (p) ret.push(p);
            return ret;
        }

        const entries: StructConn.Entry[] = [];
        for (let i = 0; i < struct_conn._rowCount; i++) {
            const partners = _ps(i);
            if (partners.length < 2) continue;

            const type = conn_type_id.value(i)! as StructConnType;
            const orderType = (pdbx_value_order.value(i) || '').toLowerCase();
            let flags = BondType.Flag.None;
            let order = 1;

            switch (orderType) {
                case 'sing': order = 1; break;
                case 'doub': order = 2; break;
                case 'trip': order = 3; break;
                case 'quad': order = 4; break;
            }

            switch (type) {
                case 'covale':
                case 'covale_base':
                case 'covale_phosphate':
                case 'covale_sugar':
                case 'modres':
                    flags = BondType.Flag.Covalent;
                    break;
                case 'disulf': flags = BondType.Flag.Covalent | BondType.Flag.Sulfide; break;
                case 'hydrog': flags = BondType.Flag.Hydrogen; break;
                case 'metalc': flags = BondType.Flag.MetallicCoordination; break;
                case 'saltbr': flags = BondType.Flag.Ion; break;
            }

            entries.push({ flags, order, distance: pdbx_dist_value.value(i), partners });
        }

        return new StructConn(entries);
    }
}

export class ComponentBondInfo implements IntraUnitBonds.ComponentBondInfo {
    entries: Map<string, ComponentBondInfo.Entry> = new Map();

    newEntry(id: string) {
        let e = new ComponentBondInfo.Entry(id);
        this.entries.set(id, e);
        return e;
    }
}

export namespace ComponentBondInfo {
    export class Entry implements IntraUnitBonds.ComponentBondInfoEntry {
        map: Map<string, Map<string, { order: number, flags: number }>> = new Map();

        add(a: string, b: string, order: number, flags: number, swap = true) {
            let e = this.map.get(a);
            if (e !== void 0) {
                let f = e.get(b);
                if (f === void 0) {
                    e.set(b, { order, flags });
                }
            } else {
                let map = new Map<string, { order: number, flags: number }>();
                map.set(b, { order, flags });
                this.map.set(a, map);
            }

            if (swap) this.add(b, a, order, flags, false);
        }

        constructor(public id: string) {
        }
    }

    export function create(model: Model): ComponentBondInfo | undefined {
        if (model.sourceData.kind !== 'mmCIF') return
        const { chem_comp_bond } = model.sourceData.data;
        if (!chem_comp_bond._rowCount) return void 0;

        let info = new ComponentBondInfo();

        const { comp_id, atom_id_1, atom_id_2, value_order, pdbx_aromatic_flag, _rowCount: rowCount } = chem_comp_bond;

        let entry = info.newEntry(comp_id.value(0)!);

        for (let i = 0; i < rowCount; i++) {

            const id = comp_id.value(i)!;
            const nameA = atom_id_1.value(i)!;
            const nameB = atom_id_2.value(i)!;
            const order = value_order.value(i)!;
            const aromatic = pdbx_aromatic_flag.value(i) === 'Y';

            if (entry.id !== id) {
                entry = info.newEntry(id);
            }

            let flags: number = BondType.Flag.Covalent;
            let ord = 1;
            if (aromatic) flags |= BondType.Flag.Aromatic;
            switch (order.toLowerCase()) {
                case 'doub':
                case 'delo':
                    ord = 2;
                    break;
                case 'trip': ord = 3; break;
                case 'quad': ord = 4; break;
            }

            entry.add(nameA, nameB, ord, flags);
        }

        return info;
    }
}