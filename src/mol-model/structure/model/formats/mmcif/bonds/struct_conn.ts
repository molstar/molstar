/**
 * Copyright (c) 2017-2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../model'
import { Element, Structure } from '../../../../structure'
import { LinkType } from '../../../types'
import { findEntityIdByAsymId, findAtomIndexByLabelName } from '../util'
import { Column } from 'mol-data/db'
import { ModelPropertyDescriptor } from '../../../properties/custom';
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';
import { SortedArray } from 'mol-data/int';
import { CifWriter } from 'mol-io/writer/cif'

export interface StructConn {
    getResidueEntries(residueAIndex: number, residueBIndex: number): ReadonlyArray<StructConn.Entry>,
    getAtomEntries(atomIndex: number): ReadonlyArray<StructConn.Entry>,
    readonly entries: ReadonlyArray<StructConn.Entry>
}

export namespace StructConn {
    export const Descriptor: ModelPropertyDescriptor = {
        isStatic: true,
        name: 'struct_conn',
        cifExport: {
            categoryNames: ['struct_conn'],
            categoryProvider(ctx) {
                const struct_conn = getStructConn(ctx.model);
                if (!struct_conn) return [];

                const strConn = get(ctx.model);
                if (!strConn || strConn.entries.length === 0) return [];

                const foundAtoms = new Set<Element>();
                const indices: number[] = [];
                for (const entry of strConn.entries) {
                    const { partners } = entry;
                    let hasAll = true;
                    for (let i = 0, _i = partners.length; i < _i; i++) {
                        const atom = partners[i].atomIndex;
                        if (foundAtoms.has(atom)) continue;
                        if (hasAtom(ctx.structure, atom)) {
                            foundAtoms.add(atom);
                        } else {
                            hasAll = false;
                            break;
                        }
                    }
                    if (hasAll) {
                        indices[indices.length] = entry.rowIndex;
                    }
                }

                return [
                    () => CifWriter.Category.ofTable('struct_conn', struct_conn, indices)
                ];
            }
        }
    }

    function hasAtom({ units }: Structure, element: Element) {
        for (let i = 0, _i = units.length; i < _i; i++) {
            if (SortedArray.indexOf(units[i].elements, element) >= 0) return true;
        }
        return false;
    }

    function _resKey(rA: number, rB: number) {
        if (rA < rB) return `${rA}-${rB}`;
        return `${rB}-${rA}`;
    }
    const _emptyEntry: Entry[] = [];

    class StructConnImpl implements StructConn {
        private _residuePairIndex: Map<string, StructConn.Entry[]> | undefined = void 0;
        private _atomIndex: Map<number, StructConn.Entry[]> | undefined = void 0;

        private getResiduePairIndex() {
            if (this._residuePairIndex) return this._residuePairIndex;
            this._residuePairIndex = new Map();
            for (const e of this.entries) {
                const ps = e.partners;
                const l = ps.length;
                for (let i = 0; i < l - 1; i++) {
                    for (let j = i + i; j < l; j++) {
                        const key = _resKey(ps[i].residueIndex, ps[j].residueIndex);
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


        getResidueEntries(residueAIndex: number, residueBIndex: number): ReadonlyArray<StructConn.Entry> {
            return this.getResiduePairIndex().get(_resKey(residueAIndex, residueBIndex)) || _emptyEntry;
        }

        getAtomEntries(atomIndex: number): ReadonlyArray<StructConn.Entry> {
            return this.getAtomIndex().get(atomIndex) || _emptyEntry;
        }

        constructor(public entries: StructConn.Entry[]) {
        }
    }

    export interface Entry {
        rowIndex: number,
        distance: number,
        order: number,
        flags: number,
        partners: { residueIndex: number, atomIndex: Element, symmetry: string }[]
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

    export function attachFromMmCif(model: Model): boolean {
        if (model.customProperties.has(Descriptor)) return true;
        if (model.sourceData.kind !== 'mmCIF') return false;
        const { struct_conn } = model.sourceData.data;
        if (struct_conn._rowCount === 0) return false;
        model.customProperties.add(Descriptor);
        model._staticPropertyData.__StructConnData__ = struct_conn;
        return true;
    }

    function getStructConn(model: Model) {
        return model._staticPropertyData.__StructConnData__ as mmCIF_Database['struct_conn'];
    }

    export const PropName = '__StructConn__';
    export function get(model: Model): StructConn | undefined {
        if (model._staticPropertyData[PropName]) return model._staticPropertyData[PropName];
        if (!model.customProperties.has(Descriptor)) return void 0;

        const struct_conn = getStructConn(model);

        const { conn_type_id, pdbx_dist_value, pdbx_value_order } = struct_conn;
        const p1 = {
            label_asym_id: struct_conn.ptnr1_label_asym_id,
            label_comp_id: struct_conn.ptnr1_label_comp_id,
            label_seq_id: struct_conn.ptnr1_label_seq_id,
            auth_seq_id: struct_conn.ptnr1_auth_seq_id,
            label_atom_id: struct_conn.ptnr1_label_atom_id,
            label_alt_id: struct_conn.pdbx_ptnr1_label_alt_id,
            ins_code: struct_conn.pdbx_ptnr1_PDB_ins_code,
            symmetry: struct_conn.ptnr1_symmetry
        };
        const p2: typeof p1 = {
            label_asym_id: struct_conn.ptnr2_label_asym_id,
            label_comp_id: struct_conn.ptnr2_label_comp_id,
            label_seq_id: struct_conn.ptnr2_label_seq_id,
            auth_seq_id: struct_conn.ptnr2_auth_seq_id,
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
                asymId,
                ps.label_comp_id.value(row),
                ps.auth_seq_id.value(row),
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
            let flags = LinkType.Flag.None;
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
                    flags = LinkType.Flag.Covalent;
                    break;
                case 'disulf': flags = LinkType.Flag.Covalent | LinkType.Flag.Sulfide; break;
                case 'hydrog': flags = LinkType.Flag.Hydrogen; break;
                case 'metalc': flags = LinkType.Flag.MetallicCoordination; break;
                case 'saltbr': flags = LinkType.Flag.Ion; break;
            }

            entries.push({ rowIndex: i, flags, order, distance: pdbx_dist_value.value(i), partners });
        }

        const ret = new StructConnImpl(entries);
        model._staticPropertyData[PropName] = ret;
        return ret;
    }
}