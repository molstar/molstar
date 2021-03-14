/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '../../../structure/element';
import { StructureProperties } from '../../../structure/properties';

export interface ResidueSetEntry {
    label_asym_id: string,
    label_comp_id: string,
    label_seq_id: number,
    alt_id: string,
    ins_code: string,
    // 1_555 by default
    operator_name?: string
}

export class ResidueSet {
    private index = new Map<string, Map<number, ResidueSetEntry[]>>();
    private list: ResidueSetEntry[] = [];

    add(entry: ResidueSetEntry) {
        let root = this.index.get(entry.label_asym_id);
        if (!root) {
            root = new Map();
            this.index.set(entry.label_asym_id, root);
        }

        let entries = root.get(entry.label_seq_id);
        if (!entries) {
            entries = [];
            root.set(entry.label_seq_id, entries);
        }

        const exists = this._find(entry, entries);
        if (!exists) {
            entries.push(entry);
            this.list.push(entry);
            return true;
        }

        return false;
    }

    private _asym_id = StructureProperties.chain.label_asym_id;
    private _seq_id = StructureProperties.residue.label_seq_id;
    private _comp_id = StructureProperties.atom.label_comp_id;
    private _alt_id = StructureProperties.atom.label_alt_id;
    private _ins_code = StructureProperties.residue.pdbx_PDB_ins_code;
    private _op_name = StructureProperties.unit.operator_name;

    has(loc: StructureElement.Location) {
        const asym_id = this._asym_id(loc);
        if (!this.index.has(asym_id)) return false;
        const root = this.index.get(asym_id)!;
        const seq_id = this._seq_id(loc);
        if (!root.has(seq_id)) return false;
        const entries = root.get(seq_id)!;

        const comp_id = this._comp_id(loc);
        const alt_id = this._alt_id(loc);
        const ins_code = this._ins_code(loc);
        const op_name = this._op_name(loc) ?? '1_555';

        for (const e of entries) {
            if (e.label_comp_id !== comp_id || e.alt_id !== alt_id || e.ins_code !== ins_code) continue;
            if (this.checkOperator && (e.operator_name ?? '1_555') !== op_name) continue;

            return true;
        }

        return false;
    }

    private _find(entry: ResidueSetEntry, xs: ResidueSetEntry[]) {
        for (const e of xs) {
            if (e.label_comp_id !== entry.label_comp_id || e.alt_id !== entry.alt_id || e.ins_code !== entry.ins_code) continue;
            if (this.checkOperator && (e.operator_name ?? '1_555') !== (entry.operator_name ?? '1_555')) continue;

            return true;
        }

        return false;
    }

    constructor(private checkOperator = false) {

    }
}