/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { StringBuilder, deepEqual } from '../../../mol-util';
import Writer from '../../../mol-io/writer/writer';
import { Encoder, Category, Field } from '../../../mol-io/writer/cif/encoder';
import { ComponentBond } from '../../../mol-model-formats/structure/property/bonds/comp';

interface Atom {
    label_atom_id: string,
    Cartn_x: number,
    Cartn_y: number,
    Cartn_z: number,
    type_symbol: string
}

function Atom(partial: any): Atom {
    return { 
        label_atom_id: partial.label_atom_id,
        Cartn_x: partial.Cartn_x,
        Cartn_y: partial.Cartn_y,
        Cartn_z: partial.Cartn_z,
        type_symbol: partial.type_symbol
    }
}

export abstract class LigandEncoder implements Encoder<string> {
    protected builder: StringBuilder;
    protected meta: StringBuilder;
    protected componentData: ComponentBond;
    protected error = false;
    protected encoded = false;
    readonly isBinary = false;
    binaryEncodingProvider = void 0;

    abstract encode(): void;

    abstract _writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx): void;

    protected abstract writeFullCategory<Ctx>(sb: StringBuilder, category: Category<Ctx>, context?: Ctx): void;

    writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx) {
        if (this.encoded) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (this.metaInformation && (category.name === 'model_server_result' || category.name === 'model_server_params' || category.name === 'model_server_stats')) {
            this.writeFullCategory(this.meta, category, context);
            return;
        }

        // if error: force writing of meta information
        if (category.name === 'model_server_error') {
            this.writeFullCategory(this.meta, category, context);
            this.error = true;
            return;
        }

        // only care about atom_site category when writing SDF
        if (category.name !== 'atom_site') {
            return;
        }

        this._writeCategory(category, context);
    }

    setComponentBondData(componentData: ComponentBond) {
        this.componentData = componentData;
    }

    writeTo(stream: Writer) {
        const chunks = StringBuilder.getChunks(this.builder);
        for (let i = 0, _i = chunks.length; i < _i; i++) {
            stream.writeString(chunks[i]);
        }
    }

    getSize() {
        return StringBuilder.getSize(this.builder);
    }

    getData() {
        return StringBuilder.getString(this.builder);
    }

    protected getAtoms<Ctx>(instance: Category.Instance<Ctx>, source: any): Atom[] {
        const sortedFields = this.getSortedFields(instance, ['Cartn_x', 'Cartn_y', 'Cartn_z', 'type_symbol']);
        const label_atom_id = this.getField(instance, 'label_atom_id');

        // all of this is used to ensure that only 1 residue is written
        const auxiliaryFields = this.getSortedFields(instance, ['label_seq_id', 'label_asym_id', 'pdbx_PDB_ins_code', 'pdbx_PDB_model_num']);

        return this._getAtoms(source, sortedFields, label_atom_id, auxiliaryFields);
    }

    private _getAtoms(source: any, fields: Field<any, any>[], label_atom_id: Field<any, any>, auxiliaryFields: Field<any, any>[]): Atom[] {
        const atoms: Atom[] = [];
        let index = 0;
        let id: (string | number)[] | undefined = void 0;

        // is outer loop even needed?
        l: for (let _c = 0; _c < source.length; _c++) {
            const src = source[_c];
            const data = src.data;

            if (src.rowCount === 0) continue;

            const it = src.keys();
            while (it.hasNext)  {
                const key = it.move();

                // ensure only a single residue is written
                // TODO fairly certain this can be done at query level
                if (!id) {
                    id = auxiliaryFields.map(af => af.value(key, data, 0));
                } else {
                    if (!deepEqual(id, auxiliaryFields.map(af => af.value(key, data, index)))) {
                        break l;
                    }
                }

                const lai = label_atom_id.value(key, data, index) as string;
                const label = this.getLabel(lai);
                if (this.skipHydrogen(label)) {
                    index++;
                    continue;
                }
                const a: { [k: string]: (string | number) } = { 'label_atom_id': lai };

                for (let _f = 0, _fl = fields.length; _f < _fl; _f++) {
                    const f: Field<any, any> = fields[_f]!;
                    a[f.name] = f.value(key, data, index);
                }

                atoms.push(Atom(a));
                index++;
            }
        }

        return atoms;
    }

    protected skipHydrogen(label: string) {
        if (this.hydrogens) {
            return false;
        }
        return label.startsWith('H');
    }

    protected getLabel(s: string) {
        return s.replace(/[^A-Z]+/g, '');
    }

    private getSortedFields<Ctx>(instance: Category.Instance<Ctx>, names: string[]) {
        return names.map(n => this.getField(instance, n));
    }

    private getField<Ctx>(instance: Category.Instance<Ctx>, name: string) {
        return instance.fields.find(f => f.name === name)!;
    }

    protected getName<Ctx>(instance: Category.Instance<Ctx>, source: any): string {
        const label_comp_id = this.getField(instance, 'label_comp_id');
        return label_comp_id.value(source[0].keys().move(), source[0].data, 0) as string;
    }

    startDataBlock() {}

    setFilter() {}

    setFormatter() {}

    isCategoryIncluded() {
        return true;
    }

    constructor(readonly encoder: string, readonly metaInformation: boolean, readonly hydrogens: boolean) {
        this.builder = StringBuilder.create();
        this.meta = StringBuilder.create();
    }
}