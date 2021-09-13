/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { StringBuilder } from '../../mol-util';
import { Writer } from './writer';
import { Encoder, Category, Field } from './cif/encoder';
import { ComponentAtom } from '../../mol-model-formats/structure/property/atoms/chem_comp';
import { ComponentBond } from '../../mol-model-formats/structure/property/bonds/chem_comp';

interface Atom {
    Cartn_x: number,
    Cartn_y: number,
    Cartn_z: number,
    type_symbol: string,
    index: number
}

function Atom(partial: any): Atom {
    return { ...partial };
}

export abstract class LigandEncoder implements Encoder<string> {
    protected builder: StringBuilder;
    protected meta: StringBuilder;
    protected componentAtomData: ComponentAtom;
    protected componentBondData: ComponentBond;
    protected error = false;
    protected encoded = false;
    readonly isBinary = false;
    binaryEncodingProvider = void 0;

    abstract encode(): void;

    protected abstract _writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx): void;

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

    setComponentAtomData(componentAtomData: ComponentAtom) {
        this.componentAtomData = componentAtomData;
    }

    setComponentBondData(componentBondData: ComponentBond) {
        this.componentBondData = componentBondData;
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

    protected getAtoms<Ctx>(instance: Category.Instance<Ctx>, source: any): Map<string, Atom> {
        const sortedFields = this.getSortedFields(instance, ['Cartn_x', 'Cartn_y', 'Cartn_z']);
        const label_atom_id = this.getField(instance, 'label_atom_id');
        const type_symbol = this.getField(instance, 'type_symbol');
        return this._getAtoms(source, sortedFields, label_atom_id, type_symbol);
    }

    private _getAtoms(source: any, fields: Field<any, any>[], label_atom_id: Field<any, any>, type_symbol: Field<any, any>): Map<string, Atom> {
        const atoms = new Map<string, Atom>();
        let index = 0;

        // is outer loop even needed?
        for (let _c = 0; _c < source.length; _c++) {
            const src = source[_c];
            const data = src.data;

            if (src.rowCount === 0) continue;

            const it = src.keys();
            while (it.hasNext) {
                const key = it.move();

                const lai = label_atom_id.value(key, data, index) as string;
                const ts = type_symbol.value(key, data, index) as string;
                if (this.skipHydrogen(ts)) {
                    index++;
                    continue;
                }
                const a: { [k: string]: (string | number) } = {};

                for (let _f = 0, _fl = fields.length; _f < _fl; _f++) {
                    const f: Field<any, any> = fields[_f]!;
                    a[f.name] = f.value(key, data, index);
                }
                a[type_symbol.name] = ts;
                a['index'] = index;

                atoms.set(lai, Atom(a));
                index++;
            }
        }

        return atoms;
    }

    protected skipHydrogen(type_symbol: string) {
        if (this.hydrogens) {
            return false;
        }
        return type_symbol === 'H';
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