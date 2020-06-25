/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { StringBuilder, deepEqual } from '../../../mol-util';
import Writer from '../../../mol-io/writer/writer';
import { Encoder, Category, Field } from '../../../mol-io/writer/cif/encoder';
import { ComponentBond } from '../../../mol-model-formats/structure/property/bonds/comp';

export interface Atom {
    id: string,
    x: number,
    y: number,
    z: number,
    type_symbol: string
}

export abstract class LigandEncoder implements Encoder<string> {
    protected builder: StringBuilder;
    protected componentData: ComponentBond;
    readonly isBinary = false;
    binaryEncodingProvider = void 0;

    abstract writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx): void;

    abstract encode(): void;

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

        return this.getAtomsInternal(source, sortedFields, label_atom_id, auxiliaryFields);
    }

    private getAtomsInternal(source: any, fields: Field<any, any>[], label_atom_id: Field<any, any>, auxiliaryFields: Field<any, any>[]): Atom[] {
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
                const d: (string | number)[] = [lai];

                for (let _f = 0, _fl = fields.length; _f < _fl; _f++) {
                    const f: Field<any, any> = fields[_f]!;
                    d.push(f.value(key, data, index));
                }

                atoms.push({ id: d[0] as string, x: d[1] as number, y: d[2] as number, z: d[3] as number, type_symbol: d[4] as string});
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

    constructor(readonly encoder: string, readonly hydrogens: boolean) {
        this.builder = StringBuilder.create();
    }
}