/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { StringBuilder, deepEqual } from '../../../mol-util';
import Writer from '../writer';
import { Encoder, Category, Field } from '../cif/encoder';
import { getCategoryInstanceData } from '../cif/encoder/util';
import { ComponentBond } from '../../../mol-model-formats/structure/property/bonds/comp';

// specification: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
export class Mol2Encoder implements Encoder<string> {
    private builder: StringBuilder;
    private meta: StringBuilder;
    private encoded = false;
    private error = false;
    private componentData: ComponentBond;
    readonly isBinary = false;
    binaryEncodingProvider = void 0;

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

    startDataBlock() {

    }

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

        // use separate builder because we still need to write Counts and Bonds line
        const ctab = StringBuilder.create();
        const bonds = StringBuilder.create();
        const charges = StringBuilder.create();

        // write Atom block and gather data for Bonds and Charges
        const { instance, source } = getCategoryInstanceData(category, context);
        const sortedFields = this.getSortedFields(instance, ['Cartn_x', 'Cartn_y', 'Cartn_z', 'type_symbol']);
        const label_comp_id = this.getField(instance, 'label_comp_id');
        const label_atom_id = this.getField(instance, 'label_atom_id');

        // all of this is used to ensure that only 1 residue is written
        // TODO potentially we could use a dedicated query mode or use the 'correct' SDF definition and allow an arbitrary number of molecules in a file
        const auxiliaryFields = this.getSortedFields(instance, ['label_seq_id', 'label_asym_id', 'pdbx_PDB_ins_code', 'pdbx_PDB_model_num']);

        // write header
        const name = label_comp_id.value(source[0].keys().move(), source[0].data, 0) as string;
        StringBuilder.write(this.builder, `# Name: ${name}\n# Created by ${this.encoder}\n\n`);

        const bondMap = this.componentData.entries.get(name)!;
        let bondCount = 0;

        // traverse once to determine all actually present atoms
        const atoms = this.getAtoms(source, sortedFields, label_atom_id, auxiliaryFields, ctab);
        for (let i1 = 0, il = atoms.length; i1 < il; i1++) {
            bondMap.map.get(atoms[i1])!.forEach((v, k) => {
                const i2 = atoms.indexOf(k);
                const label2 = this.getLabel(k);
                if (i1 < i2 && atoms.indexOf(k) > -1 && !this.skipHydrogen(label2)) {
                    const { order } = v;
                    StringBuilder.writeIntegerPadLeft(bonds, i1 + 1, 3);
                    StringBuilder.writeIntegerPadLeft(bonds, i2 + 1, 3);
                    StringBuilder.writeIntegerPadLeft(bonds, order, 3);
                    StringBuilder.writeSafe(bonds, '  0  0  0  0\n');
                    // TODO 2nd value: Single bonds: 0 = not stereo, 1 = Up, 4 = Either, 6 = Down, Double bonds: 0 = Use x-, y-, z-coords from atom block to determine cis or trans, 3 = Cis or trans (either) double bond
                    bondCount++;
                }
            });
        }

        // write counts line
        StringBuilder.writeIntegerPadLeft(this.builder, atoms.length, 3);
        StringBuilder.writeIntegerPadLeft(this.builder, bondCount, 3);
        StringBuilder.write(this.builder, '  0     0  0  0  0  0  0999 V2000\n');
        // TODO 2nd value: chiral flag: 0=not chiral, 1=chiral

        StringBuilder.writeSafe(this.builder, StringBuilder.getString(ctab));
        StringBuilder.writeSafe(this.builder, StringBuilder.getString(bonds));
        StringBuilder.writeSafe(this.builder, StringBuilder.getString(charges));
        // TODO charges?

        StringBuilder.writeSafe(this.builder, 'M  END\n');
    }

    private getAtoms(source: any, fields: Field<any, any>[], label_atom_id: Field<any, any>, auxiliaryFields: Field<any, any>[], ctab: StringBuilder): string[] {
        const atoms = [];
        let index = 0;
        let id;

        // is this loop needed?
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
                atoms.push(lai);

                for (let _f = 0, _fl = fields.length; _f < _fl; _f++) {
                    const f: Field<any, any> = fields[_f]!;
                    const v = f.value(key, data, index);
                    this.writeValue(ctab, v, f.type);
                }

                StringBuilder.writeSafe(ctab, '  0  0  0  0  0  0  0  0  0  0  0  0\n');
                index++;
            }
        }

        return atoms;
    }

    private skipHydrogen(label: string) {
        if (this.hydrogens) {
            return false;
        }
        return label.startsWith('H');
    }

    private getLabel(s: string) {
        return s.replace(/[^A-Z]+/g, '');
    }

    private writeFullCategory<Ctx>(sb: StringBuilder, category: Category<Ctx>, context?: Ctx) {
        const { instance, source } = getCategoryInstanceData(category, context);
        const fields = instance.fields;
        const src = source[0];
        const data = src.data;

        const it = src.keys();
        const key = it.move();
        for (let _f = 0; _f < fields.length; _f++) {
            const f = fields[_f]!;

            StringBuilder.writeSafe(sb, `> <${category.name}.${f.name}>\n`);
            const val = f.value(key, data, 0);
            StringBuilder.writeSafe(sb, val as string);
            StringBuilder.writeSafe(sb, '\n\n');
        }
    }

    private writeValue(sb: StringBuilder, val: string | number, t: Field.Type, floatPrecision: number = 4) {
        if (t === Field.Type.Str) {
            // type_symbol is the only string field - width 2, right-padded
            StringBuilder.whitespace1(sb);
            StringBuilder.writePadRight(sb, val as string, 2);
        } else if (t === Field.Type.Int) {
            StringBuilder.writeInteger(sb, val as number);
        } else {
            // coordinates have width 10 and are left-padded
            StringBuilder.writePadLeft(sb, (val as number).toFixed(floatPrecision), 10);
        }
    }

    private getSortedFields<Ctx>(instance: Category.Instance<Ctx>, names: string[]) {
        return names.map(n => this.getField(instance, n));
    }

    private getField<Ctx>(instance: Category.Instance<Ctx>, name: string) {
        return instance.fields.find(f => f.name === name)!;
    }

    encode() {
        // write meta-information, do so after ctab
        if (this.error || this.metaInformation) {
            StringBuilder.writeSafe(this.builder, StringBuilder.getString(this.meta));
        }

        // terminate file
        StringBuilder.writeSafe(this.builder, '$$$$\n');

        this.encoded = true;
    }

    setFilter(filter?: Category.Filter) {}

    setFormatter(formatter?: Category.Formatter) {}

    isCategoryIncluded(name: string) {
        return true;
    }

    constructor(readonly encoder: string, readonly metaInformation: boolean, readonly hydrogens: boolean) {
        this.builder = StringBuilder.create();
        this.meta = StringBuilder.create();
    }
}