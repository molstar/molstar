/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { StringBuilder } from '../../../mol-util';
import Writer from '../writer';
import { Encoder, Category, Field } from '../cif/encoder';
import { getCategoryInstanceData } from '../cif/encoder/util';

// specification: http://c4.cabrillo.edu/404/ctfile.pdf
export class SdfEncoder implements Encoder<string> {
    private builder: StringBuilder;
    private meta: StringBuilder;
    private encoded = false;
    private dataBlockCreated = false;
    private error = false;

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

    startDataBlock(name: string) {
        this.dataBlockCreated = true;
        StringBuilder.write(this.builder, `${name}\nCreated by ${this.encoder}\n\n`);
    }

    writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx) {
        if (this.encoded) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (!this.dataBlockCreated) {
            throw new Error('No data block created.');
        }

        if (!this.hideMetaInformation && (category.name === 'model_server_result' || category.name === 'model_server_params' || category.name === 'model_server_stats' || category.name === 'model_server_error')) {
            this.writeFullCategory(this.meta, category, context);
            // if error: force writing of meta information
            if (category.name === 'model_server_error') {
                this.error = true;
            }
            return;
        }

        // ignore meta, error, and misc categories when writing SDF
        if (category.name !== 'atom_site') {
            return;
        }

        // use separate builder because we still need to write Counts and Bonds line
        const ctab = StringBuilder.create();
        const bonds = StringBuilder.create();
        const charges = StringBuilder.create();

        // write Atom block and gather data for Bonds and Charges
        // 'Specifies the atomic symbol and any mass difference, charge, stereochemistry, and associated hydrogens for each atom.'
        const { instance, source, rowCount: atomCount } = getCategoryInstanceData(category, context);
        const fields = this.getSortedFields(instance);
        const fieldCount = fields.length;
        let bondCount = 0;

        let index = 0;
        for (let _c = 0; _c < source.length; _c++) {
            const src = source[_c];
            const data = src.data;

            if (src.rowCount === 0) continue;

            const it = src.keys();
            while (it.hasNext)  {
                const key = it.move();

                for (let _f = 0; _f < fieldCount; _f++) {
                    const f: Field<any, any> = fields[_f]!;
                    const val = f.value(key, data, index);
                    this.writeValue(ctab, val, f.type);
                }
                
                StringBuilder.writeSafe(ctab, '  0  0  0  0  0\n');

                index++;
            }
        }

        // write counts line
        // 'Important specifications here relate to the number of atoms, bonds, and atom lists, the chiral flag setting, and the Ctab version.'
        StringBuilder.writeIntegerPadLeft(this.builder, atomCount, 3);
        StringBuilder.writeIntegerPadLeft(this.builder, bondCount, 3);
        StringBuilder.write(this.builder, '  0  0  0  0  0  0  0  0999 V2000\n');

        StringBuilder.writeSafe(this.builder, StringBuilder.getString(ctab));
        StringBuilder.writeSafe(this.builder, StringBuilder.getString(bonds));
        StringBuilder.writeSafe(this.builder, StringBuilder.getString(charges));
        
        StringBuilder.writeSafe(this.builder, 'M  END\n');
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

    private getSortedFields<Ctx>(instance: Category.Instance<Ctx>) {
        return ['Cartn_x', 'Cartn_y', 'Cartn_z', 'type_symbol']
            .map(n => instance.fields.find(f => f.name === n));
    }

    readonly isBinary = false;
    binaryEncodingProvider = void 0;

    encode() {
        // write meta-information, do so after ctab
        if (this.error || !this.hideMetaInformation) {
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

    constructor(readonly encoder: string, readonly hideMetaInformation: boolean) {
        this.builder = StringBuilder.create();
        this.meta = StringBuilder.create();
    }
}