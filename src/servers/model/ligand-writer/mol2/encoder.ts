/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Category } from '../../../../mol-io/writer/cif/encoder';
import { LigandEncoder } from '../ligand-encoder';
import { StringBuilder } from '../../../../mol-util';
import { getCategoryInstanceData } from '../../../../mol-io/writer/cif/encoder/util';
import { BondType } from '../../../../mol-model/structure/model/types';

// specification: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
// TODO amide (and real sp/sp2/sp3) support for bonds and SYBYL atom types: see https://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
// TODO support charges
export class Mol2Encoder extends LigandEncoder {
    private out: StringBuilder;

    _writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx): void {
        const a = StringBuilder.create();
        const b = StringBuilder.create();
        const { instance, source } = getCategoryInstanceData(category, context);

        // write header
        const name = this.getName(instance, source);
        StringBuilder.writeSafe(this.builder, `# Name: ${name}\n# Created by ${this.encoder}\n\n`);

        const bondMap = this.componentData.entries.get(name)!;
        let bondCount = 0;

        const atoms = this.getAtoms(instance, source);
        StringBuilder.writeSafe(a, '@<TRIPOS>ATOM\n');
        StringBuilder.writeSafe(b, '@<TRIPOS>BOND\n');
        for (let i1 = 0, il = atoms.length; i1 < il; i1++) {
            const atom = atoms[i1];

            let aromatic = false;
            bondMap.map.get(atom.label_atom_id)!.forEach((v, k) => {
                const i2 = atoms.findIndex(e => e.label_atom_id === k);
                const label2 = this.getLabel(k);
                if (i1 < i2 && atoms.findIndex(e => e.label_atom_id === k) > -1 && !this.skipHydrogen(label2)) {
                    const { order, flags } = v;
                    const ar = flags === BondType.Flag.Aromatic;
                    if (ar) aromatic = true;
                    StringBuilder.writeSafe(b, `${++bondCount} ${i1 + 1} ${i2 + 1} ${ar ? 'ar' : order}`);
                    StringBuilder.newline(b);
                }
            });

            const sub = aromatic ? '.ar' : '';
            StringBuilder.writeSafe(a, `${i1 + 1} ${atom.type_symbol} ${atom.Cartn_x.toFixed(3)} ${atom.Cartn_y.toFixed(3)} ${atom.Cartn_z.toFixed(3)} ${atom.type_symbol}${sub} 1 ${name} 0.000\n`);
        }

        StringBuilder.writeSafe(this.out, `@<TRIPOS>MOLECULE\n${name}\n${atoms.length} ${bondCount} 0 0 0\nSMALL\nNO_CHARGES\n\n`);
        StringBuilder.writeSafe(this.out, StringBuilder.getString(a));
        StringBuilder.writeSafe(this.out, StringBuilder.getString(b));
        StringBuilder.writeSafe(this.out, `@<TRIPOS>SUBSTRUCTURE\n${name} ${name} 1\n`);
    }

    protected writeFullCategory<Ctx>(sb: StringBuilder, category: Category<Ctx>, context?: Ctx) {
        const { instance, source } = getCategoryInstanceData(category, context);
        const fields = instance.fields;
        const src = source[0];
        const data = src.data;

        const it = src.keys();
        const key = it.move();
        for (let _f = 0; _f < fields.length; _f++) {
            const f = fields[_f]!;

            StringBuilder.writeSafe(sb, `# ${category.name}.${f.name}: `);
            const val = f.value(key, data, 0);
            StringBuilder.writeSafe(sb, val as string);
            StringBuilder.newline(sb);
        }
        StringBuilder.newline(sb);
    }

    encode(): void {
        // write meta-information, do so after ctab
        if (this.error || this.metaInformation) {
            StringBuilder.writeSafe(this.builder, StringBuilder.getString(this.meta));
        }
        StringBuilder.writeSafe(this.builder, StringBuilder.getString(this.out));

        this.encoded = true;
    }

    constructor(encoder: string, metaInformation: boolean, hydrogens: boolean) {
        super(encoder, metaInformation, hydrogens);
        this.out = StringBuilder.create();
    }
}