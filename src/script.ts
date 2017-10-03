/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import * as util from 'util'
import * as fs from 'fs'

import Gro from './reader/gro/parser'
import CIF from './reader/cif/index'

// const file = '1crn.gro'
// const file = 'water.gro'
// const file = 'test.gro'
const file = 'md_1u19_trj.gro'

export function _gro() {
    fs.readFile(`./examples/${file}`, 'utf8', function (err, input) {
        if (err) {
            return console.log(err);
        }
        // console.log(data);

        console.time('parse')
        const parsed = Gro(input)
        console.timeEnd('parse')
        if (parsed.isError) {
            console.log(parsed)
            return;
        }

        const groFile = parsed.result

        console.log('structure count: ', groFile.structures.length);

        const data = groFile.structures[0];

        // const header = groFile.blocks[0].getCategory('header')
        const { header, atoms } = data;
        console.log(JSON.stringify(header, null, 2));
        console.log('number of atoms:', atoms.count);

        console.log(`'${atoms.residueNumber.value(1)}'`)
        console.log(`'${atoms.residueName.value(1)}'`)
        console.log(`'${atoms.atomName.value(1)}'`)
        console.log(atoms.z.value(1))
        console.log(`'${atoms.z.value(1)}'`)

        const n = atoms.count;
        console.log('rowCount', n)

        console.time('getFloatArray x')
        const x = atoms.x.toArray({ array: Float32Array })
        console.timeEnd('getFloatArray x')
        console.log(x.length, x[0], x[x.length - 1])

        console.time('getFloatArray y')
        const y = atoms.y.toArray({ array: Float32Array })
        console.timeEnd('getFloatArray y')
        console.log(y.length, y[0], y[y.length - 1])

        console.time('getFloatArray z')
        const z = atoms.z.toArray({ array: Float32Array })
        console.timeEnd('getFloatArray z')
        console.log(z.length, z[0], z[z.length - 1])

        console.time('getIntArray residueNumber')
        const residueNumber = atoms.residueNumber.toArray({ array: Int32Array })
        console.timeEnd('getIntArray residueNumber')
        console.log(residueNumber.length, residueNumber[0], residueNumber[residueNumber.length - 1])
    });
}

export function _cif() {
    const path = `./examples/1cbs_updated.cif`;
    //const path = 'c:/test/quick/3j3q.cif';
    fs.readFile(path, 'utf8', function (err, input) {
        if (err) {
            return console.log(err);
        }

        console.time('parseCIF');
        const parsed = CIF.parseText(input);
        console.timeEnd('parseCIF');
        if (parsed.isError) {
            console.log(parsed);
            return;
        }

        const data = parsed.result.blocks[0];

        const atom_site = data.categories._atom_site;
        console.log(atom_site.getField('Cartn_x')!.float(0));
        //console.log(atom_site.getField('label_atom_id')!.toStringArray());

        const mmcif = CIF.applySchema(CIF.schema.mmCIF, data);
        console.log(mmcif.atom_site.Cartn_x.value(0));
        console.log(mmcif.entity.type.toArray());
        console.log(mmcif.pdbx_struct_oper_list.matrix.value(0));
    });
}

_cif();