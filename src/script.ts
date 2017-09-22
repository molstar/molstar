/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import * as util from 'util'
import * as fs from 'fs'

import Gro from './reader/gro/format'

//const file = '1crn.gro'
// const file = 'water.gro'
// const file = 'test.gro'
const file = 'md_1u19_trj.gro'

fs.readFile(`./examples/${file}`, 'utf8', function (err,data) {
    if (err) {
        return console.log(err);
    }
    // console.log(data);

    console.time('parse')
    const parsed = Gro.parse(data)
    console.timeEnd('parse')
    if (parsed.isError) {
        console.log(parsed)
    } else {
        const groFile = parsed.result
        const data = Gro.schema(groFile.blocks[0])

        // const header = groFile.blocks[0].getCategory('header')
        const { header, atoms } = data;
        if (header._rowCount !== 1) {
            console.log('title', header.title.value(0))
            console.log('timeInPs', header.timeInPs.value(0))
            console.log('numberOfAtoms', header.numberOfAtoms.value(0))
            console.log('boxX', header.boxX.value(0))
            console.log('boxY', header.boxY.value(0))
            console.log('boxZ', header.boxZ.value(0))
        } else {
            console.error('no header')
        }

        if (atoms._rowCount > 0) {
            console.log(`'${atoms.residueNumber.value(1)}'`)
            console.log(`'${atoms.residueName.value(1)}'`)
            console.log(`'${atoms.atomName.value(1)}'`)
            console.log(atoms.z.value(1))
            console.log(`'${atoms.z.value(1)}'`)

            const n = atoms._rowCount
            console.log('rowCount', n)

            console.time('getFloatArray x')
            const x = atoms.x.toArray(0, n, x => new Float32Array(x))!
            console.timeEnd('getFloatArray x')
            console.log(x.length, x[0], x[x.length-1])

            console.time('getFloatArray y')
            const y = atoms.y.toArray(0, n, x => new Float32Array(x))!
            console.timeEnd('getFloatArray y')
            console.log(y.length, y[0], y[y.length-1])

            console.time('getFloatArray z')
            const z = atoms.z.toArray(0, n, x => new Float32Array(x))!
            console.timeEnd('getFloatArray z')
            console.log(z.length, z[0], z[z.length-1])

            console.time('getIntArray residueNumber')
            const residueNumber = atoms.residueNumber.toArray(0, n, x => new Int32Array(x))!
            console.timeEnd('getIntArray residueNumber')
            console.log(residueNumber.length, residueNumber[0], residueNumber[residueNumber.length-1])
        } else {
            console.error('no atoms')
        }
    }
});
