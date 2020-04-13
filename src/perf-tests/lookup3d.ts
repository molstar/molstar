import * as util from 'util';
import * as fs from 'fs';
import { CIF } from '../mol-io/reader/cif';

import { Structure } from '../mol-model/structure';

import { GridLookup3D } from '../mol-math/geometry';
// import { sortArray } from 'mol-data/util';
import { OrderedSet } from '../mol-data/int';
import { trajectoryFromMmCIF, MmcifFormat } from '../mol-model-formats/structure/mmcif';
import { getBoundary } from '../mol-math/geometry/boundary';

require('util.promisify').shim();
const readFileAsync = util.promisify(fs.readFile);

async function readData(path: string) {
    if (path.match(/\.bcif$/)) {
        const input = await readFileAsync(path);
        const data = new Uint8Array(input.byteLength);
        for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
        return data;
    } else {
        return readFileAsync(path, 'utf8');
    }
}


export async function readCIF(path: string) {
    const input = await readData(path);
    const comp = typeof input === 'string' ? CIF.parseText(input) : CIF.parseBinary(input);
    const parsed = await comp.run();
    if (parsed.isError) {
        throw parsed;
    }

    const models = await trajectoryFromMmCIF(parsed.result.blocks[0]).run();
    const structures = models.map(Structure.ofModel);

    return { mmcif: models[0].sourceData.data as MmcifFormat.Data, models, structures };
}

export async function test() {
    const { mmcif, structures } = await readCIF('e:/test/quick/1tqn_updated.cif');

    const position = { x: mmcif.db.atom_site.Cartn_x.toArray(), y: mmcif.db.atom_site.Cartn_y.toArray(), z: mmcif.db.atom_site.Cartn_z.toArray(),
        indices: OrderedSet.ofBounds(0, mmcif.db.atom_site._rowCount),
        // radius: [1, 1, 1, 1]
        // indices: [1]
    };
    const lookup = GridLookup3D(position, getBoundary(position));
    console.log(lookup.boundary.box, lookup.boundary.sphere);

    const result = lookup.find(-30.07, 8.178, -13.897, 10);
    console.log(result.count); // , sortArray(result.indices));

    // const sl = structures[0].lookup3d;
    // const result1 = sl.find(-30.07, 8.178, -13.897, 10);
    // console.log(result1.count);//, result1.indices);

    console.log(structures[0].boundary);
    console.log(lookup.boundary);
}

test();