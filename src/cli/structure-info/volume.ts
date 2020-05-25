#!/usr/bin/env node
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import * as argparse from 'argparse';
import * as util from 'util';

import { Volume } from '../../mol-model/volume';
import { downloadCif } from './helpers';
import { CIF } from '../../mol-io/reader/cif';
import { Table } from '../../mol-data/db';
import { StringBuilder } from '../../mol-util';
import { Task } from '../../mol-task';
import { createVolumeIsosurfaceMesh } from '../../mol-repr/volume/isosurface';
import { Theme } from '../../mol-theme/theme';
import { volumeFromDensityServerData, DscifFormat } from '../../mol-model-formats/volume/density-server';

require('util.promisify').shim();
const writeFileAsync = util.promisify(fs.writeFile);

async function getVolume(url: string): Promise<Volume> {
    const cif = await downloadCif(url, true);
    const data = CIF.schema.densityServer(cif.blocks[1]);
    return await volumeFromDensityServerData(data).run();
}

function print(volume: Volume) {
    if (!DscifFormat.is(volume.sourceData)) return;
    const { volume_data_3d_info } = volume.sourceData.data;
    const row = Table.getRow(volume_data_3d_info, 0);
    console.log(row);
    console.log(volume.grid.transform);
    console.log(volume.grid.stats);
}

async function doMesh(volume: Volume, filename: string) {
    const mesh = await Task.create('', runtime => createVolumeIsosurfaceMesh({ runtime }, volume, Theme.createEmpty(), { isoValue: Volume.IsoValue.absolute(1.5) } )).run();
    console.log({ vc: mesh.vertexCount, tc: mesh.triangleCount });

    // Export the mesh in OBJ format.
    const { vertexCount, triangleCount } = mesh;

    const vs = mesh.vertexBuffer.ref.value;
    const ts = mesh.indexBuffer.ref.value;

    const obj = StringBuilder.create();
    for (let i = 0; i < vertexCount; i++) {
        StringBuilder.write(obj, 'v ');
        StringBuilder.writeFloat(obj, vs[3 * i + 0], 100);
        StringBuilder.whitespace1(obj);
        StringBuilder.writeFloat(obj, vs[3 * i + 1], 100);
        StringBuilder.whitespace1(obj);
        StringBuilder.writeFloat(obj, vs[3 * i + 2], 100);
        StringBuilder.newline(obj);
    }
    for (let i = 0; i < triangleCount; i++) {
        StringBuilder.write(obj, 'f ');
        StringBuilder.writeIntegerAndSpace(obj, ts[3 * i + 0] + 1);
        StringBuilder.writeIntegerAndSpace(obj, ts[3 * i + 1] + 1);
        StringBuilder.writeInteger(obj, ts[3 * i + 2] + 1);
        StringBuilder.newline(obj);
    }

    await writeFileAsync(filename, StringBuilder.getString(obj));
}

async function run(url: string, meshFilename: string) {
    const volume = await getVolume(url);
    print(volume);
    await doMesh(volume, meshFilename);
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Info about VolumeData from mol-model module'
});
parser.addArgument([ '--emdb', '-e' ], {
    help: 'EMDB id, for example 8116',
});
parser.addArgument([ '--mesh' ], {
    help: 'Mesh filename',
    required: true
});
interface Args {
    emdb?: string,
    mesh: string
}
const args: Args = parser.parseArgs();

run(`https://ds.litemol.org/em/emd-${args.emdb}/cell?detail=4`, args.mesh);