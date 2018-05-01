
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as argparse from 'argparse'
import { VolumeData, parseDensityServerData } from 'mol-model/volume'
import { downloadCif } from './helpers'
import CIF from 'mol-io/reader/cif'
import { DensityServer_Data_Database } from 'mol-io/reader/cif/schema/density-server';
import { Run } from 'mol-task';
import { Table } from 'mol-data/db';

type Volume = { source: DensityServer_Data_Database, volume: VolumeData }

async function getVolume(url: string): Promise<Volume> {
    const cif = await downloadCif(url, false);
    const data = CIF.schema.densityServer(cif.blocks[1]);
    return { source: data, volume: await Run(parseDensityServerData(data)) };
}

function print(volume: Volume) {
    const { volume_data_3d_info } = volume.source;
    const row = Table.getRow(volume_data_3d_info, 0);
    console.log(row);
}

async function run(url: string) {
    const volume = await getVolume(url);
    print(volume);
}

const parser = new argparse.ArgumentParser({
addHelp: true,
description: 'Info about VolumeData from mol-model module'
});
parser.addArgument([ '--emdb', '-e' ], {
    help: 'EMDB id, for example 8116',
});
interface Args {
    emdb?: string
}
const args: Args = parser.parseArgs();

run(`https://webchem.ncbr.muni.cz/DensityServer/em/emd-${args.emdb}/cell?detail=1&encoding=cif`);