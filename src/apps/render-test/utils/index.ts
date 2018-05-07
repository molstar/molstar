/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF from 'mol-io/reader/cif'
import { Run, Progress } from 'mol-task'
import { Model } from 'mol-model/structure'
import { VolumeData, parseDensityServerData } from 'mol-model/volume'
import { DensityServer_Data_Database } from 'mol-io/reader/cif/schema/density-server';

export function log(progress: Progress) {
    const p = progress.root.progress
    console.log(`${p.message} ${(p.current/p.max*100).toFixed(2)}%`)
}

export async function downloadCif(url: string, isBinary: boolean) {
    const data = await fetch(url);
    return parseCif(isBinary ? new Uint8Array(await data.arrayBuffer()) : await data.text());
}

export async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data)
    const parsed = await Run(comp, log, 500);
    if (parsed.isError) throw parsed;
    return parsed.result
}

export async function getModelFromPdbId(pdbid: string) {
    const cif = await downloadCif(`https://files.rcsb.org/download/${pdbid}.cif`, false)
    return Model.create({ kind: 'mmCIF', data: CIF.schema.mmCIF(cif.blocks[0]) })
}

const readFileAsText = (file: File) => {
    const fileReader = new FileReader()
    return new Promise<string>((resolve, reject) => {
        fileReader.onerror = () => {
            fileReader.abort()
            reject(new DOMException('Error parsing input file.'))
        }
        fileReader.onload = () => resolve(fileReader.result)
        fileReader.readAsText(file)
    })
}

export async function getModelFromFile(file: File) {
    const cif = await parseCif(await readFileAsText(file))
    return Model.create({ kind: 'mmCIF', data: CIF.schema.mmCIF(cif.blocks[0]) })
}

export type Volume = { source: DensityServer_Data_Database, volume: VolumeData }

export async function getVolumeFromEmdId(emdid: string): Promise<Volume> {
    const cif = await downloadCif(`https://webchem.ncbr.muni.cz/DensityServer/em/emd-${emdid}/cell?detail=4`, true)
    const data = CIF.schema.densityServer(cif.blocks[1])
    return { source: data, volume: await Run(parseDensityServerData(data)) }
}