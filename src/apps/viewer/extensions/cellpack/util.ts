/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CIF } from '../../../../mol-io/reader/cif'
import { parsePDB } from '../../../../mol-io/reader/pdb/parser';

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

async function parsePDBfile(data: string) {
    const comp = parsePDB(data);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

async function downloadCif(url: string, isBinary: boolean) {
    const data = await fetch(url);
    return parseCif(isBinary ? new Uint8Array(await data.arrayBuffer()) : await data.text());
}

async function downloadPDB(url: string) {
    const data = await fetch(url);
    return parsePDBfile(await data.text());
}

export async function getFromPdb(id: string) {
    const parsed = await downloadCif(`https://files.rcsb.org/download/${id}.cif`, false);
    return parsed.blocks[0];
}

function getCellPackDataUrl(id: string, baseUrl: string) {
    const url = `${baseUrl}/other/${id}`
    return url.endsWith('.pdb') ? url : `${url}.pdb`
}

export async function getFromCellPackDB(id: string, baseUrl: string) {
    const parsed = await downloadPDB(getCellPackDataUrl(id, baseUrl));
    return parsed;
}