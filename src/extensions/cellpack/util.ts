/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CIF } from '../../mol-io/reader/cif';
import { parsePDB } from '../../mol-io/reader/pdb/parser';
import { AssetManager, Asset } from '../../mol-util/assets';

export async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

export async function parsePDBfile(data: string, id: string) {
    const comp = parsePDB(data, id);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

async function downloadCif(url: string, isBinary: boolean, assetManager: AssetManager) {
    const type = isBinary ? 'binary' : 'string';
    const asset = await assetManager.resolve(Asset.getUrlAsset(assetManager, url), type).run();
    return { cif: await parseCif(asset.data), asset };
}

async function downloadPDB(url: string, id: string, assetManager: AssetManager) {
    const asset = await assetManager.resolve(Asset.getUrlAsset(assetManager, url), 'string').run();
    return { pdb: await parsePDBfile(asset.data, id), asset };
}

export async function getFromPdb(pdbId: string, assetManager: AssetManager) {
    const { cif, asset } = await downloadCif(`https://models.rcsb.org/${pdbId.toUpperCase()}.bcif`, true, assetManager);
    return { mmcif: cif.blocks[0], asset };
}

export async function getFromCellPackDB(id: string, baseUrl: string, assetManager: AssetManager) {
    if (id.toLowerCase().endsWith('.cif') || id.toLowerCase().endsWith('.bcif')) {
        const isBinary = id.toLowerCase().endsWith('.bcif');
        const { cif, asset } = await downloadCif(`${baseUrl}/other/${id}`, isBinary, assetManager);
        return { mmcif: cif.blocks[0], asset };
    } else {
        const name = id.endsWith('.pdb') ? id.substring(0, id.length - 4) : id;
        return await downloadPDB(`${baseUrl}/other/${name}.pdb`, name, assetManager);
    }
}

export type IngredientFiles = { [name: string]: Asset.File }