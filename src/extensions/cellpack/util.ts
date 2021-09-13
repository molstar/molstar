/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { CIF } from '../../mol-io/reader/cif';
import { parsePDB } from '../../mol-io/reader/pdb/parser';
import { AssetManager, Asset } from '../../mol-util/assets';
import { Structure } from '../../mol-model/structure';
import { Vec3 } from '../../mol-math/linear-algebra';
import { PluginContext } from '../../mol-plugin/context';

export async function parseCif(plugin: PluginContext, data: string | Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await plugin.runTask(comp);
    if (parsed.isError) throw parsed;
    return parsed.result;
}

export async function parsePDBfile(plugin: PluginContext, data: string, id: string) {
    const comp = parsePDB(data, id);
    const parsed = await plugin.runTask(comp);
    if (parsed.isError) throw parsed;
    return parsed.result;
}

async function downloadCif(plugin: PluginContext, url: string, isBinary: boolean, assetManager: AssetManager) {
    const type = isBinary ? 'binary' : 'string';
    const asset = await plugin.runTask(assetManager.resolve(Asset.getUrlAsset(assetManager, url), type));
    return { cif: await parseCif(plugin, asset.data), asset };
}

async function downloadPDB(plugin: PluginContext, url: string, id: string, assetManager: AssetManager) {
    const asset = await assetManager.resolve(Asset.getUrlAsset(assetManager, url), 'string').run();
    return { pdb: await parsePDBfile(plugin, asset.data, id), asset };
}

export async function getFromPdb(plugin: PluginContext, pdbId: string, assetManager: AssetManager) {
    // ${pdbId.toUpperCase()}
    const { cif, asset } = await downloadCif(plugin, `https://models.rcsb.org/${pdbId}.bcif`, true, assetManager);
    return { mmcif: cif.blocks[0], asset };
}

export async function getFromOPM(plugin: PluginContext, pdbId: string, assetManager: AssetManager) {
    const asset = await plugin.runTask(assetManager.resolve(Asset.getUrlAsset(assetManager, `https://opm-assets.storage.googleapis.com/pdb/${pdbId.toLowerCase()}.pdb`), 'string'));
    return { pdb: await parsePDBfile(plugin, asset.data, pdbId), asset };
}

export async function getFromCellPackDB(plugin: PluginContext, id: string, baseUrl: string, assetManager: AssetManager) {
    if (id.toLowerCase().endsWith('.cif') || id.toLowerCase().endsWith('.bcif')) {
        const isBinary = id.toLowerCase().endsWith('.bcif');
        const { cif, asset } = await downloadCif(plugin, `${baseUrl}/other/${id}`, isBinary, assetManager);
        return { mmcif: cif.blocks[0], asset };
    } else {
        const name = id.endsWith('.pdb') ? id.substring(0, id.length - 4) : id;
        return await downloadPDB(plugin, `${baseUrl}/other/${name}.pdb`, name, assetManager);
    }
}

export type IngredientFiles = { [name: string]: Asset.File }

export function getStructureMean(structure: Structure) {
    let xSum = 0, ySum = 0, zSum = 0;
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i];
        const { elements } = unit;
        const { x, y, z } = unit.conformation;
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const eI = elements[j];
            xSum += x(eI);
            ySum += y(eI);
            zSum += z(eI);
        }
    }
    const { elementCount } = structure;
    return Vec3.create(xSum / elementCount, ySum / elementCount, zSum / elementCount);
}

export function getFloatValue(value: DataView, offset: number) {
    // if the last byte is a negative value (MSB is 1), the final
    // float should be too
    const negative = value.getInt8(offset + 2) >>> 31;

    // this is how the bytes are arranged in the byte array/DataView
    // buffer
    const [b0, b1, b2, exponent] = [
        // get first three bytes as unsigned since we only care
        // about the last 8 bits of 32-bit js number returned by
        // getUint8().
        // Should be the same as: getInt8(offset) & -1 >>> 24
        value.getUint8(offset),
        value.getUint8(offset + 1),
        value.getUint8(offset + 2),

        // get the last byte, which is the exponent, as a signed int
        // since it's already correct
        value.getInt8(offset + 3)
    ];

    let mantissa = b0 | (b1 << 8) | (b2 << 16);
    if (negative) {
        // need to set the most significant 8 bits to 1's since a js
        // number is 32 bits but our mantissa is only 24.
        mantissa |= 255 << 24;
    }

    return mantissa * Math.pow(10, exponent);
}

/*
async function loadPackingResultsBinary(plugin: PluginContext, runtime: RuntimeContext, file: Asset.File, packing: CellPack){
    const model_data = await readFromFile(file.file!, 'binary').runInContext(runtime);// async ?
    let buffer = model_data.buffer;
    let {cell, packings } = packing;
    if (!IsNativeEndianLittle) {
        // flip the byte order
        buffer = flipByteOrder(model_data, 4);
    }
    const numbers = new DataView(buffer);
    const ninst = getFloatValue(numbers, 0);
    const npoints = getFloatValue(numbers, 4);
    const ncurve = getFloatValue(numbers, 8);

    let pos = new Float32Array();
    let quat = new Float32Array();
    let ctr_pos = new Float32Array();
    let ctr_info = new Float32Array();
    let curve_ids = new Float32Array();

    let offset = 12;
    if (ninst !== 0){
        pos = new Float32Array(buffer, offset, ninst * 4);offset += ninst * 4 * 4;
        quat = new Float32Array(buffer, offset, ninst * 4);offset += ninst * 4 * 4;
    }
    if ( npoints !== 0 ) {
        ctr_pos = new Float32Array(buffer, offset, npoints * 4);offset += npoints * 4 * 4;
        offset += npoints * 4 * 4;
        ctr_info = new Float32Array(buffer, offset, npoints * 4);offset += npoints * 4 * 4;
        curve_ids = new Float32Array(buffer, offset, ncurve * 4);offset += ncurve * 4 * 4;
    }

    for (let i = 0; i < ninst; i++) {
        const x: number =  pos[i * 4 + 0];
        const y: number =  pos[i * 4 + 1];
        const z: number =  pos[i * 4 + 2];
        const ingr_id = pos[i * 4 + 3] as number;
        const pid = cell.mapping_ids![ingr_id];
        if (!packings[pid[0]].ingredients[pid[1]].results) {
            packings[pid[0]].ingredients[pid[1]].results = [];
        }
        packings[pid[0]].ingredients[pid[1]].results.push([Vec3.create(x, y, z),
            Quat.create(quat[i * 4 + 0], quat[i * 4 + 1], quat[i * 4 + 2], quat[i * 4 + 3])]);
    }
    let counter = 0;
    let ctr_points: Vec3[] = [];
    let prev_ctype = 0;
    let prev_cid = 0;

    for (let i = 0; i < npoints; i++) {
        const x: number = -ctr_pos[i * 4 + 0];
        const y: number =  ctr_pos[i * 4 + 1];
        const z: number =  ctr_pos[i * 4 + 2];
        const cid: number = ctr_info[i * 4 + 0];// curve id
        const ctype: number = curve_ids[cid * 4 + 0];// curve type
        // cid  148 165 -1 0
        // console.log("cid ",cid,ctype,prev_cid,prev_ctype);//165,148
        if (prev_ctype !== ctype){
            const pid = cell.mapping_ids![-prev_ctype - 1];
            const cname = `curve${counter}`;
            packings[pid[0]].ingredients[pid[1]].nbCurve = counter + 1;
            packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
            ctr_points = [];
            counter = 0;
        } else if (prev_cid !== cid){
            ctr_points = [];
            const pid = cell.mapping_ids![-prev_ctype - 1];
            const cname = `curve${counter}`;
            packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
            counter += 1;
        }
        ctr_points.push(Vec3.create(x, y, z));
        prev_ctype = ctype;
        prev_cid = cid;
    }
    // do the last one
    if ( npoints !== 0 ) {
        const pid = cell.mapping_ids![-prev_ctype - 1];
        const cname = `curve${counter}`;
        packings[pid[0]].ingredients[pid[1]].nbCurve = counter + 1;
        packings[pid[0]].ingredients[pid[1]][cname] = ctr_points;
    }
    return packings;
}
*/