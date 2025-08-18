/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Zip } from '../../mol-util/zip/zip';
import { MVSData } from './mvs-data';

/**
 * Creates an MVSX zip file with from the provided data and assets
 */
export async function createMVSX(data: MVSData, assets: { name: string, content: string | Uint8Array }[]) {
    const encoder = new TextEncoder();
    const files: Record<string, Uint8Array> = {
        'index.mvsj': encoder.encode(JSON.stringify(data)),
    };
    for (const asset of assets) {
        files[asset.name] = typeof asset.content === 'string'
            ? encoder.encode(asset.content)
            : asset.content;
    }

    const zip = await Zip(files).run();
    return new Uint8Array(zip);
}