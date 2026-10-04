import { zipSync } from 'fflate';
import type { MVSData } from '@molstar/mvs-builder/mvs-data';

/** Creates an MVSX archive with the JSON document and referenced assets. */
export async function createMVSX(data: MVSData, assets: { name: string, content: string | Uint8Array<ArrayBuffer> }[]): Promise<Uint8Array<ArrayBuffer>> {
    const encoder = new TextEncoder();
    const files: Record<string, Uint8Array> = {
        'index.mvsj': encoder.encode(JSON.stringify(data)),
    };
    for (const asset of assets) {
        files[asset.name] = typeof asset.content === 'string' ? encoder.encode(asset.content) : asset.content;
    }
    // `zipSync` is synchronous, but preserve the existing async MVSData API.
    // ZIP stores local date fields; use the same fields in every time zone.
    return zipSync(files, { mtime: new Date(1980, 0, 1) });
}
