import { zipSync, type ZipOptions } from 'fflate';
import type { MVSData } from '@molstar/mvs-builder/mvs-data';

/** Creates an MVSX archive with the JSON document and referenced assets.
 * Pass `options.zip.mtime` for reproducible output; timestamps default to the current time. */
export async function createMVSX(
  data: MVSData,
  assets: { name: string; content: string | Uint8Array<ArrayBuffer> }[],
  options?: { zip?: ZipOptions },
): Promise<Uint8Array<ArrayBuffer>> {
  const encoder = new TextEncoder();
  const files: Record<string, Uint8Array> = {
    'index.mvsj': encoder.encode(JSON.stringify(data)),
  };
  for (const asset of assets) {
    files[asset.name] = typeof asset.content === 'string' ? encoder.encode(asset.content) : asset.content;
  }
  // `zipSync` is synchronous, but preserve the existing async MVSData API.
  return zipSync(files, options?.zip);
}
