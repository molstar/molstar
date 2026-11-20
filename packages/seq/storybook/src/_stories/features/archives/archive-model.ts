/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { zip } from '@molstar/core/util/zip/zip';
import {
  type Archive,
  ArchiveIndexName,
  type ArchiveRead,
  ArchiveSpecName,
  type AssetReader,
  readArchive,
  writeArchive,
} from '@molstar/seqx-core/data';

// What the Archives stories show, DOM-free (the story and its headless twin print the same text): a view saved with
// its data as one archive and opened again, and archives the reader refuses (EC-91).

/** A view saved and opened: the archive's bytes, what it holds and the printed listing. */
export interface SavedView {
  readonly bytes: Uint8Array;
  readonly archive: Archive;
  readonly listing: string;
}

/** The archive as a listing: its size, its spec and index, and each file with its size. */
export function listArchive(bytes: Uint8Array, archive: Archive): string {
  const lines = [`${bytes.byteLength} bytes as zip`, `  ${ArchiveIndexName}`, `  ${ArchiveSpecName}`];
  for (const [name, data] of archive.files) lines.push(`  ${name}  ${data.byteLength} bytes`);
  return lines.join('\n');
}

/** Saves a view with the files its sources name (read through `assets`), then opens the bytes again. */
export async function saveView(
  spec: Node,
  assets: AssetReader,
  options: { readonly keep?: (uri: string) => boolean } = {},
): Promise<SavedView> {
  const written = await writeArchive(spec, { assets, ...options });
  if (!written.bytes) throw new Error(written.issues.map((i) => i.message).join('; '));
  const read = await readArchive(written.bytes);
  if (!read.archive) throw new Error(read.issues.map((i) => i.message).join('; '));
  return { bytes: written.bytes, archive: read.archive, listing: listArchive(written.bytes, read.archive) };
}

const enc = new TextEncoder();

/** An archive that is refused: what it holds that no reader may open (a path above the root, a bomb), and why. */
export interface Hostile {
  readonly id: string;
  readonly title: string;
  readonly bytes: () => Promise<Uint8Array>;
}

const index = (files: string[]) => JSON.stringify({ archive: 'seqx', version: 1, spec: ArchiveSpecName, files });
const spec = JSON.stringify({ kind: 'spec', params: { version: '0.1.0' } });

async function archiveOf(entries: Record<string, Uint8Array<ArrayBuffer>>, files: string[] = []) {
  return new Uint8Array(
    await zip(SyncRuntimeContext, {
      [ArchiveIndexName]: enc.encode(index(files)),
      [ArchiveSpecName]: enc.encode(spec),
      ...entries,
    }),
  );
}

export const HostileArchives: readonly Hostile[] = [
  {
    id: 'above-root',
    title: 'An entry above the root',
    bytes: () => archiveOf({ '../../outside.seq': enc.encode('never written anywhere') }),
  },
  {
    id: 'absolute',
    title: 'An absolute path',
    bytes: () => archiveOf({ '/etc/passwd': enc.encode('never read') }),
  },
  {
    id: 'bomb',
    title: 'A megabyte of zeros in a few hundred bytes',
    bytes: () => archiveOf({ 'zeros.seq': new Uint8Array(1 << 20) }, ['zeros.seq']),
  },
  {
    id: 'many',
    title: 'More entries than the bound',
    bytes: () => {
      const entries: Record<string, Uint8Array<ArrayBuffer>> = {};
      for (let i = 0; i < 1001; i++) entries[`e${i}`] = enc.encode('x');
      return archiveOf(entries);
    },
  },
];

/** What the reader says of a refused archive, printed. */
export async function printRefusal(h: Hostile): Promise<{ read: ArchiveRead; text: string }> {
  const bytes = await h.bytes();
  const read = await readArchive(bytes);
  const lines = [`${h.title}: ${bytes.byteLength} bytes`, read.archive ? 'opened' : 'refused'];
  for (const i of read.issues) lines.push(`${i.severity} ${i.code}: ${i.message}`);
  return { read, text: lines.join('\n') };
}

/** The assets of the page: a location is fetched relative to it (the Vite plugin serves the examples' fixtures). */
export const pageAssets: AssetReader = {
  readText: async (uri) => {
    const r = await fetch(uri);
    if (!r.ok) throw new Error(`${uri}: ${r.status}`);
    return r.text();
  },
  readBytes: async (uri) => {
    const r = await fetch(uri);
    if (!r.ok) throw new Error(`${uri}: ${r.status}`);
    return new Uint8Array(await r.arrayBuffer());
  },
};
