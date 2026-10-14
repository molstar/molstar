/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `.svsx` (seqspec:sources#archives): a zip with a document at its main file and the files it reads. Reading, the
// relative URIs inside it, and the loaders that answer them. Writing is the builder's.

import { Issue } from '@molstar/seq-core/model';
import type { Loaders, StreamSource } from './contract.js';
import type { SeqViewSpecDoc } from './doc.js';
import { KernelRegistry } from './kernel-registry.js';
import { readDocument } from './normalize.js';
import { createRegistries, register, type SvsRegistries } from './registry.js';
import { MvsRegistry } from './tables/kernel.js';
import { isRelativeUri, resolveArchivePath } from './uri-path.js';
import { validate } from './validate.js';
import { RuntimeContext } from './zip/shims.js';
import { unzip, zip } from './zip/zip.js';

/** The main files of an archive, in the order a reader looks for them (MessagePack first). */
export const SvsxMainFiles = ['index.svsm', 'index.svsj'] as const;

export interface Svsx {
  /** The document (decoded, validated with `base`); null when it could not be read. */
  doc: SeqViewSpecDoc | null;
  issues: Issue[];
  /** The main file's path in the archive. */
  main: string;
  /** The main file's directory: relative URIs resolve against it (`ValidateOptions.base`). */
  base: string;
  /** Every file of the archive by path, as stored. */
  files: ReadonlyMap<string, Uint8Array>;
}

export interface ReadSvsxOptions {
  /** The main file's path; default `index.svsm`, else `index.svsj`. */
  main?: string;
  /** Validated against these; default the kernel's and MVS's. */
  registries?: SvsRegistries;
  strict?: boolean;
  /** At most this many entries; default 4,096 (seqspec:sources#archives). */
  maxEntries?: number;
  /** The entries' declared sizes summed, at most; default 1 GiB. */
  maxBytes?: number;
}

/** The default bounds of an archive, checked from its central directory before anything is inflated. */
export const SvsxLimits = { maxEntries: 4096, maxBytes: 1 << 30 } as const;

/** Opens an archive: its files, and its main file read as a document and validated (relative URIs allowed). */
export async function readSvsx(bytes: Uint8Array, opts: ReadSvsxOptions = {}): Promise<Svsx> {
  const files = new Map<string, Uint8Array>();
  const copy = bytes.slice().buffer;
  try {
    const limits = {
      maxEntries: opts.maxEntries ?? SvsxLimits.maxEntries,
      maxBytes: opts.maxBytes ?? SvsxLimits.maxBytes,
    };
    const entries = await unzip(RuntimeContext, copy, false, limits);
    for (const [path, file] of Object.entries(entries))
      if (file instanceof Uint8Array && !path.endsWith('/')) files.set(path, file);
  } catch (e) {
    const message = `the archive cannot be read: ${e instanceof Error ? e.message : String(e)}`;
    return { doc: null, issues: [Issue.error('svs.format.invalid', '', message)], main: '', base: '', files };
  }
  const main = opts.main ?? SvsxMainFiles.find((f) => files.has(f)) ?? '';
  const file = files.get(main);
  if (file === undefined) {
    const message = opts.main
      ? `the archive has no ${opts.main}`
      : `the archive has neither ${SvsxMainFiles.join(' nor ')}`;
    return { doc: null, issues: [Issue.error('svs.format.invalid', '', message)], main, base: '', files };
  }
  const base = main.includes('/') ? main.slice(0, main.lastIndexOf('/')) : '';
  const read = readDocument(file);
  const reg = opts.registries ?? createRegistries();
  if (!opts.registries) register(reg, [KernelRegistry, MvsRegistry]);
  const issues = read.doc ? validate(read.doc, reg, { strict: opts.strict, base }) : read.issues;
  return { doc: read.doc, issues, main, base, files };
}

/** An archive's bytes from its files by path (the zip of Mol*, copied: the builder's `writeSvsx` packs with it). */
export async function zipFiles(files: Readonly<Record<string, Uint8Array>>): Promise<Uint8Array> {
  return new Uint8Array(await zip(RuntimeContext, files as Record<string, Uint8Array<ArrayBuffer>>));
}

/**
 * Loaders that answer a relative URI from the archive's files (`base`: the main file's directory) and pass everything
 * else to `fallback`. A relative URI the archive lacks, or that climbs above its root, fails to load. A streamed
 * source's relative `uri` reaches the fallback's `stream` as its path from the archive's root.
 */
export function archiveLoaders(files: ReadonlyMap<string, Uint8Array>, fallback: Loaders, base = ''): Loaders {
  const relative = (uri: string | null) => uri !== null && isRelativeUri(uri);
  return {
    base,
    async load(request, signal) {
      if (!relative(request.uri)) return fallback.load(request, signal);
      const path = resolveArchivePath(request.uri as string, base);
      const file = path === null ? undefined : files.get(path);
      if (file === undefined) {
        throw new Error(`${request.uri} is ${path === null ? 'outside the archive' : `not in the archive (${path})`}`);
      }
      return file;
    },
    parsed:
      fallback.parsed &&
      ((request, format, signal) =>
        relative(request.uri)
          ? Promise.resolve(null)
          : (fallback.parsed?.(request, format, signal) ?? Promise.resolve(null))),
    stream:
      fallback.stream &&
      ((source: StreamSource) => {
        const path = source.uri !== null && relative(source.uri) ? resolveArchivePath(source.uri, base) : null;
        return fallback.stream?.(path === null ? source : { ...source, uri: path }) ?? null;
      }),
  };
}
