/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printJson, parseJson } from '@molstar/core/repr/codecs/json';
import { type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import type { RuntimeContext } from '@molstar/core/task';
import { unzip, zip } from '@molstar/core/util/zip/zip';
import { checkSpecBounds, DefaultSpecBounds } from '@molstar/seqx-core/lang';
import { ArchiveIssues } from './issues.js';
import type { AssetReader } from './protocol.js';

// Archives over the assets port (_spec/README.md "Data", Archives): a spec with its data files and an index as one zip,
// read and written with `@molstar/core/util/zip`. An opened archive is a spec and its files; `archiveReader` is the
// assets port's reader over the files, so a source's `uri` names an entry and is read as any asset is. Bounded before
// anything is inflated (EC-91): the central directory is read first, and its names, counts, declared sizes and ratios
// are the bounds' to check; a name is a path under the archive's root, never above it.

/** What an archive may hold (EC-91): checked from the central directory before an entry is inflated. */
export interface ArchiveLimits {
  /** Entries (directories included): 1,000. */
  readonly maxEntries: number;
  /** The entries' expanded sizes together: 256 MB. */
  readonly maxBytes: number;
  /** One entry's expanded size: 64 MB, the bound of a file parsed whole (`OpenLimits.maxBytes`). */
  readonly maxEntryBytes: number;
  /** The expanded size over the compressed, of an entry larger than `ratioFloor`: 200 (deflate's own is 1,032). */
  readonly maxRatio: number;
  /** An entry no larger than this is not held to the ratio: 64 KB. */
  readonly ratioFloor: number;
  /** A name's length, in characters: 256. */
  readonly maxNameLength: number;
}

export const DefaultArchiveLimits: ArchiveLimits = {
  maxEntries: 1000,
  maxBytes: 256 * 1024 * 1024,
  maxEntryBytes: 64 * 1024 * 1024,
  maxRatio: 200,
  ratioFloor: 64 * 1024,
  maxNameLength: 256,
};

/** The entry that says what the archive is and which entries are its spec and its data. */
export const ArchiveIndexName = 'index.json';
/** The spec's entry, unless the index names another. */
export const ArchiveSpecName = 'spec.json';
/** The index's `archive` member. */
export const ArchiveKind = 'seqx';
export const ArchiveVersion = 1;

/** An opened archive: the spec (as written, not yet checked against the language) and the files its index lists. */
export interface Archive {
  readonly spec: Node;
  readonly files: ReadonlyMap<string, Uint8Array>;
}

/** An archive read, or why not (the issues are the codec's, the bounds' and `seq.archive.*`). */
export interface ArchiveRead {
  readonly archive?: Archive;
  readonly issues: readonly Issue[];
}

/** Why an entry's name is no path under the root (EC-91), or undefined when it is one. A trailing `/` is a directory. */
export function checkArchiveName(name: string, limits: ArchiveLimits = DefaultArchiveLimits): string | undefined {
  if (name.length === 0) return 'it is empty';
  if (name.length > limits.maxNameLength) return `it is longer than ${limits.maxNameLength}`;
  if (/[\u0000-\u001f\u007f]/.test(name)) return 'it has a control character';
  if (name.includes('\\')) return 'it has a backslash';
  if (name.startsWith('/')) return 'it is absolute';
  if (/^[A-Za-z]:/.test(name)) return 'it names a drive';
  const segments = (name.endsWith('/') ? name.slice(0, -1) : name).split('/');
  for (const s of segments) {
    if (s === '') return 'it has an empty segment';
    if (s === '..') return 'it goes above the root';
    if (s === '.') return 'it has a `.` segment';
  }
  return undefined;
}

const bad = (code: keyof typeof ArchiveIssues.codes, message: string): Issue => issueOf(ArchiveIssues, code, message);

function buffer(bytes: Uint8Array): ArrayBuffer {
  const b = bytes.buffer as ArrayBuffer;
  return bytes.byteOffset === 0 && bytes.byteLength === b.byteLength ? b : (bytes.slice().buffer as ArrayBuffer);
}

/** The bounds an archive's listing (names and declared sizes) must meet: the first one broken, or undefined. */
function checkListing(
  listing: Readonly<Record<string, { size: number; csize: number }>>,
  limits: ArchiveLimits,
): Issue | undefined {
  const folded = new Map<string, string>();
  for (const name of Object.keys(listing)) {
    const why = checkArchiveName(name, limits);
    if (why) return bad('seq.archive.name', `The entry ${JSON.stringify(name)} is refused: ${why}`);
    const key = name.normalize('NFC').toLowerCase();
    const other = folded.get(key);
    if (other !== undefined)
      return bad('seq.archive.name', `The entries ${JSON.stringify(other)} and ${JSON.stringify(name)} are one name`);
    folded.set(key, name);
    const { size, csize } = listing[name];
    if (size > limits.maxEntryBytes)
      return bad(
        'seq.archive.bounds',
        `The entry ${JSON.stringify(name)} expands to ${size} bytes, more than ${limits.maxEntryBytes}`,
      );
    if (size > limits.ratioFloor && size > limits.maxRatio * csize)
      return bad(
        'seq.archive.bounds',
        `The entry ${JSON.stringify(name)} expands ${size} bytes from ${csize}, a ratio over ${limits.maxRatio}`,
      );
  }
  return undefined;
}

const decoder = new TextDecoder('utf-8', { fatal: true });

/** The index's shape: its spec's entry and its data files', or why not. */
function readIndex(text: string, limits: ArchiveLimits): { spec: string; files: string[] } | string {
  let json: unknown;
  try {
    json = JSON.parse(text);
  } catch {
    return 'the index is not JSON';
  }
  const o = json as { archive?: unknown; version?: unknown; spec?: unknown; files?: unknown };
  if (typeof json !== 'object' || json === null || Array.isArray(json) || o.archive !== ArchiveKind)
    return `the index is not one of a \`${ArchiveKind}\` archive`;
  if (typeof o.version !== 'number' || o.version > ArchiveVersion)
    return `the index is of version ${String(o.version)}`;
  if (typeof o.spec !== 'string' || checkArchiveName(o.spec, limits) || o.spec.endsWith('/'))
    return 'the index names no spec entry';
  if (!Array.isArray(o.files) || o.files.length > limits.maxEntries) return 'the index lists no files';
  for (const f of o.files)
    if (typeof f !== 'string' || checkArchiveName(f, limits) || f.endsWith('/'))
      return `the index lists ${JSON.stringify(f)}`;
  return { spec: o.spec, files: o.files as string[] };
}

/**
 * Opens an archive (never throws on bad bytes: the issues say why). The listing is checked first, then the index, the
 * spec and the files it lists are inflated; entries it does not list are never inflated.
 */
export async function readArchive(
  bytes: Uint8Array,
  options: { readonly limits?: Partial<ArchiveLimits>; readonly runtime?: RuntimeContext } = {},
): Promise<ArchiveRead> {
  const limits = { ...DefaultArchiveLimits, ...options.limits };
  const rt = options.runtime ?? SyncRuntimeContext;
  const fail = (...issues: Issue[]): ArchiveRead => ({ issues });
  if (bytes.byteLength > limits.maxBytes)
    return fail(bad('seq.archive.bounds', `The archive has ${bytes.byteLength} bytes, more than ${limits.maxBytes}`));
  const buf = buffer(bytes);
  let listing: Awaited<ReturnType<typeof unzip>>;
  try {
    listing = await unzip(rt, buf, true, { maxEntries: limits.maxEntries, maxBytes: limits.maxBytes });
  } catch (e) {
    const message = e instanceof Error ? e.message : String(e);
    return fail(
      bad(/bound/.test(message) ? 'seq.archive.bounds' : 'seq.archive.unreadable', `The archive: ${message}`),
    );
  }
  const sizes = listing as Readonly<Record<string, { size: number; csize: number }>>;
  const refused = checkListing(sizes, limits);
  if (refused) return fail(refused);

  const inflate = async (want: ReadonlySet<string>) => {
    try {
      return (await unzip(rt, buf, false, undefined, (n) => want.has(n))) as Record<string, Uint8Array>;
    } catch (e) {
      throw bad('seq.archive.unreadable', `The archive: ${e instanceof Error ? e.message : String(e)}`);
    }
  };
  try {
    if (!(ArchiveIndexName in sizes)) return fail(bad('seq.archive.index', `The archive has no ${ArchiveIndexName}`));
    const indexText = decoder.decode((await inflate(new Set([ArchiveIndexName])))[ArchiveIndexName]);
    const index = readIndex(indexText, limits);
    if (typeof index === 'string') return fail(bad('seq.archive.index', index));
    for (const n of [index.spec, ...index.files])
      if (!(n in sizes))
        return fail(bad('seq.archive.index', `The index lists ${JSON.stringify(n)}, which the archive does not hold`));
    const entries = await inflate(new Set([index.spec, ...index.files]));
    const spec = entries[index.spec];
    if (spec.byteLength > DefaultSpecBounds.maxBytes)
      return fail(
        bad('seq.archive.spec', `The spec has ${spec.byteLength} bytes, more than ${DefaultSpecBounds.maxBytes}`),
      );
    const decoded = parseJson(decoder.decode(spec));
    if (!decoded.node) return { issues: [bad('seq.archive.spec', 'The spec is not a tree'), ...decoded.issues] };
    const bounded = checkSpecBounds(decoded.node);
    if (bounded.length) return { issues: [bad('seq.archive.spec', 'The spec is past the spec bounds'), ...bounded] };
    const files = new Map<string, Uint8Array>();
    for (const f of index.files) files.set(f, entries[f]);
    return { archive: { spec: decoded.node, files }, issues: decoded.issues };
  } catch (e) {
    if (typeof e === 'object' && e !== null && 'severity' in e && 'code' in e) return fail(e as Issue);
    return fail(
      bad(
        'seq.archive.spec',
        `The archive's text is not UTF-8 or not readable: ${e instanceof Error ? e.message : String(e)}`,
      ),
    );
  }
}

/** A location an archive does not hold. */
export class ArchiveEntryMissing extends Error {
  constructor(readonly uri: string) {
    super(`The location is not an entry of the archive: ${uri}`);
    this.name = 'ArchiveEntryMissing';
  }
}

/** The location as an entry's name: relative to the archive's root, a leading `./` dropped. */
function entryOf(uri: string): string {
  return uri.replace(/^(\.\/)+/, '');
}

/** The assets port's reader over an archive's files: a location is an entry's name; anything else is missing. */
export function archiveReader(archive: Archive): AssetReader {
  const read = (uri: string, signal: AbortSignal) => {
    if (signal.aborted) throw signal.reason ?? new Error('aborted');
    const bytes = archive.files.get(entryOf(uri));
    if (!bytes) throw new ArchiveEntryMissing(uri);
    return bytes;
  };
  return {
    readText: async (uri, signal) => new TextDecoder().decode(read(uri, signal)),
    readBytes: async (uri, signal) => read(uri, signal),
  };
}

/** What `writeArchive` says: the bytes, or the issues that stopped it. */
export interface ArchiveWritten {
  readonly bytes?: Uint8Array;
  readonly issues: readonly Issue[];
}

const Stored = /\.(msgpack|gz|zip|bgz|bw|bigwig|png|jpe?g)$/i;
const HasScheme = /^[a-z][a-z0-9+.-]*:/i;

/** The uris a source names (`uri`, `index`), in order of appearance, each once. */
function sourceUris(node: Node, out: Set<string>): void {
  if (node.kind === 'source')
    for (const k of ['uri', 'index'] as const) {
      const v = node.params?.[k];
      if (typeof v === 'string') out.add(v);
    }
  for (const c of node.children ?? []) sourceUris(c, out);
}

function rewrite(node: Node, names: ReadonlyMap<string, string>): Node {
  let params = node.params;
  if (node.kind === 'source' && params) {
    const next: Record<string, Value> = { ...params };
    for (const k of ['uri', 'index'] as const) {
      const v = params[k];
      const name = typeof v === 'string' ? names.get(v) : undefined;
      if (name !== undefined) next[k] = name;
    }
    if (next.uri !== params.uri || next.index !== params.index) params = next;
  }
  let children = node.children;
  if (children) {
    const next = children.map((c) => rewrite(c, names));
    if (next.some((c, i) => c !== children![i])) children = next;
  }
  return params === node.params && children === node.children
    ? node
    : { ...node, ...(params ? { params } : {}), ...(children ? { children } : {}) };
}

/** An entry's name for a location: its own path when that is one under the root, else `data/<n>-<its last segment>`. */
function nameFor(uri: string, taken: ReadonlySet<string>, limits: ArchiveLimits): string {
  const own = entryOf(uri);
  const reserved = own === ArchiveIndexName || own === ArchiveSpecName;
  if (
    !HasScheme.test(own) &&
    !/[?#]/.test(own) &&
    !own.endsWith('/') &&
    !reserved &&
    !checkArchiveName(own, limits) &&
    !taken.has(own)
  )
    return own;
  const last = own.split(/[\\/]/).filter(Boolean).pop() ?? 'file';
  const base = last.replace(/[^A-Za-z0-9._-]/g, '_').slice(-64) || 'file';
  for (let n = taken.size; ; n++) {
    const name = `data/${n}-${base}`;
    if (!taken.has(name)) return name;
  }
}

/**
 * Writes a spec with the data files its sources name (`uri`, `index`) as an archive: each is read through the assets
 * port's reader (its policy decides what may be read), saved as an entry, and the location in the written spec is the
 * entry's name. The same bounds as reading apply, so what is written opens. A file that cannot be read stops the
 * write: an archive that lacks a file its spec names is not the view saved.
 */
export async function writeArchive(
  spec: Node,
  options: {
    readonly assets: AssetReader;
    /**
     * Locations a format runtime of the host serves itself (`generated:chr17`: no file), left as written and never
     * read. The archive does not hold them: opened where that runtime is not, such a source is the host's to answer.
     */
    readonly keep?: (uri: string) => boolean;
    readonly signal?: AbortSignal;
    readonly limits?: Partial<ArchiveLimits>;
    readonly runtime?: RuntimeContext;
  },
): Promise<ArchiveWritten> {
  const limits = { ...DefaultArchiveLimits, ...options.limits };
  const signal = options.signal ?? new AbortController().signal;
  const uris = new Set<string>();
  sourceUris(spec, uris);
  const names = new Map<string, string>();
  const taken = new Set<string>();
  const files: Record<string, Uint8Array<ArrayBuffer>> = Object.create(null);
  const issues: Issue[] = [];
  let total = 0;
  for (const uri of uris) {
    if (options.keep?.(uri)) continue;
    const name = nameFor(uri, taken, limits);
    taken.add(name);
    names.set(uri, name);
    let data: Uint8Array;
    try {
      data = await options.assets.readBytes(uri, signal);
    } catch (e) {
      if (signal.aborted) throw e;
      issues.push(
        bad(
          'seq.archive.data',
          `${JSON.stringify(uri)} could not be read: ${e instanceof Error ? e.message : String(e)}`,
        ),
      );
      continue;
    }
    total += data.byteLength;
    if (data.byteLength > limits.maxEntryBytes || total > limits.maxBytes || taken.size + 2 > limits.maxEntries) {
      issues.push(bad('seq.archive.bounds', `${JSON.stringify(uri)} takes the archive past its bounds`));
      continue;
    }
    files[name] = data.slice() as Uint8Array<ArrayBuffer>;
  }
  if (issues.length) return { issues };
  const text = printJson(rewrite(spec, names));
  const encoder = new TextEncoder();
  const index = JSON.stringify({
    archive: ArchiveKind,
    version: ArchiveVersion,
    spec: ArchiveSpecName,
    files: [...taken],
  });
  const all: Record<string, Uint8Array<ArrayBuffer>> = Object.create(null);
  all[ArchiveIndexName] = encoder.encode(index);
  all[ArchiveSpecName] = encoder.encode(text);
  if (all[ArchiveSpecName].byteLength > DefaultSpecBounds.maxBytes)
    return { issues: [bad('seq.archive.spec', 'The spec is past the spec bounds')] };
  for (const n of taken) all[n] = files[n];
  const noCmpr = new Set([...taken].filter((n) => Stored.test(n)));
  return { bytes: new Uint8Array(await zip(options.runtime ?? SyncRuntimeContext, all, noCmpr)), issues: [] };
}
