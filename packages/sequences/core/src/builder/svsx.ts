/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Writing `.svsx` archives (seqspec:sources#archives). The one module of the builder that
// imports `@molstar/seq-core/format` at run time (it writes MessagePack and zips), so it is its own entry,
// `@molstar/seq-core/builder/svsx`: the builder's main entry stays free of dependencies (seqspec:format#builders).

import type { DocNode, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { toJSON, toMsgPack, zipFiles } from '@molstar/seq-core/format';

export interface WriteSvsxOptions {
  /** The main file's form: MessagePack (`index.svsm`, the default) or JSON (`index.svsj`). */
  format?: 'msgpack' | 'json';
}

/** The archive of `doc` as its main file and `files` by path (bytes as their format reads them); a path is relative, inside the archive. */
export async function writeSvsx(
  doc: SeqViewSpecDoc,
  files: ReadonlyMap<string, Uint8Array> | Readonly<Record<string, Uint8Array>> = {},
  opts: WriteSvsxOptions = {},
): Promise<Uint8Array> {
  const json = opts.format === 'json';
  const main = json ? 'index.svsj' : 'index.svsm';
  const entries: Record<string, Uint8Array> = { [main]: json ? new TextEncoder().encode(toJSON(doc)) : toMsgPack(doc) };
  const list: [string, Uint8Array][] = files instanceof Map ? [...files] : Object.entries(files);
  for (const [path, bytes] of list) {
    if (
      path === '' ||
      path.startsWith('/') ||
      path.endsWith('/') ||
      path.split('/').some((s) => s === '' || s === '.' || s === '..')
    ) {
      throw new Error(`"${path}" is not a path inside an archive`);
    }
    if (path === 'index.svsm' || path === 'index.svsj')
      throw new Error(`"${path}" is the document's name, not a data file's`);
    entries[path] = bytes;
  }
  return zipFiles(entries);
}

/** A document whose absolute-URI sources are fetched into `files`, each URI's path in the archive its new relative `uri`. */
export interface Bundled {
  doc: SeqViewSpecDoc;
  files: Map<string, Uint8Array>;
}

/** The stable path of an absolute URI in the archive: `data/` and its host and path, odd characters `_`; a clash gets `-2`, `-3`. */
function pathOf(uri: string, taken: ReadonlySet<string>): string {
  const u = new URL(uri);
  const segs = `${u.host}${u.pathname}${u.search}`.split('/').map((s) => s.replace(/[^A-Za-z0-9._-]/g, '_'));
  const clean = segs.filter((s) => s !== '' && s !== '.' && s !== '..');
  const base = `data/${clean.join('/') || 'index'}`;
  let path = base;
  for (let n = 2; taken.has(path); n++) path = `${base}-${n}`;
  return path;
}

const bundled = (n: DocNode) => {
  const uri = n.params?.uri;
  return typeof uri === 'string' && uri.startsWith('https://') && n.params?.protocol == null ? uri : null;
};

/**
 * Fetches each absolute-URI source (a node's `uri` with no `protocol`) once through `fetch`, in the order the document
 * first names them (the tree, then the nodes its snapshots insert), and rewrites the URI to its path in `files`; the document is not changed. A streamed `protocol` source
 * is a tile server, left as it is. Rejects when a fetch fails or answers a non-OK status, naming the URI.
 */
export async function bundleUris(doc: SeqViewSpecDoc, fetch: (uri: string) => Promise<Response>): Promise<Bundled> {
  const uris: string[] = [];
  const scan = (n: DocNode) => {
    const uri = bundled(n);
    if (uri !== null && !uris.includes(uri)) uris.push(uri);
    n.children?.forEach(scan);
  };
  scan(doc.root);
  for (const s of doc.snapshots ?? []) for (const op of s.ops ?? []) if (op.op === 'insert') scan(op.node);
  const bytes = await Promise.all(
    uris.map(async (uri) => {
      const r = await fetch(uri);
      if (!r.ok) throw new Error(`${uri}: HTTP ${r.status}`);
      return new Uint8Array(await r.arrayBuffer());
    }),
  );
  const files = new Map<string, Uint8Array>();
  const paths = new Map<string, string>();
  uris.forEach((uri, i) => {
    const path = pathOf(uri, new Set(files.keys()));
    paths.set(uri, path);
    files.set(path, bytes[i]);
  });
  const rewrite = (n: DocNode): DocNode => {
    const uri = bundled(n);
    const children = n.children?.map(rewrite);
    return {
      ...n,
      ...(uri === null ? {} : { params: { ...n.params, uri: paths.get(uri) } }),
      ...(children ? { children } : {}),
    };
  };
  const snapshots = doc.snapshots?.map((s) =>
    s.ops ? { ...s, ops: s.ops.map((op) => (op.op === 'insert' ? { ...op, node: rewrite(op.node) } : op)) } : s,
  );
  return { doc: { ...doc, root: rewrite(doc.root), ...(snapshots ? { snapshots } : {}) }, files };
}
