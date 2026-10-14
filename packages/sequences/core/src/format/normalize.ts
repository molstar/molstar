/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `normalize` (every default filled, tables as columns with codes decoded: seqspec:format#what-a-change-affects);
// documents to and from JSON and MessagePack (seqspec:format#values).

import {
  columnNumbers,
  decodeMsgPack,
  decodeTypedArrays,
  encodeMsgPack,
  encodeTypedArrays,
  type Issue,
} from '@molstar/seq-core/model';
import type { DocNode, SeqViewSpecDoc } from './doc.js';
import type { SvsRegistries } from './registry.js';
import { isPlainObject } from './schema/field-schema.js';
import { normalizeParams } from './schema/params-schema.js';
import type { TreeSchema } from './schema/tree-schema.js';
import { SvsTreeSchema } from './tables/svs-tree.js';
import { checkEnvelope, EnvelopeSchema } from './validate.js';

/** Every default filled; idempotent. Unknown kinds, cases and namespaces are kept as they are. */
export function normalize(doc: SeqViewSpecDoc, reg: SvsRegistries, schema: TreeSchema = SvsTreeSchema): SeqViewSpecDoc {
  const env = normalizeParams(EnvelopeSchema, doc as unknown as Record<string, unknown>, reg);
  return { ...env, root: normalizeNode(doc.root, schema, reg) } as SeqViewSpecDoc;
}

function normalizeNode(n: DocNode, schema: TreeSchema, reg: SvsRegistries): DocNode {
  const ns = isPlainObject(n) && Object.hasOwn(schema.nodes, n.kind) ? schema.nodes[n.kind] : undefined;
  if (!ns) return n;
  const out: DocNode = { kind: n.kind };
  const ref = n.ref ?? (ns.ref === 'implied' ? n.kind : undefined);
  if (ref !== undefined) out.ref = ref;
  const { params, custom } = normalizeNodeParams(n, schema, reg);
  out.params = params;
  if (custom) out.custom = custom;
  if (n.children)
    out.children = Array.isArray(n.children) ? n.children.map((c) => normalizeNode(c, schema, reg)) : n.children;
  return out;
}

/** A node's params and `custom` with every default filled (an unknown kind's as they are). */
export function normalizeNodeParams(n: DocNode, schema: TreeSchema, reg: SvsRegistries) {
  const ns = Object.hasOwn(schema.nodes, n.kind) ? schema.nodes[n.kind] : undefined;
  if (!ns) return { params: n.params ?? {}, custom: n.custom };
  const params = normalizeParams(ns.params, n.params ?? {}, reg);
  const owned = ns.owned?.schema(params, reg);
  if (ns.owned && owned) {
    const value = params[ns.owned.param];
    params[ns.owned.param] = normalizeParams(owned, isPlainObject(value) ? value : {}, reg);
  }
  const custom =
    n.custom &&
    Object.fromEntries(
      Object.entries(n.custom).map(([key, value]) => {
        const entry = reg.get('custom', `${key}/${n.kind}`);
        return [key, entry && isPlainObject(value) ? normalizeParams(entry.schema, value, reg) : value];
      }),
    );
  return { params, custom };
}

type Obj = Record<string, unknown>;
/** `v` with its typed arrays as the arrays of their numbers, float32 as its shortest decimals (JSON
 *  would write an object of indices); only the objects and arrays on the way to one copied, the caller's untouched. */
function plainArrays(v: object): object {
  if (ArrayBuffer.isView(v)) return v instanceof Float32Array ? columnNumbers(v, null) : Array.from(v as Float64Array);
  let out: Obj | undefined;
  if (Array.isArray(v)) for (let i = 0; i < v.length; i++) out = plainAt(v as never, i, out);
  else for (const k in v) if (Object.hasOwn(v, k)) out = plainAt(v as Obj, k, out);
  return out ?? v;
}
/** `out`, or a copy of `v` made once its `k` changes, with `v[k]` through `plainArrays`. */
function plainAt(v: Obj, k: string | number, out: Obj | undefined): Obj | undefined {
  const x = v[k];
  const y = x !== null && typeof x === 'object' ? plainArrays(x) : x;
  if (y === x) return out;
  const o = out ?? ((Array.isArray(v) ? v.slice() : { ...v }) as Obj);
  o[k] = y;
  return o;
}

export function toJSON(doc: SeqViewSpecDoc, opts?: { pretty?: boolean }): string {
  return JSON.stringify(plainArrays(doc), null, opts?.pretty ? 2 : undefined);
}

/**
 * MessagePack: arrays held as typed arrays written as the typed-array extension; `typedArrays: false` writes every
 * array plain, for decoders without the extension.
 */
export function toMsgPack(doc: SeqViewSpecDoc, opts?: { typedArrays?: boolean }): Uint8Array {
  if (opts?.typedArrays === false) return encodeMsgPack(JSON.parse(toJSON(doc)));
  return encodeMsgPack(doc, encodeTypedArrays);
}

export type DocumentRead = { doc: SeqViewSpecDoc | null; issues: Issue[] };
/** Parses a document and checks its envelope (format, version); `validate` checks the rest. */
export function fromJSON(text: string): DocumentRead {
  return parsed(() => JSON.parse(text), 'JSON');
}

/** As `fromJSON`; the typed-array extension read as typed arrays. */
export function fromMsgPack(bytes: Uint8Array): DocumentRead {
  return parsed(() => decodeMsgPack(bytes, decodeTypedArrays), 'MessagePack');
}

/** A MessagePack document starts with a map header, a JSON one with `{` or whitespace. */
export const isMsgPack = (bytes: Uint8Array) => (bytes[0] & 0xf0) === 0x80 || bytes[0] === 0xde || bytes[0] === 0xdf;

/** A document file's bytes, JSON or MessagePack. */
export function readDocument(bytes: Uint8Array): DocumentRead {
  return isMsgPack(bytes) ? fromMsgPack(bytes) : fromJSON(new TextDecoder().decode(bytes));
}

function parsed(parse: () => unknown, what: string): DocumentRead {
  let value: unknown;
  try {
    value = parse();
  } catch (e) {
    const message = `not ${what}: ${e instanceof Error ? e.message : String(e)}`;
    return { doc: null, issues: [{ code: 'svs.format.invalid', severity: 'error', path: '', message }] };
  }
  const issues: Issue[] = [];
  const ok = checkEnvelope(value, issues) !== null;
  return { doc: ok ? (value as SeqViewSpecDoc) : null, issues };
}
