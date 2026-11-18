/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A patch is also a tree (_spec/README.md "Patches", A patch as a tree): kind `patch`, one child per op, laid out by
// `patchToNode` the same way for every language. `at` and `parent` are plain address strings (not `ref`: the stored-
// patch rule walks them, refs.ts), embedded nodes slots of the target language. `patchLanguage(language)` is one
// object per language, its id the language's with `.patch` after it, so a slot names it; `AnyPatchLanguage` is the
// same over no particular language (embedded nodes `any`). Both print, encode and decode through the existing codecs.

import { defineLanguage, type Fields, type Language, type ParamType } from '../schema/index.js';
import {
  type Address,
  type Issue,
  isExactAddress,
  isNode,
  issueOf,
  type Node,
  own,
  parseAddress,
  printAddress,
  type Value,
} from '../tree/index.js';
import { type Op, type Patch, PatchIssues, type Place } from './patch.js';

const Revision = { type: 'string', pattern: '^[0-9a-f]{32}$' } as const;
const At = { type: 'string', required: true, doc: 'The address of the node the op changes (ids only).' } as const;
const Key = { type: 'string', required: true, doc: 'The key it sets or unsets.' } as const;
const Where = {
  type: 'string',
  enum: ['params', 'custom'],
  doc: 'Whether the key is a param (the default) or custom data.',
} as const;
const Parent = { type: 'string', required: true, doc: 'The address of the new parent (ids only).' } as const;
const PlaceFields = {
  after: { type: 'string', doc: 'The id of the sibling it goes after.' },
  first: { type: 'boolean', doc: 'It goes first among its siblings.' },
  last: { type: 'boolean', doc: 'It goes last among its siblings.' },
} as const;

/** Each patch language's target (undefined for `AnyPatchLanguage`): what the stored-patch rule types refs by. */
const Targets = new WeakMap<Language, Language | undefined>();
const ByLanguage = new WeakMap<Language, Language>();

function describe(id: string, version: string, node: ParamType, target: Language | undefined): Language {
  const embedded = { ...node, required: true, doc: 'The node it inserts or replaces with.' } as Fields[string];
  const op = (doc: string, params: Fields) => ({ doc, parents: ['patch'], ids: 'forbidden' as const, params });
  const language = defineLanguage({
    id,
    version,
    root: 'patch',
    doc: `Patches${target ? ` of ${target.id}` : ''}: ops by address of ids, applied in order, all or nothing.`,
    kinds: {
      patch: {
        doc: 'A patch: its ops in order.',
        parents: [],
        ids: 'forbidden',
        params: {
          base: { ...Revision, doc: 'The revision it was made against.' },
          result: { ...Revision, doc: 'The revision it gives.' },
        },
      },
      set: op('Sets one key.', {
        at: At,
        key: Key,
        value: { type: 'any', required: true, doc: 'The value.' },
        where: Where,
      }),
      unset: op('Unsets one key.', { at: At, key: Key, where: Where }),
      insert: op('Inserts a node under a parent, at a place.', { parent: Parent, ...PlaceFields, node: embedded }),
      remove: op('Removes a node and its subtree.', { at: At }),
      move: op('Moves a node to a parent and a place, optionally under a new id.', {
        at: At,
        parent: Parent,
        ...PlaceFields,
        rename: { type: 'string', doc: 'Its id at the new place.' },
      }),
      replace: op('Replaces a whole node, its id kept.', { at: At, node: embedded }),
    },
  });
  Targets.set(language, target);
  return language;
}

/** The patch language of a language: memoized, its id `<language id>.patch`. */
export function patchLanguage(language: Language): Language {
  let p = ByLanguage.get(language);
  if (!p) {
    p = describe(`${language.id}.patch`, language.version, { type: 'node', language: language.id }, language);
    ByLanguage.set(language, p);
  }
  return p;
}

/** Patches over no particular language: embedded nodes are `any`, addresses strings. */
export const AnyPatchLanguage: Language = describe('patch', '1.0.0', { type: 'any' }, undefined);

/** Whether a language is a patch language, and of what (undefined: `AnyPatchLanguage`). */
export function patchTarget(language: Language): { readonly target: Language | undefined } | undefined {
  return Targets.has(language) ? { target: Targets.get(language) } : undefined;
}

/** A patch as a tree: `(patch :base … (set :at /a :key k :value 1) …)`. */
export function patchToNode(patch: Patch): Node {
  const params: Record<string, Value> = {};
  if (patch.base !== undefined) params.base = patch.base;
  if (patch.result !== undefined) params.result = patch.result;
  const node: { kind: string; params?: Record<string, Value>; children: Node[] } = {
    kind: 'patch',
    children: patch.ops.map(opToNode),
  };
  if (Object.keys(params).length) node.params = params;
  return node;
}

function placeParams(place: Place): Record<string, Value> {
  if ('last' in place) return { last: true };
  return place.after === null ? { first: true } : { after: place.after };
}

function opToNode(op: Op): Node {
  const p = (params: Record<string, Value | undefined>): Node => {
    const out: Record<string, Value> = {};
    for (const [k, v] of Object.entries(params)) if (v !== undefined) out[k] = v;
    return { kind: op.op, params: out };
  };
  const a = (x: Address) => printAddress(x);
  switch (op.op) {
    case 'set':
      return p({ at: a(op.at), key: op.key, value: op.value, where: op.where });
    case 'unset':
      return p({ at: a(op.at), key: op.key, where: op.where });
    case 'insert':
      return p({ parent: a(op.parent), ...placeParams(op.place), node: op.node });
    case 'remove':
      return p({ at: a(op.at) });
    case 'move':
      return p({ at: a(op.at), parent: a(op.parent), ...placeParams(op.place), rename: op.rename });
    case 'replace':
      return p({ at: a(op.at), node: op.node });
  }
}

/** A patch read from its tree, or the issues of why not (`repr.patch.op`, with the patch tree's path). */
export function patchFromNode(node: Node): { readonly patch?: Patch; readonly issues: readonly Issue[] } {
  const issues: Issue[] = [];
  const bad = (message: string, ...path: (string | number)[]) =>
    issues.push(issueOf(PatchIssues, 'repr.patch.op', message, { path }));
  if (!isNode(node) || node.kind !== 'patch') {
    bad('A patch is a node of kind patch', 'kind');
    return { issues };
  }
  const str = (v: Value | undefined) => (typeof v === 'string' ? v : undefined);
  const base = str(node.params?.base);
  const result = str(node.params?.result);
  const ops: Op[] = [];
  (node.children ?? []).forEach((c, i) => {
    const params = c.params ?? {};
    const get = (k: string) => own(params, k);
    const address = (k: string): Address | undefined => {
      const text = get(k);
      const parsed = typeof text === 'string' ? parseAddress(text).pattern : undefined;
      if (!parsed || parsed.subtree || !isExactAddress(parsed)) {
        bad('An address of ids', 'children', i, 'params', k);
        return undefined;
      }
      return parsed.segments as Address;
    };
    const key = (): string | undefined => {
      const k = get('key');
      if (typeof k !== 'string') bad('A key is a string', 'children', i, 'params', 'key');
      return typeof k === 'string' ? k : undefined;
    };
    const where = () => {
      const w = get('where');
      if (w !== undefined && w !== 'params' && w !== 'custom')
        bad('params or custom', 'children', i, 'params', 'where');
      return w as 'params' | 'custom' | undefined;
    };
    const place = (): Place | undefined => {
      const after = get('after');
      const first = get('first') === true;
      const last = get('last') === true;
      const n = (after !== undefined ? 1 : 0) + (first ? 1 : 0) + (last ? 1 : 0);
      if (n !== 1 || (after !== undefined && typeof after !== 'string')) {
        bad('A place is one of after, first and last', 'children', i, 'params');
        return undefined;
      }
      return last ? { last: true } : { after: first ? null : (after as string) };
    };
    const embedded = (): Node | undefined => {
      const n = get('node');
      if (!isNode(n)) bad('A node', 'children', i, 'params', 'node');
      return isNode(n) ? n : undefined;
    };
    const before = issues.length;
    let op: Op | undefined;
    switch (c.kind) {
      case 'set': {
        const at = address('at');
        const k = key();
        const value = get('value');
        if (value === undefined) bad('A set has a value', 'children', i, 'params', 'value');
        const w = where();
        if (at && k !== undefined && value !== undefined)
          op = { op: 'set', at, key: k, value, ...(w ? { where: w } : {}) };
        break;
      }
      case 'unset': {
        const at = address('at');
        const k = key();
        const w = where();
        if (at && k !== undefined) op = { op: 'unset', at, key: k, ...(w ? { where: w } : {}) };
        break;
      }
      case 'insert': {
        const parent = address('parent');
        const pl = place();
        const n = embedded();
        if (parent && pl && n) op = { op: 'insert', parent, place: pl, node: n };
        break;
      }
      case 'remove': {
        const at = address('at');
        if (at) op = { op: 'remove', at };
        break;
      }
      case 'move': {
        const at = address('at');
        const parent = address('parent');
        const pl = place();
        const rename = get('rename');
        if (rename !== undefined && typeof rename !== 'string') bad('A new id', 'children', i, 'params', 'rename');
        if (at && parent && pl)
          op = { op: 'move', at, parent, place: pl, ...(typeof rename === 'string' ? { rename } : {}) };
        break;
      }
      case 'replace': {
        const at = address('at');
        const n = embedded();
        if (at && n) op = { op: 'replace', at, node: n };
        break;
      }
      default:
        bad('Not an op: set, unset, insert, remove, move or replace', 'children', i, 'kind');
    }
    if (op && issues.length === before) ops.push(op);
  });
  if (issues.length) return { issues };
  const patch: { ops: Op[]; base?: string; result?: string } = { ops };
  if (base !== undefined) patch.base = base;
  if (result !== undefined) patch.result = result;
  return { patch, issues };
}
