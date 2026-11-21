/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { validate } from '@molstar/core/repr/schema';
import type { Issue, Node, TreePath } from '@molstar/core/repr/tree';
import type { ExtensionTables, SeqExtension } from '@molstar/seqx-core/context';
import { builtInSchemes, SchemeRegistry, type SchemeEntry } from '@molstar/seqx-core/theme';
import { SeqLanguageVersion } from '@molstar/seqx-core/lang';
import type { LayerShape } from '@molstar/seqx-core/view';

// An extension's descriptor checked and completed where it is defined (_spec/README.md "Context", Extensions): here,
// in the entry an extension imports, so a composition with none bundles none of it.

/**
 * A version's rank and minor: the major (1.x on), or for 0.x its minor below every 1.x (`0.5` ranks 0.0005), so `5`
 * and `0.5` never compare equal; the minor is the next part (0.x: the patch).
 */
function rank(v: string): [number, number] {
  const [a = 0, b = 0, c = 0] = v.split('.').map(Number);
  return a ? [a, b] : [b / 1e4, c];
}

/** Whether `requires` covers a language version: the same major (0.x: minor), at least its minor (0.x: patch). */
export function requiresCovers(requires: string, version = SeqLanguageVersion): boolean {
  const [a, b] = rank(requires);
  const [la, lb] = rank(version);
  return la === a && lb >= b;
}

/**
 * How an extension at `has` reads a spec's `want` for it: `old` when the spec needs a newer major (0.x: minor) than it
 * is (drawn as the fallbacks), `newer` when a newer minor in its major (drawn, and said); else undefined.
 */
export function versionReads(has: string, want: string): 'old' | 'newer' | undefined {
  const [hm, hn] = rank(has);
  const [wm, wn] = rank(want);
  return hm < wm ? 'old' : hm === wm && hn < wn ? 'newer' : undefined;
}

const fail = (code: string, list: string[]) => {
  if (list.length) throw Object.assign(new Error(`${code}: ${list.join(', ')}`), { code });
};

/**
 * An extension's descriptor checked and completed (_spec/README.md "Context", Extensions): every entry and setting
 * named in its namespace, its id never `seq` (`seq.context.extension.namespace`, every offender in one error); its
 * `requires` covering the language (`seq.context.extension.version`); then `defined` (its layer types' kinds, its version
 * reading) and its `check` written from its `params` (a `check` it gives is kept). A composition takes only descriptors
 * made here. `viewer`: the viewer's maker calls it (a descriptor with viewer parts is refused otherwise).
 */
export function defineExtension<E extends SeqExtension>(extension: E, viewer = false): E {
  const id = extension.id;
  const bad: string[] = [];
  // a viewer's parts are added by the viewer's maker (`defineViewerExtension`): never dropped unread
  if (!viewer && ['behaviors', 'commands', 'bindings', 'headerActions'].some((k) => k in extension))
    bad.push(`${id}: viewer parts (made by defineViewerExtension)`);
  const ns = (name: string) => name.startsWith(`${id}:`) || bad.push(`${id}: ${name}`);
  if (!/^[a-z][a-z0-9_-]*$/.test(id) || id === 'seq') bad.push(`${id} (its id)`);
  const kinds: Record<string, { shape: LayerShape; fallback?: string }> = {};
  const schemes: SchemeEntry[] = [];
  for (const r of extension.registry as any[]) {
    // each entry's provider by its name (a format runtime's, a layer type's, a scheme's, a picture's; a query entry's own)
    const p = r.runtime ?? r.layout ?? r.sizing ?? r.scheme ?? r.picture ?? (r.kind === 'layer-type' ? r.type : r);
    ns(p.name);
    if (r.kind === 'layer-type' && p.shape)
      kinds[p.name] = { shape: p.shape, ...(p.fallback ? { fallback: p.fallback } : {}) };
    if (r.kind === 'scheme') schemes.push(p);
  }
  for (const c of extension.config ?? []) ns(c.key.replace('.', ':'));
  fail('seq.context.extension.namespace', bad);
  fail(
    'seq.context.extension.version',
    requiresCovers(extension.requires) ? [] : [`${id} requires ${extension.requires}`],
  );
  return {
    ...extension,
    defined: {
      compose: (t: ExtensionTables): ExtensionTables => ({
        kinds: {
          types: { ...t.kinds.types, ...kinds },
          extensions: { ...t.kinds.extensions, [id]: (want: string) => versionReads(extension.version, want) },
        },
        ...(schemes.length || t.schemes
          ? { schemes: new SchemeRegistry([...(t.schemes ?? builtInSchemes()).list(), ...schemes]) }
          : {}),
      }),
    },
    ...(!extension.check && extension.params ? { check: (spec: Node) => extensionParamIssues(spec, extension) } : {}),
  };
}

/**
 * The issues of an extension's params on a spec's nodes: a layer of its own type its plain params, a node of a kind it
 * names its keys in its namespace; their codes `repr.schema.*`, their message naming the extension, their path the
 * node's (the params' keys below it).
 */
function extensionParamIssues(spec: Node, extension: SeqExtension): Issue[] {
  const out: Issue[] = [];
  const read = Object.entries(extension.params ?? {});
  const id = extension.id;
  const walk = (n: Node, path: TreePath) => {
    for (const [by, fields] of read) {
      // its own layer type's or navigator picture's node
      const own = (n.kind === 'layer' && n.params?.type === by) || (n.kind === 'navigator' && n.params?.content === by);
      if (!own && n.kind !== by) continue;
      // its own type's params it declares, or its keys on another's node
      const params: Record<string, unknown> = {};
      for (const [k, v] of Object.entries(n.params ?? {}))
        if (k in fields || (!own && k.startsWith(`${id}:`))) params[k] = v;
      if (!own && !Object.keys(params).length) continue;
      const kinds = { p: { doc: '', parents: [], params: fields } };
      const lang = { id, version: extension.version, doc: '', root: 'p', kinds };
      for (const i of validate(lang, { kind: 'p', params: params as Node['params'] }))
        out.push({ ...i, message: `${id}: ${i.message}`, path: [...path, ...(i.path ?? [])] });
    }
    n.children?.forEach((c, k) => walk(c, [...path, 'children', k]));
  };
  walk(spec, []);
  return out;
}
