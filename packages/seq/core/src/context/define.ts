/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { validate } from '@molstar/core/repr/schema';
import type { Issue, Node, TreePath } from '@molstar/core/repr/tree';
import type { SeqExtension } from '@molstar/seqx-core/context';

// An extension's params checked where it is composed (_spec/README.md "Context", Extensions): here, in the entry an
// extension imports, so a composition with none bundles no check of them.

/**
 * An extension's descriptor with its `check` written from its `params` (a `check` it gives is kept): a layer of its own
 * type its plain params, any node of a kind it names the keys in its namespace, each against `params`.
 */
export function defineExtension<E extends SeqExtension>(extension: E): E {
  if (extension.check || !extension.params) return extension;
  return { ...extension, check: (spec: Node) => extensionParamIssues(spec, extension) };
}

/**
 * The issues of an extension's params on a spec's nodes: their codes `repr.schema.*`, their message naming the
 * extension, their path the node's (the params' keys below it).
 */
export function extensionParamIssues(spec: Node, extension: SeqExtension): Issue[] {
  const out: Issue[] = [];
  const read = Object.entries(extension.params ?? {});
  const id = extension.id;
  const walk = (n: Node, path: TreePath) => {
    for (const [by, fields] of read) {
      const own = n.kind === 'layer' && n.params?.type === by;
      if (!own && n.kind !== by) continue;
      // its own type's params it declares, or its keys on another's node
      const params: Record<string, unknown> = {};
      for (const [k, v] of Object.entries(n.params ?? {}))
        if (k in fields || (!own && k.startsWith(`${id}:`))) params[k] = v;
      if (!own && !Object.keys(params).length) continue;
      const kinds = { p: { doc: '', parents: [], params: fields } };
      for (const i of validate(
        { id, version: extension.version, doc: '', root: 'p', kinds },
        {
          kind: 'p',
          params: params as Node['params'],
        },
      ))
        out.push({ ...i, message: `${id}: ${i.message}`, path: [...path, ...(i.path ?? [])] });
    }
    n.children?.forEach((c, k) => walk(c, [...path, 'children', k]));
  };
  walk(spec, []);
  return out;
}
