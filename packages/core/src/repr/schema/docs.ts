/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What is generated from a description: defaults filled in, the reference docs (Markdown), and a signature per kind in
// the S-expression spelling, which is what an agent's skill lists (`(sequence :alphabet "protein" | "dna" | "rna"
// [:letters string])`). The docs and the skill are one source with the validator, never a second description.
// Written from meso*'s schema core.

import { isRecord, type Node, own, printKey, type Value, type ValueRecord } from '../tree/index.js';
import type { Fields, Language, ParamSpec, ParamType } from './language.js';
import { typeName } from './validate.js';

/**
 * The tree with every absent param that has a default set to it, in records too (and records in lists); kinds the
 * language lacks are left as they are, and so are slots (their language is resolved late).
 */
export function withDefaults(language: Language, node: Node): Node {
  const spec = own(language.kinds, node.kind);
  const params = spec?.params ? fill(spec.params, node.params) : node.params;
  const children = node.children?.map((c) => withDefaults(language, c));
  const changed = params !== node.params || children?.some((c, i) => c !== node.children![i]);
  return changed ? { ...node, ...(params ? { params } : {}), ...(children ? { children } : {}) } : node;
}

/** `values` with its absent defaults, its records' and its lists' records' too; the same object when nothing changes. */
function fill(specs: Fields, values: ValueRecord | undefined): ValueRecord | undefined {
  let out = values;
  const set = (key: string, v: Value) => {
    out = { ...out, [key]: v };
  };
  for (const [key, p] of Object.entries(specs)) {
    const v = values === undefined ? undefined : own(values, key);
    if (v === undefined) {
      if (p.default !== undefined) set(key, p.default);
    } else {
      const filled = deep(p, v);
      if (filled !== v) set(key, filled);
    }
  }
  return out;
}

function deep(t: ParamType, v: Value): Value {
  if (t.type === 'record' && isRecord(v)) return fill(t.fields, v) ?? v;
  if (t.type === 'map' && isRecord(v)) {
    let out: Record<string, Value> | undefined;
    for (const key of Object.keys(v)) {
      const filled = deep(t.of, v[key]);
      if (filled !== v[key]) (out ??= { ...v })[key] = filled;
    }
    return out ?? v;
  }
  if (t.type === 'list' && Array.isArray(v)) {
    const items = v.map((x) => deep(t.of, x));
    return items.some((x, i) => x !== v[i]) ? items : v;
  }
  return v;
}

const json = (v: Value) => JSON.stringify(v, (_, x) => (typeof x === 'bigint' ? `${x}n` : x));

/** The kind's form, `(kind :required type [:optional type])`. */
export function signature(language: Language, kind: string): string {
  const parts = Object.entries(own(language.kinds, kind)?.params ?? {}).map(([name, p]) =>
    p.required ? `${printKey(name)} ${typeName(p)}` : `[${printKey(name)} ${typeName(p)}]`,
  );
  return `(${[kind, ...parts].join(' ')})`;
}

function paramRows(fields: Fields, prefix = ''): string[] {
  return Object.entries(fields).flatMap(([name, p]: [string, ParamSpec]) => [
    `| \`${prefix}${name}\` | ${typeName(p).replaceAll('|', '\\|')} | ${p.required ? 'required' : p.default !== undefined ? `\`${json(p.default)}\`` : 'optional'} | ${p.doc} |`,
    ...(p.type === 'record' ? paramRows(p.fields, `${prefix}${name}.`) : []),
    ...(p.type === 'table' ? paramRows(p.columns, `${prefix}${name}.`) : []),
  ]);
}

/** The language's reference, as Markdown: one section per kind. */
export function languageToMarkdown(language: Language): string {
  const out = [`# ${language.id} ${language.version}`, '', language.doc, '', `The root is a \`${language.root}\`.`];
  for (const [kind, spec] of Object.entries(language.kinds)) {
    out.push('', `## \`${kind}\``, '', spec.doc, '', `\`${signature(language, kind)}\``, '');
    out.push(
      `- Sits under: ${spec.parents.length ? spec.parents.map((p) => `\`${p}\``).join(', ') : 'nothing (root)'}`,
    );
    out.push(`- Id: ${spec.ids ?? 'optional'}`);
    if (spec.fallback !== undefined) out.push(`- Fallback: \`${spec.fallback}\``);
    if (spec.params && Object.keys(spec.params).length) {
      out.push('', '| Param | Type | Default | What |', '| --- | --- | --- | --- |', ...paramRows(spec.params));
    }
  }
  return `${out.join('\n')}\n`;
}
