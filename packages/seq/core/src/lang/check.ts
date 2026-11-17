/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the description cannot say (language-level checks, `seq.lang.*`): one of two params, what a method needs, what
// a reference names. `checkSeqSpec` is the description's validation, then these checks on a tree that passed it.

import { type Issue, issueOf, type Node, printAddress, type TreePath, type Value } from '@molstar/core/repr/tree';
import { exactAddress, hasErrors, isColumn } from '@molstar/core/repr/tree';
import { validate } from '@molstar/core/repr/schema';
import { LangIssues } from './issues.js';
import { SeqLanguage, SeqLanguageVersion } from './language.js';

type Code = keyof typeof LangIssues.codes;

/** A node with its place: its step path from the root, and its address when it has one. */
interface Placed {
  readonly node: Node;
  readonly path: TreePath;
  readonly address: string | undefined;
}

const Letters: Readonly<Record<string, RegExp>> = {
  protein: /^[A-Za-z*-]*$/,
  dna: /^[ACGTNRYKMSWBDHVacgtnrykmswbdhv-]*$/,
  rna: /^[ACGUNRYKMSWBDHVacgunrykmswbdhv-]*$/,
};

const MaxMajor = Number(SeqLanguageVersion.split('.')[0]);

/** The issues of a spec against the description (`repr.*`), and, when it is valid, the language's own (`seq.lang.*`). */
export function checkSeqSpec(spec: Node): Issue[] {
  const issues = validate(SeqLanguage, spec);
  if (hasErrors(issues)) return issues;
  return [...issues, ...checkSeqLang(spec)];
}

/** The `seq.lang.*` checks alone, for a tree that is valid against `SeqLanguage` (a normalized one). */
export function checkSeqLang(spec: Node): Issue[] {
  const issues: Issue[] = [];
  const all: Placed[] = [];
  const byAddress = new Map<string, Placed>();
  const visit = (node: Node, path: TreePath, address: readonly string[] | undefined) => {
    const here = address && node.id !== undefined ? [...address, node.id] : undefined;
    const placed: Placed = { node, path, address: here ? printAddress(here) : undefined };
    all.push(placed);
    if (placed.address) byAddress.set(placed.address, placed);
    node.children?.forEach((c, i) => visit(c, [...path, 'children', i], here));
  };
  // the root is `/`: its children are addressed by their ids alone
  all.push({ node: spec, path: [], address: '/' });
  spec.children?.forEach((c, i) => visit(c, ['children', i], []));

  const report = (code: Code, message: string, at: Placed, param?: string) => {
    const path = param ? [...at.path, 'params', param] : at.path;
    issues.push(
      issueOf(LangIssues, code, message, {
        hint: LangIssues.codes[code].hint,
        path,
        ...(at.address && at.address !== '/' ? { address: at.address } : {}),
      }),
    );
  };

  const p = (n: Node, key: string): Value | undefined => n.params?.[key];
  const target = (n: Node, key: string): Placed | undefined => {
    const a = exactAddress(p(n, key));
    return a ? byAddress.get(a) : undefined;
  };
  /** The inline letters of a sequence a param names, when it has them. */
  const lettersOf = (n: Node, key: string): string | undefined => {
    const t = target(n, key);
    const letters = t?.node.kind === 'sequence' ? p(t.node, 'letters') : undefined;
    return typeof letters === 'string' ? letters : undefined;
  };
  const xor = (a: boolean, b: boolean) => a !== b;

  // the version
  const major = Number(String(p(spec, 'version')).split('.')[0]);
  if (major > MaxMajor)
    report('seq.lang.version.unsupported', 'The spec is written in a newer major version', all[0], 'version');

  for (const at of all) {
    const n = at.node;
    switch (n.kind) {
      case 'data':
      case 'views':
        if (n.id !== n.kind) report('seq.lang.container.id', `The id of ${n.kind} is \`${n.kind}\``, at);
        break;
      case 'source':
        if (!xor(p(n, 'uri') !== undefined, p(n, 'value') !== undefined))
          report('seq.lang.source.content', 'A source has one of `uri` and `value`', at);
        break;
      case 'sequence': {
        if (!xor(p(n, 'letters') !== undefined, p(n, 'source') !== undefined))
          report('seq.lang.sequence.content', 'A sequence has one of `letters` and `source`', at);
        const letters = p(n, 'letters');
        if (typeof letters === 'string' && !Letters[String(p(n, 'alphabet'))]?.test(letters))
          report('seq.lang.sequence.letters', 'The letters are not of the alphabet', at, 'letters');
        break;
      }
      case 'ranges': {
        const rows = p(n, 'rows') as Record<string, Value> | undefined;
        if (!xor(rows !== undefined, p(n, 'source') !== undefined))
          report('seq.lang.ranges.content', 'A ranges node has one of `rows` and `source`', at);
        if (rows) {
          const { start, end } = rows;
          if (isColumn(start) && isColumn(end)) {
            const s = start as ArrayLike<number>;
            const e = end as ArrayLike<number>;
            const length = lettersOf(n, 'target')?.length;
            let bad = false;
            let past = false;
            for (let i = 0; i < Math.min(s.length, e.length); i++) {
              if (s[i] < 0 || e[i] < s[i]) bad = true;
              else if (length !== undefined && e[i] > length) past = true;
            }
            if (bad) report('seq.lang.ranges.interval', 'A range is not an interval', at, 'rows');
            else if (past) report('seq.lang.ranges.bounds', 'A range ends past the sequence', at, 'rows');
          }
        }
        break;
      }
      case 'values': {
        const values = p(n, 'values');
        if (!xor(values !== undefined, p(n, 'source') !== undefined))
          report('seq.lang.values.content', 'A values node has one of `values` and `source`', at);
        const length = lettersOf(n, 'target')?.length;
        if (isColumn(values) && length !== undefined && (values as ArrayLike<number>).length !== length)
          report('seq.lang.values.length', 'There is not one value per position', at, 'values');
        break;
      }
      case 'map': {
        const method = p(n, 'method');
        const has = (k: string) => p(n, k) !== undefined;
        const fits =
          method === 'offset'
            ? has('by') && !has('exons') && !has('strand')
            : has('exons') && has('strand') && !has('by');
        if (!fits) report('seq.lang.map.params', 'The params do not fit the method', at);
        if (method === 'codon') {
          const from = target(n, 'from')?.node;
          const to = target(n, 'to')?.node;
          if (from && to) {
            const f = p(from, 'alphabet');
            const t = p(to, 'alphabet');
            if (f !== 'protein' || (t !== 'dna' && t !== 'rna'))
              report('seq.lang.map.alphabet', 'A codon map goes from protein to dna or rna', at);
          }
        }
        break;
      }
      case 'layer': {
        const type = p(n, 'type');
        const data = p(n, 'data');
        if (data === undefined) {
          if (type !== 'letters') report('seq.lang.layer.data', 'The layer has no data', at);
        } else {
          const kind = target(n, 'data')?.node.kind;
          const want = type === 'letters' ? 'sequence' : type;
          if (kind !== undefined && kind !== want)
            report('seq.lang.layer.data', 'The data is not what the type draws', at, 'data');
        }
        break;
      }
    }
    // a reference names a node, of a kind its param takes
    const spec_ = SeqLanguage.kinds[n.kind as keyof typeof SeqLanguage.kinds];
    const fields = (spec_ as { params?: Record<string, { type: string; kinds?: readonly string[]; in?: string }> })
      .params;
    if (!fields) continue;
    for (const key of Object.keys(n.params ?? {})) {
      const f = fields[key];
      if (f?.type !== 'ref' || f.in !== undefined) continue;
      const a = exactAddress(n.params?.[key]);
      if (a === undefined) continue; // a query form: a resolve pass's
      const t = byAddress.get(a);
      if (!t || (f.kinds && !f.kinds.includes(t.node.kind)))
        report('seq.lang.ref.target', 'The reference names no node of a kind this param takes', at, key);
    }
  }

  // follow chains end
  const follows = new Map<string, string>();
  for (const at of all) {
    if (at.node.kind !== 'section' || !at.address) continue;
    const f = exactAddress(p(at.node, 'follow'));
    if (f) follows.set(at.address, f);
  }
  for (const at of all) {
    if (at.node.kind !== 'section' || !at.address || !follows.has(at.address)) continue;
    const seen = new Set<string>([at.address]);
    for (let a = follows.get(at.address); a !== undefined; a = follows.get(a)) {
      if (a === at.address) {
        report('seq.lang.section.follow', 'The section follows itself', at, 'follow');
        break;
      }
      if (seen.has(a)) break; // a cycle elsewhere: reported at its own sections
      seen.add(a);
    }
  }
  return issues;
}
