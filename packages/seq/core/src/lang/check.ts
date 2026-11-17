/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the description cannot say (language-level checks, `seq.lang.*`; _spec/README.md, Language): one of two params,
// what a method needs, what a reference names. `checkSeqSpec` is the gates in order: the version, the description's
// validation with normalization, then these checks on the normalized tree.

import { validateAndNormalize } from '@molstar/core/repr/schema';
import {
  CategoryColumn,
  exactAddress,
  hasErrors,
  type Issue,
  isColumn,
  issueOf,
  type Node,
  printAddress,
  type TreePath,
  type Value,
} from '@molstar/core/repr/tree';
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

/**
 * The issues of a spec, in gates (_spec/README.md, Language): a newer major version is `seq.lang.version.unsupported`
 * alone (the description of this version says nothing of a newer spec); else the description's (`repr.*`, the tree
 * normalized, so a spec read from JSON has typed columns); and, when it has no error, the language's own (`seq.lang.*`).
 */
export function checkSeqSpec(spec: Node): Issue[] {
  const version = spec.kind === 'spec' ? spec.params?.version : undefined;
  const major = typeof version === 'string' ? /^(\d+)\./.exec(version) : null;
  if (major && Number(major[1]) > MaxMajor) {
    return [
      issueOf(LangIssues, 'seq.lang.version.unsupported', 'The spec is written in a newer major version', {
        hint: LangIssues.codes['seq.lang.version.unsupported'].hint,
        path: ['params', 'version'],
      }),
    ];
  }
  const read = validateAndNormalize(SeqLanguage, spec);
  if (hasErrors(read.issues) || !read.node) return [...read.issues];
  return [...read.issues, ...checkSeqLang(read.node)];
}

/**
 * The `seq.lang.*` checks alone, for a normalized tree that is valid against `SeqLanguage` and whose version is
 * readable.
 */
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

  /** `tail` is where in the node: a param, and below it a column and a row index. */
  const report = (code: Code, message: string, at: Placed, ...tail: (string | number)[]) => {
    const path = tail.length ? [...at.path, 'params', ...tail] : at.path;
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
  /** The sequence a param names, with what is known of its extent. */
  const systemOf = (n: Node, key: string) => {
    const t = target(n, key)?.node;
    if (t?.kind !== 'sequence') return undefined;
    const letters = p(t, 'letters');
    return { cyclic: p(t, 'cyclic') === true, length: typeof letters === 'string' ? letters.length : undefined };
  };
  const kindsOf = (kind: string, key: string): readonly string[] | undefined => {
    const f = (SeqLanguage.kinds as Record<string, { params?: Record<string, { kinds?: readonly string[] }> }>)[kind]
      .params?.[key];
    return f?.kinds;
  };

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
        if (p(n, 'value') !== undefined && p(n, 'index') !== undefined)
          report('seq.lang.source.index', 'An inline value has no index', at, 'index');
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
          const { start, end, strand } = rows;
          if (isColumn(start) && isColumn(end)) {
            const s = start as ArrayLike<number>;
            const e = end as ArrayLike<number>;
            const sys = systemOf(n, 'target');
            let bad = -1;
            let past = -1;
            for (let i = 0; i < Math.min(s.length, e.length); i++) {
              if (bad < 0 && (s[i] < 0 || (e[i] < s[i] && !sys?.cyclic))) bad = i;
              if (past < 0 && sys?.length !== undefined && (s[i] > sys.length || e[i] > sys.length)) past = i;
            }
            if (bad >= 0) report('seq.lang.ranges.interval', 'A range is not an interval', at, 'rows', 'start', bad);
            if (past >= 0) report('seq.lang.ranges.bounds', 'A range is past the sequence', at, 'rows', 'end', past);
          }
          if (strand instanceof CategoryColumn) {
            for (let i = 0; i < strand.length; i++) {
              const v = strand.get(i);
              if (v !== null && v !== '+' && v !== '-') {
                report('seq.lang.ranges.strand', 'A strand is + or -', at, 'rows', 'strand', i);
                break;
              }
            }
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
        const exons = p(n, 'exons') as Record<string, Value> | undefined;
        if (exons && isColumn(exons.start) && isColumn(exons.end)) {
          const s = exons.start as ArrayLike<number>;
          const e = exons.end as ArrayLike<number>;
          for (let i = 0; i < Math.min(s.length, e.length); i++) {
            if (e[i] < s[i] || (i > 0 && s[i] < e[i - 1])) {
              report(
                'seq.lang.map.exons',
                'The exons are not in order, or one ends before it starts',
                at,
                'exons',
                'start',
                i,
              );
              break;
            }
          }
        }
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
        if (p(n, 'data') === undefined) {
          if (type !== 'letters') report('seq.lang.layer.data', 'The layer has no data', at);
        } else {
          // a kind the param does not take is the reference's issue alone
          const kind = target(n, 'data')?.node.kind;
          const want = type === 'letters' ? 'sequence' : type;
          if (kind !== undefined && kindsOf('layer', 'data')?.includes(kind) && kind !== want)
            report('seq.lang.layer.data', 'The data is not what the type draws', at, 'data');
        }
        break;
      }
      case 'section': {
        const range = p(n, 'range');
        const sys = systemOf(n, 'of');
        if (Array.isArray(range) && range.length === 2) {
          const [s, e] = range as number[];
          if (e < s && !sys?.cyclic) report('seq.lang.section.range', 'The range ends before it starts', at, 'range');
          else if (sys?.length !== undefined && (s > sys.length || e > sys.length))
            report('seq.lang.section.range', 'The range is past the system', at, 'range');
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
