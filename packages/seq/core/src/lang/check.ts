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
  parseAddress,
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
import { followRef, refAddresses } from './refs.js';

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
  /**
   * What an exact address names. A layer's data may name a user-defined ref: a ref naming one exact address stands
   * for its node there (the resolved view follows it the same way); other params name nodes directly.
   */
  const target = (n: Node, key: string): Placed | undefined => {
    const a = exactAddress(p(n, key));
    if (a === undefined) return undefined;
    if (n.kind !== 'layer' || key !== 'data') return byAddress.get(a);
    const f = followRef(spec, a);
    const to = 'addresses' in f && f.addresses.length === 1 ? exactAddress(f.addresses[0]) : undefined;
    return byAddress.get(to ?? a);
  };
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

  /** The layers of template tracks: one without data draws each instance's node. */
  const templateLayers = new Set<Node>();
  for (const at of all) {
    const n = at.node;
    // every child of the spec is addressable: the description requires an id of each of its kinds, this check of an
    // extension's kind (a slot's body is no child: not walked)
    if (n.id === undefined && at.path.length && !Object.hasOwn(SeqLanguage.kinds, n.kind))
      report('seq.lang.id.required', `${n.kind} needs an id`, at);
    switch (n.kind) {
      case 'track': {
        const each = p(n, 'for_each') as Record<string, Value> | undefined;
        if (!each) break;
        for (const c of n.children ?? []) templateLayers.add(c);
        const t = byAddress.get(exactAddress(each.source) ?? '');
        if (t?.node.kind !== 'source')
          report(
            'seq.lang.ref.target',
            'The reference names no node of a kind this param takes',
            at,
            'for_each',
            'source',
          );
        break;
      }
      case 'data':
      case 'views':
      case 'refs':
        if (n.id !== n.kind) report('seq.lang.container.id', `The id of ${n.kind} is \`${n.kind}\``, at);
        break;
      case 'ref': {
        // what each address of a ref names is judged by its leading ids (forms after them are a query's, and a node
        // data makes waits for its data); a chain of refs ends (_spec/README.md, Language, Refs)
        const to = refAddresses(n);
        const named = (a: string) => {
          const segments = parseAddress(a).pattern?.segments ?? [];
          const ids: string[] = [];
          for (const s of segments) {
            if (typeof s !== 'string') break;
            ids.push(s);
          }
          return !ids.length || byAddress.has(printAddress(ids));
        };
        to.forEach((a, i) => {
          // through the refs it starts with (`/refs/p/child`, where `p` names a node with that child); a chain that
          // cannot be followed is the cycle check's
          const f = followRef(spec, a);
          const ok = 'refused' in f ? true : f.addresses.every(named);
          const tail = Array.isArray(p(n, 'to')) ? ['to', i] : ['to'];
          if (!ok) report('seq.lang.ref.unresolved', 'The ref names no node of the spec', at, ...tail);
        });
        if (at.address !== undefined) {
          const f = followRef(spec, at.address);
          if ('refused' in f && (f.refused === 'cycle' || f.refused === 'depth'))
            report('seq.lang.ref.cycle', 'The ref leads back to itself, or through too many refs', at, 'to');
        }
        break;
      }
      case 'source':
      case 'sequence':
      case 'ranges':
      case 'values':
      case 'map':
        checkDataNode(n, { target: (key) => target(n, key)?.node }, (code, message, ...tail) =>
          report(code, message, at, ...tail),
        );
        break;
      case 'layer': {
        const type = p(n, 'type');
        if (p(n, 'data') === undefined) {
          if (type !== 'letters' && !templateLayers.has(n)) report('seq.lang.layer.data', 'The layer has no data', at);
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
    const fields = (
      spec_ as { params?: Record<string, { type: string; kinds?: readonly string[]; in?: string }> } | undefined
    )?.params;
    if (!fields) continue;
    for (const key of Object.keys(n.params ?? {})) {
      const f = fields[key];
      if (f?.type !== 'ref' || f.in !== undefined) continue;
      if (exactAddress(n.params?.[key]) === undefined) continue; // a query form: a resolve pass's
      const t = target(n, key);
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

/** What a data node's checks read of a node it names: the node itself (a spec node, or what stands for a dynamic one). */
export interface DataNodeTargets {
  target(key: string): Node | undefined;
  /** The length of the sequence a param names when its letters are not inline (a dynamic node's target). */
  length?(key: string): number | undefined;
}

/**
 * The language's own checks of one data node (`source`, `sequence`, `ranges`, `values`, `map`), on a valid,
 * normalized node: one of two params, intervals, bounds, strands, a method's params, exons, a codon map's alphabets.
 * `targets` gives what a reference names (the extent and cyclic flag of a sequence it names are read from it);
 * `say` takes each issue with the step path under the node's params. Shared by the spec's checks and the overlay's,
 * so a dynamic node passes the rules a spec node does.
 */
export function checkDataNode(
  n: Node,
  ctx: DataNodeTargets,
  say: (code: Code, message: string, ...tail: (string | number)[]) => void,
): void {
  const p = (x: Node, key: string): Value | undefined => x.params?.[key];
  const xor = (a: boolean, b: boolean) => a !== b;
  const systemOf = (key: string) => {
    const t = ctx.target(key);
    if (t?.kind !== 'sequence') return undefined;
    const letters = p(t, 'letters');
    const length = typeof letters === 'string' ? letters.length : ctx.length?.(key);
    return { cyclic: p(t, 'cyclic') === true, length };
  };
  switch (n.kind) {
    case 'source':
      if (!xor(p(n, 'uri') !== undefined, p(n, 'value') !== undefined))
        say('seq.lang.source.content', 'A source has one of `uri` and `value`');
      if (p(n, 'value') !== undefined && p(n, 'index') !== undefined)
        say('seq.lang.source.index', 'An inline value has no index', 'index');
      break;
    case 'sequence': {
      if (!xor(p(n, 'letters') !== undefined, p(n, 'source') !== undefined))
        say('seq.lang.sequence.content', 'A sequence has one of `letters` and `source`');
      const letters = p(n, 'letters');
      if (typeof letters === 'string' && !Letters[String(p(n, 'alphabet'))]?.test(letters))
        say('seq.lang.sequence.letters', 'The letters are not of the alphabet', 'letters');
      break;
    }
    case 'ranges': {
      const rows = p(n, 'rows') as Record<string, Value> | undefined;
      if (!xor(rows !== undefined, p(n, 'source') !== undefined))
        say('seq.lang.ranges.content', 'A ranges node has one of `rows` and `source`');
      if (rows) {
        const { start, end, strand } = rows;
        if (isColumn(start) && isColumn(end)) {
          const s = start as ArrayLike<number>;
          const e = end as ArrayLike<number>;
          const sys = systemOf('target');
          let bad = -1;
          let past = -1;
          for (let i = 0; i < Math.min(s.length, e.length); i++) {
            if (bad < 0 && (s[i] < 0 || (e[i] < s[i] && !sys?.cyclic))) bad = i;
            if (past < 0 && sys?.length !== undefined && (s[i] > sys.length || e[i] > sys.length)) past = i;
          }
          if (bad >= 0) say('seq.lang.ranges.interval', 'A range is not an interval', 'rows', 'start', bad);
          if (past >= 0) say('seq.lang.ranges.bounds', 'A range is past the sequence', 'rows', 'end', past);
        }
        if (strand instanceof CategoryColumn) {
          for (let i = 0; i < strand.length; i++) {
            const v = strand.get(i);
            if (v !== null && v !== '+' && v !== '-') {
              say('seq.lang.ranges.strand', 'A strand is + or -', 'rows', 'strand', i);
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
        say('seq.lang.values.content', 'A values node has one of `values` and `source`');
      const length = systemOf('target')?.length;
      if (isColumn(values) && length !== undefined && (values as ArrayLike<number>).length !== length)
        say('seq.lang.values.length', 'There is not one value per position', 'values');
      break;
    }
    case 'map': {
      const method = p(n, 'method');
      const has = (k: string) => p(n, k) !== undefined;
      const fits =
        method === 'offset'
          ? has('by') && !has('exons') && !has('strand') && !has('table')
          : has('exons') && has('strand') && !has('by');
      if (!fits) say('seq.lang.map.params', 'The params do not fit the method');
      const exons = p(n, 'exons') as Record<string, Value> | undefined;
      if (exons && isColumn(exons.start) && isColumn(exons.end)) {
        const s = exons.start as ArrayLike<number>;
        const e = exons.end as ArrayLike<number>;
        const count = Math.min(s.length, e.length);
        const cyclic = systemOf('to')?.cyclic === true;
        // on a cyclic gene the last exon may cross the origin, ending at or before the first one starts (EC-11)
        const wraps = (i: number) => cyclic && i === count - 1 && e[i] <= s[0];
        for (let i = 0; i < count; i++) {
          if ((e[i] < s[i] && !wraps(i)) || (i > 0 && s[i] < e[i - 1])) {
            say('seq.lang.map.exons', 'The exons are not in order, or one ends before it starts', 'exons', 'start', i);
            break;
          }
        }
      }
      if (method === 'codon') {
        const from = ctx.target('from');
        const to = ctx.target('to');
        if (from && to) {
          const f = p(from, 'alphabet');
          const t = p(to, 'alphabet');
          if (f !== 'protein' || (t !== 'dna' && t !== 'rna'))
            say('seq.lang.map.alphabet', 'A codon map goes from protein to dna or rna');
        }
      }
      break;
    }
  }
}
