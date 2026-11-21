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
import { checkSpecBounds, DefaultSpecBounds, IssueCaps, type SpecBounds } from './bounds.js';
import {
  checkBlockForm,
  checkCompound,
  checkLinks,
  checkModified,
  checkMolecule,
  checkUnit,
  inlineLength,
  type LinkParent,
  unitLength,
} from './check-molecules.js';
import { checkMapSegments, checkNumberingSegments, mapParamsFit, numberingParamsFit } from './check-maps.js';
import { LangIssues } from './issues.js';
import { isExtensionName, SectionOfKinds, SeqLanguage, SeqLanguageVersion, SeqValidateOptions } from './language.js';
import { followRef, refAddresses } from './refs.js';
import { checkLayerThemes, rowTypes } from './theme-check.js';

export type Code = keyof typeof LangIssues.codes;

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
 * Before them all, the bounds on what the spec sizes (`checkSpecBounds`): a spec past one is refused with its issue
 * alone, nothing else read (_spec/README.md, Language, Bounds).
 */
export function checkSeqSpec(spec: Node, bounds: SpecBounds = DefaultSpecBounds): Issue[] {
  const bounded = checkSpecBounds(spec, bounds);
  if (bounded.length) return bounded;
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
  const read = validateAndNormalize(SeqLanguage, spec, SeqValidateOptions);
  if (hasErrors(read.issues) || !read.node) return [...read.issues];
  return [...read.issues, ...checkSeqLang(read.node, bounds)];
}

/**
 * The `seq.lang.*` checks alone, for a normalized tree that is valid against `SeqLanguage` and whose version is
 * readable. The issues kept are bounded (`bounds.maxPlaces` a code, `maxNodeIssues` a node, the rest counted in one
 * `seq.lang.bound.issues` a code).
 */
export function checkSeqLang(spec: Node, bounds: SpecBounds = DefaultSpecBounds): Issue[] {
  const issues: Issue[] = [];
  const caps = new IssueCaps(bounds);
  const all: Placed[] = [];
  const byAddress = new Map<string, Placed>();
  const visit = (node: Node, path: TreePath, address: readonly string[] | undefined, parent: Node) => {
    const here = address && node.id !== undefined ? [...address, node.id] : undefined;
    const placed: Placed = { node, path, address: here ? printAddress(here) : undefined, parent };
    all.push(placed);
    if (placed.address) byAddress.set(placed.address, placed);
    node.children?.forEach((c, i) => visit(c, [...path, 'children', i], here, node));
  };
  // the root is `/`: its children are addressed by their ids alone
  all.push({ node: spec, path: [], address: '/' });
  spec.children?.forEach((c, i) => visit(c, ['children', i], [], spec));
  const look: SeqLookup = { spec, get: (a) => byAddress.get(a), themeDepth: bounds.themeDepth };
  const report = reporter(issues, caps);
  for (const at of all) checkPlaced(at, look, report);
  for (const at of all) checkFollow(at, look, report);
  issues.push(...caps.summary());
  return issues;
}

/** A node with its place: its step path from the root, its address when it has one, and its parent. */
export interface Placed {
  readonly node: Node;
  readonly path: TreePath;
  readonly address: string | undefined;
  readonly parent?: Node;
}

/** What the checks read of the spec: the root (refs are followed over it) and a node by its exact address. */
export interface SeqLookup {
  readonly spec: Node;
  get(address: string): Placed | undefined;
  /** The theme forms' nesting bound (`SpecBounds.themeDepth`; the default when absent). */
  readonly themeDepth?: number;
}

/** Takes one issue of a node: `tail` is where in the node, a param and below it a column and a row index. */
export type Report = (code: Code, message: string, at: Placed, ...tail: (string | number)[]) => void;

/** A `Report` that pushes each issue, with the node's step path and address, onto `issues`. */
export function reporter(issues: Issue[], caps?: IssueCaps): Report {
  return (code, message, at, ...tail) => {
    if (caps && !caps.keep(code, at.node)) return;
    const path = tail.length ? [...at.path, 'params', ...tail] : at.path;
    issues.push(
      issueOf(LangIssues, code, message, {
        hint: LangIssues.codes[code].hint,
        path,
        ...(at.address && at.address !== '/' ? { address: at.address } : {}),
      }),
    );
  };
}

const param = (n: Node, key: string): Value | undefined => n.params?.[key];

/**
 * What an exact address in a param names. A layer's data may name a user-defined ref: a ref naming one exact address
 * stands for its node there (the resolved view follows it the same way); other params name nodes directly.
 */
function targetOf(n: Node, key: string, look: SeqLookup): Placed | undefined {
  const a = exactAddress(param(n, key));
  if (a === undefined) return undefined;
  if (n.kind !== 'layer' || key !== 'data') return look.get(a);
  const f = followRef(look.spec, a);
  const to = 'addresses' in f && f.addresses.length === 1 ? exactAddress(f.addresses[0]) : undefined;
  return look.get(to ?? a);
}

type FieldTypes = Readonly<
  Record<string, { readonly type: string; readonly kinds?: readonly string[]; readonly in?: string }>
>;

/** A kind's params as the description types them (undefined for an extension's kind). */
export function fieldsOf(kind: string): FieldTypes | undefined {
  return (SeqLanguage.kinds as Record<string, { params?: FieldTypes } | undefined>)[kind]?.params;
}

/** The leading ids of an address (its node part); none when it does not parse. */
export function idRunOf(address: string): string[] {
  const segments = parseAddress(address).pattern?.segments ?? [];
  const ids: string[] = [];
  for (const s of segments) {
    if (typeof s !== 'string') break;
    ids.push(s);
  }
  return ids;
}

/** The language's checks of one node, follow chains aside (`checkFollow`): what it alone, and what it names, decide. */
export function checkPlaced(at: Placed, look: SeqLookup, report: Report): void {
  const n = at.node;
  const p = param;
  const target = (x: Node, key: string) => targetOf(x, key, look);
  const say = (code: Code, message: string, ...tail: (string | number)[]) => report(code, message, at, ...tail);
  /** The positions of a unit a param names (its range's, else its molecule's inline length). */
  const length = (key: string) => {
    const t = target(n, key)?.node;
    return t?.kind === 'unit' ? unitLength(t, inlineLength(target(t, 'of')?.node)) : undefined;
  };
  /** The sequence a param names, with what is known of its extent. */
  /** The system a param names (a sequence, a molecule's residues, a unit's own positions), with its extent. */
  const systemOf = (x: Node, key: string) => {
    const t = target(x, key)?.node;
    if (t?.kind === 'sequence' || t?.kind === 'molecule')
      return { cyclic: p(t, 'cyclic') === true, length: inlineLength(t) };
    if (t?.kind === 'unit') return { cyclic: false, length: unitLength(t, inlineLength(target(t, 'of')?.node)) };
    return undefined;
  };
  // every child of the spec is addressable: the description requires an id of each of its kinds, this check of an
  // extension's kind (a slot's body is no child: not walked)
  if (n.id === undefined && at.path.length && !Object.hasOwn(SeqLanguage.kinds, n.kind))
    report('seq.lang.id.required', `${n.kind} needs an id`, at);
  switch (n.kind) {
    case 'track': {
      const each = p(n, 'for_each') as Record<string, Value> | undefined;
      if (!each) break;
      const t = look.get(exactAddress(each.source) ?? '');
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
    case 'snapshot':
      // a step is a full state or a patch, never both nor neither (_spec/README.md, Language, Story steps)
      if ((p(n, 'patch') === undefined) === (p(n, 'state') === undefined))
        report('seq.lang.snapshot.content', 'A story step has both a `patch` and a `state`, or neither', at);
      break;
    case 'navigator':
      // a navigator's form needs its data (_spec/README.md, Language, Navigators)
      if (p(n, 'content') === 'bands' && p(n, 'bands') === undefined)
        report('seq.lang.navigator.bands', 'A navigator draws `bands` and names no annotation of bands', at, 'bands');
      break;
    case 'data':
    case 'views':
    case 'refs':
    case 'snapshots':
    case 'navigators':
    case 'compounds':
      if (n.id !== n.kind) report('seq.lang.container.id', `The id of ${n.kind} is \`${n.kind}\``, at);
      break;
    case 'ref': {
      // what each address of a ref names is judged by its leading ids (forms after them are a query's, and a node
      // data makes waits for its data); a chain of refs ends (_spec/README.md, Language, Refs)
      const to = refAddresses(n);
      const named = (a: string) => {
        const ids = idRunOf(a);
        return !ids.length || look.get(printAddress(ids)) !== undefined;
      };
      to.forEach((a, i) => {
        // through the refs it starts with (`/refs/p/child`, where `p` names a node with that child); a chain that
        // cannot be followed is the cycle check's
        const f = followRef(look.spec, a);
        const ok = 'refused' in f ? true : f.addresses.every(named);
        const tail = Array.isArray(p(n, 'to')) ? ['to', i] : ['to'];
        if (!ok) report('seq.lang.ref.unresolved', 'The ref names no node of the spec', at, ...tail);
      });
      if (at.address !== undefined) {
        const f = followRef(look.spec, at.address);
        if ('refused' in f && (f.refused === 'cycle' || f.refused === 'depth'))
          report('seq.lang.ref.cycle', 'The ref leads back to itself, or through too many refs', at, 'to');
      }
      break;
    }
    case 'numbering': {
      // a unit's labels are one a position of its own (its range's, else its molecule's)
      const u = at.parent?.kind === 'unit' ? at.parent : undefined;
      const total = u ? unitLength(u, inlineLength(target(u, 'of')?.node)) : inlineLength(at.parent);
      checkNumbering(n, total, (code, message, ...tail) => report(code, message, at, ...tail));
      break;
    }
    case 'compound':
      checkCompound(n, say);
      break;
    case 'unit':
      checkUnit(n, inlineLength(target(n, 'of')?.node), say);
      break;
    case 'block': {
      const of = target(n, 'of')?.node;
      checkBlockForm(n, of?.kind === 'unit' ? target(of, 'of')?.node : of, say);
      const scale = p(n, 'scale');
      if (typeof scale === 'number' && !(scale > 0)) say('seq.lang.block.scale', 'The scale is not positive', 'scale');
      break;
    }
    case 'bridge':
      if (p(n, 'a') !== undefined && p(n, 'a') === p(n, 'b'))
        say('seq.lang.bridge.sides', 'The bridge joins a section to itself', 'b');
      break;
    case 'source':
    case 'sequence':
    case 'molecule':
    case 'ranges':
    case 'values':
    case 'map':
    case 'alignment':
    case 'pairs':
    case 'grid':
    case 'assembly':
    case 'matrix':
      checkDataNode(n, { target: (key) => target(n, key)?.node, length, parent: linkParent(at, look) }, say);
      break;
    case 'layer': {
      const type = p(n, 'type');
      // the channel slots (_spec/README.md "Themes", Forms): what the description cannot say of a slot's body
      const data = target(n, 'data')?.node;
      const fields = data?.kind === 'ranges' ? rowTypes(p(data, 'rows')) : () => undefined;
      checkLayerThemes(
        n,
        {
          isScale: (a) => look.get(a)?.node.kind === 'scale',
          fieldType: (f) => (f === 'value' && type === 'values' ? 'number' : f === 'letter' ? 'text' : fields(f)),
          depth: look.themeDepth ?? DefaultSpecBounds.themeDepth,
        },
        (code, message, ...tail) => report(code, message, at, ...tail),
      );
      // an extension's type (`<id>:<name>`): what it draws is its own, so its data is checked by the resolved view's
      // shape (`seq.view.layer.shape`) through its fallback; a reader without it draws the fallback (`none` when absent)
      if (isExtensionName(String(type))) {
        if (p(n, 'fallback') === undefined) report('seq.lang.layer.fallback', 'No fallback', at, 'type');
        break;
      }
      if (p(n, 'data') === undefined) {
        // a template's layer without data draws each instance's node
        const template = at.parent?.kind === 'track' && p(at.parent, 'for_each') !== undefined;
        if (type !== 'letters' && type !== 'ruler' && !template)
          report('seq.lang.layer.data', 'The layer has no data', at);
      } else {
        // a kind the param does not take is the reference's issue alone
        const kind = target(n, 'data')?.node.kind;
        const want =
          type === 'letters' || type === 'ruler'
            ? 'sequence'
            : type === 'values' || type === 'matrix'
              ? type
              : type === 'rows' || type === 'logo'
                ? 'alignment'
                : type === 'pairs' || type === 'grid'
                  ? type
                  : 'ranges';
        // a `ranges` layer takes `pairs` too: each end a mark at its site
        // a `pairs` layer takes a segments map: a run a pair (its method the resolved view's, `seq.view.layer.shape`)
        const fits = kind === want || (type === 'ranges' && kind === 'pairs') || (type === 'pairs' && kind === 'map');
        if (kind !== undefined && fieldsOf('layer')?.data?.kinds?.includes(kind) && !fits)
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
      // the references inside the axis's params (a record's field, a list's items) name nodes of the kinds they take
      const named = (a: Value | undefined, kinds: readonly string[], ...tail: (string | number)[]) => {
        const exact = exactAddress(a);
        if (exact === undefined) return;
        const t = look.get(exact);
        if (!t || !kinds.includes(t.node.kind))
          report('seq.lang.ref.target', 'The reference names no node of a kind this param takes', at, ...tail);
      };
      // `of`: one system, or a list of assemblies (a group each)
      const of = p(n, 'of');
      if (Array.isArray(of)) of.forEach((x, i) => named(x, ['assembly'], 'of', i));
      else named(of, SectionOfKinds, 'of');
      if (p(n, 'include') !== undefined && !Array.isArray(of) && target(n, 'of')?.node.kind !== 'assembly')
        say('seq.lang.section.include', 'The section is of no assembly: `include` is ignored', 'include');
      const direction = p(n, 'direction');
      if (direction && typeof direction === 'object' && !Array.isArray(direction))
        named((direction as Record<string, Value>).along, ['sequence', 'map'], 'direction', 'along');
      const via = p(n, 'via');
      if (Array.isArray(via)) via.forEach((x, i) => named(x, ['map'], 'via', i));
      else named(via, ['map'], 'via');
      for (const key of ['pieces', 'collapsed']) {
        const v = p(n, key);
        if (Array.isArray(v)) v.forEach((x, i) => named(x, ['ranges', 'map'], key, i));
        else if (typeof v === 'string' && v.startsWith('/')) named(v, ['ranges', 'map'], key);
      }
      break;
    }
  }
  // a reference names a node, of a kind its param takes
  const fields = fieldsOf(n.kind);
  if (!fields) return;
  for (const key of Object.keys(n.params ?? {})) {
    const f = fields[key];
    if (f?.type !== 'ref' || f.in !== undefined) continue;
    if (exactAddress(n.params?.[key]) === undefined) continue; // a query form: a resolve pass's
    const t = target(n, key);
    if (!t || (f.kinds && !f.kinds.includes(t.node.kind)))
      report('seq.lang.ref.target', 'The reference names no node of a kind this param takes', at, key);
  }
}

/** A section's `follow` chain ends: reported at a section that follows itself, through any chain. */
export function checkFollow(at: Placed, look: SeqLookup, report: Report): void {
  if (at.node.kind !== 'section' || !at.address) return;
  const next = (a: string) => {
    const s = look.get(a);
    return s?.node.kind === 'section' ? exactAddress(param(s.node, 'follow')) : undefined;
  };
  const seen = new Set<string>([at.address]);
  for (let a = next(at.address); a !== undefined; a = next(a)) {
    if (a === at.address) {
      report('seq.lang.section.follow', 'The section follows itself', at, 'follow');
      break;
    }
    if (seen.has(a)) break; // a cycle elsewhere: reported at its own sections
    seen.add(a);
  }
}

/**
 * A numbering's own checks (_spec/README.md, Language, Numberings): params that fit its method (`labels`: one of
 * `labels` and `source`, `field` with a source only; `offset`: `by`; `segments`: runs apart), an id no built-in
 * has (the forms read `index` and `ordinal` as theirs), inline labels one a position of inline letters, and a label
 * at two positions said once (it names the first: the queries' rule).
 */
function checkNumbering(
  n: Node,
  total: number | undefined,
  say: (code: Code, message: string, ...tail: (string | number)[]) => void,
): void {
  const labels = param(n, 'labels');
  if (!numberingParamsFit(n)) say('seq.lang.numbering.content', 'The numbering’s params do not fit its method');
  checkNumberingSegments(n, say);
  if (n.id === 'index' || n.id === 'ordinal') say('seq.lang.numbering.id', `\`${n.id}\` is a built-in numbering`);
  if (!(labels instanceof CategoryColumn)) return;
  if (total !== undefined && labels.length !== total)
    say('seq.lang.numbering.length', `There are ${labels.length} labels for ${total} positions`, 'labels');
  const seen = new Set<string>();
  for (let i = 0; i < labels.length; i++) {
    const l = labels.get(i);
    if (l === null) continue;
    if (seen.has(l)) {
      say('seq.lang.numbering.duplicate', `The label ${l} is given twice: it names its first position`, 'labels', i);
      break;
    }
    seen.add(l);
  }
}

/** What a data node's checks read of a node it names: the node itself (a spec node, or what stands for a dynamic one). */
export interface DataNodeTargets {
  target(key: string): Node | undefined;
  /** The length of the system a param names when its positions are not inline (a dynamic node's target, a unit). */
  length?(key: string): number | undefined;
  /** What a `pairs` node is under (`data` when absent): a complex's or a molecule's links. */
  parent?: LinkParent;
}

/** What a `pairs` node's checks read of its parent: a complex's units by id, or a molecule's residues. */
function linkParent(at: Placed, look: SeqLookup): LinkParent {
  const parent = at.parent;
  if (!parent) return { kind: 'data' };
  if (parent.kind === 'molecule') return { kind: 'molecule', length: inlineLength(parent) };
  if (parent.kind !== 'complex') return { kind: parent.kind };
  return {
    kind: 'complex',
    unit: (id) => {
      const u = parent.children?.find((c) => c.kind === 'unit' && c.id === id);
      if (!u) return undefined;
      const of = exactAddress(param(u, 'of'));
      return { length: unitLength(u, inlineLength(of === undefined ? undefined : look.get(of)?.node)) };
    },
  };
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
  // a sequence, a molecule's residues, a unit or an alignment's columns: what is known of its extent
  const systemOf = (key: string) => {
    const t = ctx.target(key);
    if (t?.kind === 'alignment') return { cyclic: false, length: alignmentWidth(t) };
    if (t?.kind !== 'sequence' && t?.kind !== 'molecule' && t?.kind !== 'unit') return undefined;
    const length = inlineLength(t) ?? ctx.length?.(key);
    return { cyclic: p(t, 'cyclic') === true, length };
  };
  switch (n.kind) {
    case 'source': {
      // a host's `ref` beside a `uri`, or alone (not portable, said once: EC-185); inline data is never a host's
      const uri = p(n, 'uri') !== undefined;
      const value = p(n, 'value') !== undefined;
      const ref = p(n, 'ref') !== undefined;
      if (uri === value ? !(ref && !uri) : ref && value)
        say('seq.lang.source.content', 'A source has one of `uri` and `value` (a host’s `ref` with no value)');
      else if (ref && !uri)
        say('seq.lang.source.host_only', 'The source is named by a host’s `ref` alone: not portable', 'ref');
      if (p(n, 'value') !== undefined && p(n, 'index') !== undefined)
        say('seq.lang.source.index', 'An inline value has no index', 'index');
      break;
    }
    case 'molecule':
      checkMolecule(n, say);
      break;
    case 'alignment': {
      if (!xor(p(n, 'rows') !== undefined, p(n, 'source') !== undefined))
        say('seq.lang.alignment.content', 'An alignment has one of `rows` and `source`');
      const rows = p(n, 'rows') as { label?: Value; letters?: Value } | undefined;
      const letters = rows?.letters;
      if (Array.isArray(letters) && letters.length) {
        const w = String(letters[0]).length;
        letters.forEach((l, i) => {
          if (typeof l === 'string' && l.length !== w)
            say('seq.lang.alignment.width', 'A row is not as wide as the first', 'rows', 'letters', i);
        });
        if (Array.isArray(rows?.label) && rows.label.length !== letters.length)
          say('seq.lang.alignment.width', 'There is not one label a row', 'rows', 'label');
      }
      break;
    }
    case 'sequence': {
      // under an assembly: an override of a derived sequence, its letters and alphabet the assembly's when absent
      const derived = ctx.parent?.kind === 'assembly';
      if (!derived && !xor(p(n, 'letters') !== undefined, p(n, 'source') !== undefined))
        say('seq.lang.sequence.content', 'A sequence has one of `letters` and `source`');
      if (!derived && p(n, 'alphabet') === undefined)
        say('seq.lang.sequence.alphabet', 'The sequence has no alphabet', 'alphabet');
      const letters = p(n, 'letters');
      if (typeof letters === 'string' && !Letters[String(p(n, 'alphabet'))]?.test(letters))
        say('seq.lang.sequence.letters', 'The letters are not of the alphabet', 'letters');
      checkModified(n, say);
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
    case 'pairs': {
      checkLinks(n, ctx.parent ?? { kind: 'data' }, say);
      const rows = p(n, 'rows') as Record<string, Value> | undefined;
      if (!xor(rows !== undefined, p(n, 'source') !== undefined) || (p(n, 'field') !== undefined && !p(n, 'source')))
        say('seq.lang.pairs.content', 'A pairs node has one of `rows` and `source`, `field` with a source');
      if (!rows) break;
      // each end an interval within its system, never wrapping (EC-10)
      for (const side of ['a', 'b'] as const) {
        const s = rows[`${side}_start`];
        const e = rows[`${side}_end`];
        if (!isColumn(s) || !isColumn(e)) continue;
        const sys = systemOf(side === 'b' && p(n, 'b') !== undefined ? 'b' : 'a');
        const ss = s as ArrayLike<number>;
        const ee = e as ArrayLike<number>;
        for (let i = 0; i < Math.min(ss.length, ee.length); i++) {
          if (ss[i] < 0 || ee[i] <= ss[i]) {
            say('seq.lang.pairs.interval', 'A pair’s end is not an interval', 'rows', `${side}_start`, i);
            break;
          }
          if (sys?.length !== undefined && ee[i] > sys.length) {
            say('seq.lang.pairs.bounds', 'A pair’s end is past its sequence', 'rows', `${side}_end`, i);
            break;
          }
        }
      }
      break;
    }
    case 'assembly':
      if (
        (p(n, 'sequences') === undefined && p(n, 'source') === undefined) ||
        (p(n, 'names') !== undefined && p(n, 'alias_source') === undefined)
      )
        say('seq.lang.assembly.content', 'An assembly lists no sequences, or names a naming with no alias source');
      break;
    case 'matrix': {
      const inline = Number(p(n, 'values') !== undefined) + Number(p(n, 'codes') !== undefined);
      const source = p(n, 'source') !== undefined;
      if (
        inline + Number(source) !== 1 ||
        (inline && p(n, 'states') === undefined) ||
        (p(n, 'field') !== undefined && !source)
      )
        say('seq.lang.matrix.content', 'A matrix has one of `values`, `codes` and `source`, inline cells their states');
      break;
    }
    case 'grid': {
      const values = p(n, 'values');
      if (!xor(values !== undefined, p(n, 'source') !== undefined) || (p(n, 'field') !== undefined && !p(n, 'source')))
        say('seq.lang.grid.content', 'A grid has one of `values` and `source`, `field` with a source');
      const x = systemOf('x')?.length;
      const y = p(n, 'y') === undefined ? x : systemOf('y')?.length;
      if (isColumn(values) && x !== undefined && y !== undefined && (values as ArrayLike<number>).length !== x * y)
        say('seq.lang.grid.length', `There are not ${x * y} values, one a pair of positions`, 'values');
      break;
    }
    case 'values': {
      const values = p(n, 'values');
      if (!xor(values !== undefined, p(n, 'source') !== undefined) || (p(n, 'field') !== undefined && !p(n, 'source')))
        say('seq.lang.values.content', 'A values node has one of `values` and `source`, `field` with a source');
      const length = systemOf('target')?.length;
      if (isColumn(values) && length !== undefined && (values as ArrayLike<number>).length !== length)
        say('seq.lang.values.length', 'There is not one value per position', 'values');
      break;
    }
    case 'map': {
      const method = p(n, 'method');
      if (!mapParamsFit(n)) say('seq.lang.map.params', 'The params do not fit the method');
      checkMapSegments(n, say);
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

/** An inline alignment's columns: its first row's width; undefined for a source's (known once it opens). */
function alignmentWidth(n: Node): number | undefined {
  const letters = (n.params?.rows as { letters?: Value } | undefined)?.letters;
  return Array.isArray(letters) && typeof letters[0] === 'string' ? letters[0].length : undefined;
}
