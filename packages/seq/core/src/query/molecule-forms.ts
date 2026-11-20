/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The molecule forms (_spec/README.md "Queries", Forms): `component` (positions by their chemical component, the
// non-standard ones, a class through the compounds), `molecule` (by identifier, matched as a string, or class),
// `entity` (a source's molecule by its entity id), `chain` (a complex's units by chain, operator, author chain),
// `atom` (atoms by name) and `atoms` (a molecule's or a unit's atoms system). Seq* holds pointers and labels, never
// chemistry: codes and identifiers are strings matched exactly; a class word matches a word of the compound's class.

import {
  type Atom,
  issueOf,
  type Node,
  parseAddress,
  printAddress,
  type Selected,
  type Value,
} from '@molstar/core/repr/tree';
import {
  checkLimit,
  checkParams,
  type Ctx,
  fail,
  param,
  prefixNodes,
  type Resolve,
  required,
  stage,
  systemsOf,
  tick,
} from './form-kit.js';
import { QueryIssues } from './issues.js';
import { EmptyLocation, type Location, location, rangesLocation, unionLocations } from './location.js';
import type { QueryComponents } from './model.js';

/** A param's atoms: one, or a list's. */
function atomsOf(v: Value | undefined): readonly Atom[] | undefined {
  if (v === undefined) return undefined;
  return (Array.isArray(v) ? v : [v]) as Atom[];
}

/** Strings a param names (a list or one), or a form error. */
function words(form: Node, key: string): readonly string[] | undefined {
  const v = atomsOf(param(form, key));
  if (v && !v.every((x) => typeof x === 'string' || typeof x === 'number'))
    fail('seq.query.form.param', `A form of kind ${form.kind} takes words for :${key}`);
  return v?.map(String);
}

/** Whether a class word is a word of a class (`saccharide` of `D-saccharide, beta linking`), ignoring case. */
function classHas(cls: string | undefined, word: string): boolean {
  if (!cls) return false;
  const w = word.toLowerCase();
  return cls
    .toLowerCase()
    .split(/[^a-z0-9]+/)
    .includes(w);
}

/** The ranges of a system's positions whose components satisfy `test`, merged in order. */
function positionsWhere(ctx: Ctx, c: QueryComponents, test: (codes: readonly string[]) => boolean): number[] {
  const out: number[] = [];
  for (let p = 0; p < c.length; p++) {
    tick(ctx);
    if (!test(c.at(p))) continue;
    if (out.length && out[out.length - 1] === p) out[out.length - 1] = p + 1;
    else out.push(p, p + 1);
  }
  return out;
}

/** The non-standard positions of a sequence's (or a unit's) modified rows. */
function modifiedRanges(c: QueryComponents): number[] {
  const m = c.modified;
  const out: number[] = [];
  if (!m) return out;
  for (let i = 0; i < m.position.length; i++) out.push(m.position[i], m.position[i] + m.length[i]);
  return out;
}

/** Annotation types that say a position is modified (a database's PTM known by name only). */
const ModifiedTypes = new Set(['modified', 'glycosylation', 'lipidation']);

/** The extents of the `modified`-type items of the annotations on a system (by `type`, else `category`). */
function modifiedItems(ctx: Ctx, system: string): number[] {
  const out: number[] = [];
  for (const address of ctx.model.annotationsOn?.(system) ?? []) {
    const a = ctx.model.annotation(address);
    if (!a) continue;
    const read = a.items();
    stage(ctx, read.pending);
    for (const it of read.items) {
      tick(ctx);
      const t = it.fields?.type ?? it.fields?.category;
      if (typeof t === 'string' && ModifiedTypes.has(t.toLowerCase()) && (it.system ?? a.target) === system)
        out.push(it.start, it.end);
    }
  }
  return out;
}

/**
 * `(component :id MSE)`, `(component :standard false :annotations true)`, `(component :class saccharide)`: the
 * positions of the prefix's systems (sequences, molecules, units) whose component fits every param given.
 */
const component: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['id', 'standard', 'class', 'annotations']);
  const ids = words(form, 'id');
  const cls = param(form, 'class');
  const standard = param(form, 'standard');
  const annotations = param(form, 'annotations') ?? true;
  if (standard !== undefined && typeof standard !== 'boolean')
    fail('seq.query.form.param', 'A form of kind component takes true or false for :standard');
  if (cls !== undefined && typeof cls !== 'string')
    fail('seq.query.form.param', 'A form of kind component takes a word for :class');
  if (!ids && cls === undefined && standard === undefined)
    fail('seq.query.form.param', 'A form of kind component needs :id, :class or :standard', {
      hint: '(component :id MSE), (component :standard false)',
    });
  const wanted = ids && new Set(ids);
  let out: Location = EmptyLocation;
  for (const s of systemsOf(ctx, prefix, 'component')) {
    const c = ctx.model.components?.(s.address);
    if (!c) continue;
    let ranges: number[];
    if (standard === false && !wanted && cls === undefined) {
      // the sparse rows alone: a protein of 300 residues is not scanned for its handful
      ranges = c.polymer ? modifiedRanges(c) : [];
      if (annotations === true) ranges.push(...modifiedItems(ctx, s.address));
    } else if (c.pending) {
      // the letters are a source's and not held: the standard components are unknown, so the form waits
      stage(ctx, rangesLocation(s.address, [0, c.length]));
      ranges = [];
    } else {
      const nonStandard =
        c.polymer && standard !== undefined ? new Set(rangesToPositions(modifiedRanges(c))) : undefined;
      ranges = positionsWhere(ctx, c, (codes) => {
        if (wanted && !codes.some((x) => wanted.has(x))) return false;
        if (typeof cls === 'string' && !codes.some((x) => classHas(ctx.model.compoundClass?.(x), cls))) return false;
        return true;
      });
      if (standard !== undefined) {
        // a molecule's residues have no standard component: `:standard false` holds them all, `true` none
        const keep = (p: number) => (c.polymer ? nonStandard!.has(p) !== standard : !standard);
        ranges = filterPositions(ranges, keep);
      }
      if (standard === false && annotations === true && !wanted && cls === undefined)
        ranges.push(...modifiedItems(ctx, s.address));
    }
    out = unionLocations(out, location({ parts: [{ system: s.address, ranges }] }));
  }
  return { type: 'location', location: out };
};

function rangesToPositions(ranges: readonly number[]): number[] {
  const out: number[] = [];
  for (let i = 0; i + 1 < ranges.length; i += 2) for (let p = ranges[i]; p < ranges[i + 1]; p++) out.push(p);
  return out;
}

function filterPositions(ranges: readonly number[], keep: (p: number) => boolean): number[] {
  const out: number[] = [];
  for (const p of rangesToPositions(ranges)) {
    if (!keep(p)) continue;
    if (out.length && out[out.length - 1] === p) out[out.length - 1] = p + 1;
    else out.push(p, p + 1);
  }
  return out;
}

/** An address's segments, for a selected node (an address the model gives). */
function segments(address: string): Selected['address'] {
  return parseAddress(address).pattern?.segments ?? [];
}

/**
 * `(molecule :identifier "ccd:ATP")`, `(molecule :class ligand)`: molecules under the prefix (a `data` container's
 * molecule children, a complex's units' molecules) whose identifier is one named (or a host resolver's equivalent of
 * one), and whose class is the one named.
 */
const molecule: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['identifier', 'class']);
  const ids = words(form, 'identifier');
  const cls = param(form, 'class');
  if (!ids && cls === undefined) fail('seq.query.form.param', 'A form of kind molecule needs :identifier or :class');
  const wanted = ids && new Set(ids.flatMap((i) => [i, ...(ctx.model.equivalents?.(i) ?? [])]));
  const out: Selected[] = [];
  const seen = new Set<string>();
  for (const m of prefixNodes(prefix, 'molecule')) {
    const owner = printAddress(m.address);
    const units = ctx.model.units?.(owner);
    const candidates = units
      ? units.map((u) => u.of)
      : (m.node.children ?? []).flatMap((c) =>
          c.kind === 'molecule' && c.id !== undefined ? [`${owner}/${c.id}`] : [],
        );
    for (const address of candidates) {
      tick(ctx);
      if (seen.has(address)) continue;
      const q = ctx.model.molecule?.(address);
      if (!q) continue;
      if (wanted && !q.identifiers.some((i) => wanted.has(i))) continue;
      if (cls !== undefined && q.class !== cls) continue;
      seen.add(address);
      out.push({ node: q.node, address: segments(address) });
    }
  }
  checkLimit(ctx, out.length);
  return { type: 'nodes', nodes: out };
};

/** `(entity :id 3)`: the molecule a complex's source calls entity 3 (what its units' `of` name by that key). */
const entity: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['id']);
  const id = String(required(form, 'id'));
  const out: Selected[] = [];
  for (const m of prefixNodes(prefix, 'entity')) {
    const owner = selectedAddress(m);
    const of = ctx.model.units?.(owner)?.find((u) => keyOf(u.of) === id)?.of;
    const q = of === undefined ? undefined : (ctx.model.molecule?.(of) ?? ctx.model.system(of));
    if (of === undefined || !q) {
      ctx.issues.push(missing(owner));
      continue;
    }
    out.push({ node: 'node' in q ? q.node : { kind: 'sequence' }, address: segments(of) });
  }
  checkLimit(ctx, out.length);
  return { type: 'nodes', nodes: out };
};

/** A dynamic node's key: the last segment's `(node :key k)`. */
function keyOf(address: string): string | undefined {
  const s = segments(address);
  const last = s[s.length - 1];
  return last && typeof last !== 'string' && last.kind === 'node' ? String(last.params?.key) : undefined;
}

/**
 * `(chain :id A :operator 1_555)`, `(chain :auth B)`: a complex's units by their structure's chain, operator and
 * author chain, each addressed by its id (one system, one address); over a prefix that is no complex, a chain the
 * model makes, addressed by the form.
 */
const chain: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['id', 'operator', 'auth']);
  const id = param(form, 'id');
  const op = param(form, 'operator');
  const auth = param(form, 'auth');
  for (const v of [id, op, auth])
    if (v !== undefined && (v === null || typeof v === 'object'))
      fail('seq.query.form.param', 'A form of kind chain takes an atom');
  if (id === undefined && auth === undefined) fail('seq.query.form.param', 'A form of kind chain needs :id or :auth');
  const out: Selected[] = [];
  for (const m of prefixNodes(prefix, 'chain')) {
    tick(ctx);
    const owner = selectedAddress(m);
    const units = ctx.model.units?.(owner);
    if (!units) {
      const n =
        id === undefined || op !== undefined || auth !== undefined ? undefined : ctx.model.chain?.(owner, id as Atom);
      if (n) out.push({ node: n, address: [...m.address, form] });
      else ctx.issues.push(missing(owner));
      continue;
    }
    const found = units.filter(
      (u) =>
        (id === undefined || u.chain === String(id)) &&
        (op === undefined || u.operator === String(op)) &&
        (auth === undefined || u.authChain === String(auth)),
    );
    if (!found.length) ctx.issues.push(missing(owner));
    for (const u of found) out.push({ node: u.node, address: [...m.address, u.id] });
  }
  checkLimit(ctx, out.length);
  return { type: 'nodes', nodes: out };
};

function missing(address: string) {
  return issueOf(QueryIssues, 'seq.query.node.missing', 'No node of this key', { address });
}

/** `(atoms)`: a molecule's or a unit's atoms system, addressed `<system>/(atoms)`. */
const atoms: Resolve = (ctx, form, prefix) => {
  checkParams(form, []);
  const out: Selected[] = [];
  for (const m of prefixNodes(prefix, 'atoms')) {
    const owner = selectedAddress(m);
    const a = ctx.model.atoms?.(owner);
    if (!a || a.address === owner) {
      ctx.issues.push(missing(owner));
      continue;
    }
    out.push({ node: m.node, address: [...m.address, form] });
  }
  checkLimit(ctx, out.length);
  return { type: 'nodes', nodes: out };
};

/** `(atom :name O1G)`: atoms by name on the prefix's atoms systems (a molecule's, a unit's, or one itself). */
const atom: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['name']);
  const names = new Set((atomsOf(required(form, 'name')) ?? []).map(String));
  const parts: { system: string; ranges: number[] }[] = [];
  for (const m of prefixNodes(prefix, 'atom')) {
    const owner = selectedAddress(m);
    const a = ctx.model.atoms?.(owner);
    if (!a)
      return fail('seq.query.form.prefix', 'A form of kind atom applies to a system with atoms', { address: owner });
    const ranges: number[] = [];
    for (let i = 0; i < a.length; i++) {
      tick(ctx);
      const n = a.names.get(i);
      if (n !== null && names.has(n)) ranges.push(i, i + 1);
    }
    parts.push({ system: a.address, ranges });
  }
  return { type: 'location', location: location({ parts }) };
};

/** A selected node's canonical address. */
function selectedAddress(m: Selected): string {
  return printAddress(m.address);
}

/** The molecule forms, by kind. */
export const MoleculeForms = { component, molecule, entity, chain, atom, atoms } as const;
