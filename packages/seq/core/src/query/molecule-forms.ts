/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The molecule node forms (_spec/README.md "Queries", Forms): `molecule` (by identifier, matched as a string, class or
// entity id), `chain` (a complex's units by chain, operator, author chain) and `atoms` (a molecule's or a unit's atoms
// system): Nodes, what `(section :of …)` takes. Positions are `elements`' (the author, 2026-11-20: "we should not
// really have (atom), we should have (elements ....)"): `atom`, `entity` and `component` are gone. Seq* holds pointers
// and labels, never chemistry: codes and identifiers are strings matched exactly.

import {
  type Atom,
  issueOf,
  type Node,
  parseAddress,
  printAddress,
  type Selected,
  type Value,
} from '@molstar/core/repr/tree';
import { checkLimit, checkParams, fail, param, prefixNodes, type Resolve, tick } from './form-kit.js';
import { QueryIssues } from './issues.js';

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

/** An address's segments, for a selected node (an address the model gives). */
function segments(address: string): Selected['address'] {
  return parseAddress(address).pattern?.segments ?? [];
}

/**
 * `(molecule :identifier "ccd:ATP")`, `(molecule :class ligand)`, `(molecule :entity 3)`: molecules under the prefix
 * (a `data` container's molecule children, a complex's units' molecules) whose identifier is one named (or a host
 * resolver's equivalent of one), whose class is the one named, and whose entity id in its source is one named (a
 * unit's `entity`; a spec's molecule has none).
 */
const molecule: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['identifier', 'class', 'entity']);
  const ids = words(form, 'identifier');
  const cls = param(form, 'class');
  const entities = words(form, 'entity');
  if (!ids && cls === undefined && !entities)
    fail('seq.query.form.param', 'A form of kind molecule needs :identifier, :class or :entity');
  const wanted = ids && new Set(ids.flatMap((i) => [i, ...(ctx.model.equivalents?.(i) ?? [])]));
  const out: Selected[] = [];
  const seen = new Set<string>();
  for (const m of prefixNodes(prefix, 'molecule')) {
    const owner = printAddress(m.address);
    const units = ctx.model.units?.(owner);
    const wantedEntity = entities && new Set(entities);
    const candidates = units
      ? units.filter((u) => !wantedEntity || (u.entity !== undefined && wantedEntity.has(u.entity))).map((u) => u.of)
      : wantedEntity
        ? []
        : (m.node.children ?? []).flatMap((c) =>
            c.kind === 'molecule' && c.id !== undefined ? [`${owner}/${c.id}`] : [],
          );
    for (const address of candidates) {
      tick(ctx);
      if (seen.has(address)) continue;
      const q = ctx.model.molecule?.(address);
      if (!q) {
        // an entity of a source's polymer is its sequence: a node too
        const sys = wantedEntity && !ids && cls === undefined ? ctx.model.system(address) : undefined;
        if (sys) {
          seen.add(address);
          out.push({ node: { kind: 'sequence' }, address: segments(address) });
        }
        continue;
      }
      if (wanted && !q.identifiers.some((i) => wanted.has(i))) continue;
      if (cls !== undefined && q.class !== cls) continue;
      seen.add(address);
      out.push({ node: q.node, address: segments(address) });
    }
  }
  checkLimit(ctx, out.length);
  return { type: 'nodes', nodes: out };
};

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

/** A selected node's canonical address. */
function selectedAddress(m: Selected): string {
  return printAddress(m.address);
}

/** The molecule forms, by kind. */
export const MoleculeForms = { molecule, chain, atoms } as const;
