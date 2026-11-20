/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure query tests' helpers: a spec's or a source's nodes as a query model (the overlay over the static
// model, as a viewer's), an address evaluated over it, and a form built in code evaluated over a prefix (a predicate's
// positional operands are held in `args`, which the address grammar reads once `@molstar/core` does: until then a
// test builds the form). And a small fixture recorded from 3SGJ: the IgG1 Fc's chain A around Asn297, its first
// N-glycan (entity 3, chain D: eight residues, each with its C1 and a few more atoms) and the glycosidic link between
// them, with the glycan's own links.

import type { Atom, Node } from '@molstar/core/repr/tree';
import { type NodeEntry, nodeKey } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import {
  type Evaluation,
  evaluateQuery,
  partOf,
  type QueryExtensions,
  queryResolvers,
  fromNamed,
  type QueryResult,
  toNamed,
} from '@molstar/seqx-core/query';
import { buildOverlay, overlayQueryModel, type RowsSource, resolveSpec } from '@molstar/seqx-core/resolve';

export interface Run {
  readonly spec: Node;
  readonly query: ReturnType<typeof overlayQueryModel>;
  run(address: string): Evaluation;
  /** A form built in code over what an address names (`/` alone: the root). */
  form(prefix: string, form: Node): Evaluation;
}

function runner(spec: Node, query: Run['query'], extensions?: QueryExtensions): Run {
  const run = (address: string) => evaluateQuery(spec, query, address, { limit: 100000, extensions });
  return {
    spec,
    query,
    run,
    form(prefix, form) {
      const p = run(prefix);
      if (!p.result) return p;
      const table = queryResolvers(query, undefined, extensions);
      const r = table[form.kind]({ root: spec, form, prefix: toNamed(p.result), limit: 100000 } as never);
      const result = r.result && fromNamed(r.result);
      return {
        ...(result && !r.issues?.some((i) => i.severity === 'error') ? { result } : {}),
        issues: r.issues ?? [],
        complete: r.complete !== false,
        pending: (r.pending as Evaluation['pending']) ?? p.pending,
      };
    },
  };
}

/** A spec's query model (its overlay empty). */
export function specRun(spec: Node, extensions?: QueryExtensions, src?: RowsSource): Run {
  const model = resolveSpec(spec);
  return runner(spec, overlayQueryModel(model, buildOverlay(model, reader([])), src), extensions);
}

function reader(nodes: readonly { key: string; node: Node }[], source = '') {
  const map = new Map<string, NodeEntry>();
  for (const n of nodes) map.set(nodeKey(source, n.key), { source, key: n.key, node: n.node, version: 1 });
  return { nodes: map, nodeVersion: { value: 1 } };
}

/** A source's nodes (a structure read by a runtime) under `/data/<id>`, as the viewer's overlay holds them. */
export function sourceRun(
  id: string,
  nodes: readonly { key: string; node: Node }[],
  extensions?: QueryExtensions,
): Run {
  const b = createSeqSpec();
  b.data.source(id, { format: 'mmcif', uri: `${id}.cif` });
  const spec = b.build();
  const model = resolveSpec(spec);
  const overlay = buildOverlay(model, reader(nodes, `/data/${id}`));
  return runner(spec, overlayQueryModel(model, overlay), extensions);
}

/** A Location result's ranges on a system, else the issues' codes. */
export const loc = (r: Evaluation, system: string) =>
  r.result?.type === 'location'
    ? Array.from(partOf(r.result.location, system)?.ranges ?? [])
    : r.issues.map((i) => i.code);

/** How many positions a Location result holds, else the issues' codes. */
export const count = (r: Evaluation) =>
  r.result?.type === 'location'
    ? r.result.location.parts.reduce((n, p) => {
        for (let k = 0; k < p.ranges.length; k += 2) n += p.ranges[k + 1] - p.ranges[k];
        return n;
      }, 0)
    : r.issues.map((i) => i.code);

export const codes = (r: { readonly issues: readonly { readonly code: string }[] }) => r.issues.map((i) => i.code);

export const keysOf = (r: Evaluation) =>
  r.result?.type === 'keys' ? r.result.keys.keys.map((k) => `${k.owner}#${k.key}`) : r.issues.map((i) => i.code);

/** An expression form: its kind and positional operands. */
export const ex = (kind: string, ...args: (Atom | Node | readonly Atom[])[]): Node => ({ kind, args: args as never });
export const prop = (name: string) => ex('prop', name);

export type { QueryResult };

/** 3SGJ's chain A around Asn297 (auth 293–301) and its first N-glycan (chain D), as the runtime reads them. */
export function glycanSpec(): Node {
  const b = createSeqSpec({ title: '3SGJ glycan' });
  const fc = b.data.sequence('fc', { alphabet: 'protein', letters: 'EEQYNSTYR' });
  const residues = ['NAG', 'NAG', 'BMA', 'MAN', 'NAG', 'MAN', 'NAG', 'FUC'];
  const atomNames = residues.flatMap((r) => (r === 'FUC' ? ['C1', 'C2', 'C6'] : ['C1', 'C2', 'O5']));
  const glycan = b.data.molecule(
    'glycan',
    {
      class: 'branched',
      residues,
      atoms: {
        residue: residues.flatMap((_, i) => [i, i, i]),
        name: atomNames,
        element: atomNames.map((n) => n[0]),
      },
    },
    (m) =>
      m.links('bonds', {
        rows: {
          a_start: [1, 2, 3, 4, 5, 6, 7],
          a_end: [2, 3, 4, 5, 6, 7, 8],
          b_start: [0, 1, 2, 3, 2, 5, 0],
          b_end: [1, 2, 3, 4, 3, 6, 1],
          type: Array(7).fill('glycosidic'),
        },
      }),
  );
  b.data.complex('c', { name: '3SGJ Fc and its glycan' }, (c) => {
    c.unit('A', { of: fc, chain: 'A', auth_chain: 'A', operator: '1_555' }, (u) =>
      u.numbering('auth', {
        method: 'labels',
        scheme: 'auth',
        labels: Array.from({ length: 9 }, (_, i) => `${i + 293}`),
      }),
    );
    c.unit('D', { of: glycan, chain: 'D', auth_chain: 'A', operator: '1_555' });
    c.links('links', {
      rows: { a_unit: ['A'], a_start: [4], a_end: [5], b_unit: ['D'], b_start: [0], b_end: [1], type: ['glycosidic'] },
    });
  });
  return b.build();
}
