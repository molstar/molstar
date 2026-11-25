/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Patch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, theme } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { crystal, Data, entryLigands, Lengths, ligand } from './data.js';

// 1FIN as a complex (seqx-core _spec "Resolved view", A section of a complex): CDK2 and the cyclin A construct read
// from `1FIN.seq`, a complex of the entry's units: CDK2 A and C, cyclin A B and D (two copies of the complex), each
// with its `auth` numbering (CDK2's from 1, the construct's from 173: P20248's), an ATP on each CDK2 (the ATP molecule
// with its 31 heavy atoms: `whole`, a chip, or `atoms`, one patch), and the ligands of six CC0 entries of CDK2 standing
// for AlphaFill's transplants (not taken: no CC BY 4.0 or looser licence), each a molecule placed as a unit on chain A
// with its entry's resolution as its score (a ranges row over the molecule). A transplant and its contacts come from
// its own entry (computed there against that entry's CDK2, drawn on 1FIN's chain A): its card and its arcs say so. The
// contacts (heavy atoms within 4 Å, derived at recording) are the complex's `contacts`, per pair of units: CDK2 to
// cyclin A in each copy, each ATP to its CDK2, each transplant to chain A; a `pairs` layer draws them. ATP's chip is
// colored on its own.

export const title = 'Complex: 1FIN as a complex';

/** The complex's address and the section's (what a patch names). */
export const ComplexAddress = '/data/1fin';
export const ComplexSection = '/views/main/complex';
const AtpBlock = ['views', 'main', 'complex', 'atp-a'];

/** A transplant's unit id (`t-1aq1`). */
export const transplantUnit = (entry: string) => `t-${entry.toLowerCase()}`;

/** The ATP units: each on its copy's CDK2 (the ligands' list holds 1FIN's two). */
export const AtpUnits = [
  { id: 'atp-a', chain: 'A', label: 'E' },
  { id: 'atp-c', chain: 'C', label: 'F' },
] as const;

export interface ComplexUnitsOptions {
  /** ATP as a chip (`whole`, absent) or its heavy atoms (`atoms`). */
  readonly atp?: 'whole' | 'atoms';
  readonly copies?: 'merge' | 'each';
}

/** A heavy atom's element from its name (ATP's: P, O, N, C). */
const element = (name: string) => name[0];

export function build(options: ComplexUnitsOptions = {}, d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const x = crystal(b);
  const atpData = ligand('ATP', d);
  const atp = b.data.molecule('atp', {
    name: 'ATP',
    residues: ['ATP'],
    atoms: {
      residue: atpData.atoms.map(() => 0),
      name: [...atpData.atoms],
      element: atpData.atoms.map(element),
    },
    identifiers: identifiers([['ccd', 'ATP']]),
  });
  const transplants = entryLigands(d).map((l) => ({
    l,
    unit: transplantUnit(l.entry),
    molecule: b.data.molecule(`lig-${l.entry.toLowerCase()}`, {
      // where it and its contacts come from, said on its card: another entry, not 1FIN
      name: `${l.comp} from ${l.entry} (contacts in ${l.entry}'s entry)`,
      residues: [l.comp],
      identifiers: identifiers([['ccd', l.comp]]),
    }),
  }));
  // each transplant's score: its entry's resolution, a field of its row over the molecule
  const scores = transplants.map((t) =>
    b.data.ranges(`${t.unit}-score`, {
      target: t.molecule,
      rows: {
        start: [0],
        end: [1],
        id: [t.l.entry],
        label: [`${t.l.entry}, ${t.l.resolution} Å`],
        resolution: [t.l.resolution],
      },
    }),
  );
  const chainOf = (copy: number, system: 'cdk2' | 'cyclin') =>
    d.chains.find((c) => c.copy === copy && c.system === system)!.id;
  const labels = (n: number, first: number) => Array.from({ length: n }, (_, i) => String(first + i));
  b.data.complex('1fin', { name: '1FIN' }, (c) => {
    for (const ch of d.chains)
      c.unit(ch.id, { of: ch.system === 'cdk2' ? x.cdk2 : x.cyclin, chain: ch.id, auth_chain: ch.id }, (u) =>
        u.numbering('auth', {
          method: 'labels',
          labels: ch.system === 'cdk2' ? labels(Lengths.cdk2, 1) : labels(Lengths.construct, d.cyclinFirst),
        }),
      );
    for (const a of AtpUnits) c.unit(a.id, { of: atp, chain: a.label, auth_chain: a.chain });
    for (const t of transplants) c.unit(t.unit, { of: t.molecule, chain: t.unit, auth_chain: 'A', label: t.l.comp });
    const rows: { a: string; ap: number; b: string; bp: number; label: string }[] = [];
    for (const [copy, cdk2, cyclin] of d.contacts)
      rows.push({
        a: chainOf(copy, 'cdk2'),
        ap: cdk2,
        b: chainOf(copy, 'cyclin'),
        bp: cyclin,
        label: 'contact in 1FIN',
      });
    // a transplant's contacts were computed in its own entry, against that entry's CDK2: its arcs say so
    const ligandRows = (unit: string, chain: string, l: Derived['ligands'][number]) => {
      const label = l.entry === '1FIN' ? `${l.comp} contact in 1FIN` : `${l.comp} contact in ${l.entry}'s entry`;
      for (const r of [...new Set(l.contacts.map(([res]) => res))].sort((p, q) => p - q))
        rows.push({ a: unit, ap: 0, b: chain, bp: r, label });
    };
    const atps = d.ligands.filter((l) => l.comp === 'ATP' && l.entry === '1FIN');
    for (const a of AtpUnits) ligandRows(a.id, a.chain, atps.find((l) => l.chain === a.chain)!);
    for (const t of transplants) ligandRows(t.unit, 'A', t.l);
    c.links('contacts', {
      rows: {
        a_unit: rows.map((r) => r.a),
        a_start: rows.map((r) => r.ap),
        a_end: rows.map((r) => r.ap + 1),
        b_unit: rows.map((r) => r.b),
        b_start: rows.map((r) => r.bp),
        b_end: rows.map((r) => r.bp + 1),
        type: rows.map(() => 'contact'),
        label: rows.map((r) => r.label),
      },
    });
  });
  b.views.view('main', (v) => {
    const s = v.section('complex', {
      of: ComplexAddress as never,
      ...(options.copies ? { copies: options.copies } : {}),
    });
    s.block('atp-a', { form: options.atp ?? 'whole' });
    s.track('contacts', { label: 'Contacts' }).layer('contacts', {
      type: 'pairs',
      data: `${ComplexAddress}/contacts` as never,
    });
    s.track('scores', { label: 'Resolution (Å)' }, (t) =>
      scores.forEach((data, k) => t.layer(`s${k}`, { type: 'ranges', data })),
    );
    // ATP's chip colored on its own: a rule over the molecule's component
    s.track('letters', { label: 'Sequence' }).layer('letters', {
      type: 'letters',
      color: theme.rules([theme.rule([theme.eq('component', 'ATP')], '#e15759', 'ATP')]),
    });
  });
  return b.build();
}

/** ATP by its heavy atoms. */
export const atoms = () => build({ atp: 'atoms' });

/** The editing view's toggle: the first ATP as a chip or by its atoms, one patch. */
export function atpFormPatch(form: 'whole' | 'atoms'): Patch {
  return { ops: [{ op: 'set', at: AtpBlock, key: 'form', value: form }] };
}
