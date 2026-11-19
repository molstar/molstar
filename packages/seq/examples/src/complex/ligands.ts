/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printAddress } from '@molstar/core/repr/tree';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { crystal, Data, entryLigands, ligand, ligandContactRows } from './data.js';

// The ligands as blocks after the polymers (EC-64, EC-101): CDK2 (1FIN's chain A) and each ligand a block of its own,
// 1FIN's ATP first, then the ligands of six CC0 entries of CDK2, best resolution first (AlphaFill's transplants are not
// taken: no CC BY 4.0 or looser licence is stated). A ligand is one position (a chip) or, for ATP, its heavy atoms
// (`atp: 'atoms'`, a block of 31 positions: a sequence of its own with no map to the one-position ATP, this slice). The
// contacts as arcs from each ligand to its contacting residues (≤ 4 Å); each ligand's "contacts" a ref the builder
// writes, run from the bookmarks list; its entry and resolution on hover (the contacts' `entry` and `resolution`).

export const title = 'Complex: ligands of CDK2';

export interface LigandsOptions {
  /** ATP as one position (`one`, a chip) or per heavy atom (`atoms`). */
  readonly atp?: 'one' | 'atoms';
  /** Ligands from entries worse than this resolution (Å) hidden: a filter in the contacts layer's data. */
  readonly bestResolution?: number;
}

/** A ligand's node id (`lig-stu`). */
export const ligandId = (comp: string) => `lig-${comp.toLowerCase()}`;

export function build(options: LigandsOptions = {}, d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const x = crystal(b, d);
  const atoms = options.atp === 'atoms';
  const atp = ligand('ATP', d);
  const all = [atp, ...entryLigands(d)];
  const contacts: string[] = [];
  const nodes = all.map((l) => {
    const one = !(atoms && l === atp);
    const seq = b.data.sequence(ligandId(l.comp) + (one ? '' : '-atoms'), {
      alphabet: 'protein',
      letters: one ? 'X' : 'X'.repeat(l.atoms.length),
      name: l.comp,
    });
    const rows = ligandContactRows(l, one ? 'ligand' : 'atom');
    const pairs = b.data.pairs(`${ligandId(l.comp)}-contacts`, {
      a: seq,
      b: x.cdk2,
      rows: {
        ...rows,
        label: rows.id.map(() => `${l.comp} (${l.entry}, ${l.resolution.toFixed(2)} Å)`),
      },
    });
    contacts.push(pairs);
    return { l, seq, pairs, one };
  });
  // the atoms' names on the per-atom block: a range an atom
  const names = atoms
    ? b.data.ranges('atp-atoms', {
        target: nodes[0].seq,
        rows: {
          start: atp.atoms.map((_, i) => i),
          end: atp.atoms.map((_, i) => i + 1),
          id: [...atp.atoms],
          label: [...atp.atoms],
        },
      })
    : undefined;
  // hidden by resolution: the ligands past the bound, each layer's data filtered out whole
  const hidden = new Set(
    options.bestResolution === undefined
      ? []
      : all.filter((l) => l.resolution > options.bestResolution!).map((l) => l.comp),
  );
  b.views.view('main', (v) => {
    const s = v.section('ligands', {});
    s.block('cdk2', { of: x.cdk2, label: 'CDK2 (1FIN A)' });
    for (const n of nodes) s.block(ligandId(n.l.comp), { of: n.seq, label: n.l.comp });
    // padded on top so the arcs clear the blocks' labels above the track
    s.track('contacts', { label: 'Contacts (≤ 4 Å)', padding: { top: 8 } }, (t) => {
      // a ligand's arcs in its own color of the categorical palette (by its place: the chips' order), told apart
      nodes.forEach((n, k) =>
        t.layer(ligandId(n.l.comp), {
          type: 'pairs',
          data: hidden.has(n.l.comp) ? hiddenData(n.pairs) : n.pairs,
          color: ligandColor(k),
        }),
      );
    });
    if (names) s.track('atoms', { label: 'ATP atoms' }).layer('names', { type: 'ranges', data: names });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
  });
  // each ligand's contacts: a ref the builder writes (run from the bookmarks list: its residues highlighted)
  for (const n of nodes)
    b.refs.ref(`${ligandId(n.l.comp)}-contacts`, { to: n.pairs, label: `${n.l.comp} contacts (${n.l.entry})` });
  return b.build();
}

/** The color of the `k`th ligand's arcs: the categorical palette's `k`th token (its eight in turn). */
export const ligandColor = (k: number) => `token:seq-cat-${k % 8}`;

/** A ligand's contacts hidden: every row filtered out (`(not (item :type contact))`). */
export const hiddenData = (pairs: string) =>
  `${pairs}${printAddress([{ kind: 'not', params: { of: { kind: 'item', params: { type: 'contact' } } } }])}`;

/** ATP per heavy atom. */
export const atoms = () => build({ atp: 'atoms' });
