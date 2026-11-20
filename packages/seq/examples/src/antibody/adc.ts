/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { Data } from './data.js';

// The ADC construct (molecules.md, Complexes and units: a construct is a complex with no structure behind it; Tensions,
// stochastic conjugation): trastuzumab emtansine (T-DM1), generated, not recorded. Its antibody stands as trastuzumab's
// Fab chains (PDB 1N8Z's light and heavy chains, the Fc not held), two of each; its payload, emtansine (DM1 through an
// SMCC linker), a molecule named by its components. Lysine conjugation is stochastic: the payload linked to every
// lysine of every chain, each link's occupancy a share of the average drug-to-antibody ratio (3.5) by a generated
// weight, so the occupancies sum to the DAR; the arcs colored by their occupancy (a number of the links' rows). The
// payload's CCD component ids are placeholders (`UNL`, the CCD's unknown ligand) until the author names them: no id is
// guessed.

export const title = 'T-DM1: trastuzumab emtansine (generated)';

/** The average drug-to-antibody ratio. */
export const Dar = 3.5;

/**
 * The payload's components as data, a placeholder each (`UNL`): the SMCC linker's (MCC once attached) and DM1's CCD
 * ids are not named here until checked.
 */
export const PayloadComponents = [
  { code: 'UNL', name: 'MCC linker (from SMCC), CCD id to be named' },
  { code: 'UNL', name: 'DM1 (maytansinoid), CCD id to be named' },
] as const;

/** The construct's chains: two light and two heavy, each a unit of the Fab's chain. */
export const AdcChains = [
  { id: 'L1', chain: 'light' },
  { id: 'H1', chain: 'heavy' },
  { id: 'L2', chain: 'light' },
  { id: 'H2', chain: 'heavy' },
] as const;

/** The lysines of a chain's letters (0-based). */
export const lysines = (letters: string) => [...letters].flatMap((c, i) => (c === 'K' ? [i] : []));

/**
 * Each conjugation site (a lysine of a unit) and its occupancy: a share of the DAR by a weight from a seed (generated:
 * some lysines more reactive than others), the same on both copies of a chain; the shares sum to the DAR.
 */
export function conjugationSites(d: Derived = Data, seed = 2013) {
  const r = new PCG(seed);
  const weight = new Map<string, number>();
  for (const c of ['light', 'heavy'] as const)
    for (const p of lysines(d[c].letters)) weight.set(`${c} ${p}`, 1 + (r.int() % 8));
  const sites = AdcChains.flatMap((c) =>
    lysines(d[c.chain].letters).map((p) => ({ unit: c.id, position: p, w: weight.get(`${c.chain} ${p}`)! })),
  );
  const total = sites.reduce((n, s) => n + s.w, 0);
  return sites.map(({ w, ...s }) => ({ ...s, occupancy: (Dar * w) / total }));
}

export function build(d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const light = b.data.sequence('light', { alphabet: 'protein', letters: d.light.letters, name: 'trastuzumab light' });
  const heavy = b.data.sequence('heavy', {
    alphabet: 'protein',
    letters: d.heavy.letters,
    name: 'trastuzumab heavy (Fd)',
  });
  const payload = b.data.molecule(
    'emtansine',
    { name: 'emtansine (SMCC–DM1)', residues: PayloadComponents.map((c) => c.code) },
    (m) =>
      m.links('links', {
        rows: { a_start: [0], a_end: [1], b_start: [1], b_end: [2], type: ['covalent'], label: ['MCC–DM1'] },
      }),
  );
  const sites = conjugationSites(d);
  b.data.complex('adc', { name: 'T-DM1 (generated)' }, (c) => {
    for (const ch of AdcChains) c.unit(ch.id, { of: ch.chain === 'light' ? light : heavy, label: ch.id });
    c.unit('payload', { of: payload, label: `emtansine, DAR ${Dar}` });
    c.links('conjugation', {
      rows: {
        a_unit: sites.map((s) => s.unit),
        a_start: sites.map((s) => s.position),
        a_end: sites.map((s) => s.position + 1),
        b_unit: sites.map(() => 'payload'),
        b_start: sites.map(() => 0),
        b_end: sites.map(() => 1),
        type: sites.map(() => 'covalent'),
        a_atom: sites.map(() => 'NZ'),
        occupancy: sites.map((s) => s.occupancy),
        label: sites.map((s) => `${s.unit} K${s.position + 1}`),
      },
    });
  });
  b.views.view('main', (v) => {
    const s = v.section('adc', { of: '/data/adc' as never, copies: 'each' });
    s.track('conjugation', { label: `Lysine conjugation (DAR ${Dar})` }).layer('conjugation', {
      type: 'pairs',
      data: '/data/adc/conjugation' as never,
      // the occupancies drawn: an arc's color its occupancy, from none to the most occupied site's
      color: theme.field('occupancy', { domain: [0, null], palette: theme.continuous({ colors: 'viridis' }) }),
    });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}
