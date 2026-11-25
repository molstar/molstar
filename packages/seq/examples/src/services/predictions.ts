/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { P53, seeded } from '../p53/data.js';
import { compose, type ComposedOptions, structureSection } from './composed.js';

// A lab's predictions onto a structure: a score per position of a sequence the lab indexes from 1
// and gives no numbering of (a made-up "tolerance" score, SYNTHETIC: generated from a seed, high at the hotspots R175,
// R248 and R273 and over the DNA-binding domain, never a real model's output), placed on UniProt's P04637 by an offset
// map of 0 and from there on PDB 2OCJ's chains by SIFTS's maps (`composed.ts`): the two maps chain, so a position of the
// lab's lights the residue of each chain that SIFTS pairs with it (2OCJ's author numbers are UniProt's, its entity starts
// at 94; its first two and last 23 residues have no atoms: the scores there light nothing in the structure). The track is
// colored in the view and, by the host's "Color the structure by this track", in the fake structure viewer.

export const title = 'Predictions onto a structure (synthetic scores)';

/** The lab's sequence: P04637's letters, its own positions from 1. */
export const LabAddress = '/data/lab';
export const LabScoresAddress = '/data/lab-scores';

/** The bins of the score and their colors: low (blue-gray) to high (dark red). */
export const ScoreEdges = [0, 0.2, 0.4, 0.6, 0.8, 1.0001] as const;
export const ScoreColors = ['#9ecae1', '#c7e9c0', '#fee391', '#fdae6b', '#d94801'] as const;

/** The scores, a position each, in [0, 1]: noise, the DNA-binding domain raised, the hotspots near 1. */
export function labScores(): Float32Array {
  const r = seeded(175);
  const v = Float32Array.from({ length: P53.length }, (_, i) => {
    const domain = i >= 93 && i < 292 ? 0.35 : 0;
    return Math.min(1, Math.round((r() * 0.45 + domain) * 100) / 100);
  });
  for (const label of [175, 248, 273]) v[label - 1] = 1;
  return v;
}

export const scoreColor = () =>
  theme.field('value', { palette: theme.binned({ edges: [...ScoreEdges], colors: [...ScoreColors] }) });

export function build(options: ComposedOptions = {}): Node {
  const b = createSeqSpec({ title });
  const c = compose(b, options);
  const lab = b.data.sequence('lab', { alphabet: 'protein', letters: P53, name: 'lab' });
  b.data.map('lab-on-uniprot', { from: lab, to: c.protein as never, method: 'offset', by: 0 });
  const scores = b.data.values('lab-scores', { target: lab, values: labScores() });
  b.views.view('main', (v) => {
    v.section('lab', { of: lab }, (s) => {
      s.track('ruler', { label: 'Lab (from 1)' }).layer('ruler', { type: 'ruler' });
      s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
      s.track('scores', { label: 'Predictions (synthetic)' }).layer('scores', {
        type: 'values',
        data: scores,
        domain: [0, 1],
        color: scoreColor(),
      });
    });
    structureSection(v, c, '/views/main/lab');
  });
  return b.build();
}
