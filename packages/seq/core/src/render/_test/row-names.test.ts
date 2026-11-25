/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The names of a rows system in the label column (_spec/README.md "Render", Row names): an alignment's rows and a
// matrix's states are named by the one part, in one view: the same font, ink and column, a name a row, none below the
// label minimum, and the matrix's names replace its track's label.

import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { LightTokens, recordDrawLists } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { DefaultRenderConfig } from '../config.js';
import type { MatrixLay } from '../matrix.js';
import { settle } from './helpers.js';

const Samples = ['S01', 'S02', 'a sample with a name too long for the label column to hold'];
const Rows = ['MVLSPA', 'MVHLTP', 'MGLSDG'];

function spec(matrixHeight: number) {
  const b = createSeqSpec({ title: 'names' });
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: ['A', 'B', 'C'], letters: Rows } });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
  const values = Float32Array.from({ length: 8 * Samples.length }, (_, i) => (i % 5) / 5);
  const m = b.data.matrix('m', { target: s, values, states: Samples });
  b.views.view('main', (v) => {
    v.section('aln', { of: aln }, (x) =>
      x.track('rows', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln } as never)),
    );
    v.section('s', { of: s }, (x) =>
      x.track('m', { label: 'Samples', height: matrixHeight }, (t) => t.layer('m', { type: 'matrix', data: m })),
    );
  });
  return b.build();
}

/** The label-pass texts in the item font, with their x, ink and room. */
function names(lists: Parameters<typeof recordDrawLists>[0]) {
  return recordDrawLists(lists, ['labels'])
    .split('\n')
    .filter((l) => l.startsWith('text ') && l.includes('"10px sans"'))
    .map((l) => /^text "([^"]*)" ([\d.-]+) ([\d.-]+) start "10px sans" (#\w+) width ([\d.]+) max ([\d.]+)/.exec(l))
    .filter((m): m is RegExpExecArray => m !== null)
    .map((m) => ({ text: m[1], x: Number(m[2]), ink: m[4], max: Number(m[6]) }));
}

describe('the names of a rows system', () => {
  it('an alignment’s rows and a matrix’s states, one part: one font, ink and column, the long one cut, the track label replaced', async () => {
    const s = await settle(spec(60), [SeqRuntime]);
    const f = await s.frame({ size: { width: 600, height: 420 } });
    const all = names(f.lists);
    const rows = all.filter((n) => ['A', 'B', 'C'].includes(n.text));
    const states = all.filter((n) => n.text.startsWith('S0') || n.text.startsWith('a sample'));
    expect(rows).toHaveLength(3);
    expect(states).toHaveLength(3);
    const hex = `#${LightTokens.colors['msui-text'].toString(16).padStart(8, '0')}`;
    for (const n of [...rows, ...states]) {
      expect(n.ink).toBe(hex);
      expect(n.x).toBe(DefaultRenderConfig.inset);
    }
    // the long state name is cut with an ellipsis within the column, the short ones whole
    expect(states.map((n) => n.text.endsWith('…'))).toEqual([false, false, true]);
    expect(Math.max(...states.map((n) => n.max))).toBeLessThanOrEqual(DefaultRenderConfig.labelColumn.max);
    // the matrix's track has no label of its own to draw over its names
    const track = f.layout.sections[1].tracks[0];
    expect((track.layers[0].lay as MatrixLay).nameText).toHaveLength(3);
    expect(track.labelText).toBe('');
    expect(recordDrawLists(f.lists, ['labels'])).not.toContain('text "Samples"');
  });

  it('a row under the label minimum is named by none, and the track keeps its label', async () => {
    // 3 states in 15 px: 5 px a row, under `rows.labelMinPx`
    const s = await settle(spec(15), [SeqRuntime]);
    const f = await s.frame({ size: { width: 600, height: 420 } });
    const track = f.layout.sections[1].tracks[0];
    expect((track.layers[0].lay as MatrixLay).nameText).toEqual([]);
    expect(track.labelText).toBe('Samples');
  });
});
