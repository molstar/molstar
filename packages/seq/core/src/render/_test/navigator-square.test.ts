/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import { renderConfig } from '../config.js';
import { SquareNavigator, washRuns } from '../navigator-square.js';
import { LightTokens } from '../tokens.js';
import { genomeSpec, toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';
import { withNavigator } from './navigator-fixtures.js';

// The square (_spec/README.md "Render", Navigators, The square): the wash by runs (no cell twice, at most `maxCells`
// rectangles, none past it with its own issue), the base drawn once a size and tokens (a window move redraws the window).

beforeAll(async () => {
  await loadGenome();
});

const Section = '/views/v/genome';

/** The genome of `n` chromosomes with its square (one section on both axes), the section showing `include` when given. */
async function square(n: number, include?: readonly string[], maxCells?: number) {
  let spec: Node = withNavigator(genomeSpec(toyGenome(n, 0)), {
    of: [Section, Section],
    axes: 'both',
    size: [220, 220],
  });
  if (include) {
    const views = spec.children!.find((c) => c.kind === 'views')!;
    const view = views.children![0];
    const sections = view.children!.map((c) =>
      c.id === 'genome' ? { ...c, params: { ...c.params, include: [...include] } } : c,
    );
    spec = {
      ...spec,
      children: spec.children!.map((c) =>
        c === views ? { ...views, children: [{ ...view, children: sections }] } : c,
      ),
    };
  }
  const s = await settle(spec);
  const sq = new SquareNavigator({
    model: s.model,
    view: s.view,
    navigator: s.view.navigators[0],
    tokens: LightTokens,
    config: renderConfig(maxCells ? { navigator: { square: { maxCells } } } : {}),
  });
  const frame = await s.frame({ size: { width: 900, height: 300 } });
  return { s, sq, layout: frame.layout };
}

describe('The square', () => {
  it('washes by runs: every cell a section leaves out once, the rest never (property over masks)', () => {
    let seed = 7;
    const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    for (let t = 0; t < 200; t++) {
      const n = 1 + Math.floor(rnd() * 12);
      const x = Array.from({ length: n }, () => rnd() < 0.6);
      const y = Array.from({ length: n }, () => rnd() < 0.6);
      const runs = washRuns(x, y);
      const covered = new Int32Array(n * n);
      for (let r = 0; r < runs.length; r += 4)
        for (let i = runs[r]; i < runs[r + 1]; i++)
          for (let j = runs[r + 2]; j < runs[r + 3]; j++) covered[i * n + j]++;
      for (let i = 0; i < n; i++)
        for (let j = 0; j < n; j++) expect(covered[i * n + j], `${t} ${i},${j}`).toBe(x[i] && y[j] ? 0 : 1);
      // as many rectangles as runs, never a cell each
      expect(runs.length / 4).toBeLessThanOrEqual(n + Math.ceil(n / 2) * Math.ceil(n / 2));
    }
  });

  it('a section showing some chromosomes: a few rectangles of wash, drawn once; a window move redraws the window alone', async () => {
    const { sq, layout } = await square(12, ['chr2', 'chr3', 'chr7']);
    const f = sq.frame({ size: { width: 220, height: 220 }, view: layout });
    expect(f.rebuilt).toBe(true);
    const data = f.lists.passes.data.rects.n;
    // 13 lines each way, and the wash: chr2–3 and chr7 shown make 3 column runs and 3 row runs over 2 shown runs
    expect(data).toBe(26 + 3 + 3 * 2);
    const again = sq.frame({ size: { width: 220, height: 220 }, view: layout });
    expect(again.rebuilt).toBe(false);
    expect(again.lists.passes.data.rects.n).toBe(data);
    expect(again.lists.passes.window.rects.n).toBe(1);
  });

  it('past maxCells: no lattice and no wash (never n² a frame), said with its own code', async () => {
    const { sq, layout } = await square(40, ['chr1'], 1000);
    const f = sq.frame({ size: { width: 220, height: 220 }, view: layout });
    expect(f.lists.passes.data.rects.n).toBe(0);
    expect(sq.allIssues.map((i) => i.code)).toEqual(['seq.render.navigator.square']);
  });
});
