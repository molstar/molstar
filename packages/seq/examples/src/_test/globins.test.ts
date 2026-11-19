/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { DataStore, type FormatRuntime } from '@molstar/seqx-core/data';
import {
  attachSources,
  bundledTextMetrics,
  type FrameInput,
  hitView,
  LightTokens,
  Renderer,
  type RowsLay,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { generatedLabel, generatedRow, generatedRuntime, globinFamily, HbbHuman } from '../globins/data.js';
import * as globins from '../globins/family.js';

// The Globins views drawn against their data (the examples' rule): at the fit and across a tile edge every row in the
// band has its band of cells and its label, and every cell at letters is hit as its row and its column.

async function settled(spec: Node, runtimes: readonly FormatRuntime[] = []) {
  const model = resolveSpec(spec, { runtimes });
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store);
  const renderer = new Renderer({ model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() });
  return async (input: FrameInput) => {
    let f = renderer.frame(input);
    for (let i = 0; i < 6; i++) {
      store.want(f.wants);
      await store.whenSettled();
      const next = renderer.frame(input);
      f = next;
      if (!next.rebuilt) break;
    }
    return renderer.frame(input);
  };
}

function rowsOf(f: {
  layout: { sections: readonly { lays: readonly (readonly { lay: unknown; type: { name: string } }[])[] }[] };
}) {
  for (const s of f.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows layer');
}

/** Every display row in the band: its band of cells (one y each) and its label; each cell in view hit as its row and column. */
function checkDrawn(
  f: Parameters<typeof rowsOf>[0] & { layout: Parameters<typeof hitView>[0] },
  letters: (row: number) => string,
) {
  const lay = rowsOf(f);
  const ys = new Set<number>();
  for (let k = 0; k < lay.n; k++) ys.add(lay.y[k]);
  let rows = 0;
  for (let i = lay.first; i < lay.last; i++) {
    const e = lay.display.rowAt(i);
    const y = lay.rowY(i);
    if (y < lay.box.y || y + lay.rowScale > lay.box.y + lay.box.height) continue;
    rows++;
    expect(ys.has(y), `row ${i} has its band`).toBe(true);
    if (e < 0) continue;
    const text = letters(e);
    for (let j = 0; j < Math.min(text.length, 40); j++) {
      if (!lay.axis.cell(j)) continue;
      const x = (lay.axis.sx0 + lay.axis.sx1) / 2;
      const h = hitView(f.layout, x, y + lay.rowScale / 2);
      expect(h).toMatchObject({ target: 'cell', position: j, yPosition: e });
    }
  }
  expect(rows).toBeGreaterThan(0);
  return lay;
}

describe('the Globins views against their data', () => {
  it('the generated family is deterministic, its ancestor first; the tree names every row', () => {
    const a = globinFamily();
    const b = globinFamily();
    expect(a.letters).toEqual(b.letters);
    expect(a.letters[0]).toBe(HbbHuman);
    expect(a.label[0]).toBe('HBB_HUMAN/1-147');
    for (const l of a.label) expect(a.newick).toContain(l);
    // the scale variant's rows are a function of the seed and the index
    expect(generatedRow(1234, 300, 1)).toEqual(generatedRow(1234, 300, 1));
    expect(generatedRow(1234, 300, 1)).not.toEqual(generatedRow(1235, 300, 1));
  });

  it('the family at letters: every row in the band drawn and hit as its row and column', async () => {
    const family = globinFamily();
    const frame = await settled(globins.build());
    const f = await frame({
      size: { width: 900, height: 600 },
      dpr: 1,
      settled: true,
      viewports: { '/views/main/family': { start: 0, end: 40 } },
    });
    const lay = checkDrawn(f, (r) => family.letters[r]);
    expect(lay.foldsApply).toBe(true);
    expect(lay.tree.length).toBeGreaterThan(0);
  });

  it('the 10,000 rows across a tile edge (rows 1,020 to 1,030): each row drawn, labeled and hit', async () => {
    const frame = await settled(globins.scale(), [generatedRuntime()]);
    const f0 = await frame({ size: { width: 900, height: 600 }, dpr: 1, settled: true });
    const l0 = rowsOf(f0);
    // the band's display row of data row 1,020, at the band's top, at letters
    const top = l0.display.displayRow(1020);
    const f = await frame({
      size: { width: 900, height: 600 },
      dpr: 1,
      settled: true,
      viewports: {
        '/views/main/rows': { start: 0, end: 40, bands: { '/views/main/rows/rows': { bandTop: top, rowScale: 20 } } },
      },
    });
    const text = (r: number) => String.fromCharCode(...generatedRow(r, 300, 1));
    const lay = checkDrawn(f, text);
    for (let i = lay.first; i < lay.last; i++) {
      const e = lay.display.rowAt(i);
      expect(lay.label(e)).toBe(generatedLabel(e));
    }
  });

  it('the lazy 10,000,000 rows: the rows in view read alone, drawn and hit', async () => {
    const frame = await settled(globins.lazy(), [generatedRuntime()]);
    const f = await frame({
      size: { width: 900, height: 400 },
      dpr: 1,
      settled: true,
      viewports: {
        '/views/main/rows': {
          start: 0,
          end: 40,
          bands: { '/views/main/rows/rows': { bandTop: 5_000_000, rowScale: 20 } },
        },
      },
    });
    const text = (r: number) => String.fromCharCode(...generatedRow(r, 300, 2));
    const lay = checkDrawn(f, text);
    expect(lay.held).toBeUndefined();
    expect(lay.first).toBe(5_000_000);
    expect(lay.label(5_000_000)).toBe(generatedLabel(5_000_000));
  });
});
