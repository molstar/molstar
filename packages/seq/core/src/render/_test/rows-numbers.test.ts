/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows' end numbers and the band's foot (_spec/README.md "Render", The rows layer: the header): a row's first and
// last residue numbers are those of the residues in view (a gap never counted, none for a row with no residue in view),
// the row's own at the fit view, held and lazy; the last numbers in the right gutter clipped to the band as its rows
// are, so a partly shown row ends at the band's foot on every side.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type Frame, type FrameInput, recordDrawLists } from '../index.js';
import type { RowsLay } from '../rows.js';
import { settle } from './helpers.js';
import { rowsSpec, Section, Synthetic, Track } from './rows-fixtures.js';

const Labels = ['HBA/1-6', 'HBB/2-7', 'MYG/1-4', 'GLB/3-8', 'CYG/1-6'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLS--', 'MKAVLG', 'MEKVPG'];
const S = '/views/v/s';
const T = '/views/v/s/rows';

function spec(height: number) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  b.views.view('v', (v) =>
    v.section('s', { of: aln } as never, (s) =>
      s.track('rows', { label: '', height }, (t) => t.layer('r', { type: 'rows', data: aln } as never)),
    ),
  );
  return b.build();
}

function rowsLay(f: Frame): RowsLay {
  for (const s of f.layout.sections)
    for (const l of s.lays) for (const x of l ?? []) if (x.type.name === 'rows') return x.lay as RowsLay;
  throw new Error('no rows layer');
}

const view = (start: number, end: number, extra: Record<string, unknown> = {}): FrameInput => ({
  size: { width: 600, height: 300 },
  dpr: 1,
  settled: true,
  viewports: { [S]: { start, end, bands: { [T]: { bandTop: 0, rowScale: 20, ...extra } as never } } },
});

describe('the end numbers of the residues in view', () => {
  it('held: the rows’ own at the fit view; after a pan and a zoom those of the columns in view', async () => {
    const s = await settle(spec(100));
    const fit = rowsLay(await s.frame(view(0, 6)));
    expect([fit.numberText, fit.lastText]).toEqual([
      ['1', '2', '1', '3', '1'],
      ['6', '7', '4', '8', '6'],
    ]);
    // columns 2 to 4: the residues before them counted, the gaps of MYG never
    const pan = rowsLay(await s.frame(view(2, 5)));
    expect([pan.numberText, pan.lastText]).toEqual([
      ['3', '4', '3', '5', '3'],
      ['5', '6', '4', '7', '5'],
    ]);
    // columns 4 and 5: MYG has no residue in view, so no numbers (the others keep theirs, in order)
    const zoom = rowsLay(await s.frame(view(4, 6)));
    expect([zoom.numberText, zoom.lastText]).toEqual([
      ['5', '6', '7', '5'],
      ['6', '7', '8', '6'],
    ]);
    expect(zoom.numberY.length).toBe(4);
    // back at the fit view: the rows' own again
    expect(rowsLay(await s.frame(view(0, 6))).lastText).toEqual(fit.lastText);
  });

  it('lazy: the same, counted from the rows’ tiles; the last number never past the row’s letters', async () => {
    const s = await settle(rowsSpec(200_000, false), [Synthetic]);
    const frame = (start: number, end: number) =>
      s.frame({
        size: { width: 900, height: 400 },
        dpr: 1,
        settled: true,
        viewports: { [Section]: { start, end, bands: { [Track]: { bandTop: 100_000, rowScale: 20 } } } },
      });
    const check = (lay: RowsLay, c0: number, c1: number) => {
      expect(lay.held).toBeUndefined();
      let k = 0;
      let rows = 0;
      for (let i = lay.first; i < lay.last; i++) {
        const e = lay.display.rowAt(i);
        const codes = lay.lazyCodes.get(e)!;
        let before = 0;
        let inView = 0;
        for (let j = 0; j < c1; j++) if (codes[j] !== 45) j < c0 ? before++ : inView++;
        if (!inView) continue;
        expect([lay.numberText[k], lay.lastText[k]]).toEqual([String(1 + before), String(before + inView)]);
        k++;
        rows++;
      }
      expect(rows).toBeGreaterThan(0);
      expect(lay.numberText.length).toBe(rows);
      expect(lay.lastText.length).toBe(rows);
    };
    check(rowsLay(await frame(0, 300)), 0, 300);
    check(rowsLay(await frame(50, 150)), 50, 150);
    check(rowsLay(await frame(0, 300)), 0, 300);
  });
});

describe('the band’s foot', () => {
  it('a partly shown last row: its cells reach the foot, its last number clipped to the band like its label', async () => {
    const s = await settle(spec(90));
    const f = await s.frame(view(0, 6));
    const lay = rowsLay(f);
    const foot = lay.box.y + lay.box.height;
    // the fifth row starts above the foot and ends below it
    expect(lay.rowY(4)).toBeLessThan(foot);
    expect(lay.rowY(4) + lay.rowScale).toBeGreaterThan(foot);
    const lines = recordDrawLists(f.lists, ['labels', 'data']).split('\n');
    const clips: number[][] = [];
    let cells = 0;
    let lasts = 0;
    for (const l of lines) {
      const w = l.split(' ');
      if (w[0] === 'pass') clips.length = 0;
      else if (w[0] === 'clip') clips.push(w.slice(1).map(Number));
      else if (w[0] === 'unclip') clips.pop();
      else if (w[0] === 'text' && w[4] === 'start' && Number(w[2]) > lay.box.x + lay.box.width) {
        lasts++;
        const c = clips[clips.length - 1];
        expect(c, 'the last numbers are clipped').toBeDefined();
        expect([c[1], c[1] + c[3]]).toEqual([lay.box.y, foot]);
      } else if (w[0] === 'rect' && Number(w[2]) === lay.rowY(4) && Number(w[4]) === lay.rowScale) {
        // a cell of the partly shown row: whole, clipped to the band's foot by the clip it is drawn in
        cells++;
        const c = clips[clips.length - 1];
        expect(c[1] + c[3]).toBe(foot);
      }
    }
    expect(lasts).toBe(5);
    expect(cells).toBeGreaterThan(0);
  });
});
