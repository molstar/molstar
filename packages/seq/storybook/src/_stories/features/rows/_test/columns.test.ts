/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, type RowsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { FeatureViews } from '@molstar/seqx-examples';
import { columnsTsv } from '@molstar/seqx-examples/globins/columns';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// Features › Rows › Columns, headless on a fake clock, composed as its story composes it (a viewer of the view and its
// SVG export, the joined file read through the assets port): the columns on both sides, each kind, the joined fields,
// a title's click ordering the rows, a narrow width hiding the right columns first; the canvas and the SVG alike in
// both themes. The layer's own tests are `core/src/render/_test/rows-columns.test.ts`.

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const view = FeatureViews.find((v) => v.id === 'columns')!;

async function open(): Promise<SeqViewer> {
  viewer = await openExample(view, { gene: columnsTsv(), delayMs: 0 });
  // wide enough for every column
  viewer.resize(1100);
  await vi.advanceTimersByTimeAsync(100);
  viewer.draw();
  return viewer;
}

function rowsLay(v: SeqViewer): RowsLay {
  for (const s of v.frame!.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows');
}

/** The SVG export, as the story draws it. */
async function svgOf(v: SeqViewer) {
  const p = v.renderSvg({ idPrefix: 'c-' });
  await vi.advanceTimersByTimeAsync(400);
  return p;
}

/** The shown rows' labels in display order. */
const shown = (l: RowsLay) => Array.from({ length: l.display.count }, (_, i) => l.label(l.display.rowAt(i)));

describe('Features › Rows › Columns, headless', () => {
  it('draws a column of each kind on both sides, the joined file’s too; the canvas and the SVG alike, both themes', async () => {
    const v = await open();
    expect(v.layoutIssues.value).toEqual([]);
    const l = rowsLay(v);
    expect(l.colW.every((w) => w > 0)).toBe(true);
    expect(l.rowColumns.map((c) => [c.side, c.as])).toEqual([
      ['left', 'text'],
      ['left', 'chip'],
      ['right', 'number'],
      ['right', 'bar'],
      ['right', 'number'],
    ]);
    // the joined fields by row key: every fifth row has none
    expect(l.fieldValue(0, 'class')).toBe('Mammalia');
    expect(l.fieldValue(0, 'score')).toBe(294);
    expect(l.fieldValue(4, 'class')).toBeUndefined();
    expect(l.fieldValue(0, 'identity')).toBe(1);
    expect(l.titleText).toEqual(['Species', 'Class', 'Identity', 'Residues', 'Score']);
    const light = await svgOf(v);
    expect(light.issues).toEqual([]);
    // the export lays the view out at its own size (the canvas's rows rail aside): the same columns, cells and titles
    const drawn = (svg: string) => {
      for (const t of [...l.titleText, ...l.cellText.filter(Boolean)]) expect(svg).toContain(`>${t}</text>`);
    };
    drawn(light.svg);
    for (const t of ['Homo sapiens', 'Mammalia', 'Identity', 'Score', '294']) expect(light.svg).toContain(t);
    v.setTokens(DarkTokens);
    v.draw();
    const dark = await svgOf(v);
    expect(dark.issues).toEqual([]);
    drawn(dark.svg);
    // the same marks in other colors
    expect(recordDrawLists(dark.lists, ['labels'])).not.toBe(recordDrawLists(light.lists, ['labels']));
    expect(dark.lists.passes.labels.size).toBe(light.lists.passes.labels.size);
  });

  it('a click on a title orders the rows by its field, a second click the other way', async () => {
    const v = await open();
    expect(shown(rowsLay(v))[0]).toBe('HBB_HUMAN/1-147');
    const clickScore = async () => {
      const l = rowsLay(v);
      const t = l.titleText.indexOf(l.titleText.find((x) => x.startsWith('Score'))!);
      expect(v.input.click(l.titleX[t] + 2, l.titleY - 3)).toBe(true);
      await vi.advanceTimersByTimeAsync(100);
      v.draw();
    };
    await clickScore();
    const asc = rowsLay(v);
    const scores = (l: RowsLay) =>
      Array.from({ length: l.display.count }, (_, i) => l.fieldValue(l.display.rowAt(i), 'score'));
    const a = scores(asc).filter((s) => s !== undefined) as number[];
    expect(a).toEqual([...a].sort((x, y) => x - y));
    // the rows without a score last
    expect(
      scores(asc)
        .slice(a.length)
        .every((s) => s === undefined),
    ).toBe(true);
    expect(asc.titleText).toContain('Score ↑');
    await clickScore();
    const d = scores(rowsLay(v)).filter((s) => s !== undefined) as number[];
    expect(d).toEqual([...a].reverse());
    expect(rowsLay(v).titleText).toContain('Score ↓');
    expect(shown(rowsLay(v))[0]).toBe('HBB_HUMAN/1-147');
  });

  it('narrowed, the right columns yield first, then the left; a chip counts those hidden', async () => {
    const v = await open();
    v.resize(620);
    v.draw();
    const mid = rowsLay(v);
    expect(mid.colW[4]).toBe(0);
    expect(mid.colW[0]).toBeGreaterThan(0);
    expect(mid.colChip).toMatch(/hidden$/);
    v.resize(420);
    v.draw();
    const narrow = rowsLay(v);
    expect(narrow.colW).toEqual([0, 0, 0, 0, 0]);
    expect(narrow.colChip).toBe('5 columns hidden');
    expect(narrow.titleText).toEqual([]);
  });
});
