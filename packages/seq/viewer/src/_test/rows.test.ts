/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An alignment's rows in the viewer (_spec/README.md "Behaviors", Rows), headless on a fake clock: a tree node's click
// folds its clade (session state, no undo step), the bulk folds are jumps that back restores, a reveal opens the clades
// holding a row, `clades.keep` writes the folds into `collapsed` and an edit of `collapsed` resets the session's set,
// the selection is by data rows and a fold keeps it, the wheel scrolls the rows, `]` and `[` zoom them, the export draws
// the viewer's folds, and the keys check reports one key event two commands' chords match.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import type { RowsLay } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { keyConflicts } from '../bindings.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Track = '/views/v/s/rows';
const Labels = ['HBA', 'HBB', 'MYG', 'GLB', 'CYG', 'NGB'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLSDG', 'MKAVLG', 'MEKVPG', 'MERPEP'];

function spec(collapsed: string[] = []): Node {
  const b = createSeqSpec({ title: 'rows' });
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick: '((HBA,HBB)hb,(MYG,(GLB,(CYG,NGB)cn)gc)rest)root;' });
  b.views.view('v', (v) =>
    v
      .section('s', { of: aln })
      .track('rows', { label: '', height: 60 })
      .layer('r', { type: 'rows', data: aln, tree, order: 'tree', collapsed }),
  );
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(collapsed?: string[]) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
  viewer = v;
  await v.load(spec(collapsed));
  v.resize(400);
  await v.settle();
  frames.flush();
  v.input.mac = false;
  return v;
}

const lay = (v: SeqViewer) => {
  for (const s of v.frame!.layout.sections)
    for (const band of s.lays) for (const l of band ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows');
};

describe('an alignment’s rows in the viewer', () => {
  it('a click on a tree node folds its clade (session state, no undo step), again unfolds it', async () => {
    const v = await make();
    const l = lay(v);
    expect(l.display.count).toBe(6);
    const revision = v.spec.revision.value;
    // the node of `hb`: its branch's x and y
    const k = l.nodes.findIndex((_, i) => i % 3 === 0 && l.clades!.name[l.nodes[i]] === 'hb');
    expect(v.input.click(l.nodes[k + 1], l.nodes[k + 2])).toBe(true);
    frames!.flush();
    expect(v.view.band(Track)?.folds).toEqual(['hb']);
    expect(lay(v).display.count).toBe(5);
    expect(v.spec.revision.value).toBe(revision);
    const w = lay(v).wedges;
    expect(v.input.click(w[0] + 2, (w[1] + w[3]) / 2)).toBe(true);
    frames!.flush();
    expect(v.view.band(Track)?.folds).toEqual([]);
  });

  it('collapse all and a reveal are jumps: back restores the folds; the selection is by data rows and a fold keeps it', async () => {
    const v = await make();
    v.selection.add(rangesLocation('/data/aln/(rows)', [4, 5]));
    const selected = v.selection.value.value;
    await v.execute('clades.collapseAll', { track: Track });
    frames!.flush();
    expect(lay(v).display.count).toBe(2);
    expect(v.selection.value.value).toBe(selected);
    // a reveal of NGB opens the clades holding it
    await v.execute('clades.reveal', { track: Track, rows: [5] });
    frames!.flush();
    const open = v.view.band(Track)!.folds!;
    expect(open).toEqual(['hb']);
    expect(v.view.back()).toBe(true);
    frames!.flush();
    expect(v.view.band(Track)?.folds?.length).toBe(4);
    expect(v.view.back()).toBe(true);
    expect(v.view.band(Track)?.folds).toBeUndefined();
  });

  it('`clades.keep` writes the folds shown into `collapsed`; an edit of `collapsed` resets the session’s folds', async () => {
    const v = await make(['gc']);
    expect(lay(v).display.count).toBe(4);
    await v.execute('clades.toggle', { track: Track, clade: 'hb' });
    frames!.flush();
    expect(lay(v).display.count).toBe(3);
    await v.execute('clades.keep', { track: Track });
    await v.settle();
    frames!.flush();
    const layer = v.spec.root.value.children![1].children![0].children![0].children![0].children![0];
    expect(layer.params?.collapsed).toEqual(['gc', 'hb']);
    // the opening state changed: the session's set is reset to it (shown the same)
    expect(v.view.band(Track)?.folds).toBeUndefined();
    expect(lay(v).display.count).toBe(3);
  });

  it('the wheel scrolls the rows in their band, `]` and `[` zoom them; the export draws the viewer’s folds', async () => {
    const v = await make();
    const l = lay(v);
    const x = l.box.x + 5;
    const y = l.rowY(1) + 3;
    v.input.move(x, y);
    expect(v.input.wheel(x, y, { deltaX: 0, deltaY: 30, deltaMode: 0 })).toBe(true);
    frames!.flush();
    expect(v.view.band(Track)?.bandTop).toBe(30);
    expect(lay(v).first).toBe(1);
    expect(v.input.key({ key: ']' })).toBe(true);
    frames!.flush();
    expect(v.view.band(Track)?.rowScale).toBe(40);
    expect(v.input.key({ key: '[' })).toBe(true);
    expect(v.view.band(Track)?.rowScale).toBe(20);
    await v.execute('clades.toggle', { track: Track, clade: 'hb' });
    v.view.setBand(Track, { bandTop: 0 });
    frames!.flush();
    const svg = await v.renderSvg({ dpr: 1 });
    expect(svg.svg).toContain('HBA +1');
  });

  it('the keys check reports a chord one key event matches with another command’s (Shift+= is `+`)', () => {
    expect(keyConflicts({ 'view.zoom-in': ['=', '+'], 'rows.zoom-in': [']'], 'rows.zoom-out': ['['] })).toEqual([]);
    const issues = keyConflicts({ 'view.zoom-in': ['+'], 'rows.zoom-in': ['Shift+='] });
    expect(issues.map((i) => i.code)).toEqual(['seq.viewer.binding.conflict']);
  });
});
