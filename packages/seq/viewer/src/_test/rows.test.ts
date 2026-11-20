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

import { patchToNode } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { type RowsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { headlessAnimation } from '../animation.js';
import { keyConflicts } from '../bindings.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Track = '/views/v/s/rows';
const Labels = ['HBA', 'HBB', 'MYG', 'GLB', 'CYG', 'NGB'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLSDG', 'MKAVLG', 'MEKVPG', 'MERPEP'];

const Newick = '((HBA,HBB)hb,(MYG,(GLB,(CYG,NGB)cn)gc)rest)root;';
/** Support values as inner labels: two clades labelled `95`, one `90`, each named by its first and last leaves. */
const Support = '((HBA,HBB)95,(MYG,(GLB,(CYG,NGB)95)90)rest)root;';

function spec(collapsed: string[] = [], newick = Newick): Node {
  const b = createSeqSpec({ title: 'rows' });
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick });
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

async function make(collapsed?: string[], newick?: string) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
  viewer = v;
  await v.load(spec(collapsed, newick));
  v.resize(480);
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

  it('a pointer over a tree node lights its clade through the hover model: branches and one band over its rows', async () => {
    const v = await make();
    const l = lay(v);
    const k = l.nodes.findIndex((_, i) => i % 3 === 0 && l.clades!.name[l.nodes[i]] === 'hb');
    v.input.move(l.nodes[k + 1], l.nodes[k + 2]);
    frames!.flush();
    expect(v.hover.current.target).toBe('tree-node');
    const rects = recordDrawLists(v.frame!.lists, ['hover'])
      .split('\n')
      .filter((r) => r.startsWith('rect '))
      .map((r) => r.slice(5).split(' ').map(Number));
    const bands = rects.filter((r) => r[0] === l.box.x && r[2] === l.box.width);
    expect(bands.map((r) => [r[1], r[3]])).toEqual([[l.rowY(0), 2 * l.rowScale]]);
    expect(rects.filter((r) => r[0] < l.box.x).length).toBeGreaterThan(0);
    // off the tree: nothing of the clade
    v.input.move(l.box.x + 5, l.rowY(0) + 3);
    frames!.flush();
    expect(v.hover.current.target).not.toBe('tree-node');
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

  it('clades labelled alike (support values) are named apart: a click, a hover, a toggle and `clades.keep` reach each', async () => {
    const v = await make([], Support);
    const l = lay(v);
    expect(l.clades!.name).toContain('HBA..HBB');
    expect(l.clades!.name).toContain('CYG..NGB');
    expect(l.clades!.name).toContain('GLB..NGB');
    expect(l.clades!.name).toContain('rest');
    // a pointer and a click on the first `95` node (the rows in view): that clade by its name
    const k = l.nodes.findIndex((_, i) => i % 3 === 0 && l.clades!.name[l.nodes[i]] === 'HBA..HBB');
    v.input.move(l.nodes[k + 1], l.nodes[k + 2]);
    frames!.flush();
    expect(v.hover.current.element?.key).toBe('HBA..HBB');
    expect(v.input.click(l.nodes[k + 1], l.nodes[k + 2])).toBe(true);
    frames!.flush();
    expect(v.view.band(Track)?.folds).toEqual(['HBA..HBB']);
    expect(lay(v).display.count).toBe(5);
    // the second `95`, by its own name: CYG and NGB in one summary row, the first clade still folded
    await v.execute('clades.toggle', { track: Track, clade: 'CYG..NGB' });
    frames!.flush();
    expect(v.view.band(Track)?.folds).toEqual(['HBA..HBB', 'CYG..NGB']);
    expect(lay(v).display.count).toBe(4);
    await v.execute('clades.keep', { track: Track });
    await v.settle();
    frames!.flush();
    const layer = v.spec.root.value.children![1].children![0].children![0].children![0].children![0];
    expect(layer.params?.collapsed).toEqual(['HBA..HBB', 'CYG..NGB']);
    expect(lay(v).display.count).toBe(4);
  });

  it('the wheel scrolls the rows in their band, `]` and `[` zoom them; the export draws the viewer’s folds', async () => {
    const v = await make();
    const l = lay(v);
    const x = l.box.x + 5;
    const y = l.rowY(1) + 3;
    v.input.move(x, y);
    expect(v.input.wheel(x, y, { deltaX: 0, deltaY: 30, deltaMode: 0 })).toBe(true);
    frames!.flush();
    // 30 px of wheel at 20 px a row: a row and a half
    expect(v.view.band(Track)?.bandTop).toBe(1.5);
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

  it('after a pointer click into the rows, PageDown and PageUp move the window a page (its rows less one); at an end the key is the page’s', async () => {
    const v = await make();
    const l = lay(v);
    const page = Math.floor(l.box.height / l.rowScale) - 1;
    expect(page).toBeGreaterThan(0);
    const max = l.maxTop(l.rowScale, l.display.count);
    expect(max).toBeGreaterThan(page);
    // a click into the rows, as the mount delivers it (down, up, click); the key is the viewer's with no cursor set
    const x = l.box.x + 5;
    const y = l.rowY(1) + 3;
    v.input.move(x, y);
    v.input.down(x, y, {});
    v.input.up();
    v.input.click(x, y, {});
    frames!.flush();
    expect(v.input.key({ key: 'PageUp' })).toBe(false);
    expect(v.input.key({ key: 'PageDown' })).toBe(true);
    frames!.flush();
    expect(v.view.band(Track)?.bandTop).toBe(page);
    expect(lay(v).first).toBe(page);
    expect(v.input.key({ key: 'PageDown' })).toBe(true);
    expect(v.view.band(Track)?.bandTop).toBe(max);
    frames!.flush();
    expect(v.input.key({ key: 'PageDown' })).toBe(false);
    expect(v.input.key({ key: 'PageUp' })).toBe(true);
    expect(v.view.band(Track)?.bandTop).toBe(max - page);
  });

  it('PageDown with pinned rows: the pinned band is not part of the page (200 px, 20 px rows, 2 pinned: row 8 is first)', async () => {
    vi.useFakeTimers({ toFake: ['Date', 'performance'] });
    frames = installFrames();
    const labels = Array.from({ length: 40 }, (_, i) => `R${i + 1}`);
    const b = createSeqSpec({ title: 'pinned' });
    const aln = b.data.alignment('aln', {
      alphabet: 'protein',
      rows: { label: labels, letters: labels.map(() => 'MVLSPA') },
    });
    b.views.view('v', (v) =>
      v
        .section('s', { of: aln })
        .track('rows', { label: '', height: 200 })
        .layer('r', { type: 'rows', data: aln, pinned: ['R1', 'R2'] } as never),
    );
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
    viewer = v;
    await v.load(b.build());
    v.resize(480);
    await v.settle();
    frames.flush();
    expect(lay(v).display.pinned).toBe(2);
    expect(lay(v).rowScale).toBe(20);
    expect(v.input.key({ key: 'PageDown' })).toBe(true);
    frames.flush();
    // 8 scrolled rows in the window, a page of 7: the 8th scrolled row is the first now (display row 7 after the pinned)
    expect(v.view.band(Track)?.bandTop).toBe(7);
    expect(lay(v).labelText[lay(v).display.pinned]).toBe('R8');
  });

  it('the keys check reports a chord one key event matches with another command’s (Shift+= is `+`)', () => {
    expect(keyConflicts({ 'view.zoom-in': ['=', '+'], 'rows.zoom-in': [']'], 'rows.zoom-out': ['['] })).toEqual([]);
    const issues = keyConflicts({ 'view.zoom-in': ['+'], 'rows.zoom-in': ['Shift+='] });
    expect(issues.map((i) => i.code)).toEqual(['seq.viewer.binding.conflict']);
  });
});

const Rows = '/data/aln/(rows)';
const layer = (v: SeqViewer) => v.spec.root.value.children![1].children![0].children![0].children![0].children![0];

describe('the rows’ band in display rows, kept by key', () => {
  it('a fold above the window keeps the window’s top row', async () => {
    const v = await make();
    v.view.setBand(Track, { bandTop: 3 });
    frames!.flush();
    expect(lay(v).keyAt(lay(v).top)?.key).toBe('GLB');
    await v.execute('clades.toggle', { track: Track, clade: 'hb' });
    frames!.flush();
    expect(v.view.band(Track)?.bandTop).toBe(2);
    expect(lay(v).keyAt(lay(v).top)?.key).toBe('GLB');
  });

  it('back lands on the same top row by its key, though a row above it was hidden since', async () => {
    const v = await make();
    v.view.setBand(Track, { bandTop: 3 });
    frames!.flush();
    await v.execute('clades.collapseAll', { track: Track });
    frames!.flush();
    await v.execute('rows.hide', { track: Track, rows: [0] });
    await v.settle();
    frames!.flush();
    expect(v.view.back()).toBe(true);
    frames!.flush();
    // GLB is the third row shown now (HBA hidden): its top 2, where a kept px would show CYG
    expect(v.view.band(Track)?.folds).toBeUndefined();
    expect(v.view.band(Track)?.bandTop).toBe(2);
    expect(lay(v).keyAt(lay(v).top)?.key).toBe('GLB');
  });

  it('`]` and `[` zoom about the band’s middle row', async () => {
    const v = await make();
    v.view.setBand(Track, { bandTop: 1, rowScale: 10 });
    frames!.flush();
    const l = lay(v);
    const area = l.box.height - l.pinnedHeight;
    const mid = l.top + area / 2 / l.rowScale;
    expect(v.input.key({ key: ']' })).toBe(true);
    frames!.flush();
    const z = lay(v);
    expect(z.rowScale).toBe(20);
    expect(z.top + area / 2 / z.rowScale).toBeCloseTo(mid, 9);
  });

  it('a step’s motion over a folded band keeps the session’s folds, the frame’s top and zoom over them', async () => {
    const { port, clock } = headlessAnimation();
    vi.useFakeTimers({ toFake: ['Date', 'performance'] });
    frames = installFrames();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, animation: port });
    viewer = v;
    const base = spec();
    const down = { sections: { '/views/v/s': { start: 0, end: 6, bands: { [Track]: { bandTop: 2, rowScale: 20 } } } } };
    await v.load({
      ...base,
      children: [
        ...(base.children ?? []),
        {
          kind: 'snapshots',
          id: 'snapshots',
          children: [
            {
              kind: 'snapshot',
              id: 'down',
              params: { patch: patchToNode({ ops: [] }), view: down, duration: 400, easing: 'linear' },
            },
          ],
        },
      ],
    });
    v.resize(480);
    await v.settle();
    v.draw();
    v.view.setBand(Track, { folds: ['hb'] });
    v.draw();
    expect(lay(v).display.count).toBe(5);
    const { plan } = v.playStep('down');
    expect(plan).toBeDefined();
    clock.tick(200);
    v.draw();
    const mid = lay(v);
    expect(mid.folds).toEqual(['hb']);
    expect(mid.display.count).toBe(5);
    expect(mid.top).toBeGreaterThan(0);
    expect(mid.top).toBeLessThan(2);
    clock.runAll(16);
    v.draw();
    expect(lay(v).folds).toEqual(['hb']);
  });
});

describe('the rows’ commands and gestures', () => {
  it('pins toggle; hidden rows merge and show again; an order; a collapse to a depth', async () => {
    const v = await make();
    await v.execute('rows.pin', { track: Track, rows: [2] });
    await v.settle();
    expect(layer(v).params?.pinned).toEqual(['MYG']);
    await v.execute('rows.pin', { track: Track, rows: [2] });
    await v.settle();
    expect(layer(v).params?.pinned).toEqual([]);
    await v.execute('rows.hide', { track: Track, rows: [0] });
    await v.settle();
    await v.execute('rows.hide', { track: Track, rows: [0, 1] });
    await v.settle();
    expect(layer(v).params?.hidden_rows).toEqual(['HBA', 'HBB']);
    frames!.flush();
    expect(lay(v).display.count).toBe(4);
    await v.execute('rows.show', { track: Track });
    await v.settle();
    expect(layer(v).params?.hidden_rows).toBeUndefined();
    await v.execute('rows.order', { track: Track, order: { field: 'label', direction: 'ascending' } });
    await v.settle();
    frames!.flush();
    const shown = (l: RowsLay) => Array.from({ length: l.display.count }, (_, i) => l.label(l.display.rowAt(i)));
    expect(shown(lay(v))).toEqual(['CYG', 'GLB', 'HBA', 'HBB', 'MYG', 'NGB']);
    await v.execute('rows.order', { track: Track, order: 'tree' });
    await v.settle();
    frames!.flush();
    await v.execute('clades.collapseTo', { track: Track, depth: 1 });
    frames!.flush();
    expect(v.view.band(Track)?.folds).toEqual(['hb', 'rest']);
    expect(lay(v).display.count).toBe(2);
  });

  it('Shift+click on a row header extends in display order, not data order', async () => {
    const v = await make();
    await v.execute('rows.order', { track: Track, order: { field: 'label', direction: 'ascending' } });
    await v.settle();
    frames!.flush();
    const l = lay(v);
    // CYG GLB HBA HBB MYG NGB: CYG, then Shift+HBA selects CYG, GLB and HBA (data rows 4, 3, 0)
    expect(v.input.click(l.labelX + 2, l.rowY(0) + 5)).toBe(true);
    expect(v.input.click(l.labelX + 2, l.rowY(2) + 5, { shiftKey: true })).toBe(true);
    const part = v.selection.value.value.shown.parts.find((p) => p.system === Rows)!;
    expect(Array.from(part.ranges)).toEqual([0, 1, 3, 5]);
  });

  it('a double-click on a summary row opens its clade', async () => {
    const v = await make(['hb']);
    const l = lay(v);
    expect(l.display.rowAt(0)).toBeLessThan(0);
    expect(v.input.doubleClick(l.box.x + 3, l.rowY(0) + 5)).toBe(true);
    await v.settle();
    frames!.flush();
    expect(lay(v).display.count).toBe(6);
  });

  it('the wheel over another layer’s `row` hit pans as any, never taken by the rows', async () => {
    const v = await make();
    const l = lay(v);
    const x = l.box.x + 5;
    const y = l.rowY(1) + 3;
    const input = v.input as unknown as { _hitAt: (x: number, y: number) => unknown; _hit: Record<string, unknown> };
    const own = input._hitAt.bind(input);
    vi.spyOn(input, '_hitAt').mockImplementation((px: number, py: number) => {
      own(px, py);
      Object.assign(input._hit, { layer: '/views/v/s/other/l', target: 'row', ySystem: undefined });
      return input._hit;
    });
    const before = v.view.get('/views/v/s')!.start;
    expect(v.input.wheel(x, y, { deltaX: 0, deltaY: 30, deltaMode: 0 })).toBe(true);
    expect(v.view.band(Track)?.bandTop ?? 0).toBe(0);
    void before;
  });

  it('the bindings: the matrix’s cells', async () => {
    const v = await make();
    const b = v.bindings;
    const none = {};
    expect(b.resolve('click', none, 'tree-node')).toBe('clades.toggle');
    expect(b.resolve('click', { mod: true }, 'tree-node')).toBe('selection.toggle');
    expect(b.resolve('double-click', none, 'summary-row')).toBe('clades.toggle');
    expect(b.resolve('click', none, 'summary-row')).toBe('selection.set');
    expect(b.resolve('click', { shift: true }, 'row')).toBe('selection.extend');
    expect(b.resolve('click', { mod: true }, 'row')).toBe('selection.toggle');
    expect(b.resolve('click', none, 'cell')).toBe('selection.set');
    expect(b.resolve('wheel', none, 'row')).toBe('pan');
    expect(v.config.view.keys['rows.zoom-in']).toEqual([']']);
    expect(v.config.view.keys['rows.zoom-out']).toEqual(['[']);
  });
});

describe('the rows’ columns in the viewer', () => {
  it('a click on a column’s title orders the rows by its field (one patch), a second click the other way', async () => {
    vi.useFakeTimers({ toFake: ['Date', 'performance'] });
    frames = installFrames();
    const b = createSeqSpec({ title: 'columns' });
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
    b.views.view('v', (vw) =>
      vw.section('s', { of: aln }, (s) => {
        s.track('ruler', { label: '' }, (t) => t.layer('ruler', { type: 'ruler' }));
        s.track('rows', { label: '', height: 140 }, (t) =>
          t.layer('r', { type: 'rows', data: aln, columns: [{ field: 'label', side: 'right' }] } as never),
        );
      }),
    );
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
    viewer = v;
    await v.load(b.build());
    v.resize(900);
    await v.settle();
    frames.flush();
    const shown = (l: RowsLay) => Array.from({ length: l.display.count }, (_, i) => l.label(l.display.rowAt(i)));
    const click = async () => {
      const l = lay(v);
      expect(v.input.click(l.cols.titleX[0] + 2, l.cols.titleY - 3)).toBe(true);
      await v.settle();
      frames!.flush();
    };
    const rows = () => v.spec.root.value.children![1].children![0].children![0].children![1].children![0];
    await click();
    expect(rows().params?.order).toEqual({ field: 'label', direction: 'ascending' });
    expect(shown(lay(v))).toEqual(['CYG', 'GLB', 'HBA', 'HBB', 'MYG', 'NGB']);
    expect(lay(v).cols.titleText[0]).toBe('label ↑');
    await click();
    expect(rows().params?.order).toEqual({ field: 'label', direction: 'descending' });
    expect(shown(lay(v))).toEqual(['NGB', 'MYG', 'HBB', 'HBA', 'GLB', 'CYG']);
    // a fast second click (the platform's `detail` 2) flips again, read from the spec as it is now (no frame between)
    const l = lay(v);
    expect(v.input.click(l.cols.titleX[0] + 2, l.cols.titleY - 3, { detail: 1 })).toBe(true);
    expect(v.input.click(l.cols.titleX[0] + 2, l.cols.titleY - 3, { detail: 2 })).toBe(true);
    await v.settle();
    expect(rows().params?.order).toEqual({ field: 'label', direction: 'descending' });
  });

  it('the joined fields follow an alignment read again (its rows in another order), never the old rows', async () => {
    vi.useFakeTimers({ toFake: ['Date', 'performance'] });
    frames = installFrames();
    const build = (labels: readonly string[]) => {
      const b = createSeqSpec({ title: 'joined' });
      const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: [...labels], letters: Letters } });
      b.data.source('extra', {
        format: 'tsv',
        value: 'id\tclade\nHBA\talpha\nNGB\tneuro\n',
        system: '/data/aln/(rows)',
      });
      b.views.view('v', (vw) =>
        vw.section('s', { of: aln }, (s) =>
          s.track('rows', { label: '', height: 140 }, (t) =>
            t.layer('r', { type: 'rows', data: aln, columns: [{ field: 'clade' }] } as never),
          ),
        ),
      );
      return b.build();
    };
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
    viewer = v;
    await v.load(build(Labels));
    v.resize(900);
    await v.settle();
    frames.flush();
    expect(lay(v).fieldValue(0, 'clade')).toBe('alpha');
    // the alignment's rows read again in reverse (one patch: the joined source untouched, its entries as they were)
    const o = await v.spec.patch('Reverse the rows', {
      ops: [
        {
          op: 'set',
          at: ['data', 'aln'],
          key: 'rows',
          value: { label: [...Labels].reverse(), letters: [...Letters].reverse() },
        },
      ],
    });
    expect(o.ok && o.changed).toBe(true);
    await v.settle();
    frames.flush();
    expect(lay(v).label(0)).toBe('NGB');
    expect(lay(v).fieldValue(0, 'clade')).toBe('neuro');
    expect(lay(v).fieldValue(5, 'clade')).toBe('alpha');
  });

  describe('a cut value whole only under the pointer; the hidden columns named', () => {
    const Species = [
      'Homo sapiens (human)',
      'Mus musculus (house mouse)',
      'Physeter catodon',
      'Danio rerio',
      'Gallus gallus',
      'Equus caballus',
    ];
    async function open(width: number) {
      vi.useFakeTimers({ toFake: ['Date', 'performance'] });
      frames = installFrames();
      const b = createSeqSpec({ title: 'columns' });
      const aln = b.data.alignment('aln', {
        alphabet: 'protein',
        rows: { label: Labels, letters: Letters, OS: Species },
      } as never);
      b.views.view('v', (vw) =>
        vw.section('s', { of: aln }, (s) => {
          s.track('ruler', { label: '' }, (t) => t.layer('ruler', { type: 'ruler' }));
          s.track('rows', { label: '', height: 140 }, (t) =>
            t.layer('r', {
              type: 'rows',
              data: aln,
              columns: [
                { field: 'OS', title: 'Species', width: 40 },
                { field: 'identity', side: 'right', as: 'number', title: 'Identity' },
              ],
            } as never),
          );
        }),
      );
      const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
      viewer = v;
      await v.load(b.build());
      v.resize(width);
      await v.settle();
      frames.flush();
      return v;
    }
    const tips = (v: SeqViewer) =>
      recordDrawLists(v.frame!.lists, ['tip'])
        .split('\n')
        .filter((l) => l.startsWith('text '));

    it('the pointer on a cut cell reads it whole (the hover names it); on the row’s data or label, no tip', async () => {
      const v = await open(900);
      const l = lay(v);
      const y = l.rowY(0) + l.rowScale / 2;
      v.input.move(l.cols.x[0] + 4, y);
      frames!.flush();
      expect(v.hover.current.state).toBe('Species: Homo sapiens (human)');
      expect(tips(v)).toHaveLength(1);
      expect(tips(v)[0]).toContain('"Homo sapiens (human)"');
      v.input.move(l.box.x + 4, y);
      frames!.flush();
      expect(v.hover.current.target).toBe('cell');
      expect(tips(v)).toEqual([]);
      v.input.move(l.labelX + 2, y);
      frames!.flush();
      expect(v.hover.current.target).toBe('row');
      expect(tips(v)).toEqual([]);
    });

    it('narrowed, the chip names the hidden columns in the title band, off the data; its card lists each and why', async () => {
      const v = await open(330);
      const l = lay(v);
      expect(l.cols.width).toEqual([0, 0]);
      expect(l.cols.chip).toBe('Species, Identity hidden');
      expect(l.cols.chipX + l.cols.chipW).toBeLessThanOrEqual(l.box.x);
      v.input.move(l.cols.chipX + 4, l.cols.titleY - 3);
      frames!.flush();
      expect(v.hover.current.target).toBe('columns-hidden');
      const words = new Describer({
        spec: v.spec,
        store: v.store,
        addresses: v.addresses,
        view: v.view,
        layout: () => v.frame?.layout,
      });
      const card = words.describe(v.hover.current);
      expect(card?.title).toBe('Hidden columns');
      expect(card?.rows).toEqual([
        { label: 'Species', value: 'no room at this width' },
        { label: 'Identity', value: 'no room at this width' },
      ]);
      expect(card?.text).toBe('Species, Identity hidden, no room at this width');
      // a click on it does nothing (no order, no selection)
      const revision = v.spec.revision.value;
      expect(v.input.click(l.cols.chipX + 4, l.cols.titleY - 3)).toBe(true);
      expect(v.spec.revision.value).toBe(revision);
    });
  });

  it('a move within a logo column, from one letter to another, keeps the hover (its letter is no hover state)', async () => {
    vi.useFakeTimers({ toFake: ['Date', 'performance'] });
    frames = installFrames();
    const b = createSeqSpec({ title: 'logo' });
    // column 0 half M, half A: two letters stacked
    const letters = ['MVLS', 'MVHL', 'MGLS', 'AKAV', 'AEKV', 'AERP'];
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters } });
    b.views.view('v', (vw) =>
      vw.section('s', { of: aln }, (s) => {
        s.track('logo', { label: '', height: 64 }, (t) => t.layer('logo', { type: 'logo', data: aln }));
        s.track('rows', { label: '', height: 140 }, (t) => t.layer('r', { type: 'rows', data: aln }));
      }),
    );
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
    viewer = v;
    await v.load(b.build());
    v.resize(480);
    await v.settle();
    frames.flush();
    const logo = v.frame!.layout.sections[0].tracks[0].box;
    const l = lay(v);
    l.axis.cell(0);
    const x = (l.axis.sx0 + l.axis.sx1) / 2;
    // two points in the column over different letters
    const states = new Map<string, number>();
    for (let y = logo.y + 1; y < logo.y + logo.height; y += 1) {
      const s = v.hit(x, y)?.state;
      if (s && !states.has(s)) states.set(s, y);
    }
    expect(states.size).toBe(2);
    const [y1, y2] = [...states.values()];
    v.input.move(x, y1);
    frames!.flush();
    const version = v.hover.version.value;
    expect(v.hover.current.position).toBe(0);
    v.input.move(x, y2);
    frames!.flush();
    expect(v.hover.version.value).toBe(version);
  });
});
