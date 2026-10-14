/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// UX1, the context's UI state (seqspec:viewer#state, seqspec:viewer#the-document): `setActive` (it cancels the press,
// finishes the viewports' animation, pauses the streams, the player and the input), tasks and issues, the runtime
// actions that join `DefaultActions` with their managers, the track menu asked for by the keyboard, `Ui.keep`, and the
// label providers' new reach (field providers, a renderer's entry). The `set:<ref>` channels from a section's
// `highlights` are LS1's.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { Loaders, Patch, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { tilesStream } from '@molstar/seq-core/stream/tiles';
import { tileServer } from '@molstar/seq-core/stream/tiles/server';
import { describe, expect, it, vi } from 'vitest';
import { case3, case10 } from '../../core/_test/format/cases.js';
import type { Documents } from '../src/document.js';
import {
  type ActivateEvent,
  createSeqViewer,
  DefaultLabelProviders,
  DefaultSeqViewerSpec,
  DragPhase,
  HoverKind,
  KernelRenderers,
  type LayerRenderer,
  manualScheduler,
  type SeqIssue,
} from '../src/index.js';
import type { Snapshots } from '../src/snapshots.js';
import { StarsFormat, StarsRenderer } from '../src/testing/index.js';
import { Ui } from '../src/ui.js';
import { starsDoc } from './conformance-fixtures.js';
import { manualLoaders } from './graph-fixtures.js';
import { driven } from './input-fixtures.js';
import { opened, preparedRenderers } from './layout-fixtures.js';

const S = 'family.main';
const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};
const patch = (ops: Patch['ops']): Patch => ({ format: 'seq_view_spec_patch', ops });
const rowY = (t: Awaited<ReturnType<typeof driven>>, d: number) => t.layout(S).scrolledTop + d * 18 + 9;

describe('setActive: a covered or hidden viewer', () => {
  it('a paused viewer emits nothing: keys and pointers are ignored until it is active again (ported)', async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    let keys = 0;
    let pointers = 0;
    t.ctx.interaction.key.subscribe(() => keys++);
    t.ctx.interaction.pointer.subscribe(() => pointers++);
    t.key(S, 'Escape', { typing: true });
    t.key(S, 'Escape');
    expect(t.ctx.ui.active.value).toBe(true);
    t.ctx.setActive(false);
    expect(t.ctx.ui.active.value).toBe(false);
    t.key(S, 'Escape');
    t.pointer('move', S, 400, rowY(t, 0));
    t.ctx.setActive(true);
    t.pointer('move', S, 400, rowY(t, 0));
    expect([keys, pointers]).toEqual([1, 1]);
    t.ctx.dispose();
  });

  it("cancels a press in progress: the drag ends as a cancel and the next move is a hover (IN2's wiring)", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const seen: string[] = [];
    t.ctx.interaction.drag.subscribe((e) => seen.push(`${e.phase}:${e.owner}`));
    const y = rowY(t, 0);
    t.pointer('down', S, 500, y);
    t.pointer('move', S, 480, y);
    expect(t.ctx.ui.drag.value?.owner).toBe('dragPan');
    const scroll = t.vp(S).live.current.scrollX;
    t.ctx.setActive(false);
    expect([t.ctx.ui.drag.value, seen.at(-1)]).toEqual([null, `${DragPhase.cancel}:dragPan`]);
    t.ctx.setActive(true);
    t.pointer('move', S, 470, y);
    expect([t.ctx.ui.drag.value, t.vp(S).live.current.scrollX]).toEqual([null, scroll]);
    t.ctx.dispose();
  });

  it("finishes the viewports' animations at once, and later moves jump instead of animating (VC3's item)", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const main = t.vp(S);
    main.fit();
    main.zoomCentered(28, 80, { animate: true, origin: 'behavior', durationMs: 200 });
    expect(main.live.current.animating).toBe(true);
    t.ctx.setActive(false);
    expect([main.live.current.animating, main.live.current.interacting, main.live.current.cell]).toEqual([
      false,
      false,
      28,
    ]);
    main.zoomCentered(25, 80, { animate: true, origin: 'behavior', durationMs: 200 });
    expect([main.live.current.animating, main.live.current.cell]).toEqual([false, 25]);
    t.ctx.dispose();
  });

  it('reaches the snapshot player; the streams ask for nothing while inactive, and for the view once active again', async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const spy = vi.spyOn(t.ctx.snapshots as unknown as Snapshots, 'setActive');
    t.ctx.setActive(false);
    t.ctx.setActive(true);
    expect(spy.mock.calls).toEqual([[false], [true]]);
    t.ctx.dispose();

    const URI = 'https://example.org/lazy';
    const calls: string[] = [];
    const server = tileServer({ rowCount: 2000, columnCount: 400 });
    const loaders: Loaders = {
      load: () => Promise.reject(new Error('no fixture')),
      stream: tilesStream((u) =>
        u === URI
          ? (path, signal) => {
              calls.push(path);
              return server.transport(path, signal);
            }
          : null,
      ),
    };
    const b = createSeqViewSpec({ title: 'lazy' });
    const aln = b.data.alignment({
      ref: 'lazy',
      row_count: 2000,
      column_count: 400,
      source: { uri: URI, protocol: 'seqstar-tiles/1' },
    });
    b.view({ ref: 'v', label: 'A', style: 'bare', height: 'box' })
      .block({ alignment: aln })
      .track({ label: 'Rows' })
      .rows();
    const l = await opened(b.getDocument(), {
      renderers: preparedRenderers(),
      loaders,
      size: { 'v.main': [800, 400] },
    });
    const tiles = () => calls.filter((c) => c.startsWith('tiles/')).length;
    l.ctx.setActive(false);
    l.scheduler.advance(16);
    await flush();
    expect(tiles()).toBe(0);
    // Active again, the next frame asks for the view's regions without a viewport change.
    l.ctx.setActive(true);
    l.scheduler.advance(16);
    await flush();
    expect(tiles()).toBeGreaterThan(0);
    l.ctx.dispose();
  });

  it('is refused after dispose, as the other verbs are', async () => {
    const ctx = createSeqViewer();
    ctx.dispose();
    expect(() => ctx.setActive(false)).toThrow(/after dispose/);
  });
});

describe('tasks and issues', () => {
  it('busy only after its delay, until the work ends (ported)', async () => {
    const io = manualLoaders();
    const scheduler = manualScheduler();
    const ctx = createSeqViewer(undefined, { loaders: io.loaders, scheduler });
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 's', source: { uri: 'https://x.test/slow.fa', format: 'fasta' } });
    b.view({ ref: 'v', label: 'V' }).track().letters();
    const done = ctx.setDocument(b.getDocument());
    await flush();
    const p = ctx.tasks.progress.value;
    expect(p.pending).toBeGreaterThan(p.done);
    scheduler.advance(749);
    expect(ctx.tasks.busy.value).toBe(false);
    scheduler.advance(1);
    expect(ctx.tasks.busy.value).toBe(true);
    io.resolve('https://x.test/slow.fa', '>s\nMKV\n');
    await done;
    expect([ctx.tasks.busy.value, ctx.tasks.progress.value]).toEqual([false, { pending: 0, done: 0 }]);
    expect(scheduler.pendingTimers).toBe(0);
    ctx.dispose();
  });

  it("issues.current: the graph's by node first, the others after, each once; a new document clears the others", async () => {
    const t = await opened(case10(), { renderers: preparedRenderers() });
    const { ctx } = t;
    const codes = () => ctx.issues.current.value.map((i) => `${i.ref}:${i.code}`);
    const before = codes();
    ctx.report('seqv.test', 'one');
    ctx.report('seqv.test', 'one');
    ctx.report('seqv.test', 'two');
    expect(codes().slice(before.length)).toEqual(['null:seqv.test', 'null:seqv.test']);
    expect(ctx.issues.current.value.filter((i) => i.ref === null).map((i) => i.message)).toEqual(['one', 'two']);
    await ctx.setDocument(case3());
    expect(ctx.issues.current.value.filter((i) => i.code === 'seqv.test')).toEqual([]);
    ctx.dispose();
  });

  it("a preview's `svs.snapshot.invalid` warnings are returned and listed in issues.current while it is shown (VC2's item)", async () => {
    const doc: SeqViewSpecDoc = case3();
    doc.snapshots = [{ ref: 'drops', ops: [{ op: 'set', ref: 'family.logo', params: { label: 'x' } }] }];
    const { ctx } = await opened(doc, { renderers: preparedRenderers() });
    const listed = () => ctx.issues.current.value.filter((i) => i.code === 'svs.snapshot.invalid').length;
    expect(listed()).toBe(0);
    const r = await ctx.applyPatch(patch([{ op: 'remove', ref: 'family.logo' }]), { layer: 'preview' });
    expect(r.ok).toBe(true);
    expect(r.issues.map((i) => i.code)).toContain('svs.snapshot.invalid');
    expect(listed()).toBe(1);
    // The same preview again lists it once; `setPreview` goes the same way; an invalid preview changes nothing.
    await ctx.document.setPreview(patch([{ op: 'remove', ref: 'family.logo' }]));
    expect(listed()).toBe(1);
    const bad = await ctx.applyPatch(patch([{ op: 'remove', ref: 'nothing' }]), { layer: 'preview' });
    expect([bad.ok, listed()]).toEqual([false, 1]);
    // Cleared, or replaced by a preview that breaks nothing, the warning leaves the list.
    await ctx.document.setPreview(null);
    expect(listed()).toBe(0);
    await ctx.applyPatch(patch([{ op: 'remove', ref: 'family.logo' }]), { layer: 'preview' });
    await ctx.applyPatch(patch([{ op: 'set', ref: 'family.main', params: { label: 'M' } }]), { layer: 'preview' });
    expect(listed()).toBe(0);
    ctx.dispose();
  });
});

describe("the runtime actions join DefaultActions with their managers (VC2's item)", () => {
  it('section.style, group.toggle, clade.toggle, viewport.fit and stream.retry; `when` gates them (ported)', async () => {
    const t = await driven(case10());
    expect(t.ctx.actions.can('track.hide', 'chains.nothing')).toBe(false);
    expect(t.ctx.actions.run('lab:none', null)).toBe(false);
    expect(t.ctx.actions.run('section.style', { section: 'chains.main', value: 'bare' })).toBe(true);
    expect(t.ctx.ui.style.value.get('chains.main')).toBe('bare');
    expect(t.ctx.actions.run('section.style', { section: 'chains.main', value: null })).toBe(true);
    expect(t.ctx.ui.style.value.has('chains.main')).toBe(false);
    expect(t.ctx.actions.can('section.style', { section: 'nope', value: 'bare' })).toBe(false);
    expect(t.ctx.actions.run('clade.toggle', { alignment: 'sh2', node: 1 })).toBe(true);
    expect(Array.from(t.ctx.ui.clades.value.get('sh2') ?? [])).toEqual([1]);
    expect(t.ctx.actions.can('clade.toggle', { alignment: 'sh2', node: -1 })).toBe(false);
    expect(t.ctx.actions.can('group.toggle', 'chains.main')).toBe(false);
    expect(t.ctx.actions.can('stream.retry', 'chains.main')).toBe(false);
    t.ctx.dispose();
  });

  it('group.toggle toggles a group by ref; viewport.fit fits one section, or all', async () => {
    const f = await driven(case3(), { size: { [S]: [1000, 150] } });
    const main = f.vp(S);
    main.zoomTo(0, 5, 15);
    const zoomed = main.live.current.cell;
    expect(f.ctx.actions.run('viewport.fit', S)).toBe(true);
    const fitted = main.live.current.cell;
    expect(fitted).not.toBe(zoomed);
    main.zoomTo(0, 5, 15);
    expect(f.ctx.actions.run('viewport.fit', undefined)).toBe(true);
    expect(main.live.current.cell).toBe(fitted);
    f.ctx.dispose();

    const b = createSeqViewSpec({ title: 'Group' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR' });
    const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
    v.group({ ref: 'v.g', label: 'G', collapsed: true }).track({ ref: 'v.t', label: 'T' }).letters();
    const { ctx } = await opened(b.getDocument(), { renderers: preparedRenderers() });
    expect([ctx.actions.can('group.toggle', 'v.g'), ctx.actions.can('group.toggle', 'v.t')]).toEqual([true, false]);
    expect(ctx.actions.run('group.toggle', 'v.g')).toBe(true);
    expect([...ctx.ui.expanded.value]).toEqual(['v.g']);
    expect(ctx.actions.run('group.toggle', 'v.g')).toBe(true);
    expect(ctx.ui.expanded.value.size).toBe(0);
    ctx.dispose();
  });
});

describe('the track menu asked for by the keyboard (events.menu)', () => {
  async function focused() {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    t.focus(S);
    const asked: ActivateEvent[] = [];
    t.ctx.events.menu.subscribe((e) => asked.push({ ...e }));
    return { t, asked };
  }

  it("Shift+F10 and ContextMenu on the cursor's track ask for its menu, with the keyboard as the origin", async () => {
    const { t, asked } = await focused();
    t.key(S, 'ArrowDown');
    const track = t.ctx.ui.keyCursor.value?.track;
    expect(track).toBeTruthy();
    expect(t.key(S, 'F10', { shift: true })).toBe(true);
    expect(t.key(S, 'ContextMenu')).toBe(true);
    expect(asked).toEqual([
      { section: S, track, origin: 'keyboard', x: -1, y: -1 },
      { section: S, track, origin: 'keyboard', x: -1, y: -1 },
    ]);
    t.ctx.dispose();
  });

  it("F10 alone is the page's; nothing is asked of a cursor on a group's header", async () => {
    const { t, asked } = await focused();
    t.key(S, 'ArrowDown');
    expect(t.key(S, 'F10')).toBe(false);
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });
});

describe('Ui.keep: per-section state follows the shown document', () => {
  it('a cursor state of a section that leaves completes after the update; one that returns within it keeps its object', async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const c = t.ctx.ui.cursor(S);
    let done = false;
    c.subscribe({ complete: () => (done = true) });
    // A document without the section, then back, within one update: the same state.
    const documents = t.ctx.document as unknown as Documents;
    documents.reset();
    documents.replace(case3(), undefined, []);
    await flush();
    expect([t.ctx.ui.cursor(S) === c, done]).toEqual([true, false]);
    await t.ctx.setDocument(case10());
    await flush();
    expect(done).toBe(true);
    t.ctx.dispose();
  });
});

describe("labels: field providers and a renderer's own entry (IN3's, XR1's items)", () => {
  const STARS = { renderers: [...KernelRenderers, StarsRenderer as LayerRenderer], registry: [{ svs: StarsFormat }] };

  it("the stars renderer's label entry writes its affinity line for the star under the pointer", async () => {
    const t = await opened(starsDoc(), STARS);
    const m = t.ctx.graph.assembled.value?.view;
    const k = m?.prepared.findIndex((p) => p.renderer === 'org.example:stars') ?? -1;
    const ax = t.ctx.layout.axis('v.main');
    const l = t.ctx.layout.get('v.main')?.current;
    const e = l ? Array.from(l.entryLayer.subarray(0, l.entryCount)).indexOf(k) : -1;
    const y = l ? l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2 : 0;
    const at = (col: number) => {
      const s = ax?.shownOf(0, col) ?? 0;
      return ax ? (ax.x(0, s as never) + ax.x(0, (s + 1) as never)) / 2 : 0;
    };
    // Pointer events through the dispatcher need the behaviors; the hit itself is enough for the label.
    const hit = (await import('../src/index.js')).Hit.create();
    t.ctx.layout.hitTest('v.main', at(12), y, hit);
    t.ctx.hover.setFromHit(hit, 'pointer', 'element');
    const label = t.ctx.labels.compute(t.ctx.hover.live.current);
    expect(label?.title).toBe('s2');
    expect(label?.lines.find((x) => x.key === 'affinity')).toEqual({ key: 'affinity', value: '8.0 pKd' });
    t.ctx.dispose();
  });

  it("the item provider's data column reads a field provider's values, as the data's own fields", async () => {
    const t = await opened(starsDoc(), STARS);
    const view = t.ctx.graph.assembled.value?.view;
    const item = DefaultLabelProviders.find((p) => p.name === 'item');
    if (!view || !item) throw new Error('no view');
    const column = (field: string) => [{ data_columns: [{ data: 'sites', field, label: 'pKd' }] }];
    const lines: string[][] = [];
    const out = { title: () => {}, line: (k: string, v: string) => void lines.push([k, v]), replaceTitle: () => {} };
    const label = (field: string) => {
      const layers = { ...view.resolved.layers, params: column(field) };
      const lc = {
        set: view.set,
        view: { ...view, resolved: { ...view.resolved, layers } },
        section: 0,
        registry: t.ctx.registry,
      };
      const h = { kind: HoverKind.row, dataColumn: 0, layer: 0, row: { row: 2, alignment: 0, index: 0 } };
      item.label(h as never, lc as never, out);
    };
    label('org.example:affinity');
    label('kd_nm');
    label('org.example:unknown');
    // Row 2 holds a K_d of 100 nM: pKd 7 from the provider, 100 from the column itself, nothing for an unknown field.
    expect(lines).toEqual([
      ['pKd', '7'],
      ['pKd', '100'],
    ]);
    t.ctx.dispose();
  });

  it('the `layer` provider (priority 75) runs after `item` and before `value`', async () => {
    const { ctx } = await opened(starsDoc(), STARS);
    expect(ctx.labels.providers.value).toEqual([
      'position',
      'copies',
      'row',
      'codon',
      'item',
      'layer',
      'value',
      'pair',
      'column',
      'stream',
    ]);
    ctx.dispose();
  });
});

describe('issues from the engine reach the list through report()', () => {
  it('an engine report is an issue with no ref, in the list and on the event', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    const seen: SeqIssue[] = [];
    ctx.events.issue.subscribe((i) => seen.push(i));
    ctx.report('seqv.renderer.error', 'x failed');
    expect([seen.map((i) => i.code), ctx.issues.current.value.map((i) => i.code)]).toEqual([
      ['seqv.renderer.error'],
      ['seqv.renderer.error'],
    ]);
    ctx.dispose();
  });
});

describe('open legend groups across a sync (VB-05)', () => {
  it("a template instance's group (`<section>/<track>/<node ref>`) stays while its section is shown", async () => {
    const ui = new Ui();
    ui.toggleLegend('main/genes/hg38/T1');
    ui.toggleLegend('main/plain');
    ui.toggleLegend('gone/plain');
    ui.keep(['main']);
    await Promise.resolve();
    await Promise.resolve();
    expect([...ui.legendOpen.value].sort()).toEqual(['main/genes/hg38/T1', 'main/plain']);
    ui.dispose();
  });
});
