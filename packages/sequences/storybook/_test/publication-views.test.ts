/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Seq*/Text view › Publication tracks story's two views (support/seqstar-annotations.tsx): Text and Tracks are two
// presets of one section (`ctx.ui.setStyle`), so switching rebuilds nothing: the document object is the same, every
// track is laid out and drawn in both, and an annotation added in one is in the other. The Annotations panel's two
// sets (Above and Below) each have their own ranges track, one before the letters and one after: an annotation of
// either is drawn in both views, above or below the letters' line, one undo step each.

import { type SeqSet } from '@molstar/seq-core/model';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  manualScheduler,
  type SeqViewerContext,
} from '@molstar/seq-viewer';
import { renderSvg } from '@molstar/seq-viewer/svg';
import { fixtureMeasurer } from '@molstar/seq-viewer/testing';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { seqLoci } from '../../viewer/_test/overlay-docs.js';
import {
  annotationRows,
  PUBLICATION_ANNOTATIONS,
  PUBLICATION_ANNOTATIONS_ABOVE,
  PUBLICATION_PLACES,
  PUBLICATION_PRESETS,
  PUBLICATION_SECTION,
  publicationTracksDoc,
} from '../src/_stories/support/seqstar-publication.js';

const TRACKS = ['annotations-above', 'letters', 'annotations-track', 'hyd', 'conf', 'zn', 'regions-track', 'ss'];
const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};
const open: SeqViewerContext[] = [];
afterEach(() => {
  for (const c of open.splice(0)) c.dispose();
});

async function story() {
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry] },
    { scheduler: manualScheduler(), text: fixtureMeasurer() },
  );
  open.push(ctx);
  const done = ctx.setDocument(publicationTracksDoc());
  ctx.viewports.get(PUBLICATION_SECTION)?.setSize(1000, 600, 1);
  expect((await done).ok).toBe(true);
  await ctx.graph.update();
  const layout = () => ctx.layout.get(PUBLICATION_SECTION)?.current;
  /** The refs of the tracks the section lays out (a track repeats on each line of the text). */
  const laid = () => {
    const view = ctx.graph.assembled.value?.view;
    const l = layout();
    const out = new Set<string>();
    for (let e = 0; view && l && e < l.entryCount; e++)
      if (l.entryTrack[e] >= 0) out.add(view.set.strings.get(view.resolved.tracks.ref[l.entryTrack[e]]));
    return [...out];
  };
  const show = async (view: 'text' | 'tracks') => {
    ctx.ui.setStyle(PUBLICATION_SECTION, PUBLICATION_PRESETS[view]);
    await ctx.graph.update();
    ctx.viewports.get(PUBLICATION_SECTION)?.setSize(1000, 600, 1);
  };
  return { ctx, layout, laid, show };
}

describe('the publication story: Text and Tracks are one document', () => {
  it('switching keeps the document and draws every track in both views', async () => {
    const { ctx, layout, laid, show } = await story();
    const set = vi.spyOn(ctx, 'setDocument');
    const doc = ctx.document.shown.value;
    // The text: wrapped lines, every track on them.
    expect(layout()?.geometry).toBe(1);
    expect(laid().sort()).toEqual([...TRACKS].sort());
    const text = await renderSvg(ctx, { width: 1000 });

    await show('tracks');
    expect(layout()?.geometry).toBe(0);
    expect(laid().sort()).toEqual([...TRACKS].sort());
    const tracks = await renderSvg(ctx, { width: 1000 });
    expect(tracks).not.toBe(text);
    // The tracks keep their heights: the values and the area are taller than the one-row tracks.
    const l = layout();
    const heights = Array.from({ length: l?.entryCount ?? 0 }, (_, e) => l?.entryHeight[e] ?? 0);
    expect(Math.max(...heights)).toBeGreaterThan(Math.min(...heights));

    await show('text');
    expect(layout()?.geometry).toBe(1);
    expect(laid().sort()).toEqual([...TRACKS].sort());
    // Nothing was rebuilt: the same document object, no document set again.
    expect(ctx.document.shown.value).toBe(doc);
    expect(set).not.toHaveBeenCalled();
  });

  it('an annotation added in one view is in the other, one undo step', async () => {
    const { ctx, show } = await story();
    const loci = seqLoci(ctx.graph.assembled.value?.set as SeqSet, 'p53', 100, 130);
    const rows = () => annotationRows(ctx.document.shown.value);
    expect(rows()).toEqual([]);
    expect(
      ctx.actions.run('bookmark.add', { section: PUBLICATION_SECTION, loci, label: 'Core', color: '#e76f51' }),
    ).toBe(true);
    await flush();
    const colored = (svg: string) => /e76f51|rgb\(231,\s*111,\s*81\)/i.test(svg);
    expect(colored(await renderSvg(ctx, { width: 1000 }))).toBe(true);
    await show('tracks');
    expect(colored(await renderSvg(ctx, { width: 1000 }))).toBe(true);
    expect(rows()).toHaveLength(1);
    expect(await ctx.history.undo()).toBe(true);
    expect(rows()).toEqual([]);
    expect(await ctx.history.undo()).toBe(false);
  });

  it('the panel offers Above (before the letters) and Below (after), each a set with its own track in both views', async () => {
    const { ctx, layout, show } = await story();
    expect(PUBLICATION_PLACES.map((p) => [p.label, p.data])).toEqual([
      ['Above', PUBLICATION_ANNOTATIONS_ABOVE],
      ['Below', PUBLICATION_ANNOTATIONS],
    ]);
    const set = () => ctx.graph.assembled.value?.set as SeqSet;
    const rowsOf = (data: string) => annotationRows(ctx.document.shown.value, data);
    /** The first top of a track's entries, and the letters', in the layout. */
    const tops = () => {
      const view = ctx.graph.assembled.value?.view;
      const l = layout();
      const out: Record<string, number[]> = {};
      for (let e = 0; view && l && e < l.entryCount; e++)
        if (l.entryTrack[e] >= 0)
          (out[view.set.strings.get(view.resolved.tracks.ref[l.entryTrack[e]])] ??= []).push(l.entryTop[e]);
      return out;
    };
    const colored = (svg: string, hex: string, rgb: string) => new RegExp(`${hex}|${rgb}`, 'i').test(svg);
    const [up, down] = PUBLICATION_PLACES;
    const base = await renderSvg(ctx, { width: 1000 });
    expect(colored(base, '2a9d8f', 'rgb\\(42,\\s*157,\\s*143\\)')).toBe(false);
    expect(
      ctx.actions.run('bookmark.add', {
        section: PUBLICATION_SECTION,
        loci: seqLoci(set(), 'p53', 100, 130),
        label: 'Up',
        color: '#e76f51',
        data: up.data,
      }),
    ).toBe(true);
    await flush();
    expect(
      ctx.actions.run('bookmark.add', {
        section: PUBLICATION_SECTION,
        loci: seqLoci(set(), 'p53', 200, 230),
        label: 'Down',
        color: '#2a9d8f',
        data: down.data,
      }),
    ).toBe(true);
    await flush();
    expect([rowsOf(up.data).length, rowsOf(down.data).length]).toEqual([1, 1]);
    for (const view of ['text', 'tracks'] as const) {
      await show(view);
      const svg = await renderSvg(ctx, { width: 1000 });
      expect(colored(svg, 'e76f51', 'rgb\\(231,\\s*111,\\s*81\\)'), `${view}: Above`).toBe(true);
      expect(colored(svg, '2a9d8f', 'rgb\\(42,\\s*157,\\s*143\\)'), `${view}: Below`).toBe(true);
      // Above's track is laid out over the letters', Below's under.
      const t = tops();
      expect(t['annotations-above'][0], view).toBeLessThan(t.letters[0]);
      expect(t['annotations-track'][0], view).toBeGreaterThan(t.letters[0]);
    }
    // A delete from one set leaves the other; one undo step each (two adds, one remove).
    expect(ctx.actions.run('bookmark.remove', { section: PUBLICATION_SECTION, item: 'b1', data: up.data })).toBe(true);
    await flush();
    expect([rowsOf(up.data).length, rowsOf(down.data).length]).toEqual([0, 1]);
    for (const expected of [
      [1, 1],
      [1, 0],
      [0, 0],
    ]) {
      expect(await ctx.history.undo()).toBe(true);
      expect([rowsOf(up.data).length, rowsOf(down.data).length]).toEqual(expected);
    }
    expect(await ctx.history.undo()).toBe(false);
  });

  it('a wavy annotation is drawn in both views, unlike the same one solid; a solid after it is fine', async () => {
    const { ctx, show } = await story();
    // The loci are made again after an undo (the set was built again).
    const add = (line_style?: 'wavy') =>
      ctx.actions.run('bookmark.add', {
        section: PUBLICATION_SECTION,
        loci: seqLoci(ctx.graph.assembled.value?.set as SeqSet, 'p53', 100, 130),
        color: '#e76f51',
        data: PUBLICATION_ANNOTATIONS_ABOVE,
        ...(line_style ? { line_style } : {}),
      });
    const both = async () => {
      const out: string[] = [];
      for (const view of ['text', 'tracks'] as const) {
        await show(view);
        out.push(await renderSvg(ctx, { width: 1000 }));
      }
      return out;
    };
    expect(add()).toBe(true);
    await flush();
    const solid = await both();
    expect(await ctx.history.undo()).toBe(true);
    await ctx.graph.update();
    expect(add('wavy')).toBe(true);
    await flush();
    expect(annotationRows(ctx.document.shown.value, PUBLICATION_ANNOTATIONS_ABOVE)).toMatchObject([
      { id: 'b1', line_style: 'wavy' },
    ]);
    const wavy = await both();
    for (const [k, view] of ['text', 'tracks'].entries()) {
      expect(wavy[k], `${view}: drawn`).toMatch(/e76f51|rgb\(231,\s*111,\s*81\)/i);
      expect(wavy[k], `${view}: not the solid one`).not.toBe(solid[k]);
    }
    expect(add()).toBe(true);
    await flush();
    expect(ctx.issues.current.value.filter((i) => i.severity === 'error' || i.severity === 'warning')).toEqual([]);
  });
});
