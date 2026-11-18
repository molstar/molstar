/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { colorHex, type Frame, type RangesLay, recordDrawLists } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';

// The Gene by region story, headless on a fake clock: the viewer composed as `example-root.ts` composes it, with the
// story's SVG exports reading through the viewer's store (after the first settle, and again each time the view rests
// or the layout's issues move). The canvas shows each transcript's domains (a nested source the overlay opens, one
// level down) once their answer has arrived, never only after an export has let go of their source.

/** Each answer's delay: the region's, then the domains' (opened when their transcript's region answers). */
const DelayMs = 300;
const FrameMs = 16;

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** The domain tracks of a frame: their pending placeholders and marks as laid out, and as drawn (the data pass). */
function domains(frame: Frame, pending: string) {
  const tracks = frame.layout.sections.flatMap((s) => s.tracks).filter((t) => t.track.address.includes('/domains'));
  const rects = recordDrawLists(frame.lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => l.split(' '));
  return tracks.map((t) => {
    const lay = t.layers[0].lay as RangesLay;
    const inBox = rects.filter((r) => {
      const y = Number(r[2]);
      return y >= t.box.y && y < t.box.y + t.box.height;
    });
    return {
      track: t.track.label,
      pending: lay.pending.length,
      marks: lay.m,
      drawnPending: inBox.filter((r) => r[5] === pending).length,
      drawnMarks: inBox.filter((r) => r[5] !== pending).length,
    };
  });
}

describe('Gene by region, on the canvas', () => {
  it('draws the transcripts’ domains one frame after their answer, while the story’s exports share its store', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
    const view = P53Views.find((v) => v.id === 'gene-by-region')!;
    const text = regionSeqText();
    const base = DefaultSeqViewerSpec();
    viewer = createSeqViewer({
      ...base,
      ui: undefined,
      registry: [...(base.registry ?? []), ...(view.registry?.({ delayMs: DelayMs }) ?? [])],
      ports: { assets: () => ({ reader: { readText: async () => text } }) },
    });
    const v = viewer;
    v.resize(720);
    // the story's SVG: exports through the viewer's store, each with its own want set and overlay sources
    const exports: Promise<unknown>[] = [];
    const svg = () => exports.push(v.renderSvg());
    v.rested.listen(svg);
    v.layoutIssues.listen(svg);
    void v
      .load(view.build())
      .then(() => v.settle())
      .then(svg);

    // the region answers, its transcripts' domains are opened and asked for, and answer: then one frame
    await vi.advanceTimersByTimeAsync(2 * DelayMs + 4 * FrameMs);
    const pending = colorHex(v.tokens.colors['seq-pending']);
    const now = domains(v.frame!, pending);
    expect(now.map((d) => d.track)).toEqual(['Domains T1-domains', 'Domains T2-domains', 'Domains T3-domains']);
    for (const d of now) {
      expect(d, d.track).toMatchObject({ pending: 0, drawnPending: 0 });
      expect(d.marks, d.track).toBeGreaterThan(0);
      expect(d.drawnMarks, d.track).toBeGreaterThan(0);
    }
    expect(exports.length).toBeGreaterThan(0);

    // and they stay once every export has ended
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    await Promise.all(exports);
    await vi.advanceTimersByTimeAsync(FrameMs);
    for (const d of domains(v.frame!, pending)) expect(d, d.track).toMatchObject({ pending: 0, drawnPending: 0 });
  });

  it('draws the story’s SVG again once the picture rests after data arrives, at the canvas’s height', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
    const view = P53Views.find((v) => v.id === 'gene-by-region')!;
    const text = regionSeqText();
    const base = DefaultSeqViewerSpec();
    viewer = createSeqViewer({
      ...base,
      ui: undefined,
      registry: [...(base.registry ?? []), ...(view.registry?.({ delayMs: DelayMs }) ?? [])],
      ports: { assets: () => ({ reader: { readText: async () => text } }) },
    });
    const v = viewer;
    v.resize(720);
    // the story's SVG: the last export drawn wins
    let last: { svg: string; height: number } | undefined;
    let latest = 0;
    const svg = async () => {
      const mine = ++latest;
      const height = v.height;
      const r = await v.renderSvg();
      if (mine === latest) last = { svg: r.svg, height };
    };
    v.rested.listen(() => void svg());
    // the first export before any frame has laid the view out (a page whose frames lag): nothing to draw yet
    void v.load(view.build()).then(svg);

    await vi.advanceTimersByTimeAsync(20 * DelayMs);
    expect(v.height).toBeGreaterThan(32);
    expect(last?.height).toBe(v.height);
    expect(last?.svg).toContain(`height="${v.height}"`);
    expect(last?.svg).toContain('T3-domains');
  });
});
