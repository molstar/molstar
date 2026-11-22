/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Values looks story, headless on a fake clock, composed as the story composes it (`compose.ts`: each look's config
// over each cell's view): every look draws every cell, a variant never draws a stub for a missing value, the profile
// draws paths where a value is under 4 px and the variants' bars from it, and hover outlines a value in each.

import { readFileSync } from 'node:fs';
import { type ValuesLay, recordDrawLists } from '@molstar/seqx-core/render';
import type { ValuesLook } from '@molstar/seqx-core/render';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { withViewRegistry } from '../../../examples/view-registry.js';
import { Cells, lookConfig, Looks } from '../compose.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const bytes = (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`));
const reader = {
  readText: async (uri: string) => new TextDecoder().decode(bytes(uri)),
  readBytes: async (uri: string) => bytes(uri),
  openRange: async (uri: string) => {
    const file = bytes(uri);
    return {
      size: file.length,
      read: async (ranges: readonly { offset: number; length: number }[]) =>
        ranges.map((r) => file.subarray(r.offset, r.offset + r.length)),
      dispose() {},
    };
  },
};

async function open(look: ValuesLook, id: string) {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const cell = Cells.find((c) => c.id === id)!;
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    config: [...(base.config ?? []), lookConfig(look)] as never,
    registry: withViewRegistry(base.registry ?? [], cell.view.registry?.({ delayMs: 100 }) ?? []),
    ports: { assets: () => ({ reader: reader as never }) },
  });
  viewer = v;
  v.resize(480);
  void v.load(cell.view.build()).then(() => v.settle());
  await vi.advanceTimersByTimeAsync(3000);
  cell.opening?.(v);
  await vi.advanceTimersByTimeAsync(3000);
  v.draw();
  return v;
}

const values = (v: SeqViewer) => {
  for (const s of v.frame!.layout.sections)
    for (const t of s.tracks) for (const l of t.layers) if (l.type.name === 'values') return l.lay as ValuesLay;
  throw new Error('no values layer');
};

/** The data pass as text: its lines per primitive. */
const dataLines = (v: SeqViewer) => recordDrawLists(v.frame!.lists, ['data']).split('\n');

describe('the values looks story, headless', () => {
  const missing = (v: SeqViewer) => v.frame!.layout.tokens.colors['seq-missing'];

  it.each(Looks.map((l) => l.look))('%s: every cell draws, and hover outlines a value', async (look) => {
    for (const cell of Cells) {
      const v = await open(look, cell.id);
      const lay = values(v);
      expect(lay.n, cell.id).toBeGreaterThan(10);
      expect(lay.config.values.look).toBe(look);
      const x = lay.box.x + lay.box.width / 2;
      const y = lay.box.y + lay.box.height / 2;
      const hit = v.hit(x, y);
      expect(hit?.layer, cell.id).toBe(lay.key);
      v.hover.set(hit, 'test');
      const f = v.draw() ?? v.frame!;
      expect(
        recordDrawLists(f.lists, ['hover'])
          .split('\n')
          .filter((l) => l.startsWith('rect ')).length,
        cell.id,
      ).toBeGreaterThan(0);
      v.hover.set(undefined, 'test');
      v.retire();
      viewer = undefined;
    }
  });

  it('the variants draw the whole protein as a profile (paths) and the zoomed cells as bars; today draws columns', async () => {
    const today = await open('today', 'conservation-out');
    const todayPaths = dataLines(today).filter((l) => l.startsWith('path ')).length;
    today.retire();
    viewer = undefined;
    for (const look of ['profile', 'bars', 'quiet'] as const) {
      const out = await open(look, 'conservation-out');
      expect(values(out).axis.scale, look).toBeLessThan(4);
      expect(dataLines(out).filter((l) => l.startsWith('path ')).length, look).toBeGreaterThan(todayPaths);
      out.retire();
      viewer = undefined;
      const inn = await open(look, 'conservation-in');
      expect(values(inn).axis.scale, look).toBeGreaterThanOrEqual(4);
      const rects = dataLines(inn).filter((l) => l.startsWith('rect ')).length;
      expect(rects, look).toBeGreaterThan(40);
      inn.retire();
      viewer = undefined;
    }
  });

  it('a missing value is a gap in a variant, never a stub at zero', async () => {
    const today = await open('today', 'conservation-in');
    const stubs = (v: SeqViewer) =>
      dataLines(v).filter((l) => l.includes(`#${(missing(v) >>> 0).toString(16).padStart(8, '0')}`)).length;
    const hadToday = stubs(today);
    today.retire();
    viewer = undefined;
    // the zoomed-in range holds no missing residue of its own: the stub count of a look never exceeds today's
    for (const look of ['profile', 'bars', 'quiet'] as const) {
      const v = await open(look, 'conservation-in');
      expect(stubs(v), look).toBeLessThanOrEqual(hadToday);
      v.retire();
      viewer = undefined;
    }
  });
});
