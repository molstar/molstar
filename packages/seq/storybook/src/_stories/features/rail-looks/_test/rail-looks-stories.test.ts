/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Rail looks story, headless on a fake clock, composed as the story composes it (`compose.ts`: each look's config
// over each cell's view): every look draws the one rail path (a capsule's rounded rects in the variants, none in
// today's), insets its ends by the well's corner radius, and keeps hover and drag working; the rows' and the lanes'
// rails take the same look.

import { readFileSync } from 'node:fs';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import type { RailTint } from '@molstar/seqx-core/render';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { withViewRegistry } from '../../../examples/view-registry.js';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { Cells, Looks, lookConfig } from '../compose.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/**
 * The GFF3 lanes cell as the GFF stories' twin opens it: the runtime's opens awaited (a gzip file inflates by no clock),
 * then the twin's delays on the fake clock.
 */
async function openLanes(look: RailTint): Promise<SeqViewer> {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const view = cell('lanes').view;
  const opens: Promise<unknown>[] = [];
  const watched = (view.registry?.({ delayMs: 150 }) ?? []).map(
    (e): SeqRegistryEntry =>
      e.kind === 'format-runtime'
        ? {
            ...e,
            runtime: {
              ...e.runtime,
              open: (source, ctx) => {
                const p = e.runtime.open(source, ctx);
                opens.push(p);
                return p;
              },
            } satisfies FormatRuntime,
          }
        : e,
  );
  const base = DefaultSeqViewerSpec();
  const bytes = async (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`));
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    config: [...(base.config ?? []), lookConfig(look)] as never,
    registry: withViewRegistry(base.registry ?? [], watched),
    ports: {
      assets: () => ({ reader: { readText: async (u) => new TextDecoder().decode(await bytes(u)), readBytes: bytes } }),
    },
  });
  viewer = v;
  v.resize(900);
  void v.load(view.build()).then(() => v.settle());
  await vi.advanceTimersByTimeAsync(200);
  for (let seen = -1; seen !== opens.length; ) {
    seen = opens.length;
    await Promise.allSettled(opens);
    await vi.advanceTimersByTimeAsync(200);
  }
  await vi.advanceTimersByTimeAsync(5000);
  v.draw();
  return v;
}

const cell = (id: string) => Cells.find((c) => c.id === id)!;

async function open(look: RailTint, id: string) {
  const c = cell(id);
  viewer = await openExample(c.view, { config: [lookConfig(look)] as never });
  c.opening?.(viewer);
  viewer.draw();
  return viewer;
}

/** The window pass's paths: the knobs' discs and their rings. */
function discs(v: SeqViewer): number {
  return v.frame!.lists.passes.window.paths.n;
}

describe('the rail looks story, headless', () => {
  it.each(Looks.map((l) => l.look))('%s: the ruler rail, zoomed in; hover grows it and a drag pans', async (look) => {
    const v = await open(look, 'ruler-zoomed');
    const r = v.frame!.layout.sections[0].rail!;
    expect(r.active).toBe(true);
    // ends inset by the well's corner radius in the variants, not in today's look
    expect(r.inset).toBe(v.context.renderConfig.squircleMax);
    expect(discs(v)).toBeGreaterThan(0);
    const k = r.width / r.total;
    const y = r.bottom - 2;
    const mid = r.x + ((r.start + r.end) / 2) * k;
    v.input.move(mid, y);
    expect(v.input.windowHit?.part).toBe('body');
    const before = v.view.get('/views/main/protein')!.start;
    v.input.down(mid, y);
    v.input.move(mid + 20 * k, y + 60);
    v.input.up();
    v.draw();
    expect(v.view.get('/views/main/protein')!.start).toBeGreaterThan(before);
  });

  it.each(Looks.map((l) => l.look))('%s: the ruler rail at the fit stays a target at its ends', async (look) => {
    const v = await open(look, 'ruler-fit');
    const r = v.frame!.layout.sections[0].rail!;
    expect(r.start).toBe(0);
    expect(r.end).toBe(r.total);
    // the whole axis is the window: the grips still resize it
    v.input.move(r.x + 0.5, r.bottom - 2, { pointerType: 'touch' });
    expect(v.input.windowHit?.part).toBe('start');
  });

  it.each(Looks.map((l) => l.look))('%s: the rows rail and the lanes band take the same look', async (look) => {
    const rows = await open(look, 'rows');
    const rr = rows.frame!.layout.sections[0].rowsRails![0];
    expect(rr.vertical).toBe(true);
    expect(rr.inset).toBe(rows.context.renderConfig.squircleMax);
    expect(discs(rows)).toBeGreaterThan(0);
    rows.retire();
    const lanes = await openLanes(look);
    const lr = lanes.frame!.layout.sections.flatMap((x) => x.lanesRails ?? [])[0];
    expect(lr).toBeDefined();
    expect(lr!.inset).toBe(lanes.context.renderConfig.squircleMax);
  });
});
