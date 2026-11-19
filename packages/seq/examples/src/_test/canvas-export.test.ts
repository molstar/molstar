/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInRender } from '@molstar/seqx-core/render/built-ins';
import { DataStore } from '@molstar/seqx-core/data';
import {
  attachSources,
  bundledTextMetrics,
  LightTokens,
  Renderer,
  recordDrawLists,
  renderConfig,
  renderSvg,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { build } from '../p53/tracks.js';

// The canvas path and the export path draw the same picture. The viewer lays out before its data arrives and then once
// with it (into the buffers the first layout made); the export lays out again before it draws. p53's 393 conservation
// values are more bars than a values layer starts with: its buffers grow while it fills, and once dropped the bars
// pushed before the growth (the canvas showed only the last 137).

const Size = { width: 700, height: 160 };

/** The data pass as the viewer's frames make it: a frame with no data, then the first (settled) one with it. */
async function canvasPath(dpr: number): Promise<string> {
  const model = resolveSpec(build(), {});
  const view = resolveViews(model)[0];
  const store = new DataStore();
  const renderer = new Renderer({
    ...BuiltInRender,
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    config: renderConfig(),
  });
  renderer.frame({ size: Size, dpr, settled: false });
  await attachSources(model, store);
  let frame = renderer.frame({ size: Size, dpr, settled: true });
  for (let i = 0; i < 4; i++) {
    store.want(frame.wants);
    await store.whenSettled();
    frame = renderer.frame({ size: Size, dpr, settled: true });
    if (!frame.rebuilt) break;
  }
  store.retire();
  return recordDrawLists(frame.lists, ['data']);
}

async function exportPath(dpr: number): Promise<string> {
  const r = await renderSvg(build(), { size: Size, dpr });
  expect(r.settled).toBe(true);
  return recordDrawLists(r.lists, ['data']);
}

/** The values track's bars: the last clip block's rects. */
function valueBars(record: string): number {
  const lines = record.trimEnd().split('\n');
  let clip = lines.length - 1;
  while (clip >= 0 && !lines[clip].startsWith('clip ')) clip--;
  return lines.slice(clip + 1).filter((l) => l.startsWith('rect ')).length;
}

describe('p53 tracks: the canvas path and the export path', () => {
  it.each([1, 2])('draw the same data pass at dpr %i, a conservation bar a residue', async (dpr) => {
    const canvas = await canvasPath(dpr);
    const svg = await exportPath(dpr);
    expect(valueBars(svg)).toBe(393);
    expect(valueBars(canvas)).toBe(393);
    expect(canvas).toBe(svg);
  });
});
