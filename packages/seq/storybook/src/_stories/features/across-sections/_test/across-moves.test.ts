/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LightTokens } from '@molstar/seqx-core/render';
import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { AcrossViews } from '@molstar/seqx-examples/across/views';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { openZoomed, openZoomedGrape } from '../openings.js';
import { bridges } from './across-checks.js';

// Features/Across sections after a side moves (FX28): a bridge's ribbons on the canvas against the SVG export's. The
// export lays out afresh from the store, the canvas from what it kept: a count the canvas keeps from an earlier window
// shows (an assembly read kept its rows of the window it was first built for).

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
  await import('@molstar/seqx-viewer/bridge-view');
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const canvas = (v: SeqViewer) => bridges(v).map((b) => b.count);

/** The export's ribbons a bridge, its clock run until it settles (-1: it never did). */
async function exported(v: SeqViewer): Promise<number[]> {
  const pending = v.renderSvg({ idPrefix: 'x-', tokens: LightTokens });
  let done = false;
  void pending.then(() => (done = true));
  for (let i = 0; i < 30 && !done; i++) await vi.advanceTimersByTimeAsync(2000);
  if (!done) return [-1];
  const r = await pending;
  return [...r.svg.matchAll(/data-count="(\d+)"/g)].map((m) => Number(m[1]));
}

async function settle(v: SeqViewer) {
  await vi.advanceTimersByTimeAsync(4000);
  v.draw();
}

/**
 * A deterministic run of pans and zooms of the bridge's two sides (`panPx`, the wheel's `zoom` about an offset), the
 * canvas's ribbon count against the export's after each; the steps that differ.
 */
async function sweep(id: string, steps: number, opening?: (v: SeqViewer) => unknown, seed = 7): Promise<string[]> {
  viewer = await openExample(AcrossViews.find((x) => x.id === id)!, { delayMs: 300, wait: 30_000 });
  const v = viewer;
  v.draw();
  await opening?.(v);
  await settle(v);
  const [bridge] = bridges(v);
  const sides = [bridge.a, bridge.b];
  let state = seed;
  const rnd = () => (state = (state * 1664525 + 1013904223) >>> 0) / 2 ** 32;
  const wrong: string[] = [];
  for (let i = 0; i < steps; i++) {
    const section = sides[rnd() < 0.7 ? 0 : 1];
    const w = v.view.get(section)!;
    let name: string;
    if (rnd() < 0.4) {
      const dx = (rnd() - 0.5) * 800;
      name = `${section} panPx ${Math.round(dx)}`;
      v.view.panPx(section, dx);
    } else {
      const f = rnd() < 0.5 ? 1.5 + rnd() * 3 : 1 / (1.5 + rnd() * 3);
      const at = w.start + rnd() * (w.end - w.start);
      name = `${section} zoom ${f.toFixed(2)} at ${Math.round(at)}`;
      v.view.zoom(section, f, at);
    }
    await settle(v);
    await settle(v);
    const svg = await exported(v);
    const c = canvas(v);
    if (JSON.stringify(c) !== JSON.stringify(svg)) wrong.push(`${i} ${name}: canvas ${c} svg ${svg}`);
  }
  return wrong;
}

/** Every side of every bridge in turn: a pan, a zoom in, a zoom out past where it began, a pan back; the steps that differ. */
async function everySide(id: string, opening?: (v: SeqViewer) => unknown): Promise<string[]> {
  viewer = await openExample(AcrossViews.find((x) => x.id === id)!, {
    delayMs: 300,
    wait: id === 'across-pending' ? 400 : 30_000,
  });
  const v = viewer;
  await vi.advanceTimersByTimeAsync(30_000);
  v.draw();
  await opening?.(v);
  await settle(v);
  const wrong: string[] = [];
  const sides = [...new Set(bridges(v).flatMap((b) => [b.a, b.b]))];
  for (const section of sides) {
    const moves: [string, () => void][] = [
      ['panPx 300', () => v.view.panPx(section, 300)],
      ['zoom 3', () => v.view.zoom(section, 3, mid(v, section))],
      ['zoom 1/6', () => v.view.zoom(section, 1 / 6, mid(v, section))],
      ['panPx -300', () => v.view.panPx(section, -300)],
      ['zoom out to the whole', () => v.view.zoom(section, 1e-3, mid(v, section))],
    ];
    for (const [name, move] of moves) {
      move();
      await settle(v);
      await settle(v);
      const svg = await exported(v);
      const c = canvas(v);
      if (JSON.stringify(c) !== JSON.stringify(svg)) wrong.push(`${id} ${section} ${name}: canvas ${c} svg ${svg}`);
    }
  }
  return wrong;
}

const mid = (v: SeqViewer, section: string) => {
  const w = v.view.get(section)!;
  return (w.start + w.end) / 2;
};

describe('Features/Across sections after a side moves (FX28)', () => {
  it('one side zoomed in: a zoom out and a pan leave the ribbons the export draws', async () => {
    expect(await sweep('across-zoomed', 3, openZoomed)).toEqual([]);
  });

  for (const id of [
    'across-fit',
    'across-zoomed',
    'across-bins',
    'across-stub',
    'across-pending',
    'across-self',
    'across-chain',
  ])
    it(`${id}: every side panned and zoomed, the canvas's ribbons the export's`, { tags: ['medium'] }, async () => {
      expect(await everySide(id)).toEqual([]);
    });

  // the zoomed openings: a side zoomed in, then zoomed out past the fit (the gate's other state)
  it('across-zoomed opened zoomed: every side moved', { tags: ['medium'] }, async () => {
    expect(await everySide('across-zoomed', openZoomed)).toEqual([]);
  });

  it('across-bins opened on grape zoomed: every side moved', { tags: ['medium'] }, async () => {
    expect(await everySide('across-bins', openZoomedGrape)).toEqual([]);
  });
});
