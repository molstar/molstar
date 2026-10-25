/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Canvas2D backend's paths, by the calls it makes (seqspec:viewer#the-render-backend): Chrome's closePath costs as
// much as its path holds (3J3Q's unmerged secondary structure, 20,000 arrows a frame: seconds a drag), so a fill leaves
// out every close its fill makes itself, a stroke keeps its closes and is drawn in chunks, and a close followed by a
// line within one subpath stays. The pixels are the Chromium suite's (`canvas/backend.test.ts`).

import { beforeEach, describe, expect, it, vi } from 'vitest';
import { Canvas2DBackend, Canvas2DShared, type Paint2D } from '../src/canvas/index.js';
import { type DrawResources, Pass, PathBatch, PathMode, SkinToken } from '../src/draw/api.js';

const resources: DrawResources = {
  font: () => '20px sans-serif',
  fontPx: () => 20,
  token: (t) => (t === SkinToken.ink ? '#112233' : '#ff00ff'),
  glyph: () => null,
};

/** A 2D context that records its calls into `calls` by name (properties set are taken and ignored). */
function context(calls: string[]) {
  return new Proxy({} as Record<string, unknown>, {
    get: (t, k) =>
      k in t ? t[k as string] : (...args: unknown[]) => (calls.push(String(k)), args.length ? undefined : 0),
    set: (t, k, v) => ((t[k as string] = v), true),
  });
}
// The shared state's scratch canvas, in Node.
beforeEach(() => {
  vi.stubGlobal(
    'OffscreenCanvas',
    class {
      constructor(
        public width: number,
        public height: number,
      ) {}
      getContext() {
        return context([]);
      }
    },
  );
});

function recorder() {
  const calls: string[] = [];
  const g = context(calls);
  const be = new Canvas2DBackend([g as unknown as Paint2D, null], new Canvas2DShared(resources, { atlas: false }), 0);
  return { be, calls };
}

/** `n` closed triangles in one style. */
function triangles(n: number, mode: number = PathMode.fill, alpha = 255): PathBatch {
  const p = new PathBatch();
  for (let i = 0; i < n; i++) {
    p.begin(0xff0000, alpha, mode, 2);
    p.moveTo(i, 0);
    p.lineTo(i + 5, 0);
    p.lineTo(i, 5);
    p.close();
  }
  return p;
}

const count = (calls: string[], name: string) => calls.filter((c) => c === name).length;

function drawn(p: PathBatch) {
  const { be, calls } = recorder();
  be.beginPass(Pass.main, 100, 100, 1);
  const from = calls.length;
  be.paths(p);
  return calls.slice(from);
}

describe('canvas paths: no closePath a draw does not need', () => {
  it('a fill of 1,000 closed subpaths is one path, closed by its fill alone', () => {
    const calls = drawn(triangles(1000));
    expect([
      count(calls, 'beginPath'),
      count(calls, 'moveTo'),
      count(calls, 'closePath'),
      count(calls, 'fill'),
    ]).toEqual([1, 1000, 0, 1]);
  });

  it('a close followed by a line in its subpath stays (the line starts from the closed start)', () => {
    const p = new PathBatch();
    p.begin(0xff0000);
    p.moveTo(0, 0);
    p.lineTo(5, 0);
    p.lineTo(0, 5);
    p.close();
    p.lineTo(5, 5);
    p.lineTo(0, 9);
    p.close();
    expect(count(drawn(p), 'closePath')).toBe(1);
  });

  it('an opaque stroke keeps its closes, drawn every 256 subpaths; a translucent one in one path', () => {
    const opaque = drawn(triangles(1000, PathMode.stroke));
    expect([count(opaque, 'closePath'), count(opaque, 'beginPath'), count(opaque, 'stroke')]).toEqual([1000, 4, 4]);
    const translucent = drawn(triangles(1000, PathMode.stroke, 128));
    expect([count(translucent, 'closePath'), count(translucent, 'stroke')]).toEqual([1000, 1]);
  });
});
