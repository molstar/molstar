/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tracks of one kind share the tallest lanes only when close (_spec/README.md "Render", Track heights; layouts.md
// "Tracks of one kind"): within 2 lanes or 2x of the cluster's smallest, clusters by sorted lanes, so never by order.

import { describe, expect, it } from 'vitest';
import { DefaultRenderConfig } from '../index.js';
import { ContentSizing, lanesHeight, type SizedTrack } from '../sizing.js';

const C = DefaultRenderConfig;

/** A track of `lanes` lanes (and a label lane when `label`), of a kind. */
function track(lanes: number, kind = 'ranges', label = 0, height?: number): SizedTrack {
  return {
    track: {} as never,
    kind,
    content: lanesHeight(lanes, C) + label,
    label,
    fixed: 22,
    ...(height !== undefined ? { height } : {}),
  };
}

const sized = (tracks: SizedTrack[]) => ContentSizing.heights(tracks, C);

describe('tracks of one kind share their lanes only when close', () => {
  it('one tall track and many short ones: the short ones fit their own lanes', () => {
    const h = sized([track(40), ...Array.from({ length: 30 }, () => track(1))]);
    expect(h[0]).toBe(lanesHeight(40, C));
    for (let i = 1; i < h.length; i++) expect(h[i]).toBe(lanesHeight(1, C));
  });

  it('two close tracks share the taller (within 2 lanes, or 2x)', () => {
    expect(sized([track(3), track(5)])).toEqual([lanesHeight(5, C), lanesHeight(5, C)]);
    expect(sized([track(10), track(20)])).toEqual([lanesHeight(20, C), lanesHeight(20, C)]);
  });

  it('far tracks fit their own (1 and 4 lanes: 3 lanes and 4x apart)', () => {
    expect(sized([track(1), track(4)])).toEqual([lanesHeight(1, C), lanesHeight(4, C)]);
  });

  it('tracks join a cluster within 2 lanes or 2x of its smallest: 1, 2 and 4 lanes are {1, 2} and {4}', () => {
    expect(sized([track(1), track(2), track(4)])).toEqual([lanesHeight(2, C), lanesHeight(2, C), lanesHeight(4, C)]);
  });

  it('1, 3 and 5 lanes are {1, 3} and {5}: 5 is not close to the cluster’s smallest, though close to 3', () => {
    expect(sized([track(1), track(3), track(5)])).toEqual([lanesHeight(3, C), lanesHeight(3, C), lanesHeight(5, C)]);
  });

  it('a fractional lane height does not flip a boundary: 2 and 4 lanes of 4.2 px share', () => {
    const cfg = { ...C, laneHeight: 4.2 };
    const t = (n: number): SizedTrack => ({ ...track(1), content: lanesHeight(n, cfg) });
    expect(ContentSizing.heights([t(2), t(4)], cfg)).toEqual([lanesHeight(4, cfg), lanesHeight(4, cfg)]);
  });

  it('a throwing sort does not keep the tracks alive: the next sizing is right', () => {
    const bad = {
      get kind(): string {
        throw new Error('x');
      },
      content: 30,
      fixed: 22,
      track: {} as never,
    } as SizedTrack;
    expect(() => sized([bad, track(2)])).toThrow();
    expect(sized([track(2), track(3)])).toEqual([lanesHeight(3, C), lanesHeight(3, C)]);
  });

  it('the result does not depend on the tracks’ order', () => {
    const lanes = [1, 1, 2, 3, 4, 7, 9, 14, 40, 3, 2];
    const base = lanes.map((n) => track(n));
    const expected = new Map<SizedTrack, number>();
    sized(base).forEach((h, i) => expected.set(base[i], h));
    let seed = 7;
    const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    for (let r = 0; r < 50; r++) {
      const shuffled = base.slice();
      for (let i = shuffled.length - 1; i > 0; i--) {
        const j = Math.floor(rnd() * (i + 1));
        [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
      }
      sized(shuffled).forEach((h, i) => expect(h).toBe(expected.get(shuffled[i])));
    }
  });

  it('kinds never share, and each track keeps its own label lane', () => {
    const h = sized([track(3, 'ranges', 14), track(5, 'ranges', 0), track(3, 'letters+ranges')]);
    expect(h).toEqual([lanesHeight(5, C) + 14, lanesHeight(5, C), lanesHeight(3, C)]);
  });

  it('a track without content keeps its fixed height', () => {
    const none: SizedTrack = { track: {} as never, kind: 'values', fixed: 40 };
    expect(sized([none, track(2)])).toEqual([40, lanesHeight(2, C)]);
  });
});
