/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { recordDrawLists } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeViewer } from '../../selection/compose.js';
import {
  Channels,
  clearChannel,
  cycleStyle,
  describeChannels,
  HighlightConfig,
  setChannel,
  toggleChannel,
} from '../compose.js';

// The Highlights story, headless on a fake clock, composed as the story composes it: three overlapping channels set and
// cleared as the buttons do; the stripes' geometry and colors in the canvas's draw lists and in the SVG export beside it,
// and each style drawn.

const tracks = P53Views.find((v) => v.id === 'tracks')!;
const [Domain, Motif, Host] = Channels;

let viewer: SeqViewer | undefined;
const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function make() {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const v = composeViewer(tracks, { noRails: true, config: HighlightConfig });
  viewer = v;
  v.resize(900);
  await v.load(tracks.build());
  await vi.advanceTimersByTimeAsync(400);
  v.draw();
  return v;
}

const redraw = async (v: SeqViewer) => {
  await vi.advanceTimersByTimeAsync(100);
  v.draw();
};

/** The highlight pass's rects: `[x, y, w, h, paint]`. */
const rects = (lists: Parameters<typeof recordDrawLists>[0]) =>
  recordDrawLists(lists, ['highlight'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => {
      const p = l.split(' ');
      return { x: Number(p[1]), y: Number(p[2]), w: Number(p[3]), h: Number(p[4]), paint: p.slice(5).join(' ') };
    });

describe('Highlights, headless', () => {
  it('the three channels set one by one: stripes where they overlap, the whole band where one lights alone', async () => {
    const v = await make();
    await setChannel(v, Domain, 'tint');
    await redraw(v);
    const alone = rects(v.frame!.lists);
    expect(alone.length).toBe(1);
    const H = alone[0].h;
    await setChannel(v, Motif, 'tint');
    await redraw(v);
    const two = rects(v.frame!.lists);
    // domain alone, domain over motif, motif under domain, (the motif lies inside the domain: nothing more)
    expect(two.map((r) => r.h).sort((a, b) => a - b)).toEqual([H / 2, H / 2, H, H].filter((_, i) => i < two.length));
    expect(new Set(two.map((r) => r.y)).size).toBe(2);
    await setChannel(v, Host, 'tint');
    await redraw(v);
    const three = rects(v.frame!.lists);
    // three stripes where all of them light: the domain's on top, then the motif's, then the host's
    const triple = three.filter((r) => Math.abs(r.h - H / 3) < 0.01);
    expect(triple.length).toBe(3);
    expect(triple.map((r) => r.y)).toEqual([...triple.map((r) => r.y)].sort((a, b) => a - b));
    expect(new Set(triple.map((r) => r.paint)).size).toBe(3);
    expect(describeChannels(v)).toBe('channels: domain (tint), motif (tint), host (tint)');
    // the SVG draws the same stripes
    const svg = await v.renderSvg({ idPrefix: 'h-' });
    expect(rects(svg.lists)).toEqual(three);
    // cleared from the buttons: back to the whole band
    await toggleChannel(v, Motif);
    await toggleChannel(v, Host);
    await redraw(v);
    expect(rects(v.frame!.lists).length).toBe(1);
    await clearChannel(v, Domain);
    await redraw(v);
    expect(rects(v.frame!.lists)).toEqual([]);
  });

  it('each style is drawn: a tint, a bar along the foot, an outline; the next style is one button away', async () => {
    const v = await make();
    // the channels are apart on p53 for this: the domain and the host's own, set on their own
    await setChannel(v, Domain, 'tint');
    await setChannel(v, Host, 'outline');
    await redraw(v);
    const lists = rects(v.frame!.lists);
    expect(lists.some((r) => /^#[0-9a-f]{6}4d( squircle [\d.]+)?$/.test(r.paint))).toBe(true);
    expect(lists.some((r) => r.paint.startsWith('none stroke'))).toBe(true);
    await cycleStyle(v, Host);
    await redraw(v);
    expect(v.highlight.get(Host.name)?.style).toBe('tint');
    await cycleStyle(v, Host);
    await setChannel(v, Motif, 'band');
    await redraw(v);
    expect(v.highlight.get(Host.name)?.style).toBe('band');
    expect(rects(v.frame!.lists).some((r) => r.h === 3)).toBe(true);
    const svg = await v.renderSvg({ idPrefix: 's-' });
    expect(rects(svg.lists)).toEqual(rects(v.frame!.lists));
  });
});
