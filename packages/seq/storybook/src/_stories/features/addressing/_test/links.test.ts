/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import { Passages } from '@molstar/seqx-examples/p53/addressing';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeViewer, describeRun, followLink, hoverLink, leaveLink } from '../compose.js';

// The Links story, headless on a fake clock, composed as the story composes it: each link of the text run as a click
// runs it; a link's hover on its own channel; one residue in two numberings focused alike; copy link of the selection.

const view = FeatureViews.find((v) => v.id === 'addressing')!;
const links = Passages.flat().filter((r): r is { text: string; href: string } => typeof r !== 'string');
const href = (text: string) => links.find((l) => l.text === text)!.href;

let viewer: SeqViewer | undefined;

// A frame tick asked for by an earlier file on timers since dropped would never run, and the frame lane would wait for
// it: a host function of this file's own asks for a fresh one (16 ms on the faked clock).
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
  const v = composeViewer(view);
  viewer = v;
  v.resize(720);
  await v.load(view.build());
  await vi.advanceTimersByTimeAsync(200);
  return v;
}

const protein = '/views/main/protein';

describe('Links, headless', () => {
  it('runs every link of the text with no issue', async () => {
    const v = await make();
    expect(links.length).toBe(6);
    for (const l of links) {
      const run = followLink(v, l.href);
      expect(run?.issues, l.text).toEqual([]);
      expect(
        run?.outcomes.every((o) => o.done),
        l.text,
      ).toBe(true);
    }
  });

  it('focuses, highlights on a channel and selects', async () => {
    const v = await make();
    expect(describeRun(followLink(v, href('DNA-binding domain')))).toBe('Focused');
    expect([v.view.get(protein)!.start, v.view.get(protein)!.end]).toEqual([93, 292]);
    followLink(v, href('transactivation domain'));
    expect(v.highlight.get('link')?.address).toBe('/refs/tad');
    followLink(v, href('hotspot residues'));
    expect(v.addresses.selectionAddress()).toBe('/refs/hotspots');
    expect(v.selection.value.value.shown.parts[0].ranges).toEqual(new Int32Array([174, 175, 247, 248, 272, 273]));
    // the selection's link, as the story's Copy link button gives it
    expect(v.addresses.copyLink()).toBe('!seq-focus=%2Frefs%2Fhotspots');
  });

  it('focuses one residue named in two numberings alike', async () => {
    const v = await make();
    followLink(v, href('R175 (1-based)'));
    const a = [v.view.get(protein)!.start, v.view.get(protein)!.end];
    v.view.fit(protein);
    followLink(v, href('174 (0-based)'));
    expect([v.view.get(protein)!.start, v.view.get(protein)!.end]).toEqual(a);
  });

  it('hovers a link on its own channel, leaving the pointer’s hover alone', async () => {
    const v = await make();
    hoverLink(v, href('residues 245 to 249'));
    expect(v.highlight.get('link-hover')?.address).toBe('/data/p53/(range :from 245 :to 249)');
    expect(v.hover.current.owner).toBeUndefined();
    leaveLink(v);
    expect(v.highlight.get('link-hover')).toBeUndefined();
  });
});
