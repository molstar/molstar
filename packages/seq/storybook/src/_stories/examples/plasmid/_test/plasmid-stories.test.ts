/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolveAddress } from '@molstar/core/repr/tree';
import { PlasmidViews } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { PlasmidLength } from '@molstar/seqx-examples/plasmid/derive';
import type { RangesLay } from '@molstar/seqx-core/render';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession, historyRows } from '../../../features/editing/compose.js';
import { openExample } from '../../_test/draw-checks.js';
import { bookmarkRange, rangeOf, setStart, showEnzymes, sitesDataOf, strandHidden, toggleStrand } from '../editing.js';

// The Plasmid's stories, headless on a fake clock, composed as the stories compose them: the crossing view's bookmarks
// focus the feature that crosses the origin by its shortest arc (EC-12; a Shift+drag across the origin is pending);
// the editing views (a strand's track hidden, the
// enzymes shown, bla bookmarked from its start to its stop, the linear view's start set) each one patch and one step of
// the host's stack, undone and redone.

const view = (id: string) => PlasmidViews.find((v) => v.id === id)!;
let session: EditingSession | undefined;

const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  session?.dispose();
  session = undefined;
  vi.useRealTimers();
});

async function make(id: string) {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const v = view(id);
  // the record's `seq` files read as committed, at once (a read of the disk is not on the fake clock)
  const read = (uri: string) => readFileSync(`${PackageFolder}${uri}`);
  const s = composeEditing(v, {
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => read(uri).toString('utf8'),
          readBytes: async (uri) => new Uint8Array(read(uri)),
        },
      }),
    },
  });
  session = s;
  s.viewer.resize(720);
  await s.viewer.load(v.build());
  await vi.advanceTimersByTimeAsync(64);
  return s;
}

describe('Plasmid, features across the origin, headless', () => {
  it('a bookmark on bla or the primer site focuses it: the address is the focus, the view takes a window over it', async () => {
    const s = await make('crossing');
    const v = s.viewer;
    for (const id of ['bla', 'primer']) {
      expect(v.runCommand('address.focus', { address: `/refs/${id}` })).toBe('Focused');
      await vi.advanceTimersByTimeAsync(64);
      v.draw();
      expect(v.addresses.focusAddress()).toBe(`/refs/${id}`);
    }
  });

  // the focus takes the arc the item covers, never the hull of its runs (0 to the length across the origin)
  it('a bookmark on bla focuses it by its shortest arc: its two runs, not the whole ring (EC-12)', async () => {
    const s = await make('crossing');
    const v = s.viewer;
    expect(v.runCommand('address.focus', { address: '/refs/bla' })).toBe('Focused');
    await vi.advanceTimersByTimeAsync(64);
    v.draw();
    const w = v.viewport!;
    expect(w.end - w.start).toBeGreaterThanOrEqual(861);
    expect(w.end - w.start).toBeLessThan(PlasmidLength / 2);
  });
});

describe('Plasmid, editing, headless', () => {
  it('a strand’s track hidden, the enzymes shown, bla bookmarked: each one step, undone and redone in order', async () => {
    const s = await make('ring');
    const r0 = s.viewer.spec.root.value;
    expect(await toggleStrand(s, '-')).toBe('Hidden');
    expect(strandHidden(s, '-')).toBe(true);
    expect(strandHidden(s, '+')).toBe(false);
    expect(await showEnzymes(s, ['EcoRI', 'HindIII'])).toBe(true);
    expect(sitesDataOf(s)).toBe('/data/sites/(item :id [EcoRI HindIII])');
    expect(await bookmarkRange(s, 2311, 486, 'bla-range', 'bla, start to stop')).toBe('Bookmarked');
    expect(historyRows(s).done.length).toBe(3);
    // the bookmark is a range over the origin: it resolves and focuses by its shortest arc
    const ref = s.viewer.addresses.refs.value.find((r) => r.id === 'bla-range');
    expect(ref?.to).toEqual(['/data/plasmid/(range :from 2312 :to 486)']);
    expect(s.viewer.runCommand('address.focus', { address: '/refs/bla-range' })).toBe('Focused');
    // Undo walks them back, Redo forward
    await s.stack.undo();
    expect(s.viewer.addresses.refs.value.find((r) => r.id === 'bla-range')).toBeUndefined();
    await s.stack.undo();
    expect(sitesDataOf(s)).toBe('/data/sites');
    await s.stack.undo();
    expect(strandHidden(s, '-')).toBe(false);
    expect(s.viewer.spec.root.value).toEqual(r0);
    await s.stack.redo();
    expect(strandHidden(s, '-')).toBe(true);
  });

  it('the rotated linear view’s start is one patch: set, undone, redone', async () => {
    const s = await make('rotated-linear');
    expect(rangeOf(s)).toEqual([2000, 700]);
    expect(await setStart(s, 2400, 1386)).toBe(true);
    expect(rangeOf(s)).toEqual([2400, 1100]);
    await vi.advanceTimersByTimeAsync(64);
    expect(resolveAddress(s.viewer.spec.root.value, ['views', 'main', 'linear'])?.params?.range).toEqual([2400, 1100]);
    expect(s.viewer.frame!.layout.sections[0].axis!.pieces).toEqual([2400, 2686, 0, 1100]);
    await s.stack.undo();
    expect(rangeOf(s)).toEqual([2000, 700]);
    await s.stack.redo();
    expect(rangeOf(s)).toEqual([2400, 1100]);
  });
});

describe('Plasmid ring, drawn on the canvas', () => {
  it('bla’s arrow on the ring at fit, and on the detail at its start, on the piece holding it', async () => {
    const v = await openExample(view('ring'), {});
    try {
      const ring = v.frame!.layout.sections[0];
      expect(ring.tracks.length).toBe(4);
      // the ring at 720 px draws bars, and each feature its arrowhead by its own width: bla's on the minus strand
      const ringLay = ring.tracks[1].layers[0].lay as RangesLay;
      expect(ringLay.mode).toBe('bars');
      const ringTips = Array.from({ length: ringLay.m }, (_, j) => j).filter(
        (j) => ringLay.mitem[j] >= 0 && ringLay.keys[ringLay.mitem[j]] === 'bla' && ringLay.mtip[j] > 0,
      );
      expect(ringTips).toHaveLength(1);
      // the detail zoomed to bla's start draws the minus strand's arrow on the piece holding its start
      v.view.setRange('/views/main/ring', 2250, 2450);
      await vi.advanceTimersByTimeAsync(50);
      v.draw();
      const lay = v.frame!.layout.sections[1].tracks[1].layers[0].lay as RangesLay;
      const tips: number[] = [];
      for (let j = 0; j < lay.m; j++) if (lay.keys[lay.mitem[j]] === 'bla') tips.push(lay.mtip[j]);
      expect(tips).toEqual([1]);
      expect(v.frame!.layout.issues).toEqual([]);
    } finally {
      v.retire();
      vi.useRealTimers();
    }
  });
});
