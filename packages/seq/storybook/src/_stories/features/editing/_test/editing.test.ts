/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashNode, resolveAddress } from '@molstar/core/repr/tree';
import { P53Views } from '@molstar/seqx-examples';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { DomHost } from '@molstar/ui/dom';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  bookmarkView,
  cancelDrag,
  composeEditing,
  dragHeight,
  type EditingSession,
  endDrag,
  historyRows,
  moveFirst,
  printedSpec,
  recolorStructure,
  setDomainsColor,
  toggleTrack,
  trackHidden,
} from '../compose.js';

// The Editing story, headless, composed as the story composes it: a stand-in host's one stack with the viewer's spec
// and a stand-in structure as participants; edits as commands, a drag as a preview ending in one patch or Esc, undone
// and redone in order; the printed spec following. The performance representative: an edit (patch, resolve, redraw)
// on the p53 spec and its undo.

const view = P53Views.find((v) => v.id === 'tracks')!;
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

async function make(fake = true) {
  if (fake) vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const s = composeEditing(view);
  session = s;
  s.viewer.resize(720);
  await s.viewer.load(view.build());
  if (fake) await vi.advanceTimersByTimeAsync(32);
  else await s.viewer.settle();
  return s;
}

const order = (s: EditingSession) =>
  resolveAddress(s.viewer.spec.root.value, ['views', 'main', 'protein'])!.children!.map((c) => c.id);
const height = (root: Parameters<typeof resolveAddress>[0], id: string) =>
  resolveAddress(root, ['views', 'main', 'protein', id])?.params?.height;

describe('Editing, headless', () => {
  it('edits as commands, each one step between the host’s own, undone and redone in order', async () => {
    const s = await make();
    const r0 = s.viewer.spec.root.value;
    expect(await toggleTrack(s, 'conservation')).toBe('Hidden');
    expect(trackHidden(s, 'conservation')).toBe(true);
    await recolorStructure(s, 'blue');
    expect(await moveFirst(s, 'domains')).toBe('Moved');
    expect(await setDomainsColor(s, '#d62728')).toBe('Colored');
    expect(await bookmarkView(s)).toBe('Bookmarked');
    expect(order(s)).toEqual(['domains', 'letters', 'conservation']);
    expect(historyRows(s).done).toEqual([
      'Add a bookmark',
      'Color a layer',
      'Move a track',
      'Recolor the structure',
      'Hide a track',
    ]);
    expect(printedSpec(s.viewer.spec.shown.value)).toContain('#d62728');

    for (let i = 0; i < 3; i++) await s.stack.undo();
    expect(order(s)).toEqual(['letters', 'domains', 'conservation']);
    expect(s.structure.value).toBe('blue');
    await s.stack.undo();
    expect(s.structure.value).toBe('grey');
    expect(trackHidden(s, 'conservation')).toBe(true);
    await s.stack.undo();
    expect(s.viewer.spec.root.value).toBe(r0);
    expect(historyRows(s).redo).toHaveLength(5);
    await s.stack.redo();
    await s.stack.redo();
    expect([trackHidden(s, 'conservation'), s.structure.value]).toEqual([true, 'blue']);
  });

  it('a drag previews each move and ends in one step; Esc ends it in nothing', async () => {
    const s = await make();
    const root = s.viewer.spec.root.value;
    for (const h of [30, 40, 52]) expect(dragHeight(s, 'domains', h)).toBe(true);
    expect(height(s.viewer.spec.shown.value, 'domains')).toBe(52);
    expect(printedSpec(s.viewer.spec.shown.value)).toContain(':height 52');
    expect(s.viewer.spec.root.value).toBe(root);
    await endDrag(s);
    expect(height(s.viewer.spec.root.value, 'domains')).toBe(52);
    expect(historyRows(s).done).toEqual(['Resize a track']);

    dragHeight(s, 'letters', 70);
    expect(cancelDrag(s)).toBe('Cancelled');
    expect(height(s.viewer.spec.shown.value, 'letters')).toBeUndefined();
    expect(historyRows(s).done).toEqual(['Resize a track']);
  });

  it('a drag survives an edit landing mid-gesture, and is dropped when the edit removes what it resizes', async () => {
    const s = await make();
    dragHeight(s, 'domains', 44);
    await moveFirst(s, 'conservation');
    expect(height(s.viewer.spec.shown.value, 'domains')).toBe(44);
    await endDrag(s);
    expect(historyRows(s).done).toEqual(['Resize a track', 'Move a track']);

    dragHeight(s, 'letters', 20);
    await s.viewer.spec.patch('Remove the letters', {
      ops: [{ op: 'remove', at: ['views', 'main', 'protein', 'letters'] }],
    });
    expect(s.viewer.spec.preview.value).toBeUndefined();
    expect(s.viewer.spec.issues.value.map((i) => i.code)).toEqual(['seq.viewer.preview.dropped']);
  });

  it('an arrival never changes the spec’s hash, nor enters the history', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
    const region = P53Views.find((v) => v.id === 'gene-by-region')!;
    const text = regionSeqText();
    const s = composeEditing(region, {
      delayMs: 300,
      ports: { assets: () => ({ reader: { readText: async () => text } }) },
    });
    session = s;
    s.viewer.resize(720);
    await s.viewer.load(region.build());
    await s.viewer.spec.patch('Retitle', { ops: [{ op: 'set', at: [], key: 'title', value: 'edited' }] });
    const revision = s.viewer.spec.revision.value;
    const nodes = s.viewer.store.store.nodeVersion.value;
    await vi.advanceTimersByTimeAsync(4 * 300);
    // the region's transcripts arrived: the store moved, the spec did not
    expect(s.viewer.store.store.nodeVersion.value).not.toBe(nodes);
    expect(s.viewer.spec.revision.value).toBe(revision);
    expect(hashNode(s.viewer.spec.root.value)).toBe(revision);
    expect(historyRows(s).done).toEqual(['Retitle']);
  });
});

describe('Editing: the performance representative', () => {
  itWithin('an edit on the p53 spec (patch, resolve, redraw) and its undo, under 16 ms each', async (within) => {
    const s = await make(false);
    const v = s.viewer;
    const edit = async (h: number) => {
      await v.spec.patch('Resize', {
        ops: [{ op: 'set', at: ['views', 'main', 'protein', 'domains'], key: 'height', value: h }],
      });
      v.view.resolved.value;
      v.draw();
    };
    const undo = async () => {
      await s.stack.undo();
      v.view.resolved.value;
      v.draw();
    };
    // warm-ups
    for (let i = 0; i < 3; i++) await edit(20 + i);
    for (let i = 0; i < 3; i++) await undo();
    const n = 5;
    await within(16 * n, async () => {
      for (let i = 0; i < n; i++) await edit(30 + i);
    });
    await within(16 * n, async () => {
      for (let i = 0; i < n; i++) await undo();
    });
    expect(height(v.spec.root.value, 'domains')).toBeUndefined();
  });
});
