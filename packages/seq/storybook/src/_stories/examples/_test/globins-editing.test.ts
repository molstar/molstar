/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { hashNode, resolveAddress } from '@molstar/core/repr/tree';
import { theme } from '@molstar/seqx-core/lang';
import type { BlockPlacement, RowsLay } from '@molstar/seqx-core/render';
import { GlobinsViews } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession } from '../../features/editing/compose.js';
import {
  copiesOf,
  keepFolds,
  keptFolds,
  logoHidden,
  pinnedRows,
  pinRows,
  RowsLayer,
  setScheme,
  toggleClade,
  toggleCopies,
  toggleLogo,
} from '../globins/editing.js';

// The Globins' editing views (examples.md, Globins: Views), headless on a fake clock, composed as the story composes
// them: each edit one patch, one step of a host's stack, undone to the spec it started from. The family: the residue
// scheme changed (a theme patch) and the logo hidden; the folded clades: `clades.keep` writes the folds shown, a row
// pinned above the tree's order; hemoglobin: the copies split and merged again.

let session: EditingSession | undefined;
const raf = DomHost.requestAnimationFrame;

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  session?.dispose();
  session = undefined;
  vi.useRealTimers();
});

/** The fixtures read at once (the file system would not settle on the fake clock). */
const fixtures = () => ({
  reader: {
    readText: async (uri: string) => readFileSync(`${PackageFolder}${uri}`, 'utf8'),
    readBytes: async (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`)),
  },
});

async function open(id: string) {
  const view = GlobinsViews.find((v) => v.id === id)!;
  const s = composeEditing(view, { ports: { assets: fixtures }, delayMs: 0 });
  session = s;
  s.viewer.resize(900);
  await s.viewer.load(view.build());
  await vi.advanceTimersByTimeAsync(200);
  s.viewer.draw();
  return s;
}

async function drawn(s: EditingSession) {
  await vi.advanceTimersByTimeAsync(50);
  s.viewer.draw();
  return s.viewer.frame!.layout.sections[0];
}

const hash = (s: EditingSession) => hashNode(s.viewer.spec.root.value);

describe('Globins editing views, headless', () => {
  it('family: the residue scheme changed and the logo hidden, a step each, undone to the start', async () => {
    const s = await open('family');
    const start = hash(s);
    expect(await setScheme(s, 'percent-identity')).toBe('Colored');
    expect(resolveAddress(s.viewer.spec.root.value, RowsLayer.split('/').filter(Boolean))?.params?.color).toEqual(
      theme.scheme('percent-identity'),
    );
    expect(await toggleLogo(s)).toBe('Hidden');
    expect(logoHidden(s)).toBe(true);
    expect((await drawn(s)).tracks.some((t) => t.track.address.endsWith('/logo'))).toBe(false);
    await s.stack.undo();
    await s.stack.undo();
    expect(hash(s)).toBe(start);
    expect((await drawn(s)).tracks.some((t) => t.track.address.endsWith('/logo'))).toBe(true);
  });

  it('folded: a clade unfolded in the session, the folds kept (one patch), a row pinned; undone', async () => {
    const s = await open('folded');
    const start = hash(s);
    const opened = keptFolds(s);
    expect(opened.length).toBe(2);
    await toggleClade(s, opened[0]);
    // session state: no step, the spec unchanged
    expect(hash(s)).toBe(start);
    expect(await keepFolds(s)).toBe('Kept');
    expect(keptFolds(s)).toEqual([opened[1]]);
    const sec = await drawn(s);
    let lay: RowsLay | undefined;
    for (const b of sec.lays) for (const l of b ?? []) if (l.type.name === 'rows') lay = l.lay as RowsLay;
    const last = lay!.display.rowAt(lay!.display.count - 1);
    await pinRows(s, [last]);
    expect(pinnedRows(s)).toEqual([lay!.label(last)]);
    await s.stack.undo();
    await s.stack.undo();
    expect(hash(s)).toBe(start);
  });

  it('hemoglobin: the copies split (one patch) and merged again by undo', async () => {
    const s = await open('hemoglobin');
    const start = hash(s);
    const labels = async () => ((await drawn(s)).placement as BlockPlacement).labels;
    expect(await labels()).toEqual(['α (A) ×2', 'β (B) ×2', 'heme ×4']);
    expect(await toggleCopies(s)).toBe('Split');
    expect(copiesOf(s)).toBe('each');
    expect(await labels()).toEqual(['α (A)', 'β (B)', 'α (C)', 'β (D)', 'heme', 'heme', 'heme', 'heme']);
    await s.stack.undo();
    expect(hash(s)).toBe(start);
    expect(await labels()).toEqual(['α (A) ×2', 'β (B) ×2', 'heme ×4']);
  });
});
