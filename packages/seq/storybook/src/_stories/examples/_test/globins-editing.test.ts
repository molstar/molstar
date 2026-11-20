/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { hashNode, resolveAddress } from '@molstar/core/repr/tree';
import { theme } from '@molstar/seqx-core/lang';
import { type BlockPlacement, type PairsLay, PlacePoint, type RowsLay } from '@molstar/seqx-core/render';
import { GlobinsViews } from '@molstar/seqx-examples';
import { HemeContacts } from '@molstar/seqx-examples/globins/hemes';
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

const Pt = new PlacePoint();

/**
 * Every heme contact's arc at its data: its feet on the first block showing the heme and its chain (split, the data
 * names no copy; merged, the ×4 and ×2 blocks), the chain's at the contact's residue.
 */
function contactsAtData(sec: Awaited<ReturnType<typeof drawn>>) {
  const blocks = sec.section.blocks!;
  const x = (system: string, p: number) => {
    const b = blocks.find((q) => q.axis.system === system)!;
    expect(sec.placement.place(b.offset! + (b.form === 'whole' ? 0.5 : p + 0.5), Pt)).toBe(true);
    return Pt.x;
  };
  let n = 0;
  for (const l of sec.sectionLays) {
    if (l.type.name !== 'pairs') continue;
    const lay = l.lay as PairsLay;
    const d = lay.prep!.data!;
    for (let i = 0; i < lay.n; i++, n++) {
      const k = lay.pair[i];
      const want = [x(lay.prep!.aSystem, d.aStart[k]), x(lay.prep!.bSystem, d.bStart[k])].sort((a, b) => a - b);
      const got = [lay.curve[8 * i], lay.curve[8 * i + 6]].sort((a, b) => a - b);
      expect(got[0]).toBeCloseTo(want[0], 6);
      expect(got[1]).toBeCloseTo(want[1], 6);
    }
  }
  expect(n).toBe(HemeContacts.hba.length + HemeContacts.hbb.length);
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

  it('hemoglobin: the copies split (one patch) and merged again by undo, every heme contact an arc at its data', async () => {
    const s = await open('hemoglobin');
    const start = hash(s);
    const labels = async () => ((await drawn(s)).placement as BlockPlacement).labels;
    expect(await labels()).toEqual(['α (A) ×2', 'β (B) ×2', 'heme ×4']);
    contactsAtData(await drawn(s));
    expect(await toggleCopies(s)).toBe('Split');
    expect(copiesOf(s)).toBe('each');
    expect(await labels()).toEqual(['α (A)', 'β (B)', 'α (C)', 'β (D)', 'heme', 'heme', 'heme', 'heme']);
    // the split hemoglobin keeps every contact (on the first copy: its data names none)
    contactsAtData(await drawn(s));
    expect(s.viewer.layoutIssues.value).toEqual([]);
    await s.stack.undo();
    expect(hash(s)).toBe(start);
    expect(await labels()).toEqual(['α (A) ×2', 'β (B) ×2', 'heme ×4']);
    contactsAtData(await drawn(s));
  });
});
