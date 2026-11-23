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
import { HemeUnitContacts } from '@molstar/seqx-examples/globins/hemes';
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

// The Globins' editing views, headless on a fake clock, composed as the story composes
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
 * Each heme block's arcs reach exactly one chain block, the heme's own chain's (merged, the heme ×4's reach α ×2 and
 * β ×2); every heme contact's arc at its data (the complex's pairs, derived per pair of units): its feet on the blocks
 * showing its heme's unit and its chain's, the chain's at its residue.
 */
function contactsAtData(sec: Awaited<ReturnType<typeof drawn>>) {
  const blocks = sec.section.blocks!;
  const x = (b: number, p: number) => {
    const q = blocks[b];
    expect(sec.placement.place(q.offset! + (q.form === 'whole' ? 0.5 : p + 0.5), Pt)).toBe(true);
    return Pt.x;
  };
  const last = (q: (typeof blocks)[number]) => (q.form === 'whole' ? 0 : q.axis.ranges![1] - q.axis.ranges![0] - 1);
  const blockAt = (px: number) => blocks.findIndex((q, b) => x(b, 0) - 1 <= px && px <= x(b, last(q)) + 1);
  const lays = sec.sectionLays.filter((l) => l.type.name === 'pairs').map((l) => l.lay as PairsLay);
  // the drawing: each heme block's arcs reach exactly one chain block (split: the heme's own chain's, in order)
  const reached = new Map<number, Set<number>>();
  for (const lay of lays)
    for (let i = 0; i < lay.n; i++) {
      const feet = [blockAt(lay.curve[8 * i]), blockAt(lay.curve[8 * i + 6])];
      const heme = feet.find((b) => blocks[b].form === 'whole' || blocks[b].label.startsWith('heme'))!;
      let r = reached.get(heme);
      if (!r) reached.set(heme, (r = new Set()));
      r.add(feet[0] === heme ? feet[1] : feet[0]);
    }
  const hemes = blocks.map((q, b) => (q.label.startsWith('heme') ? b : -1)).filter((b) => b >= 0);
  expect([...reached.keys()].sort()).toEqual(hemes);
  for (const [heme, chains] of reached) {
    const own = blocks.length === 8 ? [hemes.indexOf(heme)] : [0, 1];
    expect([...chains].sort()).toEqual(own);
  }
  // the data: each arc's feet at its heme's unit and its chain's residue
  let n = 0;
  const blockOf = (unit: string) => blocks.findIndex((q) => q.systems.includes(unit));
  for (const lay of lays) {
    const d = lay.prep!.data!;
    expect(lay.prep!.endUnits, 'the contacts name their units').toBeDefined();
    const units = lay.prep!.endUnits!;
    for (let i = 0; i < lay.n; i++, n++) {
      const k = lay.pair[i];
      const want = [x(blockOf(units[d.aUnit![k]]), d.aStart[k]), x(blockOf(units[d.bUnit![k]]), d.bStart[k])];
      const got = [lay.curve[8 * i], lay.curve[8 * i + 6]].sort((a, b) => a - b);
      want.sort((a, b) => a - b);
      expect(got[0]).toBeCloseTo(want[0], 6);
      expect(got[1]).toBeCloseTo(want[1], 6);
    }
  }
  expect(n).toBe(HemeUnitContacts.reduce((m, c) => m + c.residues.length, 0));
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
    // split: each heme's arcs reach the chain it sits in, and only it
    contactsAtData(await drawn(s));
    expect(s.viewer.layoutIssues.value).toEqual([]);
    await s.stack.undo();
    expect(hash(s)).toBe(start);
    expect(await labels()).toEqual(['α (A) ×2', 'β (B) ×2', 'heme ×4']);
    contactsAtData(await drawn(s));
  });
});
