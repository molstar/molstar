/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure link (_spec/README.md "Structure link") on a recording side, headless on a fake clock: hover,
// selection and focus each way, echoes never sent back, a stale null, one selection message a frame, focus with the
// jump's duration, the clearing at each of the three ends, a binding's `via`, the gate, the directions; Seq*'s frame
// and the side as two consumers of one link.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { EmptyLocation, type Location, location, locationsEqual } from '@molstar/seqx-core/query';
import { seqToValue } from '@molstar/seqx-core/data';
import { atomsSystemOf, nodeAddress } from '@molstar/seqx-core/resolve';
import type { RangesLay } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { type LinkStructureOptions, linkStructure } from '../link-structure.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';
import {
  Aln,
  Complex,
  ComplexSection,
  EntStart,
  onChains,
  onUp,
  p53Structure,
  type RecordingSide,
  recordingSide,
  Up,
  UpSection,
} from './structure-fixture.js';

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(options: LinkStructureOptions = {}, root: Node = p53Structure()) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewer = v;
  await v.load(root);
  v.resize(800);
  await v.settle();
  frames.flush();
  const side = recordingSide();
  const dispose = linkStructure(v, side, options);
  return { v, side, dispose, frames };
}

const pointer = { pointer: true };
const hoverUp = (v: SeqViewer, position: number) => v.hover.set({ section: UpSection, system: Up, position }, pointer);
const same = (a: Location | null | undefined, b: Location) => !!a && locationsEqual(a, b);
const kinds = (side: RecordingSide) => side.told.map((t) => t.kind);

/**
 * A sequence `x` and a spec map onto unit A of a complex a source describes (SIFTS per chain): the description
 * arrives with the source's region (`describe`), not before.
 */
function described(): Node {
  const b = createSeqSpec({ title: 'a map onto a described unit' });
  const ent = { kind: 'sequence', params: { alphabet: 'protein', letters: 'CPECGKSF' } } as Node;
  const cx = {
    kind: 'complex',
    id: 'cx',
    children: [{ kind: 'unit', id: 'A', params: { of: nodeAddress(Zf, 'ent'), chain: 'A' } }],
  } as Node;
  const data = {
    coverage: { systems: [{ name: 'zf', length: 8 }], shapes: [{ system: 'zf', shape: 'letters' }], levels: [1] },
    pairs: [
      {
        request: { system: 'zf', start: 0, end: 8, level: 1, shapes: ['letters'] },
        answer: {
          ok: true,
          tables: { letters: { form: 'letters', start: 0, letters: 'CPECGKSF' } },
          nodes: [
            { key: 'ent', node: ent },
            { key: 'cx', node: cx },
          ],
        },
      },
    ],
  };
  b.data.source('zf', { format: 'seq', value: seqToValue(data as never) });
  const xs = b.data.sequence('x', { alphabet: 'protein', letters: 'CPECGKS' });
  b.data.map('xa', { from: xs, to: UnitA as never, method: 'offset', by: 1 });
  b.views.view('main', (v) => v.section('x', { of: xs }).track('letters').layer('letters', { type: 'letters' }));
  return b.build();
}
const Zf = '/data/zf';
const UnitA = `${nodeAddress(Zf, 'cx')}/A`;

describe('linkStructure over a map that waits for its description (EC-169)', () => {
  it('holds a hover and a selection until the description arrives, never sends them as nothing; then sends them', async () => {
    const { v, side, frames } = await make({}, described());
    const loc = location({ parts: [{ system: '/data/x', ranges: [0, 2] }] });
    v.hover.setLocation(loc, pointer);
    v.selection.set({ elements: v.selection.value.value.elements, location: loc, shown: loc });
    frames.flush();
    // pending: nothing told (not a null, not an empty target)
    expect(side.told).toEqual([]);
    // the description arrives with the source's region
    v.store.store.want([{ source: Zf, system: 'zf', start: 0, end: 8, level: 1 }]);
    await v.store.store.whenSettled();
    frames.flush();
    expect(kinds(side).sort()).toEqual(['highlight', 'select']);
    for (const t of side.told)
      expect(same(t.target, location({ parts: [{ system: UnitA, ranges: [1, 3] }] }))).toBe(true);
  });
});

describe('linkStructure', () => {
  it('hover out: a residue of P04637 lit on every chain of the entity (EC-48), unobserved lights nothing (EC-118)', async () => {
    const { v, side } = await make();
    hoverUp(v, EntStart + 10);
    expect(kinds(side)).toEqual(['highlight']);
    expect(same(side.told[0].target, onChains(10, 11))).toBe(true);
    // the same image sends nothing (the zinc and the other residues unchanged)
    hoverUp(v, EntStart + 10);
    expect(side.told).toHaveLength(1);
    // an unobserved end: nothing to light, so the hover is cleared
    hoverUp(v, EntStart + 0);
    expect(side.told[1]).toEqual({ kind: 'highlight', target: null });
    // the pointer leaving after that sends nothing more
    v.hover.clear(pointer);
    expect(side.told).toHaveLength(2);
  });

  it('hover in: drawn as the hover on P04637 and the complex, never sent back; the pointer takes it (EC-141)', async () => {
    const { v, side, frames } = await make();
    side.hover(onChains(10, 11, ['A']));
    expect(v.hover.active).toBe(true);
    expect(v.hover.current.layer).toBeUndefined();
    frames.flush();
    // Seq*'s frame, the second consumer: the hover projected onto P04637
    const drawn = v.emphasis.projected?.hover;
    expect(drawn?.parts.find((p) => p.system === Up)?.ranges).toEqual(Int32Array.of(EntStart + 10, EntStart + 11));
    expect(side.told).toEqual([]);
    // the pointer takes the hover; the side's late leave is stale and clears nothing
    hoverUp(v, EntStart + 20);
    expect(kinds(side)).toEqual(['highlight']);
    side.hover(null);
    expect(v.hover.owner).toBe(pointer);
    // the side's hover again: it ends the pointer's sent before it (cleared there), and is itself never sent; its
    // leave clears only the link's own and sends nothing
    side.hover(onChains(30, 31, ['B']));
    expect(side.told[1]).toEqual({ kind: 'highlight', target: null });
    side.hover(null);
    expect(v.hover.active).toBe(false);
    expect(side.told).toHaveLength(2);
  });

  it('hover in: a hover that happened before the pointer’s last own is dropped, a later one taken (EC-141)', async () => {
    const { v, side } = await make();
    hoverUp(v, EntStart + 20);
    const own = performance.now();
    side.hover(onChains(30, 31, ['B']), own - 1);
    expect(v.hover.owner).toBe(pointer);
    side.hover(onChains(30, 31, ['B']), own + 1);
    expect(v.hover.owner).not.toBe(pointer);
    expect(v.hover.active).toBe(true);
  });

  it('selection out: one message a frame, the last value; an unchanged image sends nothing', async () => {
    const { v, side, frames } = await make();
    for (const end of [20, 30, 40]) {
      const loc = onUp(EntStart + 10, EntStart + end);
      v.selection.set({ elements: v.selection.value.value.elements, location: loc, shown: loc });
    }
    expect(side.told).toEqual([]);
    frames.flushOnce();
    expect(kinds(side)).toEqual(['select']);
    expect(same(side.told[0].target, onChains(10, 40))).toBe(true);
    v.selection.clear();
    frames.flush();
    expect(side.told[1]).toEqual({ kind: 'select', target: null });
    v.selection.clear();
    frames.flush();
    expect(side.told).toHaveLength(2);
  });

  it('selection in: set and revealed, never sent back; a null clears it', async () => {
    const { v, side, frames } = await make();
    v.addresses.focusLocation(onUp(0, 20));
    const jumps = v.view.jumps.value;
    side.click(onChains(150, 160, ['A']));
    frames.flush();
    expect(same(v.selection.value.value.shown, onChains(150, 160, ['A']))).toBe(true);
    // revealed: a jump, the zoom kept, and no focus sent for it
    expect(v.view.jumps.value).toBe(jumps + 1);
    expect(v.view.lastJump?.kind).toBe('reveal');
    expect(side.told).toEqual([]);
    side.click(null);
    frames.flush();
    expect(v.selection.value.value.shown).toBe(EmptyLocation);
    expect(side.told).toEqual([]);
  });

  it('focus out with the jump’s duration on the fake clock; back brings the camera; a fit names nothing', async () => {
    const { v, side } = await make();
    const jumpMs = v.config.view.jumpMs;
    v.addresses.focusLocation(onUp(EntStart + 50, EntStart + 60));
    expect(side.told).toEqual([{ kind: 'focus', target: expect.anything(), durationMs: jumpMs }]);
    expect(same(side.told[0].target, onChains(50, 60))).toBe(true);
    v.addresses.focusLocation(onUp(EntStart + 100, EntStart + 120));
    expect(v.view.back()).toBe(true);
    expect(kinds(side)).toEqual(['focus', 'focus', 'focus']);
    expect(same(side.told[2].target, onChains(50, 60))).toBe(true);
    // a move by hand since a jump: the step it pushes names no Location, so its back sends nothing
    v.view.pan(UpSection, 5);
    v.addresses.focusLocation(onUp(EntStart + 150, EntStart + 160));
    const n = side.told.length;
    expect(v.view.back()).toBe(true);
    expect(side.told).toHaveLength(n);
    // the view's fit: nothing sent
    v.runCommand('view.fit', { section: UpSection });
    expect(v.view.lastJump?.kind).toBe('fit');
    expect(side.told).toHaveLength(n);
  });

  it('a double-click focuses the domain: told to the side with the jump’s duration, the motion run on a fake clock', async () => {
    // every timer faked (the frames' fallback a 16 ms timer): the clock drives the jump's motion
    vi.useFakeTimers();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewer = v;
    v.resize(800);
    void v.load(p53Structure());
    await vi.advanceTimersByTimeAsync(2000);
    const side = recordingSide();
    linkStructure(v, side);
    v.draw();
    const lay = v.frame!.layout.sections[0].tracks[1].layers[0].lay as RangesLay;
    const [x, y] = [(lay.mx0[0] + lay.mx1[0]) / 2, lay.laneY + (lay.mlane[0] + 0.5) * lay.laneHeight];
    expect(v.input.doubleClick(x, y)).toBe(true);
    expect(v.view.lastJump).toMatchObject({ kind: 'focus' });
    const jumpMs = v.config.view.jumpMs;
    expect(side.told).toEqual([{ kind: 'focus', target: expect.anything(), durationMs: jumpMs }]);
    expect(same(side.told[0].target, onChains(50, 60))).toBe(true);
    // the sections move over the duration told: under way at its middle, at rest past it
    expect(v.moving).toBe(true);
    vi.advanceTimersByTime(jumpMs / 2);
    expect(v.moving).toBe(true);
    vi.advanceTimersByTime(jumpMs / 2 + 40);
    expect(v.moving).toBe(false);
    // a step zoom names no Location (`zoom`): nothing sent
    v.runCommand('view.zoom-out', { section: UpSection });
    expect(v.view.lastJump?.kind).toBe('zoom');
    expect(side.told).toHaveLength(1);
  });

  it('focus in: one jump, not sent back; followJumps false keeps back and forward here', async () => {
    const { v, side } = await make({ followJumps: false });
    const jumps = v.view.jumps.value;
    side.doubleClick(onChains(50, 60, ['B']));
    expect(v.view.jumps.value).toBe(jumps + 1);
    expect(v.view.lastJump?.kind).toBe('focus');
    expect(side.told).toEqual([]);
    v.addresses.focusLocation(onUp(EntStart + 100, EntStart + 120));
    expect(kinds(side)).toEqual(['focus']);
    v.view.back();
    expect(kinds(side)).toEqual(['focus']);
  });

  describe('clearing at each end: highlight and select null once, the link’s hover dropped', () => {
    const ends: [string, (x: { v: SeqViewer; side: RecordingSide; dispose: () => void; frames: Frames }) => void][] = [
      ['the disposer', (x) => x.dispose()],
      ['the side’s signal', (x) => x.side.abort.abort()],
      ['the viewer retiring', (x) => x.v.retire()],
    ];
    for (const [name, end] of ends)
      it(name, async () => {
        const x = await make();
        x.side.hover(onChains(10, 11, ['A']));
        // a selection message waiting for its frame is dropped with the link
        const loc = onUp(EntStart + 10, EntStart + 20);
        x.v.selection.set({ elements: x.v.selection.value.value.elements, location: loc, shown: loc });
        end(x);
        x.frames.flush();
        expect(x.side.told).toEqual([
          { kind: 'highlight', target: null },
          { kind: 'select', target: null },
        ]);
        expect(x.side.listeners).toBe(0);
        expect(x.v.hover.active).toBe(false);
        x.dispose();
        x.side.abort.abort();
        expect(x.side.told).toHaveLength(2);
      });
  });

  it('a binding’s via: the alignment’s path one residue off, the ranked SIFTS without it (EC-23)', async () => {
    const ranked = await make({}, p53Structure({ aln: true }));
    hoverUp(ranked.v, EntStart + 10);
    expect(same(ranked.side.told[0].target, onChains(10, 11))).toBe(true);
    ranked.v.retire();
    frames?.restore();
    const via = await make({ complexes: [{ complex: Complex, via: [Aln] }] }, p53Structure({ aln: true }));
    hoverUp(via.v, EntStart + 10);
    expect(same(via.side.told[0].target, onChains(9, 10))).toBe(true);
  });

  it('the gate (EC-167): past the hover’s bounds whole units, or nothing', async () => {
    const units = await make({ limits: { hover: { maxRuns: 2, maxLength: 100, over: 'units' } } });
    hoverUp(units.v, EntStart + 10);
    expect(same(units.side.told[0].target, onChains(2, 196))).toBe(true);
    units.v.retire();
    frames?.restore();
    const skip = await make({ limits: { hover: { maxRuns: 2, maxLength: 100, over: 'skip' } } });
    hoverUp(skip.v, EntStart + 10);
    expect(skip.side.told).toEqual([]);
  });

  it('directions: out only ignores the side; in only tells it nothing', async () => {
    const o = await make({ hover: 'out', select: 'out', focus: 'out' });
    o.side.hover(onChains(10, 11));
    o.side.click(onChains(10, 11));
    o.side.doubleClick(onChains(10, 11));
    expect(o.v.hover.active).toBe(false);
    expect(o.v.selection.value.value.shown).toBe(EmptyLocation);
    expect(o.side.listeners).toBe(0);
    o.v.retire();
    frames?.restore();
    const i = await make({ hover: 'in', select: 'in', focus: 'in' });
    hoverUp(i.v, EntStart + 10);
    i.v.addresses.focusLocation(onUp(EntStart + 50, EntStart + 60));
    i.frames.flush();
    i.v.retire();
    expect(i.side.told).toEqual([]);
  });

  it('a section of the complex: a residue hovered on chain C is sent on chain C alone, the copy rule (EC-25)', async () => {
    const { v, side } = await make();
    v.hover.set({ section: ComplexSection, system: `${Complex}/C`, position: 40 }, pointer);
    expect(side.told).toHaveLength(1);
    expect(same(side.told[0].target, onChains(40, 41, ['C']))).toBe(true);
  });

  it('atoms from the side: kept on a unit whose atoms form is shown, collapsed onto its residue on another', async () => {
    const b = createSeqSpec({ title: 'atoms' });
    const atp = b.data.molecule('atp', {
      residues: ['ATP'],
      atoms: { residue: [0, 0, 0], name: ['PG', 'O1G', 'N1'], element: ['P', 'O', 'N'] },
    });
    b.data.complex('c', {}, (c) => {
      c.unit('l', { of: atp });
      c.unit('l2', { of: atp });
    });
    b.views.view('v', (v) =>
      v.section('complex', { of: '/data/c', copies: 'each' }, (s) => {
        s.block('l', { form: 'atoms' });
        s.track('letters').layer('letters', { type: 'letters' });
      }),
    );
    const { v, side } = await make({ revealOnSelect: 'none' }, b.build());
    const [l, l2] = [atomsSystemOf('/data/c/l'), atomsSystemOf('/data/c/l2')];
    expect(v.view.resolved.value?.sections[0].blocks?.map((x) => x.axis.system)).toContain(l);
    side.click(
      location({
        parts: [
          { system: l, ranges: [1, 2] },
          { system: l2, ranges: [1, 2] },
        ],
      }),
    );
    const shown = v.selection.value.value.shown;
    expect(shown.parts.map((p) => [p.system, Array.from(p.ranges)])).toEqual([
      [l, [1, 2]],
      ['/data/c/l2', [0, 1]],
    ]);
  });
});
