/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { location, locationsEqual } from '@molstar/seqx-core/query';
import { LightTokens } from '@molstar/seqx-core/render';
import { colorRuns } from '@molstar/seqx-core/structure';
import { P53 } from '@molstar/seqx-examples/p53/data';
import { P53Views } from '@molstar/seqx-examples';
import { ChainTableFile, EntityStart, EntityLength, StructureAddresses } from '@molstar/seqx-examples/p53/structure';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterAll, afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import type { Cell, ChainTable, Span } from '../../../support/fake-structure/chain-table.js';
import { type HeadersBehavior, HeadersProvider } from '@molstar/seqx-viewer/behaviors/headers';
import { BaseColor } from '../../../support/fake-structure/mvs-state.js';
import { ColorByEntry, colorByBehavior, viewerModel } from '../../../support/host-adapter.js';
import { colorRunner, type FakeLink, linkFake, litLine } from '../compose.js';
import { fakeClock, openView, tableOf } from './open.js';

// The structure link's stories, headless on a fake clock, composed as the stories compose them (compose.ts): p53 and
// PDB 2OCJ in Seq*, the fake structure viewer on its own recorded chain table, `linkStructure` through the host's
// adapter. Residue numbers are the independent witness: 2OCJ's author numbers are UniProt's, so a P04637 position `p`
// is the author number `p + 1` of the table, a fact of the table and of nothing in Seq*. The look of the fake and the
// stories' hints are not tested; what crosses is.

const Protein = StructureAddresses.protein;
const Section = '/views/main/protein';
const Complex = StructureAddresses.complex;
const pointer = { pointer: true };

const table: ChainTable = tableOf(ChainTableFile);
const view = (id: string) => P53Views.find((v) => v.id === id)!;

let viewer: SeqViewer | undefined;
let link: FakeLink | undefined;
let restore: (() => void) | undefined;

beforeAll(async () => {
  viewer = await openView(view('structure'));
});
afterAll(() => viewer?.retire());
afterEach(() => {
  link?.dispose();
  link = undefined;
  viewer?.hover.clear(pointer);
  restore?.();
  restore = undefined;
});

/** The link of a test: on the fake clock, the fake's latency as given. */
function open(latency = 0, options: Parameters<typeof linkFake>[3] = {}) {
  restore = fakeClock();
  link = linkFake(viewer!, table, Complex, { latency, ...options });
  return { v: viewer!, fake: link.fake };
}

const cellsOf = (fake: FakeLink['fake'], spans: readonly Span[]) =>
  spans.flatMap((s) =>
    Array.from({ length: s.end - s.start }, (_, k) => ({
      row: s.row,
      index: s.start + k,
      asym: fake.rows[s.row].asym,
      res: fake.rows[s.row].residues[s.start + k],
    })),
  );

/** The cells of the table's polymer rows that are residue `n` (author number) and modeled. */
const modeled = (fake: FakeLink['fake'], n: number): Cell[] =>
  fake.rows.flatMap((r) =>
    r.polymer
      ? r.residues.flatMap((x, index) => (x.observed && x.auth === String(n) ? [{ row: r.index, index }] : []))
      : [],
  );

const hoverUp = (v: SeqViewer, position: number) =>
  v.hover.set({ section: Section, system: Protein, position }, pointer);

describe('Hover both ways', () => {
  it('every lit cell of the fake is in the image of Seq*’s hover, on every chain, the gaps lit never (EC-48, EC-118)', () => {
    const { v, fake } = open();
    for (let p = EntityStart; p < EntityStart + EntityLength; p++) {
      hoverUp(v, p);
      const lit = cellsOf(fake, fake.highlight.value).map((c) => `${c.row}:${c.index}`);
      const expected = modeled(fake, p + 1).map((c) => `${c.row}:${c.index}`);
      expect(lit, `P04637 ${p + 1}`).toEqual(expected);
      for (const c of cellsOf(fake, fake.highlight.value)) expect(c.res.observed).toBe(true);
    }
    hoverUp(v, 174);
    expect(fake.highlight.value).toHaveLength(4);
    // an unobserved end lights nothing, a zinc's chip is not a P04637 residue
    hoverUp(v, EntityStart);
    expect(fake.highlight.value).toEqual([]);
    expect(fake.rows.filter((r) => !r.polymer).every((r) => !fake.highlight.value.some((s) => s.row === r.index))).toBe(
      true,
    );
  });

  it('and back: a cell of the fake is Seq*’s hover on its own chain and, projected, on P04637 at its author number, never a sibling copy (EC-65)', async () => {
    const { v, fake } = open();
    for (const row of fake.rows.filter((r) => r.polymer))
      for (let index = 0; index < row.residues.length; index += 9) {
        if (!row.residues[index].observed) continue;
        fake.pointerOver({ row: row.index, index });
        const loc = v.hover.location();
        const unit = loc!.parts.filter((p) => p.system.startsWith(`${Complex}/`));
        // the hovered residue: author number `auth`, position `auth - 1` on P04637, the chain's unit position `index`
        expect(unit.map((p) => [p.system, ...p.ranges])).toEqual([[`${Complex}/${row.asym}`, index, index + 1]]);
        expect(v.hover.active).toBe(true);
        if (index % 27 === 0) {
          v.draw();
          const on = v.emphasis.projected?.hover?.parts.find((p) => p.system === Protein);
          expect(Array.from(on?.ranges ?? []), `${row.asym} ${row.residues[index].auth}`).toEqual([
            Number(row.residues[index].auth) - 1,
            Number(row.residues[index].auth),
          ]);
        }
      }
    fake.pointerOver(null);
  });

  it('a gap of the fake is no hover and sends nothing: it has no atoms (EC-118)', () => {
    const { v, fake } = open();
    const gap = fake.rows[0].residues.findIndex((r) => !r.observed);
    fake.pointerOver({ row: 0, index: gap });
    expect(fake.sent).toEqual([]);
    expect(v.hover.active).toBe(false);
  });

  it('latency: the fake’s late messages never take the hover from the pointer, and its late leave clears nothing (EC-141)', async () => {
    const { v, fake } = open(300);
    const cell = { row: 0, index: 80 };
    fake.pointerOver(cell);
    expect(v.hover.active).toBe(false);
    await vi.advanceTimersByTimeAsync(300);
    // the message arrived: the link’s hover is drawn on P04637 (no card)
    expect(v.hover.active).toBe(true);
    // the pointer moves into Seq*; the fake’s pointer left a moment ago: its null is on the way
    fake.pointerOver(null);
    hoverUp(v, EntityStart + 120);
    expect(v.hover.owner).toBe(pointer);
    const pointerAt = v.hover.location();
    await vi.advanceTimersByTimeAsync(400);
    // the stale null arrived: the pointer’s hover is untouched
    expect(v.hover.active).toBe(true);
    expect(locationsEqual(v.hover.location()!, pointerAt!)).toBe(true);
    // and the fake lit what the pointer hovers, apart from its own hover
    expect(cellsOf(fake, fake.highlight.value).map((c) => c.res.auth)).toEqual(
      Array(4).fill(String(EntityStart + 121)),
    );
  });
});

describe('Hover with latency, the other way', () => {
  it('a late hover of the fake, older than the pointer’s own, is dropped; a newer one is taken (EC-141)', async () => {
    const { v, fake } = open(600);
    // the pointer sweeps across the fake (its hover is on the way) and back onto P04637 within the latency
    fake.pointerOver({ row: 0, index: 80 });
    await vi.advanceTimersByTimeAsync(200);
    hoverUp(v, EntityStart + 120);
    const pointerAt = v.hover.location()!;
    await vi.advanceTimersByTimeAsync(500);
    // the stale message arrived: it took nothing
    expect(fake.sent.map((m) => m.kind)).toEqual(['hover']);
    expect(v.hover.owner).toBe(pointer);
    expect(locationsEqual(v.hover.location()!, pointerAt)).toBe(true);
    // a hover of the fake after the pointer's last is taken, when it arrives
    fake.pointerOver({ row: 1, index: 90 });
    await vi.advanceTimersByTimeAsync(600);
    expect(v.hover.owner).not.toBe(pointer);
    expect(v.hover.location()!.parts.map((p) => [p.system, ...p.ranges])).toEqual([[`${Complex}/B`, 90, 91]]);
  });
});

describe('The adapter’s words', () => {
  it('a selector it places is a target; one with a field it does not place is none, never widened to every unit', () => {
    const { fake } = open();
    const { adapter } = link!;
    const t = adapter.targetOfSelectors([{ label_asym_id: 'A', beg_label_seq_id: 3, end_label_seq_id: 5 }]);
    expect(t!.parts.map((p) => [p.system, ...p.ranges])).toEqual([[`${Complex}/A`, 2, 5]]);
    for (const sel of [{ label_entity_id: '1' }, { auth_asym_id: 'A' }, { label_asym_id: 'A', auth_seq_id: 100 }])
      expect(adapter.targetOfSelectors([sel]), JSON.stringify(sel)).toBeNull();
    expect(fake.rows.length).toBeGreaterThan(0);
  });
});

describe('Selection both ways', () => {
  it('Seq*’s selection is the fake’s on every chain; a click, a toggle, Esc and a structure-only selection come back whole (EC-48, EC-118)', async () => {
    const { v, fake } = open();
    const at = (a: number, b: number) => location({ parts: [{ system: Protein, ranges: [a, b] }] });
    // Seq* → fake: one message a frame; the unobserved end is cut
    v.selection.set({
      elements: v.selection.value.value.elements,
      location: at(EntityStart, EntityStart + 6),
      shown: at(EntityStart, EntityStart + 6),
    });
    await vi.advanceTimersByTimeAsync(20);
    const sel = cellsOf(fake, fake.selection.value);
    expect(sel.map((c) => c.res.auth)).toEqual(
      fake.rows.filter((r) => r.polymer).flatMap(() => ['96', '97', '98', '99']),
    );
    expect(new Set(sel.map((c) => c.asym))).toEqual(new Set(['A', 'B', 'C', 'D']));

    // fake → Seq*: a click selects its cell alone, replacing; Seq* reveals it and does not send it back
    const told = fake.sent.length;
    const unitB = (...r: number[]) => [[`${Complex}/B`, ...r]];
    const shown = () => v.selection.value.value.shown.parts.map((p) => [p.system, ...p.ranges]);
    fake.click({ row: 1, index: 100 });
    expect(shown()).toEqual(unitB(100, 101));
    await vi.advanceTimersByTimeAsync(20);
    expect(link!.log.value.filter((l) => l.startsWith('→ select'))).toHaveLength(1);
    expect(fake.sent.length).toBe(told + 1);
    // Cmd/Ctrl+click toggles a second cell in, again out
    fake.click({ row: 1, index: 110 }, true);
    expect(shown()).toEqual(unitB(100, 101, 110, 111));
    fake.click({ row: 1, index: 110 }, true);
    expect(shown()).toEqual(unitB(100, 101));
    // "Select chain": a selection only the structure had, whole on the unit (its observed residues)
    fake.selectChain();
    expect(shown()).toEqual(unitB(2, 196));
    // Esc clears both
    fake.clear();
    expect(v.selection.value.value.shown.parts).toEqual([]);
  });
});

describe('Focus and back', () => {
  it('a double-click focuses Seq* and a focus in Seq* moves the fake’s camera over the jump’s duration; back brings it', async () => {
    const { v, fake } = open();
    fake.doubleClick({ row: 2, index: 60 });
    expect(v.view.lastJump).toMatchObject({ kind: 'focus' });
    const first = v.view.lastJump!;
    // Seq* focuses a residue: the fake’s frame is on its cells on every chain, in the jump’s time
    v.addresses.focusLocation(
      location({ parts: [{ system: Protein, ranges: [EntityStart + 150, EntityStart + 151] }] }),
    );
    await vi.advanceTimersByTimeAsync(20);
    const cam = fake.camera.value!;
    expect(cam.durationMs).toBe(v.view.lastJump!.durationMs);
    expect(cellsOf(fake, cam.spans).map((c) => c.res.auth)).toEqual(Array(4).fill(String(EntityStart + 151)));
    expect(v.view.lastJump!.n).toBeGreaterThan(first.n);
    // back: the step before named a location (the double-click’s): the camera goes there
    v.view.back();
    await vi.advanceTimersByTimeAsync(20);
    expect(v.view.lastJump!.kind).toBe('back');
    // the step before named the double-clicked cell (a unit's residue): the camera is on it, row 2, cell 60, again
    expect(fake.camera.value!.spans).toEqual([{ row: 2, start: 60, end: 61 }]);
    v.view.forward();
    await vi.advanceTimersByTimeAsync(20);
    expect(v.view.lastJump!.kind).toBe('forward');
    expect(cellsOf(fake, fake.camera.value!.spans).map((c) => c.res.auth)).toEqual(
      Array(4).fill(String(EntityStart + 151)),
    );
  });
});

describe('Chain maps', () => {
  let compare: SeqViewer;
  beforeAll(async () => {
    compare = await openView(view('structure-compare'));
  });
  afterAll(() => compare.retire());

  it('the ranked map (SIFTS’s) lights the residue by number; the link via the shifted alignment lights the one before it (EC-24)', () => {
    restore = fakeClock();
    const p = EntityStart + 107;
    const lit = (via?: string[]) => {
      link = linkFake(compare, table, Complex, via ? { link: { complexes: [{ complex: Complex, via }] } } : {});
      compare.hover.set({ section: Section, system: Protein, position: p }, pointer);
      const out = cellsOf(link.fake, link.fake.highlight.value).map((c) => `${c.asym}${c.res.auth}`);
      // the story's readout names the same residues (F8: a one-cell shift read as numbers)
      expect(litLine(link.fake)).toBe(`Lit in the fake: ${out.map((o) => `${o[0]} ${o.slice(1)}`).join(', ')}`);
      link.dispose();
      link = undefined;
      compare.hover.clear(pointer);
      return out;
    };
    // SIFTS’s `structure` map is ranked: no `via` needed, nothing ambiguous
    expect(lit()).toEqual(['A', 'B', 'C', 'D'].map((c) => `${c}${p + 1}`));
    // the alignment, shifted by one residue, puts P04637 `p` on the entity's residue before
    expect(lit([StructureAddresses.aln])).toEqual(['A', 'B', 'C', 'D'].map((c) => `${c}${p}`));
  });
});

describe('Color by this track', () => {
  let colored: SeqViewer;
  let colorLink: FakeLink | undefined;
  beforeAll(async () => {
    colored = await openView(view('structure'), {
      registry: [ColorByEntry],
      behaviors: [
        colorByBehavior((layer) =>
          colorRunner(
            colored,
            () => colorLink,
            () => 'light',
          )(layer),
        ),
      ],
    });
  });
  afterAll(() => colored.retire());

  it('the header entry’s command colors the fake: its colors equal colorRuns at each observed residue’s preimage (EC-65, EC-118)', () => {
    restore = fakeClock();
    colored.draw();
    colorLink = link = linkFake(colored, table, Complex);
    const { fake } = link;
    const track = colored.view.resolved.value!.sections[0].tracks.find((t) => t.label === 'P04637')!;
    const layer = track.layers[0];
    // through the header entry: its `when` and `args` make the command's context, as a click on its slot does
    const headers = colored.behavior<HeadersBehavior>(HeadersProvider)!;
    const target = { section: Section, track: track.address };
    expect(
      headers
        .menuOf(target)!
        .groups.flat()
        .map((i) => i.name),
    ).toContain('fake-host:color-by');
    expect(headers.run('fake-host:color-by', target)).toMatch(/^Colored \d+ runs$/);
    expect(link.log.value.at(-1)).toMatch(/^→ color .*: Colored \d+ runs$/);
    const colors = colored.renderer!.layerColors();
    const runs = colorRuns(viewerModel(colored), colored.store.store, layer, undefined, { colors }).runs;
    const channel = colors.channel(layer, LightTokens)!;
    const hex = (c: number) => `#${(c >>> 8).toString(16).padStart(6, '0')}`;
    let checked = 0;
    for (const row of fake.rows) {
      if (!row.polymer) continue;
      row.residues.forEach((r, index) => {
        const unit = `${Complex}/${row.asym}`;
        const run = runs.find((x) =>
          x.target.parts.some(
            (q) =>
              q.system === unit &&
              Array.from(q.ranges).some((e, i) => i % 2 === 0 && e <= index && index < q.ranges[i + 1]),
          ),
        );
        if (!r.observed) {
          // a gap: no atoms, in no run, the base color
          expect(run).toBeUndefined();
          return;
        }
        // the run is the color of the letter at the residue's preimage on P04637 (its author number minus one)
        const p = Number(r.auth) - 1;
        const datum = { value: Number.NaN, letter: P53[p], start: p, end: p + 1, field: () => undefined };
        expect(hex(run!.color), `${row.asym} ${r.auth}`).toBe(fake.colorOf(row.index, index));
        expect(hex(channel.at(datum) as number)).toBe(fake.colorOf(row.index, index));
        checked++;
      });
    }
    expect(checked).toBe(4 * 194);
    // the same track again: back to the base state
    expect(headers.run('fake-host:color-by', target)).toBe('Restored the base colors');
    expect(fake.colorOf(0, 20)).toBe(BaseColor);
    colorLink = undefined;
  });

  it('the track the structure is colored by has its icon pressed, one at a time; its tooltip names it (F9)', () => {
    restore = fakeClock();
    colored.draw();
    colorLink = link = linkFake(colored, table, Complex);
    const headers = colored.behavior<HeadersBehavior>(HeadersProvider)!;
    const target = (id: string) => ({ section: Section, track: `${Section}/${id}` });
    const item = (id: string) =>
      headers
        .menuOf(target(id))!
        .groups.flat()
        .find((i) => i.name === 'fake-host:color-by')!;
    const pressed = () => ['domains', 'letters'].map((id) => item(id).pressed);
    expect(pressed()).toEqual([false, false]);
    expect(item('letters').label).toBe('Color the structure by this track');
    headers.run('fake-host:color-by', target('letters'));
    expect(pressed()).toEqual([false, true]);
    // the canvas header's slot is pressed too (drawn at rest by the chrome), the other track's not
    colored.draw();
    const canvas = (id: string) =>
      colored
        .frame!.layout.sections[0].tracks.find((t) => t.track.address === `${Section}/${id}`)!
        .header!.slots.find((x) => x.part.name === 'fake-host:color-by')?.part.pressed;
    expect([canvas('domains'), canvas('letters')]).toEqual([undefined, true]);
    expect(item('letters').label).toBe('The structure is colored by this track (click to restore)');
    // another track takes the color: the press moves, never two at once
    headers.run('fake-host:color-by', target('domains'));
    expect(pressed()).toEqual([true, false]);
    // the same again restores the base colors: none pressed
    headers.run('fake-host:color-by', target('domains'));
    expect(pressed()).toEqual([false, false]);
    colorLink = undefined;
  });
});
