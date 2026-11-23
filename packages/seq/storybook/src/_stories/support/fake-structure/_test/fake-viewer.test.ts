/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  type ChainTable,
  normalizeSpans,
  observedSpans,
  rowsOf,
  selectorsOfSpans,
  spansOfSelector,
} from '../chain-table.js';
import { frameBox } from '../frame-box.js';
import { FakeStructureViewer } from '../fake-viewer.js';
import { colorNodesOf, stateIssues, structureState } from '../mvs-state.js';

// The fake structure viewer alone (it imports nothing of Seq*): its selectors round-trip over its own table, its
// MolViewSpec states are validated before it takes them, its messages are delayed by `latency`, and the grammar's
// selections are what its words say.

const res = (label: number | null, observed = true) => ({
  label,
  auth: String(label ?? 1),
  ins: '',
  comp: 'ALA',
  observed,
});
const table: ChainTable = {
  entry: 'TEST',
  chains: [
    {
      label_asym_id: 'A',
      auth_asym_id: 'A',
      entity: '1',
      operators: ['1_555', '2_555'],
      residues: [1, 2, 3, 4, 5, 6].map((n) => res(n, n !== 1 && n !== 6)),
    },
    { label_asym_id: 'B', auth_asym_id: 'A', entity: '2', operators: ['1_555'], residues: [res(null)] },
  ],
};

afterEach(() => vi.useRealTimers());

describe('the fake structure viewer’s hover name', () => {
  it('names the hovered residue, a hovered range by its ends, and nothing without a hover', () => {
    const fake = new FakeStructureViewer(table);
    expect(fake.hoverName()).toBeUndefined();
    fake.pointerOver({ row: 0, index: 2 });
    expect(fake.hoverName()).toBe('A 3 ALA');
    fake.hover.set([{ row: 0, start: 1, end: 4 }]);
    expect(fake.hoverName()).toBe('A 2 ALA to 4 ALA');
    fake.pointerOver(null);
    expect(fake.hoverName()).toBeUndefined();
  });
});

describe('the fake structure viewer’s selectors', () => {
  it('spans in selectors and back are the same observed cells; a copy is a row of its own', () => {
    const rows = rowsOf(table);
    expect(rows.map((r) => `${r.asym} ${r.operator}`)).toEqual(['A 1_555', 'A 2_555', 'B 1_555']);
    const spans = observedSpans(rows, [
      { row: 0, start: 0, end: 6 },
      { row: 1, start: 2, end: 4 },
      { row: 2, start: 0, end: 1 },
    ]);
    // the unobserved ends are cut
    expect(spans).toEqual([
      { row: 0, start: 1, end: 5 },
      { row: 1, start: 2, end: 4 },
      { row: 2, start: 0, end: 1 },
    ]);
    const selectors = selectorsOfSpans(rows, spans);
    expect(selectors[0]).toEqual({
      label_asym_id: 'A',
      instance_id: '1_555',
      beg_label_seq_id: 2,
      end_label_seq_id: 5,
    });
    expect(selectors[2]).toEqual({ label_asym_id: 'B', instance_id: '1_555' });
    expect(normalizeSpans(selectors.flatMap((s) => spansOfSelector(rows, s)))).toEqual(spans);
    // no instance: every copy
    expect(spansOfSelector(rows, { label_asym_id: 'A', label_seq_id: 3 }).map((s) => s.row)).toEqual([0, 1]);
  });
});

describe('the fake structure viewer’s state', () => {
  it('takes a valid MolViewSpec state and reads its colors; refuses an invalid one and keeps what it had', () => {
    const fake = new FakeStructureViewer(table);
    const state = structureState('TEST', [
      {
        selector: [{ label_asym_id: 'A', instance_id: '2_555', beg_label_seq_id: 2, end_label_seq_id: 3 }],
        color: '#ff0000',
      },
    ]);
    expect(stateIssues(state)).toBeUndefined();
    expect(colorNodesOf(state)).toHaveLength(2);
    expect(fake.setState(state)).toBeUndefined();
    expect(fake.colorOf(1, 1)).toBe('#ff0000');
    expect(fake.colorOf(0, 1)).not.toBe('#ff0000');
    const bad = structureState('TEST', [{ selector: [{ label_asym_id: 'A' }], color: 'notacolor' as never }]);
    expect(fake.setState(bad)?.length).toBeGreaterThan(0);
    expect(fake.colorOf(1, 1)).toBe('#ff0000');
  });
});

describe('the fake structure viewer’s grammar and latency', () => {
  it('click selects, Cmd/Ctrl+click toggles, Esc clears; a gap is nothing; its messages come `latency` late and not after it is gone', async () => {
    vi.useFakeTimers();
    const fake = new FakeStructureViewer(table, { latency: 50 });
    const told: string[] = [];
    fake.onSelect((s) =>
      told.push(JSON.stringify(s.map((x) => x.beg_label_seq_id ?? x.label_seq_id ?? x.label_asym_id))),
    );
    fake.click({ row: 0, index: 0 });
    expect(fake.selection.value).toEqual([]);
    fake.click({ row: 0, index: 2 });
    fake.click({ row: 0, index: 3 }, true);
    expect(fake.selection.value).toEqual([{ row: 0, start: 2, end: 4 }]);
    expect(told).toEqual([]);
    await vi.advanceTimersByTimeAsync(50);
    expect(told).toEqual(['[3]', '[3]']);
    fake.click({ row: 0, index: 2 }, true);
    fake.clear();
    expect(fake.selection.value).toEqual([]);
    fake.dispose();
    await vi.advanceTimersByTimeAsync(100);
    expect(told).toHaveLength(2);
  });

  it('"select chain" is the selection’s chains’ observed cells and "within 5 Å" grows by its neighbors and the chain’s ligand', () => {
    const fake = new FakeStructureViewer(table);
    fake.click({ row: 0, index: 2 });
    fake.within();
    expect(fake.selection.value).toEqual([
      { row: 0, start: 1, end: 5 },
      { row: 2, start: 0, end: 1 },
    ]);
    // the selection's chains: A, and the ligand chain B it reached
    fake.selectChain();
    expect(fake.selection.value).toEqual([
      { row: 0, start: 1, end: 5 },
      { row: 2, start: 0, end: 1 },
    ]);
  });
});

describe('the fake structure viewer’s "within 5 Å" on a whole chain', () => {
  it('selects the residues of the other chains within 5 Å, never the selection again, and its note says how many', () => {
    const fake = new FakeStructureViewer(table);
    fake.click({ row: 0, index: 2 });
    fake.selectChain();
    const chain = fake.selection.value;
    expect(chain).toEqual([{ row: 0, start: 1, end: 5 }]);
    fake.within();
    // the copy 2_555 of the chain: its residues in the interface, the unobserved cell never taken
    expect(fake.selection.value).toEqual([{ row: 1, start: 1, end: 3 }]);
    expect(fake.selection.value).not.toEqual(chain);
    expect(fake.note.value).toBe('2 residues within 5 Å of chain A');
    // the note goes with the next selection
    fake.clear();
    expect(fake.note.value).toBeUndefined();
  });
});

describe('the fake structure viewer’s camera frame', () => {
  it('is placed against its container, cell rectangles scrolled, so a cell of the third row is on the third row', () => {
    const cells = [
      { left: 130, top: 90, right: 134, bottom: 104 },
      { left: 150, top: 90, right: 154, bottom: 104 },
    ];
    const container = { left: 20, top: 40, right: 800, bottom: 200 };
    expect(frameBox(cells, container, { left: 0, top: 0 }, 6)).toEqual({ left: 104, top: 44, width: 36, height: 26 });
    expect(frameBox(cells, container, { left: 10, top: 5 }, 6).left).toBe(114);
  });
});

describe('the fake structure viewer’s other selections', () => {
  it('"select chain" takes every chain of the selection, "within 5 Å" the ligand of the selection’s copy alone', () => {
    const fake = new FakeStructureViewer(table);
    fake.click({ row: 0, index: 2 });
    fake.click({ row: 1, index: 2 }, true);
    fake.selectChain();
    expect(fake.selection.value).toEqual([
      { row: 0, start: 1, end: 5 },
      { row: 1, start: 1, end: 5 },
    ]);
    // the ligand is of copy 1_555: a selection on copy 2_555 alone does not take it
    fake.click({ row: 1, index: 2 });
    fake.within();
    expect(fake.selection.value).toEqual([{ row: 1, start: 1, end: 5 }]);
  });
});
