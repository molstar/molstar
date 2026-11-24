/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's track groups in the viewer (_spec/README.md "Behaviors", Groups), headless on a fake clock: folds are
// session state in the section's viewport, the groups' `collapsed` the opening set; a toggle is no step of the history,
// collapse all and expand all are jumps that back restores; `groups.keep` writes the folds into `collapsed` and an
// edit of `collapsed` resets the session's set; a fold keeps the selection.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { HiddenTrackMessage } from '../models/address.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

export const GroupsSection = '/views/v/p53';
const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWF';

/** Letters, then Features (domains, motifs), Variants (folded at opening), Signals (two values tracks). */
export function groupsViewerSpec(variants: Record<string, unknown> = { collapsed: true }): Node {
  const b = createSeqSpec({ title: 'groups' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: Letters });
  const ranges = (id: string, start: number[], end: number[]) =>
    b.data.ranges(id, {
      target: p53,
      rows: { start: new Int32Array(start), end: new Int32Array(end), label: start.map((s) => `${id}${s}`) },
    });
  const domains = ranges('domains', [0, 10, 30], [8, 28, 50]);
  const motifs = ranges('motifs', [4, 20], [6, 22]);
  const clinvar = ranges('clinvar', [3, 17, 40], [4, 18, 41]);
  const gnomad = ranges('gnomad', [7, 44], [8, 45]);
  const wave = (k: number) => Float32Array.from(Letters, (_, i) => 0.5 + 0.5 * Math.sin((i + k) / 4));
  const a = b.data.values('a', { target: p53, values: wave(0) });
  const c = b.data.values('c', { target: p53, values: wave(5) });
  b.views.view('v', (v) => {
    const s = v.section('p53', { of: p53 });
    s.track('letters', {}, (t) => t.layer('l', { type: 'letters' }));
    const f = s.group('features', { label: 'Features' });
    f.track('domains', { label: 'Domains' }, (t) => t.layer('l', { type: 'ranges', data: domains }));
    f.track('motifs', { label: 'Motifs' }, (t) => t.layer('l', { type: 'ranges', data: motifs }));
    const va = s.group('variants', { label: 'Variants', ...variants });
    va.track('clinvar', { label: 'ClinVar' }, (t) => t.layer('l', { type: 'variants', data: clinvar }));
    va.track('gnomad', { label: 'gnomAD' }, (t) => t.layer('l', { type: 'variants', data: gnomad }));
    const sg = s.group('signals', { label: 'Signals' });
    sg.track('a', { label: 'Signal a' }, (t) => t.layer('l', { type: 'values', data: a }));
    sg.track('c', { label: 'Signal c' }, (t) => t.layer('l', { type: 'values', data: c }));
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(spec = groupsViewerSpec()) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    history: uiHistory(),
    config: [[viewerItem('view.jumpMs'), 0]],
  });
  viewer = v;
  await v.load(spec);
  v.resize(480);
  await v.settle();
  frames.flush();
  return v;
}

const stackOf = (v: SeqViewer) => (v.spec.history as UiHistory).stack;
const folded = (v: SeqViewer) => v.view.foldedGroups(GroupsSection);
const backSteps = (v: SeqViewer) => v.view.history.value.back;

describe('track groups in the viewer', () => {
  it('opens with each group’s `collapsed`; a toggle is session state in the viewport and no step', async () => {
    const v = await make();
    expect(folded(v)).toEqual(['variants']);
    expect(v.view.get(GroupsSection)?.groups).toBeUndefined();
    const steps = backSteps(v);
    await v.execute('groups.toggle', { address: `${GroupsSection}/features` });
    expect(folded(v)).toEqual(['variants', 'features']);
    expect(v.view.get(GroupsSection)?.groups).toEqual(['variants', 'features']);
    expect(v.view.isFolded(`${GroupsSection}/features`)).toBe(true);
    await v.execute('groups.toggle', { address: `${GroupsSection}/variants` });
    expect(folded(v)).toEqual(['features']);
    // a hand's move: no step of the view's history, nothing on the spec's undo
    expect(backSteps(v)).toBe(steps);
    expect(stackOf(v).steps.value.length).toBe(0);
    expect(await v.execute('groups.toggle', { address: `${GroupsSection}/none` })).toBe('No group there');
  });

  it('collapse all and expand all are jumps: back restores the folds; a fold keeps the selection', async () => {
    const v = await make();
    v.selection.add(rangesLocation('/data/p53', [3, 4]));
    const selected = v.selection.value.value;
    const steps = backSteps(v);
    await v.execute('groups.collapseAll', { section: GroupsSection });
    expect(folded(v)).toEqual(['features', 'variants', 'signals']);
    expect(backSteps(v)).toBe(steps + 1);
    expect(await v.execute('groups.collapseAll', { section: GroupsSection })).toBe('Already so');
    await v.execute('groups.expandAll', {});
    expect(folded(v)).toEqual([]);
    expect(backSteps(v)).toBe(steps + 2);
    expect(v.selection.value.value).toBe(selected);
    expect(v.view.back()).toBe(true);
    expect(folded(v)).toEqual(['features', 'variants', 'signals']);
    expect(v.view.back()).toBe(true);
    // the opening state again: no session set
    expect(v.view.get(GroupsSection)?.groups).toBeUndefined();
    expect(folded(v)).toEqual(['variants']);
    expect(v.view.forward()).toBe(true);
    expect(folded(v)).toEqual(['features', 'variants', 'signals']);
  });

  it('`track.hide` and `track.show` on a group’s address: one entry in the host’s list; `track.move` moves the group', async () => {
    const v = await make();
    const signals = `${GroupsSection}/signals`;
    await v.execute('track.hide', { address: signals });
    await v.settle();
    frames!.flush();
    expect(v.hiddenTracks.value).toEqual([
      { address: signals, label: 'Signals (2 tracks)', section: GroupsSection, kind: 'group' },
    ]);
    expect(v.view.resolved.value!.sections[0].trackGroups!.map((g) => g.address)).not.toContain(signals);
    await v.execute('track.show', { address: signals });
    await v.settle();
    expect(v.hiddenTracks.value).toEqual([]);
    // a hidden track in a group names its group
    await v.execute('track.hide', { address: `${signals}/a` });
    await v.settle();
    expect(v.hiddenTracks.value[0]).toMatchObject({
      address: `${signals}/a`,
      group: { address: signals, label: 'Signals' },
    });
    // a group moved among the section's children: after Signals
    await v.execute('track.move', { address: `${GroupsSection}/features`, after: 'signals' });
    await v.settle();
    const children = v.spec.root.value.children![1].children![0].children![0].children!.map((c) => c.id);
    expect(children).toEqual(['letters', 'variants', 'signals', 'features']);
  });

  it('the header cursor stops on a group’s header (a folded group one stop); Space folds and opens it there', async () => {
    const v = await make();
    frames!.flush();
    v.selection.setCursor({
      section: GroupsSection,
      system: '/data/p53',
      position: 0,
      part: 'header',
      track: `${GroupsSection}/letters`,
      keyboard: true,
    });
    const key = (k: string) => v.input.key({ key: k });
    const at = () => v.selection.cursor.value?.track;
    key('ArrowDown');
    expect(at()).toBe(`${GroupsSection}/features`);
    key('ArrowDown');
    expect(at()).toBe(`${GroupsSection}/features/domains`);
    key('ArrowDown');
    key('ArrowDown');
    // Variants is folded: one stop, its tracks skipped
    expect(at()).toBe(`${GroupsSection}/variants`);
    key('ArrowDown');
    expect(at()).toBe(`${GroupsSection}/signals`);
    key(' ');
    expect(folded(v)).toEqual(['variants', 'signals']);
    frames!.flush();
    key(' ');
    expect(folded(v)).toEqual(['variants']);
    frames!.flush();
    // a track's header whose group folds: the cursor moves to the group's header
    key('ArrowDown');
    expect(at()).toBe(`${GroupsSection}/signals/a`);
    await v.execute('groups.toggle', { address: `${GroupsSection}/signals` });
    frames!.flush();
    expect(at()).toBe(`${GroupsSection}/signals`);
  });

  it('a reveal of an item only a folded group draws opens it in the same jump; back folds it again', async () => {
    const v = await make();
    const steps = backSteps(v);
    // ClinVar's item: drawn only by the folded Variants group
    expect(v.addresses.reveal('/data/clinvar/(item :label clinvar17)').done).toBe(true);
    expect(folded(v)).toEqual([]);
    expect(backSteps(v)).toBe(steps + 1);
    expect(v.view.back()).toBe(true);
    expect(folded(v)).toEqual(['variants']);
    // an item an open track draws: nothing opened
    await v.execute('groups.collapseAll', { section: GroupsSection });
    await v.execute('groups.toggle', { address: `${GroupsSection}/features` });
    v.addresses.reveal('/data/domains/(item :label domains10)');
    expect(folded(v)).toEqual(['variants', 'signals']);
  });

  it('`groups.keep` writes the folds into `collapsed` (one undo step); an edit of `collapsed` resets the session’s set', async () => {
    const v = await make();
    await v.execute('groups.toggle', { address: `${GroupsSection}/signals` });
    await v.execute('groups.keep', { section: GroupsSection });
    await v.settle();
    const groups = v.spec.root.value.children![1].children![0].children![0].children!.filter((c) => c.kind === 'group');
    expect(groups.map((g) => [g.id, g.params?.collapsed ?? false])).toEqual([
      ['features', false],
      ['variants', true],
      ['signals', true],
    ]);
    // the spec now opens with the set shown: the session's set kept, nothing moves (its own patch changes no session)
    expect(folded(v)).toEqual(['variants', 'signals']);
    // a toggle, then the undo of the kept folds: back to the spec's opening set
    await v.execute('groups.toggle', { address: `${GroupsSection}/features` });
    expect((await stackOf(v).undo()).ok).toBe(true);
    await v.settle();
    expect(v.view.get(GroupsSection)?.groups).toBeUndefined();
    expect(folded(v)).toEqual(['variants']);
  });

  it('a group a patch adds with `collapsed: true` opens folded; a renamed folded group stays folded, its old id gone', async () => {
    const v = await make();
    await v.execute('groups.toggle', { address: `${GroupsSection}/signals` });
    expect(folded(v)).toEqual(['variants', 'signals']);
    // renamed: its fold carried to the new id
    const at = ['views', 'v', 'p53'];
    let o = await v.spec.patch('Rename', {
      ops: [{ op: 'move', at: [...at, 'signals'], parent: at, place: { last: true }, rename: 'sig' }],
    });
    expect(o.ok && o.changed).toBe(true);
    await v.settle();
    expect(v.view.get(GroupsSection)?.groups).toEqual(['variants', 'sig']);
    // added folded: the section opens with the spec's folds, the new group among them
    o = await v.spec.patch('Add', {
      ops: [
        {
          op: 'insert',
          parent: at,
          place: { last: true },
          node: {
            kind: 'group',
            id: 'extra',
            params: { label: 'Extra', collapsed: true },
            children: [
              {
                kind: 'track',
                id: 't',
                children: [{ kind: 'layer', id: 'l', params: { type: 'values', data: '/data/a' } }],
              },
            ],
          },
        },
      ],
    });
    expect(o.ok && o.changed).toBe(true);
    await v.settle();
    expect(folded(v)).toContain('extra');
  });

  it('a reveal of an item only a hidden group draws says it is on a hidden track', async () => {
    const v = await make();
    await v.execute('track.hide', { address: `${GroupsSection}/variants` });
    await v.settle();
    const r = v.addresses.reveal('/data/clinvar/(item :label clinvar17)');
    expect([r.done, r.message]).toEqual([false, HiddenTrackMessage]);
  });

  it('a folded row’s header slots act on the group: its hide slot hides the group', async () => {
    const v = await make();
    await v.execute('groups.toggle', { address: `${GroupsSection}/features` });
    await (v.behaviors.get('headers') as { loaded(): Promise<void> }).loaded();
    await v.settle();
    frames!.flush();
    const geo = () => v.headers!.value.find((g) => g.target.track === `${GroupsSection}/features`);
    // the row's line is the group's: no line on its summary band
    expect(v.headers!.value.some((g) => g.target.track?.includes('#'))).toBe(false);
    const g = geo()!;
    v.input.move(g.label.x + 2, g.label.y + g.label.height / 2);
    frames!.flush();
    const slot = geo()!.slots.find((s) => s.name === 'seq:hide')!;
    v.input.move(slot.box.x + slot.box.width / 2, slot.box.y + slot.box.height / 2);
    frames!.flush();
    v.input.click(slot.box.x + slot.box.width / 2, slot.box.y + slot.box.height / 2);
    await v.settle();
    expect(v.hiddenTracks.value.map((t) => [t.address, t.kind])).toEqual([[`${GroupsSection}/features`, 'group']]);
  });

  it('a click on a group’s header folds and opens it (no selection, no step)', async () => {
    const v = await make();
    const steps = backSteps(v);
    const selected = v.selection.value.value;
    const header = () => v.frame!.layout.sections[0].groups!.find((g) => g.group.address.endsWith('/features'))!.box;
    let b = header();
    v.input.click(b.x + 20, b.y + b.height / 2);
    expect(folded(v)).toEqual(['variants', 'features']);
    await v.settle();
    frames!.flush();
    b = header();
    v.input.click(b.x + 20, b.y + b.height / 2);
    expect(folded(v)).toEqual(['variants']);
    expect(v.selection.value.value).toBe(selected);
    expect(backSteps(v)).toBe(steps);
  });

  it('`layer.lanes` sets a ranges layer’s `lanes` (one undo step); `packed` unsets it; another layer is refused', async () => {
    const v = await make();
    const layer = `${GroupsSection}/features/domains/l`;
    const domains = () =>
      v.spec.root.value.children![1].children![0].children![0].children![1].children![0].children![0];
    await v.execute('layer.lanes', { address: layer, lanes: 'one' });
    await v.settle();
    expect(domains().params?.lanes).toBe('one');
    await v.execute('layer.lanes', { address: layer, lanes: 'packed' });
    await v.settle();
    expect(domains().params?.lanes).toBeUndefined();
    expect(await v.execute('layer.lanes', { address: `${GroupsSection}/signals/a/l`, lanes: 'one' })).toBe(
      'No ranges layer there',
    );
  });
});
