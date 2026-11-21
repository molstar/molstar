/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The "n filtered" chip in a viewer without the `headers` behavior (_spec/README.md "Input"; "Behaviors", Headers:
// Entries): a minimal composition listing the one entry `FilteredHeaderEntry` (and the `layer.filter` command's
// behavior): a click on the chip clears every filtered layer's filter in one patch, undone in one step; the section's
// hidden track draws nothing (the application lists it).

import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { LettersType, RangesType, TracksLayout } from '@molstar/seqx-core/render';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { HoverProvider, LegendProvider, NavigateProvider } from '../behaviors/index.js';
import { FilteredHeaderEntry } from '../filter-chip.js';
import { defaultPorts } from '../ports/index.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Domains = '/views/main/protein/domains';

function spec() {
  const b = createSeqSpec({ title: 'filter chip' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
  const d = b.data.ranges('domains', {
    target: p,
    rows: {
      start: new Int32Array([0, 12, 20]),
      end: new Int32Array([10, 18, 26]),
      id: ['A', 'B', 'C'],
      category: ['x', 'y', 'x'],
    },
  });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters', { label: 'Letters' }).layer('l', { type: 'letters' });
    const t = s.track('domains', { label: 'Domains' });
    // two layers, each thinned by a filter of its own
    t.layer('d', { type: 'ranges', data: `${d}/(except :of (item :category x))` });
    t.layer('e', { type: 'ranges', data: `${d}/(except :of (item :category y))` });
    s.track('gone', { label: 'Gone', hidden: true }).layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames;
const viewers: SeqViewer[] = [];
beforeEach(() => {
  frames = installFrames();
});
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
});

/** The minimal composition (the size script's), with the legends' behavior for `layer.filter` and a history. */
async function make(): Promise<SeqViewer> {
  const v = createSeqViewer({
    registry: [
      { kind: 'format-runtime', runtime: SeqRuntime },
      { kind: 'layer-type', type: LettersType },
      { kind: 'layer-type', type: RangesType },
      { kind: 'layout', layout: TracksLayout },
      { kind: 'header-action', entry: FilteredHeaderEntry },
    ],
    behaviors: [{ provider: HoverProvider }, { provider: NavigateProvider }, { provider: LegendProvider }],
    ports: defaultPorts(),
    history: uiHistory(),
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const layersData = (v: SeqViewer) =>
  v.spec.root.value
    .children!.find((c) => c.id === 'views')!
    .children![0].children![0].children!.find((c) => c.id === 'domains')!
    .children!.map((l) => l.params!.data);
const mid = (b: { x: number; y: number; width: number; height: number }) => [b.x + b.width / 2, b.y + b.height / 2];

describe('the "n filtered" chip without the `headers` behavior', () => {
  it('its click clears every filtered layer in one patch, undone in one step; the other chips are facts', async () => {
    const v = await make();
    expect(v.behaviors.get('headers')).toBeUndefined();
    const s = v.frame!.layout.sections[0];
    const chip = s.tracks.find((t) => t.track.address === Domains)!.header!.states[0];
    expect(chip.text).toBe('3 filtered');
    const before = layersData(v);
    const stack = (v.spec.history as UiHistory).stack;
    const steps = stack.steps.value.length;
    // the hidden track: the view says nothing of it (no section line, no chip); the application's list restores it
    expect(s.header).toBeUndefined();
    expect(v.hiddenTracks.value.map((t) => t.label)).toHaveLength(1);
    expect(v.view.resolved.value!.sections[0].tracks).toHaveLength(2);
    // the chip: both layers' filters cleared, one step
    v.input.click(...(mid(chip.box) as [number, number]));
    await v.settle();
    expect(layersData(v)).toEqual(['/data/domains', '/data/domains']);
    expect(stack.steps.value.length).toBe(steps + 1);
    await stack.undo();
    await v.settle();
    expect(layersData(v)).toEqual(before);
  });
});
