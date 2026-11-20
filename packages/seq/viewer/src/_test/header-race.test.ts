/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `headers` behavior before its late parts land (_spec/README.md "Behaviors", Headers): a click on a slot or a
// chip is not lost, it runs once they land; and the track's "n filtered" chip counts every thinned layer's rows and
// clears every layer's filter in one patch. The late parts' chunk is held back here (a mock whose factory waits for a
// gate): this file runs isolated.

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

/** Held until opened: the late parts' module evaluates (once, the first viewer's import) only then. */
const gate = vi.hoisted(() => {
  let open = () => {};
  const wait = new Promise<void>((done) => {
    open = done;
  });
  return { wait, open: () => open() };
});
vi.mock('../behaviors/header-extras.js', async (importOriginal) => {
  await gate.wait;
  return importOriginal();
});

const Section = '/views/main/protein';
const Domains = `${Section}/domains`;
const Gone = `${Section}/gone`;

function spec() {
  const b = createSeqSpec({ title: 'header race' });
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
    s.track('more', { label: 'More' }).layer('l', { type: 'letters' });
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
  gate.open();
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
});

async function make(): Promise<SeqViewer> {
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[SeqCoreConfig.render, { labelWidth: 140 }]],
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  return v;
}

/** A click at the middle of a part (a slot or a chip) of a header line as laid out (a track's, else the section's). */
function clickPart(v: SeqViewer, track: string | undefined, name: string) {
  frames.flush();
  const s = v.frame!.layout.sections[0];
  const line = track === undefined ? s.header! : s.tracks.find((x) => x.track.address === track)!.header!;
  const box = line.slots.find((s) => s.part.name === name) ?? line.states.find((s) => s.name === name)?.box;
  if (!box) throw new Error(`no ${name} on ${track}`);
  v.input.click(box.x + box.width / 2, box.y + box.height / 2);
}

describe('a click before the late parts land', () => {
  it('runs once they land: the section’s "1 hidden" chip (core’s built-in states meanwhile) shows its track', async () => {
    const v = await make();
    clickPart(v, undefined, 'seq:hidden-tracks');
    // not landed: nothing ran yet
    expect(v.view.resolved.value!.sections[0].tracks.map((t) => t.address)).not.toContain(Gone);
    gate.open();
    await v.settle();
    expect(v.view.resolved.value!.sections[0].tracks.map((t) => t.address)).toContain(Gone);
  });
});

describe('the track’s "n filtered" chip', () => {
  it('counts every thinned layer’s rows and clears every layer’s filter in one patch', async () => {
    const v = await make();
    await v.settle();
    frames.flush();
    const line = v.frame!.layout.sections[0].tracks.find((x) => x.track.address === Domains)!.header!;
    // two rows of x out of `d`, one of y out of `e`
    expect(line.states.find((s) => s.name === 'seq:filtered')!.text).toBe('3 filtered');
    let patches = 0;
    const stop = v.spec.revision.listen(() => patches++);
    clickPart(v, Domains, 'seq:filtered');
    await v.settle();
    stop();
    const track = v.spec.root.value
      .children!.find((c) => c.id === 'views')!
      .children![0].children![0].children!.find((c) => c.id === 'domains')!;
    expect(track.children!.map((l) => l.params!.data)).toEqual(['/data/domains', '/data/domains']);
    expect(patches).toBe(1);
  });
});
