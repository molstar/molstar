/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ring's grouped label in the viewer (core _spec "Render", Rings layout: a dense fan is one label, its fan placed with
// the layout and opened by a hover): the hover on the grouped label is a `group` hit, its card lists the members, a
// click selects them all, and the hover's emphasis draws the fan with no layout.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { RingPlacement } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Ring = '/views/v/ring';

function spec() {
  const n = 2686;
  const rows = { start: [] as number[], end: [] as number[], id: [] as string[], label: [] as string[] };
  for (let k = 0; k < 13; k++) {
    rows.start.push(1081 + 4 * k);
    rows.end.push(1085 + 4 * k);
    rows.id.push(`s${k}`);
    rows.label.push(`Site${k}`);
  }
  const b = createSeqSpec({ title: 'fan' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'A'.repeat(n), cyclic: true });
  const f = b.data.ranges('sites', { target: c, rows });
  b.views.view('v', (v) => {
    v.section('ring', { of: c, layout: 'rings' }).track('sites').layer('sites', { type: 'ranges', data: f });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
});

async function make() {
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewer = v;
  await v.load(spec());
  v.resize(640);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

const ringOf = (v: SeqViewer) =>
  v.frame!.layout.sections.find((s) => s.section.address === Ring)!.placement as RingPlacement;

/** The page point inside the grouped label. */
function groupPoint(v: SeqViewer) {
  const p = ringOf(v);
  const l = p.labels;
  for (let k = 0; k < l.m; k++)
    if (l.of[l.placed[k]] >= 0) return { x: l.right[k] ? l.x[k] + 3 : l.x[k] - 3, y: l.y[k] - 0.35 * p.labelSize };
  throw new Error('no grouped label');
}

describe('a ring’s grouped label', () => {
  it('is a group hit; its card lists the members, a click selects them all', async () => {
    const v = await make();
    const at = groupPoint(v);
    const hit = v.hit(at.x, at.y)!;
    expect(hit.target).toBe('group');
    expect(hit.count).toBe(13);
    // the card
    v.input.move(at.x, at.y);
    const words = new Describer({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    });
    const card = words.describe(v.hover.current)!;
    expect(card.title).toBe('13 items');
    expect(card.rows.length).toBe(13);
    expect(card.rows[0].label).toBe('Site0');
    expect(card.rows[12]).toEqual({ label: 'More', value: '1 more' });
    expect(card.text).toContain('13 items: Site0');
    v.input.leave();
    // a click selects every member
    expect(v.input.click(at.x, at.y)).toBe(true);
    expect(v.selection.value.value.elements.keys.map((k) => k.key).sort()).toEqual(
      Array.from({ length: 13 }, (_, k) => `s${k}`).sort(),
    );
  });

  it('a pointer on it opens the fan, laying nothing out; the hit reads the hover, with no draw between', async () => {
    const v = await make();
    const p = ringOf(v);
    const at = groupPoint(v);
    const layout = v.frame!.layout;
    v.input.move(at.x, at.y);
    // no draw yet: a hit on a fan label already answers as its item
    const slot = p.fan.slot[p.fan.members[3]];
    const pt = { x: p.fan.right[slot] ? p.fan.x[slot] + 3 : p.fan.x[slot] - 3, y: p.fan.y[slot] - 0.35 * p.labelSize };
    expect(v.hit(pt.x, pt.y)?.element?.key).toBe(p.fan.hits[p.fan.members[3]]!.element!.key);
    frames!.flush();
    v.draw();
    expect(v.frame!.layout).toBe(layout);
    expect(v.hit(pt.x, pt.y)?.element?.key).toBe(p.fan.hits[p.fan.members[3]]!.element!.key);
    // away from it: the grouped form again, the fan's place no hit
    v.input.leave();
    expect(v.hit(pt.x, pt.y)?.element?.key).not.toBe(p.fan.hits[p.fan.members[3]]!.element!.key);
  });
});
