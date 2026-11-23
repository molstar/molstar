/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a reader without the extension draws (EC-235, _spec/README.md, The fallback): the same spec, the ribbons as
// core's `pairs` (the same bins, chords), canvas and SVG alike, the issue once with the version said; a viewer that
// registers it late (`viewer.extend`) draws the ribbons, and its undo the chords again. The viewer's parts: the
// double-click scoped to the ribbons opens their chromosomes (`rondo:open`, one patch), `rondo:partner` walks the
// cursor along the links by count.

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { afterEach, describe, expect, it } from 'vitest';
import { rondo } from '../index.js';
import type { RondoLay } from '../lay.js';
import { rondoViewer } from '../viewer.js';
import { ringSpec, Size, settle } from './helpers.js';

const Bed = { kind: 'format-runtime' as const, runtime: BedRuntime };
const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
});

const typesOf = (f: { layout: { sections: readonly { sectionLays: readonly { type: { name: string } }[] }[] } }) =>
  f.layout.sections.flatMap((s) => s.sectionLays.map((l) => l.type.name));

describe('without the extension (EC-235)', () => {
  it('draws the ribbons as core’s pairs, canvas and SVG alike, the issue once with the version said', async () => {
    const spec = ringSpec();
    const f = await (await settle(spec, [])).frame({ size: Size });
    expect(typesOf(f)).toEqual(['pairs']);
    const issues = f.layout.issues.filter((i) => i.code === 'seq.view.layer.extension');
    const ctx = createSeqContext({ registry: [Bed] });
    const view = ctx.view(ctx.resolve(spec))!;
    const said = view.issues.filter((i) => i.code === 'seq.view.layer.extension');
    expect(said.map((i) => i.message)).toEqual(['1 layer(s) need `rondo` 0.1: drawn as `pairs`']);
    expect(issues.length).toBeLessThanOrEqual(1);
    // the SVG export draws the same lists as the canvas's frame
    const svg = await ctx.renderSvg(spec, { size: Size });
    expect(svg.lists.passes.data.paths.n).toBe(f.lists.passes.data.paths.n);
    expect(svg.issues.filter((i) => i.code === 'seq.view.layer.extension').length).toBe(1);
    // with it: no such issue, the ribbons drawn
    const withIt = await createSeqContext({ extensions: [rondo], registry: [Bed] }).renderSvg(spec, { size: Size });
    expect(withIt.issues.filter((i) => i.code === 'seq.view.layer.extension')).toEqual([]);
    expect(withIt.svg).not.toEqual(svg.svg);
  });
});

async function viewerOf(spec: ReturnType<typeof ringSpec>, extensions: readonly never[] = []) {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    registry: [...(base.registry ?? []), Bed],
    extensions,
    ui: undefined,
  });
  viewers.push(v);
  await v.load(spec);
  v.resize(720, 760);
  await v.settle();
  v.draw();
  return v;
}

const layOf = (v: SeqViewer) =>
  v.frame?.layout.sections.flatMap((s) => s.sectionLays).find((l) => l.type.name === 'rondo:ribbons')?.lay as
    | RondoLay
    | undefined;

describe('the viewer', () => {
  it('`viewer.extend` draws the ribbons; its undo the chords again', async () => {
    const v = await viewerOf(ringSpec());
    expect(v.frame && typesOf(v.frame)).toEqual(['pairs']);
    const undo = await v.extend(rondoViewer);
    await v.settle();
    v.draw();
    expect(v.frame && typesOf(v.frame)).toEqual(['rondo:ribbons']);
    expect(layOf(v)!.nr).toBeGreaterThan(0);
    undo();
    await v.settle();
    v.draw();
    expect(v.frame && typesOf(v.frame)).toEqual(['pairs']);
  });

  it('a double-click on a ribbon opens its chromosomes in the section `opens` names, one patch; not a focus', async () => {
    const v = await viewerOf(ringSpec({ opens: true }), [rondoViewer as never]);
    const l = layOf(v)!;
    const r = l.strongest;
    const P = l.polys;
    const i = l.ribPoly[r];
    // a point inside the ribbon: its polygon's middle point toward its centroid until it hits
    let x = 0;
    let y = 0;
    for (let k = 0; k < P.len[i]; k++) {
      x += P.xy[2 * (P.start[i] + k)];
      y += P.xy[2 * (P.start[i] + k) + 1];
    }
    x /= P.len[i];
    y /= P.len[i];
    let hit = v.hit(x, y);
    for (let k = 0; k < P.len[i] && hit?.element?.key !== l.ribKey[r]; k++) {
      const px = P.xy[2 * (P.start[i] + k)];
      const py = P.xy[2 * (P.start[i] + k) + 1];
      hit = v.hit((px * 3 + x) / 4, (py * 3 + y) / 4);
      if (hit?.element?.key === l.ribKey[r]) {
        x = (px * 3 + x) / 4;
        y = (py * 3 + y) / 4;
      }
    }
    expect(hit?.element?.key).toBe(l.ribKey[r]);
    const names = [l.names[l.stemSystem[l.ribA[r]]], l.names[l.stemSystem[l.ribB[r]]]];
    expect(v.input.doubleClick(x, y)).toBe(true);
    await v.settle();
    const section = v.view.resolved.value!.sections[0];
    expect(section.blocks!.map((b) => b.label).sort()).toEqual([...new Set(names)].sort());
  });

  it('`rondo:partner` moves the cursor to the strongest partner, then the next by count', async () => {
    const v = await viewerOf(ringSpec(), [rondoViewer as never]);
    const l = layOf(v)!;
    const r = l.strongest;
    const a = l.ribA[r];
    v.selection.setCursor({ section: l.section, system: l.systems[l.stemSystem[a]], position: l.stemStart[a] });
    const said = await v.execute('rondo:partner', {});
    const c = v.selection.cursor.value!;
    const b = l.ribB[r];
    expect([c.system, c.position]).toEqual([l.systems[l.stemSystem[b]], l.stemStart[b]]);
    expect(said).toMatch(/^chr\d, [\d,]+ · [\d,]+ links/);
  });
});
