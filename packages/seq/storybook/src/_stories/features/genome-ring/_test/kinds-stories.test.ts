/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { ColorMode, describeRondo, Kind, type RondoLay } from '@molstar/seqx-ext-rondo';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { hideKind, kindsOf, ringStoryView } from '../views.js';

// Features/Genome Ring's kinds of links, headless on a fake clock as the stories compose them (a viewer with
// `rondoViewer`, its SVG export light and dark): the contacts faint in the hues and the made-up translocations in the
// accent over them, at the fit, on chr9 and on chr9 × chr22; each kind's legend entry; each hidden in turn.

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open(id: string) {
  viewer?.retire();
  viewer = await openExample(ringStoryView(id), { delayMs: 0, wait: 3000, fixtures: true, extensions: [rondoViewer] });
  viewer.draw();
  return viewer;
}

/** The view's SVG export light and dark: settled, no error, two pictures. */
async function exported(v: SeqViewer) {
  const svg = async (dark: boolean) => {
    const r = v.renderSvg({ idPrefix: 'k-', tokens: dark ? DarkTokens : LightTokens });
    await vi.advanceTimersByTimeAsync(400);
    return r;
  };
  const light = await svg(false);
  const dark = await svg(true);
  expect(light.settled).toBe(true);
  expect(light.issues.filter((i) => i.severity !== 'info')).toEqual([]);
  expect(dark.svg).not.toBe(light.svg);
  return { light, dark };
}

const errors = (v: SeqViewer) => [
  ...v.store.issues.value.filter((i) => i.severity === 'error'),
  ...v.layoutIssues.value.filter((i) => i.severity !== 'info'),
];

/** The ribbons (or stubs) of a lay between two chromosomes (a stub: from the first toward the second). */
function between(l: RondoLay, a: string, b: string) {
  const name = (s: number) => l.names[l.stemSystem[s]];
  const ribbons: number[] = [];
  for (let r = 0; r < l.nr; r++) {
    const x = name(l.ribA[r]);
    const y = name(l.ribB[r]);
    if ((x === a && y === b) || (x === b && y === a)) ribbons.push(r);
  }
  const stubs: number[] = [];
  for (let s = 0; s < l.nb; s++)
    if (name(l.stubStem[s]) === a && l.stubPartner[s] >= 0 && l.names[l.stubPartner[s]] === b) stubs.push(s);
  return { ribbons, stubs };
}

describe('Kinds of links, headless', () => {
  it('at the fit: faint contacts in the hues, the translocations in the accent over them; a legend entry each', async () => {
    const v = await open('ring-kinds');
    expect(errors(v)).toEqual([]);
    const [contacts, tx] = kindsOf(v);
    expect([contacts.label, tx.label]).toEqual(['Contacts', 'Translocations']);
    expect(contacts.mode).toBe(ColorMode.system);
    expect(contacts.opacity).toBe(0.45);
    expect(contacts.items).toBe(false);
    expect(tx.mode).toBe(ColorMode.fixed);
    expect(tx.one >>> 8).toBe(LightTokens.colors['msui-accent'] >>> 8);
    expect(tx.items).toBe(true);
    // the five made-up translocations, each a ribbon between its two chromosomes
    for (const [a, b] of [
      ['chr9', 'chr22'],
      ['chr1', 'chr17'],
      ['chr4', 'chr7'],
      ['chr12', 'chr18'],
      ['chr2', 'chr11'],
    ])
      expect(between(tx, a, b).ribbons.length, `${a}–${b}`).toBeGreaterThan(0);
    const titles = v.legends.legends.value.map((l) => l.title);
    expect(titles).toContain('Contacts');
    expect(titles).toContain('Translocations');
    // the BCR-ABL1-like junction's card: its kind, strands and genes (made up)
    const r = between(tx, 'chr9', 'chr22').ribbons[0];
    const lines = describeRondo(tx.hits[tx.ribPoly[r]] as never, tx)!;
    expect(lines).toContain('Type: translocation');
    expect(lines.some((x) => /^gene_a: XABL$|^gene_b: XABL$/.test(x))).toBe(true);
    await exported(v);
  });

  it('on chr9: its contacts as ribbons, its translocation to chr22 a stub in the accent', async () => {
    const v = await open('ring-kinds-chromosome');
    expect(errors(v)).toEqual([]);
    const [contacts, tx] = kindsOf(v);
    expect(contacts.nr).toBeGreaterThan(0);
    const { stubs } = between(tx, 'chr9', 'chr22');
    expect(stubs.length).toBeGreaterThan(0);
    expect(tx.polyColor[tx.stubPoly[stubs[0]]] >>> 8).toBe(LightTokens.colors['msui-accent'] >>> 8);
    await exported(v);
  });

  it('on chr9 × chr22: the contacts between them, the junction a ribbon over them', async () => {
    const v = await open('ring-kinds-pair');
    expect(errors(v)).toEqual([]);
    const [contacts, tx] = kindsOf(v);
    expect(between(contacts, 'chr9', 'chr22').ribbons.length).toBeGreaterThan(0);
    expect(between(tx, 'chr9', 'chr22').ribbons.length).toBeGreaterThan(0);
    // drawn after (on top of) the contacts: the translocations' layer is the later section layer
    const order = v.frame!.layout.sections.flatMap((s) => s.sectionLays.map((l) => l.lay));
    expect(order.indexOf(tx)).toBeGreaterThan(order.indexOf(contacts));
    expect(tx.polys.n).toBeGreaterThan(0);
    expect(tx.polyKind.subarray(0, tx.polys.n).includes(Kind.ribbon)).toBe(true);
    await exported(v);
  });

  it('each kind hidden in turn: the other alone, drawn as before', async () => {
    let v = await open('ring-kinds');
    const before = kindsOf(v).map((l) => [l.label, l.nr]);
    await hideKind('contact')(v);
    await v.settle();
    v.draw();
    expect(kindsOf(v).map((l) => [l.label, l.nr])).toEqual([before[1]]);
    v = await open('ring-kinds');
    await hideKind('translocation')(v);
    await v.settle();
    v.draw();
    expect(kindsOf(v).map((l) => [l.label, l.nr])).toEqual([before[0]]);
    await exported(v);
  });
});
