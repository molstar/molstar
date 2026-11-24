/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { resolveAddress } from '@molstar/core/repr/tree';
import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { describerFor } from '@molstar/seqx-viewer/describe';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { ColorMode, describeRondo, Kind, type RondoLay } from '@molstar/seqx-ext-rondo';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { legendRows } from '@molstar/seqx-viewer/legends';
import { hideKind, kindLayer, kindsOf, ringStoryView } from '../views.js';

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
    // the showcase's contacts are unstranded: read as `.`, never `+`
    const held = contacts.prep!.source as unknown as {
      data: { n: number };
      datum: { set(k: number): { field(n: string): unknown } };
      shown: Uint8Array;
    };
    const k = held.shown.indexOf(1);
    expect([held.datum.set(k).field('a_strand'), held.datum.set(k).field('b_strand')]).toEqual(['.', '.']);
    expect(tx.mode).toBe(ColorMode.fixed);
    expect(tx.one >>> 8).toBe(LightTokens.colors['msui-accent'] >>> 8);
    expect(tx.items).toBe(true);
    // one kind reads as one colour: its ribbons opaque, never a tint of the hues under them
    for (let r = 0; r < tx.nr; r++) if (tx.ribDrawn[r]) expect(tx.polyColor[tx.ribPoly[r]] & 0xff).toBe(255);
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

  it('the stem card: its first row names what it counts; each layer’s links on its chromosome under a heading', async () => {
    const v = await open('ring-kinds');
    const [contacts, tx] = kindsOf(v);
    // the translocations' chr22 stem: its 24 links, all to chr9
    const s = Array.from({ length: tx.ns }, (_, k) => k).find(
      (k) => tx.stemPlaced[k] && tx.stemPoly[k] >= 0 && tx.names[tx.stemSystem[k]] === 'chr22',
    )!;
    const lines = describeRondo({ element: { owner: tx.owner, key: tx.stemKey[s] }, position: 0 } as never, tx)!;
    expect(lines[1]).toBe('This stem: 24 links to chr9');
    // the chromosome's links on each layer, counted from the rows (each once, whichever end is on chr22)
    const onChr22 = (l: RondoLay) => {
      const d = l.read!.data!;
      const w = l.read!.weight;
      let n = 0;
      for (let k = 0; k < d.n; k++)
        if (l.read!.shown?.[k] !== 0 && (l.names[d.aUnit![k]] === 'chr22' || l.names[d.bUnit![k]] === 'chr22'))
          n += w ? w[k] : 1;
      return n;
    };
    const heading = lines.indexOf('chr22, all links');
    expect(heading).toBeGreaterThan(1);
    expect(lines.slice(heading + 1)).toEqual([
      `Contacts: ${onChr22(contacts).toLocaleString('en-US')} links`,
      `Translocations: ${onChr22(tx)} links`,
    ]);
    // on the card: the heading a row of its own, its value empty
    const card = describerFor({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    }).describe({
      section: tx.section,
      layer: tx.key,
      system: tx.systems[tx.stemSystem[s]],
      position: tx.stemStart[s],
      element: { owner: tx.owner, key: tx.stemKey[s] } as never,
    })!;
    expect(card.rows.slice(0, 2)).toEqual([{ label: 'This stem', value: '24 links to chr9' }, expect.anything()]);
    expect(card.rows).toContainEqual({ label: 'chr22, all links', value: '' });
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

  /** The stub labels of the kinds' layers as boxes, and the wedges (triangles) of their stubs, in page px. */
  function stubScene(v: SeqViewer) {
    const labels: { box: number[]; text: string; own: number[]; ring: number }[] = [];
    const wedges: number[][] = [];
    for (const l of kindsOf(v)) {
      const size = Number.parseFloat(l.font);
      for (let s = 0; s < l.nb; s++) {
        const p = l.stubPoly[s];
        if (p < 0) continue;
        const tri: number[] = [];
        for (let k = 0; k < l.polys.len[p]; k++)
          tri.push(l.polys.xy[2 * (l.polys.start[p] + k)], l.polys.xy[2 * (l.polys.start[p] + k) + 1]);
        wedges.push(tri);
        if (!l.stubLabel[s]) continue;
        const y = l.stubLabelY[s] + 4 - 0.35 * size;
        const hw = l.stubLabelW[s] / 2;
        const x0 = l.stubLabelX[s] - hw;
        const x1 = l.stubLabelX[s] + hw;
        const y0 = y - 0.65 * size;
        const y1 = y + 0.65 * size;
        const far = Math.max(
          ...[
            [x0, y0],
            [x1, y0],
            [x0, y1],
            [x1, y1],
          ].map(([x, yy]) => Math.hypot(x - l.frame.cx, yy - l.frame.cy) - l.frame.r),
        );
        labels.push({ box: [x0, y0, x1, y1], text: l.stubLabel[s], own: tri, ring: far });
      }
    }
    return { labels, wedges };
  }

  /** Whether a box and a triangle meet (separating axes). */
  function meets(box: number[], t: number[]): boolean {
    const [x0, y0, x1, y1] = box;
    const xs = [t[0], t[2], t[4]];
    const ys = [t[1], t[3], t[5]];
    if (Math.max(...xs) <= x0 || Math.min(...xs) >= x1 || Math.max(...ys) <= y0 || Math.min(...ys) >= y1) return false;
    for (let e = 0; e < 3; e++) {
      const nx = ys[(e + 1) % 3] - ys[e];
      const ny = xs[e] - xs[(e + 1) % 3];
      const d = [0, 1, 2].map((i) => nx * xs[i] + ny * ys[i]);
      const c = [nx * x0 + ny * y0, nx * x1 + ny * y0, nx * x0 + ny * y1, nx * x1 + ny * y1];
      if (Math.max(...c) <= Math.min(...d) || Math.min(...c) >= Math.max(...d)) return false;
    }
    return true;
  }

  const scenes: [string, string, string[] | undefined][] = [
    ['one chromosome (chr9)', 'ring-kinds-chromosome', undefined],
    ['the author’s scene: Kinds Of Links zoomed to chr2 from the square', 'ring-kinds', ['chr2']],
  ];
  for (const [what, id, include] of scenes)
    it(`stub labels, ${what}: none over another, none over a wedge, all inside the ring`, async () => {
      const v = await open(id);
      if (include) {
        v.runCommand('section.include', { address: '/views/main/ring', include });
        await vi.advanceTimersByTimeAsync(2000);
        v.draw();
      }
      const { labels, wedges } = stubScene(v);
      expect(labels.length).toBeGreaterThan(0);
      for (let i = 0; i < labels.length; i++) {
        const a = labels[i];
        // inside the ring's stems: the chromosomes' own labels are outside it, never met
        expect(a.ring, a.text).toBeLessThan(0);
        // never over its own wedge
        expect(meets(a.box, a.own), `${a.text} over its wedge`).toBe(false);
        for (let j = i + 1; j < labels.length; j++) {
          const b = labels[j];
          const overlap = a.box[0] < b.box[2] && a.box[2] > b.box[0] && a.box[1] < b.box[3] && a.box[3] > b.box[1];
          expect(overlap, `${a.text} over ${b.text}`).toBe(false);
        }
        // covering no wedge at all, but the strongest label's (it is kept whatever is under it)
        if (i > 0) for (const w of wedges) expect(meets(a.box, w), `${a.text} over a wedge`).toBe(false);
      }
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

  it('read by region (BEDPE, `type` mapped): each kind’s layer its own links, at the item level', async () => {
    const v = await open('ring-kinds-region');
    // the bins land first, then each layer's items (a want of their own): the clock turned until they have
    for (let i = 0; i < 40 && !kindsOf(v).every((l) => l.items); i++) {
      await vi.advanceTimersByTimeAsync(100);
      v.draw();
    }
    expect(errors(v)).toEqual([]);
    const [contacts, tx] = kindsOf(v);
    for (const [l, type] of [
      [contacts, 'contact'],
      [tx, 'translocation'],
    ] as const) {
      expect(l.items, type).toBe(true);
      const d = l.read!.data!;
      let n = 0;
      for (let k = 0; k < d.n; k++) if (l.read!.shown[k]) n++, expect(d.kinds[d.kind[k]]).toBe(type);
      expect(n, type).toBeGreaterThan(0);
    }
    expect(between(tx, 'chr22', 'chr9').stubs.length).toBeGreaterThan(0);
    await exported(v);
  });

  it('the card on a stem names the block, position and band under the pointer (a band once, never the first block’s)', async () => {
    const v = await open('ring-kinds');
    // the cytobands' rows land on real I/O (their pages inflate on a thread pool)
    const landed = () =>
      v
        .frame!.layout.sections[0].lays.flat()
        .some((l) => l.key.endsWith('/bands/bands') && (l.lay as { rows?: { n: number } }).rows!.n > 0);
    for (let i = 0; i < 50 && !landed(); i++) {
      await new Promise((r) => setImmediate(r));
      await vi.advanceTimersByTimeAsync(200);
      v.draw();
    }
    expect(landed()).toBe(true);
    const [, tx] = kindsOf(v);
    const d = describerFor({ ...v, layout: () => v.frame?.layout });
    const seen = new Set<string>();
    for (let s = 0; s < tx.ns; s++) {
      const name = tx.names[tx.stemSystem[s]];
      if (!tx.stemPlaced[s] || seen.has(name)) continue;
      tx.frame.point((tx.stemT0[s] + tx.stemT1[s]) / 2, 2);
      v.input.move(tx.frame.px, tx.frame.py);
      v.draw();
      const rows = d.describe(v.hover.current)!.rows;
      const at = rows.find((r) => r.label === name);
      if (!at) continue;
      seen.add(name);
      const position = Number(at.value.replace(/,/g, ''));
      // the bands the card lists: one row, on this chromosome, holding the position
      const bands = rows.filter((r) => r.item);
      expect(bands.length, name).toBe(1);
      const m = /^(\S+) ([\d,]+)–([\d,]+)$/.exec(bands[0].value)!;
      expect(m[1], name).toBe(name);
      expect(position, name).toBeGreaterThanOrEqual(Number(m[2].replace(/,/g, '')));
      expect(position, name).toBeLessThanOrEqual(Number(m[3].replace(/,/g, '')));
    }
    // several blocks, each its own band
    expect(seen.size).toBeGreaterThan(3);
  });

  it('each kind hidden in turn: the other alone, drawn as before; the hidden kind’s legend row kept in its place, dimmed', async () => {
    let v = await open('ring-kinds');
    const before = kindsOf(v).map((l) => [l.label, l.nr]);
    const titles = legendRows(v).map((r) => r.title);
    await hideKind('contact')(v);
    await v.settle();
    v.draw();
    expect(kindsOf(v).map((l) => [l.label, l.nr])).toEqual([before[1]]);
    // the Contacts row stays where it was (the story "Kinds Of Links Contacts Hidden"): every layer it keys hidden, its
    // eye "Show Contacts", its entries no action; the translocations' row as before
    expect(legendRows(v).map((r) => r.title)).toEqual(titles);
    const contacts = legendRows(v).find((r) => r.title === 'Contacts')!;
    expect(contacts.hidden).toBe(true);
    expect(contacts.layers).toEqual([{ address: kindLayer('contact'), name: 'Contacts', hidden: true }]);
    expect(contacts.entries.every((e) => !e.action)).toBe(true);
    expect(legendRows(v).find((r) => r.title === 'Translocations')!.hidden).toBe(false);
    // the hidden kind is no peer: a stem's card counts no layer apart
    const [alone] = kindsOf(v);
    expect(alone.peers.filter((p) => p !== alone && !p.empty)).toEqual([]);
    const stem = alone.polyKind.subarray(0, alone.polys.n).indexOf(Kind.trumpet);
    expect(describeRondo(alone.hits[stem] as never, alone)!.some((x) => x.startsWith('Contacts:'))).toBe(false);
    v = await open('ring-kinds');
    await hideKind('translocation')(v);
    await v.settle();
    v.draw();
    expect(kindsOf(v).map((l) => [l.label, l.nr])).toEqual([before[0]]);
    await exported(v);
  });

  it('the contacts hidden from the start (the story "Kinds Of Links Contacts Hidden" before its first frame): a Contacts row of its name alone, dimmed, its eye "Show Contacts", in its place', async () => {
    const shown = await open('ring-kinds');
    const titles = legendRows(shown).map((r) => r.title);
    shown.retire();
    viewer = undefined;
    const view = ringStoryView('ring-kinds');
    const build = () => {
      const root = view.build();
      const layer = resolveAddress(root, kindLayer('contact').split('/').filter(Boolean))!;
      (layer as { params: Record<string, unknown> }).params = { ...layer.params, hidden: true };
      return root;
    };
    viewer = await openExample({ ...view, build } as never, {
      delayMs: 0,
      wait: 3000,
      fixtures: true,
      extensions: [rondoViewer],
    });
    const v = viewer;
    v.draw();
    expect(legendRows(v).map((r) => r.title)).toEqual(titles);
    const contacts = legendRows(v).find((r) => r.title === 'Contacts')!;
    expect(contacts.hidden).toBe(true);
    expect(contacts.entries).toEqual([]);
    expect(contacts.layers).toEqual([{ address: kindLayer('contact'), name: 'Contacts', hidden: true }]);
  });
});
