/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { loadGenome } from '@molstar/seqx-core/resolve';
import type { RondoLay } from '@molstar/seqx-ext-rondo';
import { partnerChromosomes, partnersOf, rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { Announcer } from '@molstar/seqx-viewer/announcer';
import type { HeadersBehavior } from '@molstar/seqx-viewer/behaviors/headers';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { kindsOf, ribbonsOf, ringStoryView, strongestStem } from '../views.js';
import { bandsLanded } from './ring-helpers.js';

// The stem menu and the partner keys (`rondoViewer`): a secondary click on a stem opens its element menu (the label
// menu's shape, scoped to `rondo:ribbons`' stems): its partners by count shown as a pair, the stem and its partners
// selected, its chromosome alone; `.` and `,` walk the cursor's stem's partners by count, said.

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open() {
  viewer = await openExample(ringStoryView('ring'), {
    delayMs: 0,
    wait: 3000,
    fixtures: true,
    extensions: [rondoViewer],
  });
  viewer.draw();
  const headers = viewer.behavior<HeadersBehavior>('headers')!;
  await headers.loaded();
  return { v: viewer, headers };
}

const ringInclude = (v: SeqViewer) => {
  const views = v.spec.root.value.children!.find((c) => c.kind === 'views')!;
  return views.children![0].children!.find((c) => c.id === 'ring')!.params?.include;
};

describe('The stem menu', () => {
  it('a secondary click on a stem: its partners by count, its partners selected, its chromosome alone', async () => {
    const { v, headers } = await open();
    const l = ribbonsOf(v)!;
    const { stem, x, y } = strongestStem(l);
    expect(headers.context(x, y)).toBe(true);
    const menu = headers.menu.value!.menu;
    expect(menu.target.element?.type).toBe('rondo:ribbons');
    const element = menu.groups[0];
    expect(element.map((i) => i.name)).toEqual([
      'rondo:show-partner',
      'rondo:select-stem-partners',
      'rondo:show-chromosome',
    ]);
    // the partners by chromosome, most links first, each a pair; the stem's own chromosome last, within
    const partners = partnersOf(l, stem);
    const rows = partnerChromosomes(l, stem);
    const items = element[0].items!;
    expect(items.length).toBe(rows.length);
    const name = l.names[l.stemSystem[stem]];
    const first = l.names[rows[0].system];
    expect(items[0].label).toMatch(new RegExp(`^${first} · [\\d,]+ links`));
    // the partners selected: the stem's range and each partner's
    expect(await headers.runItem(element[1])).toBe(
      `${partners.length} partner${partners.length === 1 ? '' : 's'} selected`,
    );
    expect(v.selection.value.value.shown.parts.length).toBeGreaterThan(0);
    // its chromosome alone is offered by name; the strongest partner shown with it as a pair
    expect(element[2].label).toBe(`Show only ${name}`);
    await headers.runItem(items[0]);
    expect(new Set(ringInclude(v) as string[])).toEqual(new Set([name, first]));
  });

  it('only the stem’s actions: the track’s (Move up, Reverse, Hide track) stay in its header’s menu', async () => {
    const { v, headers } = await open();
    const l = ribbonsOf(v)!;
    const { x, y } = strongestStem(l);
    expect(headers.context(x, y)).toBe(true);
    const menu = headers.menu.value!.menu;
    expect(menu.groups).toHaveLength(1);
    expect(menu.groups[0].map((i) => i.name)).toEqual([
      'rondo:show-partner',
      'rondo:select-stem-partners',
      'rondo:show-chromosome',
    ]);
    headers.closeMenu();
    // the track's header menu keeps them
    const track = l.key.slice(0, l.key.lastIndexOf('/'));
    const names = headers
      .menuOf({ section: l.section, track })!
      .groups.flat()
      .map((i) => i.name);
    expect(names).toEqual(expect.arrayContaining(['seq:hide', 'seq:move-up', 'seq:reverse']));
    // an element with no entries of its own (a ribbon) keeps its track's menu
    const r = l.strongest;
    const ribbon = headers.menuOf({
      section: l.section,
      track,
      element: { layer: l.key, key: l.ribKey[r], type: 'rondo:ribbons' },
    })!;
    expect(ribbon.groups.flat().map((i) => i.name)).toContain('seq:hide');
  });

  it('“Show only” shows the stem’s chromosome alone', async () => {
    const { v, headers } = await open();
    const l = ribbonsOf(v)!;
    const { stem, x, y } = strongestStem(l);
    expect(headers.context(x, y)).toBe(true);
    const only = headers.menu.value!.menu.groups[0][2];
    await headers.runItem(only);
    expect(ringInclude(v)).toEqual([l.names[l.stemSystem[stem]]]);
  });

  it('elsewhere (a ribbon, a band) no stem entries', async () => {
    const { v, headers } = await open();
    const l = ribbonsOf(v)!;
    // the ring's centre: a ribbon or nothing, never a stem
    headers.context(l.frame.cx, l.frame.cy);
    const names = headers.menu.value?.menu.groups.flat().map((i) => i.name) ?? [];
    expect(names.filter((n) => n.startsWith('rondo:'))).toEqual([]);
  });
});

/** The stem menu of a lay's stem (a layer under another: its menu by its own element), its partner rows. */
function partnerRowsOf(headers: HeadersBehavior, l: RondoLay, stem: number) {
  const track = l.key.slice(0, l.key.lastIndexOf('/'));
  const element = { layer: l.key, key: l.stemKey[stem], type: 'rondo:ribbons' };
  expect(headers.openMenu({ section: l.section, track, element }, { x: 0, y: 0 })).toBe(true);
  const rows = headers.menu.value!.menu.groups[0].find((i) => i.name === 'rondo:show-partner')!.items!;
  headers.closeMenu();
  return rows;
}

describe('The stem menu’s partners, by chromosome (Kinds of links)', () => {
  it('a row a partner chromosome with its links; its pieces by band range and strand; within last', async () => {
    viewer = await openExample(ringStoryView('ring-kinds'), {
      delayMs: 0,
      wait: 3000,
      fixtures: true,
      extensions: [rondoViewer],
    });
    const v = viewer;
    await bandsLanded(v);
    const headers = v.behavior<HeadersBehavior>('headers')!;
    await headers.loaded();
    // the contacts' strongest stem: chr22's first (under the translocations' stem there, so opened by its element)
    const contacts = kindsOf(v)[0];
    const { stem } = strongestStem(contacts);
    expect(contacts.names[contacts.stemSystem[stem]]).toBe('chr22');
    const rows = partnerRowsOf(headers, contacts, stem);
    const labels = rows.map((r) => r.label);
    // one row a chromosome: chr9 once, its links the sum of its pieces, most first
    const chr9 = rows.filter((r) => /^chr9 /.test(r.label));
    expect(chr9).toHaveLength(1);
    expect(labels[0]).toBe(chr9[0].label);
    const pieces = partnersOf(contacts, stem).filter((p) => contacts.names[contacts.stemSystem[p.stem]] === 'chr9');
    expect(pieces.length).toBeGreaterThan(1);
    const total = pieces.reduce((n, p) => n + p.links, 0);
    expect(chr9[0].label).toBe(`chr9 · ${total.toLocaleString('en-US')} links`);
    // its pieces in a submenu, named by band range and strand, by count: q34.1-like names, never a bare chr9
    const sub = chr9[0].items!.map((i) => i.label);
    expect(sub).toHaveLength(pieces.length);
    for (const s of sub) expect(s).toMatch(/^[pq][\d.]+(–[pq][\d.]+)? · [\d,]+ links?(, reversed)?$/);
    expect(sub[0]).toMatch(new RegExp(` · ${pieces[0].links.toLocaleString('en-US')} links$`));
    // each piece shows the same pair as its row
    for (const i of chr9[0].items!) expect(i.args).toEqual(chr9[0].args);
    // chr22 is never its own partner: "Within chr22", listed last, its two pieces (its other stem, this stem's own)
    expect(labels.filter((t) => /^chr22 /.test(t))).toEqual([]);
    const within = rows.at(-1)!;
    expect(within.label).toMatch(/^Within chr22 · [\d,]+ links$/);
    expect(within.items!.map((i) => i.label)).toHaveLength(2);
    expect(within.items!.some((i) => /, this stem$/.test(i.label))).toBe(true);
    expect(within.args).toMatchObject({ include: ['chr22'] });
    // a chromosome reached by one ribbon is a plain row (no submenu)
    for (const r of rows) if (!r.items) expect(r.label).toMatch(/^chr\w+ · [\d,]+ links?(, reversed)?$/);
  });

  it('the translocations’ chr22 stem: chr9 a plain row, reversed', async () => {
    viewer = await openExample(ringStoryView('ring-kinds'), {
      delayMs: 0,
      wait: 3000,
      fixtures: true,
      extensions: [rondoViewer],
    });
    const v = viewer;
    const headers = v.behavior<HeadersBehavior>('headers')!;
    await headers.loaded();
    const tx = kindsOf(v)[1];
    const stem = Array.from({ length: tx.ns }, (_, s) => s).find(
      (s) => tx.stemPlaced[s] && tx.names[tx.stemSystem[s]] === 'chr22' && tx.stemPoly[s] >= 0,
    )!;
    const rows = partnerRowsOf(headers, tx, stem);
    expect(rows.map((r) => [r.label, !!r.items])).toEqual([['chr9 · 24 links, reversed', false]]);
  });

  it('bounded on a crowded ring: at most 8 partner rows and 8 pieces a row, the within row besides', async () => {
    viewer = await openExample(ringStoryView('ring-crowded'), {
      delayMs: 0,
      wait: 3000,
      fixtures: true,
      extensions: [rondoViewer],
    });
    const l = ribbonsOf(viewer)!;
    let most = 0;
    for (let s = 0; s < l.ns; s++) {
      if (!l.stemPlaced[s]) continue;
      const rows = partnerChromosomes(l, s);
      most = Math.max(most, partnersOf(l, s).length);
      expect(rows.filter((r) => !r.within).length).toBeLessThanOrEqual(8);
      expect(rows.filter((r) => r.within).length).toBeLessThanOrEqual(1);
      for (const r of rows) expect(r.pieces.length).toBeLessThanOrEqual(8);
    }
    // the bound bites: a stem there has more partners than the menu lists
    expect(most).toBeGreaterThan(8);
  });
});

describe('The partner keys', () => {
  it('`.` and `,` walk the stem menu’s partners in its order (the menu and the keys read one list), said', async () => {
    const { v, headers } = await open();
    const announcer = new Announcer(v, { minGapMs: 0, restMs: 0 });
    const said: string[] = [];
    announcer.message.listen((m) => said.push(m.text));
    const l = ribbonsOf(v)!;
    const { stem } = strongestStem(l);
    // the menu's rows and their pieces, in order: the keys' steps
    const rows = partnerRowsOf(headers, l, stem);
    const steps = partnerChromosomes(l, stem).flatMap((c) => c.pieces.filter((p) => p.stem !== stem));
    expect(rows.map((r) => r.label.split(' · ')[0].replace(/^Within /, ''))).toEqual(
      partnerChromosomes(l, stem).map((c) => l.names[c.system]),
    );
    v.selection.setCursor({ section: l.section, system: l.systems[l.stemSystem[stem]], position: l.stemStart[stem] });
    const partners = steps;
    expect(v.input.key({ key: '.' })).toBe(true);
    const at = (k: number) => l.systems[l.stemSystem[partners[k].stem]];
    expect(v.selection.cursor.value?.system).toBe(at(0));
    expect(rows[0].label.startsWith(l.names[l.stemSystem[partners[0].stem]])).toBe(true);
    expect(said.at(-1)).toMatch(new RegExp(`^${l.names[l.stemSystem[partners[0].stem]]}, [\\d,]+ · [\\d,]+ links`));
    expect(v.input.key({ key: '.' })).toBe(true);
    expect(v.selection.cursor.value?.system).toBe(at(1 % partners.length));
    expect(v.input.key({ key: ',' })).toBe(true);
    expect(v.selection.cursor.value?.system).toBe(at(0));
    // the keys are listed with their commands
    const keys = v.commands.filter((c) => c.id.startsWith('rondo:partner')).map((c) => c.shortcut?.key);
    expect(keys).toEqual(['.', ',']);
    announcer.dispose();
  });
});

describe('The stem menu in two viewers of one spec', () => {
  it('each lists its own viewer’s partners (the menu reads the viewer from its facts, no state of its own)', async () => {
    const { v: a } = await open();
    viewer = undefined;
    const { v: b, headers } = await open();
    try {
      const la = ribbonsOf(a)!;
      const { stem } = strongestStem(la);
      const name = la.names[la.stemSystem[stem]];
      const partnersA = partnersOf(la, stem);
      // b shows the stem's chromosome and its strongest partner only: its partners on those two, a's all over
      const first = la.names[la.stemSystem[partnersA[0].stem]];
      await b.execute('section.include', { address: '/views/main/ring', include: [name, first] });
      await b.settle();
      b.draw();
      const lb = ribbonsOf(b)!;
      const sb = strongestStem(lb);
      // its strongest stem on one of the two it shows (which one is the data's)
      expect([name, first]).toContain(lb.names[lb.stemSystem[sb.stem]]);
      const partnersB = partnersOf(lb, sb.stem);
      expect(partnersB.length).toBeLessThan(partnersA.length);
      expect(headers.context(sb.x, sb.y)).toBe(true);
      const items = headers.menu.value!.menu.groups[0][0].items!;
      expect(items.length).toBe(partnerChromosomes(lb, sb.stem).length);
      // b's partners only: on the two chromosomes it shows (a's are all over the genome)
      for (const i of items)
        expect([name, first]).toContain(
          String(i.label)
            .split(' · ')[0]
            .replace(/^Within /, ''),
        );
    } finally {
      a.retire();
    }
  });

  it('the same stem in both: its menu in the second viewer selects there, never in the first', async () => {
    const { v: a } = await open();
    viewer = undefined;
    const { v: b, headers } = await open();
    try {
      const l = ribbonsOf(b)!;
      const { stem, x, y } = strongestStem(l);
      expect(headers.context(x, y)).toBe(true);
      const select = headers.menu.value!.menu.groups[0][1];
      const n = partnersOf(l, stem).length;
      expect(await headers.runItem(select)).toBe(`${n} partner${n === 1 ? '' : 's'} selected`);
      expect(b.selection.value.value.shown.parts.length).toBeGreaterThan(0);
      expect(a.selection.value.value.shown.parts.length).toBe(0);
    } finally {
      a.retire();
    }
  });
});
