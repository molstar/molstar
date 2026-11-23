/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { loadGenome } from '@molstar/seqx-core/resolve';
import { partnersOf, rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { Announcer } from '@molstar/seqx-viewer/announcer';
import type { HeadersBehavior } from '@molstar/seqx-viewer/behaviors/headers';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { ribbonsOf, ringStoryView, strongestStem } from '../views.js';

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
    // the partners by count, strongest first, each a pair
    const partners = partnersOf(l, stem);
    const items = element[0].items!;
    expect(items.length).toBe(Math.min(8, partners.length));
    const name = l.names[l.stemSystem[stem]];
    const first = l.names[l.stemSystem[partners[0].stem]];
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

describe('The partner keys', () => {
  it('`.` and `,` walk the cursor’s stem’s partners by count, said', async () => {
    const { v } = await open();
    const announcer = new Announcer(v, { minGapMs: 0, restMs: 0 });
    const said: string[] = [];
    announcer.message.listen((m) => said.push(m.text));
    const l = ribbonsOf(v)!;
    const { stem } = strongestStem(l);
    v.selection.setCursor({ section: l.section, system: l.systems[l.stemSystem[stem]], position: l.stemStart[stem] });
    const partners = partnersOf(l, stem);
    expect(v.input.key({ key: '.' })).toBe(true);
    const at = (k: number) => l.systems[l.stemSystem[partners[k].stem]];
    expect(v.selection.cursor.value?.system).toBe(at(0));
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
      expect(lb.names[lb.stemSystem[sb.stem]]).toBe(name);
      const partnersB = partnersOf(lb, sb.stem);
      expect(partnersB.length).toBeLessThan(partnersA.length);
      expect(headers.context(sb.x, sb.y)).toBe(true);
      const items = headers.menu.value!.menu.groups[0][0].items!;
      expect(items.length).toBe(Math.min(8, partnersB.length));
      // b's partners only: on the two chromosomes it shows (a's are all over the genome)
      for (const i of items) expect([name, first]).toContain(String(i.label).split(' · ')[0]);
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
