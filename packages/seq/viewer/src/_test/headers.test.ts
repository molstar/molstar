/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `headers` behavior (_spec/README.md "Behaviors", Headers) in Node: each built-in through its slot, its menu item
// and its command alike (one patch, one undo); a host's entry added and removed with its undo, a conflict refused
// whole; a disabled entry says why; the menu's data; the secondary click under each `contextClick` and the menu key
// once per ask (EC-156); the cursor onto a header and back, and off a hidden track; the geometry published on change
// only (a hover publishes nothing); the parts' items the menu's; the export alike under both renderings.

import type { Node } from '@molstar/core/repr/tree';
import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import type { HeadersBehavior } from '../behaviors/headers.js';
import { HeadersProvider } from '../behaviors/headers.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import {
  type HeaderEntry,
  type HeaderMenuItem,
  type HeadersParams,
  headerMenu,
  headerParts,
} from '../header-entries.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/main/protein';
const Letters = `${Section}/letters`;
const Domains = `${Section}/domains`;
const More = `${Section}/more`;

/** A protein with three tracks: letters, domains (two layers, colored by id) and more letters. */
function spec(edit?: (b: ReturnType<typeof createSeqSpec>) => void): Node {
  const b = createSeqSpec({ title: 'headers' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
  const d = b.data.ranges('domains', {
    target: p,
    rows: { start: new Int32Array([0, 15]), end: new Int32Array([10, 25]), id: ['A', 'B'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters', { label: 'Letters' }).layer('l', { type: 'letters' });
    const t = s.track('domains', { label: 'Domains' });
    t.layer('d', { type: 'ranges', data: d });
    t.layer('e', { type: 'ranges', data: d });
    s.track('more', { label: 'More' }).layer('l', { type: 'letters' });
  });
  edit?.(b);
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

async function make(params?: HeadersParams, root = spec()): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    history: uiHistory(),
    // a label column of 140 px: two slots and the "⋯"
    config: [[SeqCoreConfig.render, { labelWidth: 140 }]],
    behaviors: base.behaviors.map((b) => (b.provider === HeadersProvider ? { provider: HeadersProvider, params } : b)),
  });
  viewers.push(v);
  await v.load(root);
  v.resize(600);
  // the headers' late parts (a lazy chunk): `settle` waits for them
  await v.settle();
  frames.flush();
  return v;
}

const headers = (v: SeqViewer) => v.behaviors.get('headers') as HeadersBehavior;
const stack = (v: SeqViewer) => (v.spec.history as UiHistory).stack;
const geo = (v: SeqViewer, track: string | undefined) =>
  v.headers!.value.find((g) => g.target.track === track && (track !== undefined || g.target.section === Section))!;
const mid = (b: { x: number; y: number; width: number; height: number }) => [b.x + b.width / 2, b.y + b.height / 2];
/** The pointer onto a track's label (its actions show), then onto a part of it, and a click there. */
function clickPart(v: SeqViewer, track: string | undefined, name: string, kind: 'slots' | 'states' = 'slots') {
  const g = geo(v, track);
  v.input.move(g.label.x + 2, g.label.y + g.label.height / 2);
  frames.flush();
  const box = geo(v, track)[kind].find((s) => s.name === name)!.box;
  const [x, y] = mid(box);
  v.input.move(x, y);
  frames.flush();
  return v.input.click(x, y);
}
const items = (v: SeqViewer, track?: string): HeaderMenuItem[] =>
  headerMenu(v, { section: Section, ...(track ? { track } : {}) })!.groups.flat();
const item = (v: SeqViewer, track: string | undefined, name: string) => items(v, track).find((i) => i.name === name)!;
const tracks = (v: SeqViewer) => v.view.resolved.value!.sections[0].tracks.map((t) => t.address);
/** A run's patch landed (the history runs it in turn) and its frame drawn. */
const settled = async (r: unknown) => {
  await r;
  await new Promise((done) => setTimeout(done, 0));
  frames.flush();
};

describe('the built-in entries: a slot, a menu item and the command alike', () => {
  it('seq:hide: the slot, the menu item and track.hide each one patch, one undo', async () => {
    const v = await make();
    const root = v.spec.revision.value;
    const ways: (() => unknown)[] = [
      () => clickPart(v, Domains, 'seq:hide'),
      () => headers(v).runItem(item(v, Domains, 'seq:hide')),
      () => v.execute('track.hide', { address: Domains }),
    ];
    let hidden: string | undefined;
    for (const way of ways) {
      const steps = stack(v).steps.value.length;
      await settled(way());
      expect(tracks(v)).toEqual([Letters, More]);
      expect(stack(v).steps.value.length).toBe(steps + 1);
      hidden ??= v.spec.revision.value;
      expect(v.spec.revision.value).toBe(hidden);
      await settled(stack(v).undo());
      expect(v.spec.revision.value).toBe(root);
    }
  });

  it('seq:color: its slot opens the schemes, a scheme item runs layer.channel (one patch), previewed while hovered', async () => {
    const v = await make();
    expect(clickPart(v, Letters, 'seq:color')).toBe(true);
    const shown = headers(v).menu.value!;
    const schemes = shown.menu.groups.flat();
    expect(schemes.length).toBeGreaterThan(1);
    expect(schemes.every((s) => s.command === 'layer.channel')).toBe(true);
    // hovered: previewed, no step; closed without a choice: nothing
    expect(headers(v).preview(schemes[1])).toBe(true);
    expect(v.spec.preview.value).toBeDefined();
    headers(v).closeMenu();
    expect(v.spec.preview.value).toBeUndefined();
    const steps = stack(v).steps.value.length;
    await settled(headers(v).runItem(schemes[1]));
    expect(stack(v).steps.value.length).toBe(steps + 1);
    const color = v.spec.root.value.children!.find((c) => c.kind === 'views')!.children![0].children![0].children![0]
      .children![0].params?.color as { params: { name: string } };
    expect(color.params.name).toBe((schemes[1].args!.value as { params: { name: string } }).params.name);
    // the full menu's "Color ▸" holds the same items, the scheme now checked
    const sub = item(v, Letters, 'seq:color').items!;
    expect(sub.map((s) => s.label)).toEqual(schemes.map((s) => s.label));
    expect(sub[1].checked).toBe(true);
  });

  it('seq:layers: a checked item a layer, hidden and shown by layer.hide and layer.show; the chip shows all (one patch)', async () => {
    const v = await make();
    const layers = item(v, Domains, 'seq:layers').items!;
    expect(layers.map((l) => [l.label, l.checked])).toEqual([
      ['d', true],
      ['e', true],
    ]);
    await settled(headers(v).runItem(layers[1]));
    await settled(v.execute('layer.hide', { address: `${Domains}/d` }));
    expect(item(v, Domains, 'seq:layers').items!.map((l) => l.checked)).toEqual([false, false]);
    // the state chip "2 layers hidden": its click shows both in one patch
    expect(geo(v, Domains).states.map((s) => s.text)).toEqual(['2 layers hidden']);
    const steps = stack(v).steps.value.length;
    expect(clickPart(v, Domains, 'seq:hidden-layers', 'states')).toBe(true);
    await settled(undefined);
    expect(stack(v).steps.value.length).toBe(steps + 1);
    expect(geo(v, Domains).states).toEqual([]);
  });

  it('seq:move-up and seq:move-down: track.move among the siblings, absent at the ends (EC-237)', async () => {
    const v = await make();
    expect(items(v, Letters).map((i) => i.name)).not.toContain('seq:move-up');
    expect(items(v, More).map((i) => i.name)).not.toContain('seq:move-down');
    await settled(headers(v).runItem(item(v, More, 'seq:move-up')));
    expect(tracks(v)).toEqual([Letters, More, Domains]);
    await settled(headers(v).runItem(item(v, More, 'seq:move-up')));
    expect(tracks(v)).toEqual([More, Letters, Domains]);
    await settled(headers(v).runItem(item(v, More, 'seq:move-down')));
    expect(tracks(v)).toEqual([Letters, More, Domains]);
  });

  it('seq:hidden-tracks: the section’s chip shows every hidden track in one patch; its list shows one (EC-237)', async () => {
    const v = await make();
    await settled(v.execute('track.hide', { address: Domains }));
    await settled(v.execute('track.hide', { address: More }));
    expect(geo(v, undefined).states.map((s) => s.text)).toEqual(['2 hidden']);
    // its secondary click lists each, then all
    const chip = geo(v, undefined).states[0].box;
    expect(headers(v).context(chip.x + 2, chip.y + chip.height / 2)).toBe(true);
    expect(
      headers(v)
        .menu.value!.menu.groups.flat()
        .map((i) => i.label),
    ).toEqual(['Show Domains', 'Show More', 'Show all']);
    headers(v).closeMenu();
    const steps = stack(v).steps.value.length;
    expect(clickPart(v, undefined, 'seq:hidden-tracks', 'states')).toBe(true);
    await settled(undefined);
    expect(stack(v).steps.value.length).toBe(steps + 1);
    expect(tracks(v)).toEqual([Letters, Domains, More]);
  });

  it('seq:filtered: the chip clears the filter as layer.filter does', async () => {
    const v = await make();
    await settled(
      v.execute('layer.channel', { address: `${Domains}/d`, value: { kind: 'field', params: { name: 'id' } } }),
    );
    await settled(v.execute('layer.filter', { address: `${Domains}/d`, entry: 'B' }));
    expect(geo(v, Domains).states.map((s) => s.name)).toContain('seq:filtered');
    expect(item(v, Domains, 'seq:filtered').text).toBe('1 filtered');
    expect(clickPart(v, Domains, 'seq:filtered', 'states')).toBe(true);
    await settled(undefined);
    expect(geo(v, Domains).states.map((s) => s.name)).not.toContain('seq:filtered');
  });

  it('seq:reverse and seq:collapse: the section’s menu, section.reverse and section.collapse', async () => {
    const v = await make();
    const section = items(v).map((i) => i.name);
    expect(section).toEqual(expect.arrayContaining(['seq:reverse', 'seq:collapse']));
    await settled(headers(v).runItem(item(v, undefined, 'seq:reverse')));
    expect(v.view.resolved.value!.sections[0].axis.dir).toBe(-1);
  });
});

describe('a host’s entries', () => {
  const host: HeaderEntry = {
    name: 'myapp:pin',
    on: 'track',
    icon: 'pin',
    label: 'Pin',
    command: 'track.move',
    args: (f) => ({ address: f.address, after: null }),
    order: 0,
  };

  it('registered at run time: a slot and a menu item; its undo takes both away', async () => {
    const v = await make();
    const layouts = v.counts.layouts;
    const undo = headers(v).register(host);
    frames.flush();
    expect(v.counts.layouts).toBe(layouts + 1);
    expect(item(v, More, 'myapp:pin').slot).toBe(true);
    expect(geo(v, More).slots.map((s) => s.name)).toEqual(['myapp:pin', 'seq:hide', 'seq:menu']);
    expect(clickPart(v, More, 'myapp:pin')).toBe(true);
    await settled(undefined);
    expect(tracks(v)[0]).toBe(More);
    undo();
    frames.flush();
    expect(items(v, More).map((i) => i.name)).not.toContain('myapp:pin');
    expect(geo(v, More).slots.map((s) => s.name)).toEqual(['seq:hide', 'seq:color', 'seq:menu']);
  });

  it('a conflict is refused whole: a name twice, no namespace, a command the viewer lacks', async () => {
    const v = await make();
    const before = headers(v).entries;
    expect(() => headers(v).register(host, { ...host, label: 'Again' })).toThrow(/myapp:pin: registered twice/);
    expect(() => headers(v).register({ ...host, name: 'pin' })).toThrow(/not namespaced/);
    expect(() => headers(v).register({ ...host, name: 'myapp:x', command: 'no.such' })).toThrow(
      /seq.viewer.header.command/,
    );
    expect(headers(v).entries).toBe(before);
    // in the spec's registry: the viewer is not made
    expect(() =>
      createSeqViewer({
        ...DefaultSeqViewerSpec(),
        ui: undefined,
        registry: [
          ...DefaultSeqViewerSpec().registry!,
          { kind: 'header-action', entry: { ...host, name: 'seq:hide' } },
        ],
      }),
    ).toThrow(/seq:hide: registered twice/);
  });

  it('a built-in left out of the registry is gone; a disabled entry is drawn muted and says why', async () => {
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({
      ...base,
      ui: undefined,
      config: [[SeqCoreConfig.render, { labelWidth: 140 }]],
      registry: base.registry!.filter((e) => !(e.kind === 'header-action' && e.entry.name === 'seq:hide')),
    });
    viewers.push(v);
    await v.load(spec());
    v.resize(600);
    await v.settle();
    frames.flush();
    expect(items(v, More).map((i) => i.name)).not.toContain('seq:hide');
    // block.form has no block here: its own `when` says why
    headers(v).register({
      name: 'myapp:form',
      on: 'track',
      icon: 'box',
      label: 'As blocks',
      command: 'block.form',
      order: 0,
    });
    frames.flush();
    const it = item(v, More, 'myapp:form');
    expect([it.enabled, it.why]).toEqual([false, 'No block or form given']);
    const slot = v.frame!.layout.sections[0].tracks.find((t) => t.track.address === More)!.header!.slots[0];
    expect([slot.part.name, slot.part.enabled]).toEqual(['myapp:form', false]);
  });
});

describe('the menu', () => {
  it('its data: the element’s entries, the track’s, the section’s; each item’s key, enabled or why', async () => {
    const v = await make();
    headers(v).register({
      name: 'myapp:show-in-3d',
      on: 'element',
      label: 'Show in 3D',
      command: 'address.focus',
      args: (f) => ({ address: `${f.element!.layer}` }),
    });
    const m = headerMenu(v, { section: Section, track: Domains, element: { layer: `${Domains}/d`, key: 'A' } })!;
    expect(m.groups.map((g) => g.map((i) => i.name))).toEqual([
      ['myapp:show-in-3d'],
      ['seq:hide', 'seq:color', 'seq:layers', 'seq:move-up', 'seq:move-down'],
      ['seq:reverse', 'seq:collapse'],
    ]);
    expect(item(v, undefined, 'seq:reverse')).toMatchObject({ key: expect.stringMatching(/R/), enabled: true });
  });

  it('a secondary click opens it under each contextClick; Shift+secondary click is never taken', async () => {
    for (const contextClick of ['everywhere', 'header', 'none'] as const) {
      const v = await make({ contextClick });
      const g = geo(v, Domains);
      const data = [v.dataX + 40, g.line.y + g.line.height / 2] as const;
      const label = [g.label.x + 2, g.line.y + g.line.height / 2] as const;
      expect(headers(v).context(...label, { shiftKey: true })).toBe(false);
      expect(headers(v).context(...data), contextClick).toBe(contextClick === 'everywhere');
      headers(v).closeMenu();
      expect(headers(v).context(...label), contextClick).toBe(contextClick !== 'none');
      if (contextClick !== 'none') expect(headers(v).menu.value!.menu.target.track).toBe(Domains);
      headers(v).closeMenu();
    }
  });

  it('the menu key once per ask: the browser’s secondary click after it opens nothing again (EC-156)', async () => {
    const v = await make();
    const shows: unknown[] = [];
    headers(v).menu.listen((m) => m && shows.push(m));
    v.input.key({ key: 'ArrowRight' });
    expect(v.input.key({ key: 'ContextMenu' })).toBe(true);
    expect(shows).toHaveLength(1);
    expect(headers(v).menu.value!.menu.target).toEqual({ section: Section, track: Letters });
    // the platform's own `contextmenu` that follows the key: taken (the browser's menu prevented), no second menu
    expect(headers(v).context(10, 10)).toBe(true);
    expect(shows).toHaveLength(1);
    // the next ask opens again
    headers(v).closeMenu();
    expect(v.input.key({ key: 'F10', shiftKey: true })).toBe(true);
    expect(shows).toHaveLength(2);
  });
});

describe('the cursor on headers', () => {
  it('← at the first position onto the track’s header, ↑ ↓ across headers, → back; Enter opens its menu', async () => {
    const v = await make();
    const key = (k: string) => v.input.key({ key: k });
    key('ArrowRight');
    key('ArrowLeft');
    expect(v.selection.cursor.value).toMatchObject({ part: 'header', track: Letters });
    // its actions drawn while it is there (the emphasis's cursor names the track)
    frames.flush();
    expect(v.emphasis.cursor?.track).toBe(Letters);
    key('ArrowDown');
    expect(v.selection.cursor.value?.track).toBe(Domains);
    key('ArrowDown');
    expect(v.selection.cursor.value?.track).toBe(More);
    key('ArrowUp');
    expect(v.selection.cursor.value?.track).toBe(Domains);
    expect(key('Enter')).toBe(true);
    expect(headers(v).menu.value!.menu.target.track).toBe(Domains);
    headers(v).closeMenu();
    // its words: the label, its states, the key
    expect(v.behaviors.get('hoverCard')).toBeDefined();
    key('ArrowRight');
    expect(v.selection.cursor.value).toEqual({ section: Section, system: '/data/p', position: 0, keyboard: true });
  });

  it('a hide by the keyboard moves the cursor to the next track’s header, else the previous one’s', async () => {
    const v = await make();
    v.selection.setCursor({
      section: Section,
      system: '/data/p',
      position: 0,
      part: 'header',
      track: Domains,
      keyboard: true,
    });
    frames.flush();
    await settled(v.execute('track.hide', { address: Domains }));
    expect(v.selection.cursor.value?.track).toBe(More);
    await settled(v.execute('track.hide', { address: More }));
    expect(v.selection.cursor.value?.track).toBe(Letters);
  });
});

describe('the geometry and the parts', () => {
  it('published only on change: a hover publishes nothing, a layout that moves a header does', async () => {
    const v = await make();
    const sets: unknown[] = [];
    v.headers!.listen((g) => sets.push(g));
    const g = geo(v, Domains);
    for (let i = 0; i < 4; i++) {
      v.input.move(g.label.x + 2 + i, g.label.y + 4);
      frames.flush();
    }
    // the header hovered lights none of the track's items (EC-146: refused by default)
    expect([v.hover.current.target, v.hover.current.element, v.hover.current.extent]).toEqual([
      'header-label',
      undefined,
      undefined,
    ]);
    expect(sets).toEqual([]);
    await settled(v.execute('track.hide', { address: Letters }));
    expect(sets).toHaveLength(1);
  });

  it('headerParts gives each shown header the items headerMenu gives it; actions only where they show', async () => {
    const v = await make();
    const parts = headerParts(v)!;
    expect(parts.headers.map((h) => h.geometry.target.track)).toEqual([Letters, Domains, More]);
    for (const h of parts.headers)
      expect(h.items).toEqual(headerMenu(v, h.geometry.target)!.groups.at(h.geometry.target.track ? 0 : -1));
    expect(parts.headers.map((h) => h.actions)).toEqual([false, false, false]);
    const g = geo(v, Domains);
    v.input.move(g.label.x + 2, g.label.y + 4);
    frames.flush();
    expect(headerParts(v)!.headers.map((h) => h.actions)).toEqual([false, true, false]);
  });

  it('the export draws the states alike whether the parts are drawn or a host’s DOM', async () => {
    const root = spec((b) => void b);
    const a = await make(undefined, root);
    await settled(a.execute('track.hide', { address: More }));
    const b = await make({ parts: () => ({ update() {} }) }, root);
    await settled(b.execute('track.hide', { address: More }));
    // with DOM parts the same boxes are laid out (the parts are placed over them)
    expect(b.frame!.layout.sections[0].tracks[0].header!.slots.map((s) => s.x)).toEqual(
      a.frame!.layout.sections[0].tracks[0].header!.slots.map((s) => s.x),
    );
    expect((await b.renderSvg()).svg).toBe((await a.renderSvg()).svg);
  });

  it('the words: a slot its label and key, a chip its text and what a click does, a header its states', async () => {
    const v = await make();
    await settled(v.execute('track.hide', { address: More }));
    const card = v.behaviors.get('hoverCard') as { shown: { value: unknown } };
    expect(card).toBeDefined();
    const d = (await import('../describe.js')).describerFor({ ...v, layout: () => v.frame?.layout });
    expect(
      d.text({
        section: Section,
        layer: Domains,
        system: undefined,
        position: Number.NaN,
        target: 'header-slot',
        state: 'seq:hide',
      }),
    ).toBe('Hide track');
    expect(
      d.text({
        section: Section,
        layer: Section,
        system: undefined,
        position: Number.NaN,
        target: 'header-state',
        state: 'seq:hidden-tracks',
      }),
    ).toBe('1 hidden: click show all');
    expect(
      d.text({ cursor: { section: Section, system: '/data/p', position: 0, part: 'header', track: Domains } }),
    ).toBe('Domains, Enter for actions');
  });
});

describe('the hover path (the quick representative; the allocation case is the storybook tier’s)', () => {
  it('a pointer down the label column of 200 tracks lays nothing out and publishes nothing', async () => {
    const b = createSeqSpec({ title: '200 tracks' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
    b.views.view('main', (v) => {
      const s = v.section('protein', { of: p });
      for (let i = 0; i < 200; i++) s.track(`t${i}`, { label: `Track ${i}` }).layer('l', { type: 'letters' });
    });
    const v = await make(undefined, b.build());
    const layouts = v.counts.layouts;
    const sets: unknown[] = [];
    v.headers!.listen((g) => sets.push(g));
    const lines = v.headers!.value.filter((g) => g.target.track !== undefined);
    expect(lines).toHaveLength(200);
    for (const g of lines) {
      v.input.move(g.label.x + 2, g.line.y + g.line.height / 2);
      v.draw();
    }
    expect(v.counts.layouts).toBe(layouts);
    expect(sets).toEqual([]);
  });
});
