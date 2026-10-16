/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

// The DOM types, for this file under the headless test config (`tsconfig.test.json` has no DOM lib).
/// <reference lib="dom" />

// The Annotations part (seqspec:interaction#the-default-ui) and what it is a client of: `bookmark.add` with a `color`
// (a `color` column on the row, shown by the ranges layer's `field_color`), `bookmark.remove`, one undo step each; the
// DOM part and its React twin over a section whose bookmarks set is declared with a track over it.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { SeqSet } from '@molstar/seq-core/model';
import { act, createElement } from 'react';
import { createRoot } from 'react-dom/client';
import { afterEach, describe, expect, it, vi } from 'vitest';
import * as D from '../src/dom/index.js';
import { createSeqViewer, manualScheduler, type SeqViewerContext } from '../src/index.js';
import * as R from '../src/react/index.js';
import { renderSvg } from '../src/svg/index.js';
import { fixtureMeasurer } from '../src/testing/index.js';
import { headlessSpec } from './dom/fixtures.js';
import { seqLoci } from './overlay-docs.js';

(globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }).IS_REACT_ACT_ENVIRONMENT = true;
// jsdom has no `CSS.escape`, which the parts' focus keeping reads.
(globalThis as { CSS?: unknown }).CSS ??= { escape: (s: string) => s.replace(/[^\w-]/g, (c) => `\\${c}`) };

const SECTION = 's';
const SET = `${SECTION}.bookmarks`;
const ABOVE = `${SECTION}.bookmarks.above`;
const STYLE = { theme: 'field_line_style', params: { field: 'line_style' } };
const PLACES = [
  { data: ABOVE, label: 'Above' },
  { data: SET, label: 'Below' },
] as const;
const PROTEIN = 'MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDD';
const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};

/** A section with its bookmarks set declared empty and a ranges track over it, colored by the rows' `color`. */
function doc(declared = true, above = false): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Annotations' });
  const a = b.data.sequence({ ref: 'a', value: PROTEIN });
  const notes = declared
    ? b.data.annotation({ ref: SET, label: 'Annotations', geometry: 'ranges', value: [] } as never)
    : null;
  const up = above
    ? b.data.annotation({ ref: ABOVE, label: 'Annotations above', geometry: 'ranges', value: [] } as never)
    : null;
  const sec = b.view({ ref: 'v', label: 'V' }).section({ ref: SECTION });
  sec.block({ sequence: a });
  if (up) {
    sec
      .track({ ref: 'notes-above', label: 'Annotations (above)' })
      .ranges({ data: up, labels: true, color: { theme: 'field_color' } } as never);
  }
  sec.track({ ref: 'letters', label: 'Sequence' }).letters();
  if (notes) {
    sec
      .track({ ref: 'notes', label: 'Annotations' })
      .ranges({ data: notes, labels: true, color: { theme: 'field_color' }, line_style: STYLE } as never);
  }
  return b.getDocument();
}

const cleanups: (() => void | Promise<void>)[] = [];
afterEach(async () => {
  for (const c of cleanups.splice(0).reverse()) await c();
});

async function open(d = doc()) {
  const ctx = createSeqViewer(headlessSpec(), { scheduler: manualScheduler(), text: fixtureMeasurer() });
  cleanups.push(() => ctx.dispose());
  const done = ctx.setDocument(d);
  ctx.viewports.get(SECTION)?.setSize(600, 200, 1);
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  await ctx.graph.update();
  const set = () => ctx.graph.assembled.value?.set as SeqSet;
  const select = (s: number, e: number) => ctx.selection.set(seqLoci(set(), 'a', s, e), 'api');
  const rows = (ref = SET) => {
    const node = (ctx.document.shown.value?.root.children ?? [])
      .find((c) => c.kind === 'data')
      ?.children?.find((c) => c.ref === ref);
    return node?.params?.value as Record<string, unknown>[] | undefined;
  };
  return { ctx, set, select, rows };
}

describe('bookmark.add with a color', () => {
  it('writes the `color` column; without one, no column (an old set is as it was)', async () => {
    const { ctx, set, rows } = await open();
    const at = (s: number, e: number) => seqLoci(set(), 'a', s, e);
    expect(ctx.actions.can('bookmark.add', { section: SECTION, loci: at(2, 6), color: 7 })).toBe(false);
    expect(ctx.actions.run('bookmark.add', { section: SECTION, loci: at(2, 6), label: 'Loop', color: '#e76f51' })).toBe(
      true,
    );
    await flush();
    expect(rows()).toEqual([{ target: 'a/unit', start: 2, stop: 6, id: 'b1', label: 'Loop', color: '#e76f51' }]);
    ctx.actions.run('bookmark.add', { section: SECTION, loci: at(10, 12) });
    await flush();
    expect(rows()?.[1]).toEqual({ target: 'a/unit', start: 10, stop: 12, id: 'b2' });
    // The set is still a valid document: no issue came of rows with and without the column.
    expect(ctx.issues.current.value.filter((i) => i.severity === 'error' || i.severity === 'warning')).toEqual([]);
  });

  it('is one patch and one undo step each; a declared set takes no `highlights` entry', async () => {
    const { ctx, set, rows } = await open();
    const at = (s: number, e: number) => seqLoci(set(), 'a', s, e);
    ctx.actions.run('bookmark.add', { section: SECTION, loci: at(2, 6), color: '#2a9d8f' });
    await flush();
    ctx.actions.run('bookmark.add', { section: SECTION, loci: at(10, 12), color: '#264653' });
    await flush();
    expect(rows()?.map((r) => r.id)).toEqual(['b1', 'b2']);
    expect(ctx.highlight.channels.value.size).toBe(0);
    ctx.actions.run('bookmark.remove', { section: SECTION, item: 'b1' });
    await flush();
    expect(rows()?.map((r) => [r.id, r.color])).toEqual([['b2', '#264653']]);
    // Undo: the remove, then the second add, then the first.
    expect(await ctx.history.undo()).toBe(true);
    expect(rows()?.map((r) => r.id)).toEqual(['b1', 'b2']);
    expect(await ctx.history.undo()).toBe(true);
    expect(rows()?.map((r) => r.id)).toEqual(['b1']);
    expect(await ctx.history.undo()).toBe(true);
    expect(rows()).toEqual([]);
    expect(await ctx.history.undo()).toBe(false);
  });

  it('without a declared set it makes one and its highlights entry, as before, and keeps the color', async () => {
    const { ctx, set, rows } = await open(doc(false));
    ctx.actions.run('bookmark.add', { section: SECTION, loci: seqLoci(set(), 'a', 2, 6), color: '#e76f51' });
    await flush();
    expect(rows()?.map((r) => r.color)).toEqual(['#e76f51']);
    expect([...ctx.highlight.channels.value.keys()].filter((c) => c.startsWith('set:'))).toEqual([`set:${SECTION}/0`]);
  });

  it("is drawn in the row's color by the track's `field_color`", async () => {
    const { ctx, set } = await open();
    ctx.actions.run('bookmark.add', { section: SECTION, loci: seqLoci(set(), 'a', 2, 6), color: '#e76f51' });
    await flush();
    await ctx.graph.update();
    const svg = await renderSvg(ctx, { width: 600 });
    expect(svg.toLowerCase()).toMatch(/e76f51|rgb\(231,\s*111,\s*81\)/);
  });
});

/** The part's text and structure (tag, class, the attributes compared), for the two UIs to be compared. */
function shape(n: Element, out: string[] = []): string[] {
  const attrs = [...n.attributes].filter((a) => /^(class|type|aria-.+|data-seqv-.+)$/.test(a.name));
  out.push(
    `${n.tagName.toLowerCase()} ${attrs
      .map((a) => `${a.name}=${a.value}`)
      .sort()
      .join(' ')} ${[...n.childNodes]
      .filter((c) => c.nodeType === 3)
      .map((c) => c.textContent)
      .join('')}`,
  );
  for (const c of n.children) shape(c, out);
  return out;
}

async function domPart(
  ctx: SeqViewerContext,
  places?: readonly { data: string; label: string; kind?: 'ranges' | 'points' | 'frame' }[],
) {
  const host = document.createElement('div');
  document.body.appendChild(host);
  const part = D.AnnotationsPanel(host, ctx, { section: SECTION, places });
  cleanups.push(() => {
    part.dispose();
    host.remove();
  });
  return host;
}
async function reactPart(
  ctx: SeqViewerContext,
  places?: readonly { data: string; label: string; kind?: 'ranges' | 'points' | 'frame' }[],
) {
  const host = document.createElement('div');
  document.body.appendChild(host);
  const root = createRoot(host);
  await act(async () => root.render(createElement(R.AnnotationsPanel, { context: ctx, section: SECTION, places })));
  cleanups.push(async () => {
    await act(async () => root.unmount());
    host.remove();
  });
  return host;
}
const q = (host: HTMLElement, key: string) =>
  host.querySelector<HTMLElement>(`[data-seqv-key="${key}"]`) as HTMLElement;

describe('the Annotations part (DOM)', () => {
  it('adds the selection with its label and color, lists, reveals and deletes it; Undo brings it back', async () => {
    const { ctx, select, rows } = await open();
    const host = await domPart(ctx);
    const root = host.querySelector('.seqv-annotations') as HTMLElement;
    expect(root.getAttribute('data-seqv-part')).toBe('annotations');
    expect(host.textContent).toContain('No annotations');
    // Nothing selected: Add cannot act (focusable, not disabled).
    expect(q(host, 'add').getAttribute('aria-disabled')).toBe('true');
    q(host, 'add').click();
    await flush();
    expect(rows()).toEqual([]);

    select(4, 9);
    expect(q(host, 'add').getAttribute('aria-disabled')).toBe('false');
    (q(host, 'label') as HTMLInputElement).value = 'Helix';
    (q(host, 'color') as HTMLInputElement).value = '#2a9d8f';
    q(host, 'add').click();
    await flush();
    expect(rows()).toEqual([{ target: 'a/unit', start: 4, stop: 9, id: 'b1', label: 'Helix', color: '#2a9d8f' }]);
    expect((q(host, 'label') as HTMLInputElement).value).toBe('');
    expect(host.querySelector('.seqv-annotation-label')?.textContent).toBe('Helix');
    expect(host.textContent).toContain('5–9');
    expect(host.textContent).toContain('1 annotation');
    expect((host.querySelector('.seqv-swatch') as HTMLElement).style.background).toBe('rgb(42, 157, 143)');

    // A click on the row reveals it, changing nothing else.
    const reveal = vi.spyOn(ctx, 'reveal');
    q(host, 'go:b1').click();
    expect(reveal).toHaveBeenCalledTimes(1);

    q(host, 'del:b1').click();
    await flush();
    expect(rows()).toEqual([]);
    expect(host.querySelector('.seqv-annotation')).toBeNull();
    expect(await ctx.history.undo()).toBe(true);
    expect(rows()?.map((r) => r.id)).toEqual(['b1']);
    expect(host.querySelector('.seqv-annotation-label')?.textContent).toBe('Helix');
  });

  it('Enter in the label field adds', async () => {
    const { ctx, select, rows } = await open();
    const host = await domPart(ctx);
    select(0, 3);
    const label = q(host, 'label') as HTMLInputElement;
    label.value = 'N-terminus';
    label.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
    await flush();
    expect(rows()?.map((r) => r.label)).toEqual(['N-terminus']);
  });

  it('stands down on dispose: its nodes and subscriptions go', async () => {
    const { ctx } = await open();
    const host = document.createElement('div');
    document.body.appendChild(host);
    const part = D.AnnotationsPanel(host, ctx, { section: SECTION });
    expect(host.querySelector('.seqv-annotations')).not.toBeNull();
    part.dispose();
    expect(host.querySelector('.seqv-annotations')).toBeNull();
    host.remove();
  });
});

describe('the Annotations part (React)', () => {
  it("has the DOM part's structure and text, with and without annotations", async () => {
    const { ctx, select } = await open();
    const [dom, react] = [await domPart(ctx), await reactPart(ctx)];
    expect(shape(react.firstElementChild as Element)).toEqual(shape(dom.querySelector('.seqv-annotations') as Element));
    select(4, 9);
    (q(dom, 'label') as HTMLInputElement).value = 'Helix';
    q(dom, 'add').click();
    await act(async () => flush());
    expect(shape(react.querySelector('.seqv-annotations') as Element)).toEqual(
      shape(dom.querySelector('.seqv-annotations') as Element),
    );
    expect(react.querySelector('.seqv-annotation-label')?.textContent).toBe('Helix');
  });

  it('adds, reveals and deletes through the same actions', async () => {
    const { ctx, select, rows } = await open();
    const host = await reactPart(ctx);
    select(4, 9);
    await act(async () => flush());
    const label = q(host, 'label') as HTMLInputElement;
    await act(async () => {
      // React tracks the value it set: the native setter, then the event it listens for.
      Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')?.set?.call(label, 'Helix');
      label.dispatchEvent(new Event('input', { bubbles: true }));
    });
    await act(async () => {
      q(host, 'add').click();
      await flush();
    });
    expect(rows()?.map((r) => r.label)).toEqual(['Helix']);
    const reveal = vi.spyOn(ctx, 'reveal');
    await act(async () => q(host, 'go:b1').click());
    expect(reveal).toHaveBeenCalledTimes(1);
    await act(async () => {
      q(host, 'del:b1').click();
      await flush();
    });
    expect(rows()).toEqual([]);
    expect(host.querySelector('.seqv-annotation')).toBeNull();
  });
});

describe('the Annotations part with several sets (Above and Below)', () => {
  const pressed = (host: HTMLElement) =>
    [...host.querySelectorAll('.seqv-annotations-place button')].map((b) => [
      b.textContent,
      b.getAttribute('aria-pressed'),
    ]);

  it('offers the choice only with several places; the first is chosen to start with', async () => {
    const { ctx } = await open(doc(true, true));
    const [one, many] = [await domPart(ctx), await domPart(ctx, PLACES)];
    expect(one.querySelector('.seqv-annotations-place')).toBeNull();
    expect(pressed(many)).toEqual([
      ['Above', 'true'],
      ['Below', 'false'],
    ]);
  });

  it('the choice picks the set `bookmark.add` writes to, one patch each; the list names both sets', async () => {
    const { ctx, select, rows } = await open(doc(true, true));
    const host = await domPart(ctx, PLACES);
    select(4, 9);
    (q(host, 'label') as HTMLInputElement).value = 'Helix';
    q(host, 'add').click();
    await flush();
    expect(rows(ABOVE)?.map((r) => r.label)).toEqual(['Helix']);
    expect(rows()).toEqual([]);

    q(host, 'place:1').click();
    expect(pressed(host)).toEqual([
      ['Above', 'false'],
      ['Below', 'true'],
    ]);
    select(20, 25);
    (q(host, 'label') as HTMLInputElement).value = 'Loop';
    q(host, 'add').click();
    await flush();
    expect(rows()?.map((r) => r.label)).toEqual(['Loop']);
    expect(rows(ABOVE)?.map((r) => r.label)).toEqual(['Helix']);
    // Both ids are `b1` (one in each set): the rows are told apart by their set.
    const labels = [...host.querySelectorAll('.seqv-annotation')].map((li) => [
      li.querySelector('.seqv-annotation-label')?.textContent,
      li.querySelector('.seqv-annotation-place')?.textContent,
    ]);
    expect(labels).toEqual([
      ['Helix', 'Above'],
      ['Loop', 'Below'],
    ]);
    expect(host.textContent).toContain('2 annotations');

    // A delete acts on its row's set; one undo step each (the two adds, then the delete).
    q(host, `del:${ABOVE}:b1`).click();
    await flush();
    expect(rows(ABOVE)).toEqual([]);
    expect(rows()?.map((r) => r.label)).toEqual(['Loop']);
    expect(await ctx.history.undo()).toBe(true);
    expect(rows(ABOVE)?.map((r) => r.label)).toEqual(['Helix']);
    expect(await ctx.history.undo()).toBe(true);
    expect(rows()).toEqual([]);
    expect(rows(ABOVE)?.map((r) => r.label)).toEqual(['Helix']);
    expect(await ctx.history.undo()).toBe(true);
    expect(rows(ABOVE)).toEqual([]);
    expect(await ctx.history.undo()).toBe(false);
  });

  it('a row reveals its own set; each row has its own accessible names', async () => {
    const { ctx, set, rows } = await open(doc(true, true));
    ctx.actions.run('bookmark.add', { section: SECTION, loci: seqLoci(set(), 'a', 4, 9), label: 'Up', data: ABOVE });
    ctx.actions.run('bookmark.add', { section: SECTION, loci: seqLoci(set(), 'a', 4, 9), label: 'Up' });
    await flush();
    expect([rows(ABOVE)?.length, rows()?.length]).toEqual([1, 1]);
    const host = await domPart(ctx, PLACES);
    expect(q(host, `go:${ABOVE}:b1`).getAttribute('aria-label')).toBe('Show Up (5–9, Above)');
    expect(q(host, `del:${SET}:b1`).getAttribute('aria-label')).toBe('Delete Up (Below)');
    const toLoci = vi.spyOn(ctx.refs, 'toLoci');
    q(host, `go:${SET}:b1`).click();
    expect(toLoci).toHaveBeenCalledWith(SET, { items: ['b1'] });
  });

  it("React has the DOM part's structure and does the same", async () => {
    const { ctx, select, rows } = await open(doc(true, true));
    const [dom, react] = [await domPart(ctx, PLACES), await reactPart(ctx, PLACES)];
    expect(shape(react.firstElementChild as Element)).toEqual(shape(dom.querySelector('.seqv-annotations') as Element));
    select(4, 9);
    await act(async () => flush());
    await act(async () => q(react, 'place:1').click());
    expect(pressed(react)).toEqual([
      ['Above', 'false'],
      ['Below', 'true'],
    ]);
    await act(async () => {
      q(react, 'add').click();
      await flush();
    });
    expect(rows()?.map((r) => r.id)).toEqual(['b1']);
    expect(rows(ABOVE)).toEqual([]);
    await act(async () => q(react, `place:0`).click());
    await act(async () => {
      q(react, 'add').click();
      await flush();
    });
    expect(rows(ABOVE)?.map((r) => r.id)).toEqual(['b1']);
    expect(react.querySelectorAll('.seqv-annotation')).toHaveLength(2);
    await act(async () => {
      q(react, `del:${SET}:b1`).click();
      await flush();
    });
    expect(rows()).toEqual([]);
    expect(rows(ABOVE)?.map((r) => r.id)).toEqual(['b1']);
  });
});

describe('the Annotations part: the line style', () => {
  const setStyle = (host: HTMLElement, v: string) => {
    const select = q(host, 'style') as HTMLSelectElement;
    select.value = v;
    select.dispatchEvent(new Event('change', { bubbles: true }));
  };

  it("DOM: a picker beside the color, solid to start; the pick is the row's `line_style`, listed; solid writes no column", async () => {
    const { ctx, select, rows } = await open();
    const host = await domPart(ctx);
    const picker = q(host, 'style') as HTMLSelectElement;
    expect(picker.value).toBe('solid');
    expect([...picker.options].map((o) => o.value)).toEqual([
      'solid',
      'dashed',
      'dotted',
      'wavy',
      'double',
      'box',
      'highlight',
    ]);
    select(4, 9);
    setStyle(host, 'wavy');
    q(host, 'add').click();
    await flush();
    expect(rows()).toEqual([{ target: 'a/unit', start: 4, stop: 9, id: 'b1', color: '#e76f51', line_style: 'wavy' }]);
    expect(host.querySelector('.seqv-annotation-style')?.textContent).toBe('wavy');
    // The choice stays for the next add; solid again writes none, and lists none.
    setStyle(host, 'solid');
    select(12, 15);
    q(host, 'add').click();
    await flush();
    expect(rows()?.[1]).toEqual({ target: 'a/unit', start: 12, stop: 15, id: 'b2', color: '#e76f51' });
    expect(host.querySelectorAll('.seqv-annotation-style')).toHaveLength(1);
    expect(ctx.issues.current.value.filter((i) => i.severity === 'error' || i.severity === 'warning')).toEqual([]);
  });

  it("React: the DOM part's structure, and the pick is written the same", async () => {
    const { ctx, select, rows } = await open();
    const [dom, react] = [await domPart(ctx), await reactPart(ctx)];
    expect(shape(react.firstElementChild as Element)).toEqual(shape(dom.querySelector('.seqv-annotations') as Element));
    select(4, 9);
    await act(async () => flush());
    const picker = q(react, 'style') as HTMLSelectElement;
    await act(async () => {
      Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')?.set?.call(picker, 'highlight');
      picker.dispatchEvent(new Event('change', { bubbles: true }));
    });
    await act(async () => {
      q(react, 'add').click();
      await flush();
    });
    expect(rows()?.map((r) => r.line_style)).toEqual(['highlight']);
    expect(react.querySelector('.seqv-annotation-style')?.textContent).toBe('highlight');
  });
});

describe('the Annotations part: markers (a place of `kind` points)', () => {
  /** Below's set as markers: a glyph each, no line style. */
  const MARKERS = [
    { data: ABOVE, label: 'Above' },
    { data: SET, label: 'Markers', kind: 'points' as const },
  ];
  it("DOM and React: the glyph's picker in the line style's place; the pick is the row's `glyph`, listed", async () => {
    const { ctx, select, rows } = await open(doc(true, true));
    const host = await domPart(ctx, MARKERS);
    expect(q(host, 'glyph')).toBeNull();
    q(host, 'place:1').click();
    expect(q(host, 'style')).toBeNull();
    const picker = q(host, 'glyph') as HTMLSelectElement;
    expect([...picker.options].map((o) => o.value)).toEqual(['circle', 'triangle', 'star', 'diamond', 'square']);
    picker.value = 'star';
    picker.dispatchEvent(new Event('change', { bubbles: true }));
    select(4, 5);
    q(host, 'add').click();
    await flush();
    expect(rows()).toEqual([{ target: 'a/unit', start: 4, stop: 5, id: 'b1', color: '#e76f51', glyph: 'star' }]);
    expect(host.querySelector('.seqv-annotation-style')?.textContent).toBe('star');
    // Back to a place of ranges: the line style's picker again.
    q(host, 'place:0').click();
    expect(q(host, 'glyph')).toBeNull();
    expect(q(host, 'style')).not.toBeNull();
    const react = await reactPart(ctx, MARKERS);
    await act(async () => q(react, 'place:1').click());
    expect(q(react, 'style')).toBeNull();
    expect((q(react, 'glyph') as HTMLSelectElement).value).toBe('circle');
    select(8, 9);
    await act(async () => flush());
    await act(async () => {
      q(react, 'add').click();
      await flush();
    });
    expect(rows()?.[1]).toEqual({ target: 'a/unit', start: 8, stop: 9, id: 'b2', color: '#e76f51', glyph: 'circle' });
  });
});
