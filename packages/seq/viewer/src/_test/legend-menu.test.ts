// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The legend's layers menu is the track header's Layers menu (_spec/README.md "Legends", A layer's eye; "Behaviors",
// Headers): the same component and rows, in the part's own element (inside the host's UI root, whose reset and tokens
// style the rows) and never on `document.body`; both menus end with a separator, "Show all" and "Hide all", each one
// patch, one undo step, disabled with its reason when it would change nothing.

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { type HeaderMenuItem, headerMenu } from '../header-entries.js';
import { contextMenuPort } from '../header-menu-dom.js';
import { loadLegendEyes, mountLegend } from '../mount-legend.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { stubCanvasContext } from './canvas-stub.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/main/p';
const Pair = `${Section}/pair`;

/** One track of three layers over one data and field (one legend over all three), and a track of one. */
function spec() {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
  const d = b.data.ranges('d', {
    target: p,
    rows: {
      start: new Int32Array([0, 6, 14]),
      end: new Int32Array([5, 12, 20]),
      id: ['A', 'B', 'C'],
      family: ['x', 'y', 'x'],
    },
  });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    const t = s.track('pair', { label: 'Pair' });
    for (const id of ['a', 'b', 'c']) t.layer(id, { type: 'ranges', data: d, color: theme.field('family') });
    s.track('solo', { label: 'Solo' }).layer('d', { type: 'ranges', data: d });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
const disposers: (() => void)[] = [];
afterEach(() => {
  for (const d of disposers.splice(0)) d();
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
  document.body.replaceChildren();
});

async function mounted() {
  stubCanvasContext();
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, history: uiHistory() });
  viewer = v;
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  await loadLegendEyes();
  // the host's UI root: where the header's menu and the legend's both go
  const el = document.createElement('div');
  el.className = 'msui-root';
  document.body.append(el);
  disposers.push(mountLegend(el, v));
  await loadLegendEyes();
  return { v, el };
}

/** The header's Layers menu opened as the default adapter opens it (its flyout), in `container`. */
async function headerLayers(v: SeqViewer, container: HTMLElement) {
  const port = contextMenuPort()(container);
  const menu = headerMenu(v, { section: Section, track: Pair })!;
  port.show(menu, { x: 0, y: 0 }, {
    preview() {},
    cancelPreview() {},
    closeMenu() {},
    submenuRestMs: 0,
    runItem: async (it: HeaderMenuItem) => {
      if (it.command) await v.execute(it.command, it.args);
    },
  } as never);
  disposers.push(() => port.hide());
  const layers = [...container.querySelectorAll<HTMLElement>('[role="menuitem"]')].find((r) =>
    r.textContent?.startsWith('Layers'),
  )!;
  layers.click();
  return [...container.querySelectorAll<HTMLElement>('[role="menu"]')][1];
}

async function legendMenu(el: HTMLElement) {
  el.querySelector<HTMLButtonElement>('[data-legend-eye] button')!.click();
  return vi.waitFor(() => {
    const m = document.querySelector<HTMLElement>('[role="menu"]');
    if (!m) throw new Error('no menu yet');
    return m;
  });
}

describe('the legend’s layers menu looks as the header’s', () => {
  it('opens in the part’s element (inside the UI root), never on document.body', async () => {
    const { el } = await mounted();
    const menu = await legendMenu(el);
    expect(menu.parentElement).toBe(el);
    expect(menu.closest('.msui-root')).toBe(el);
    expect(menu.parentElement).not.toBe(document.body);
  });

  it('its menu and its rows carry the classes of the header’s Layers menu, row for row', async () => {
    const { v, el } = await mounted();
    const flyout = await headerLayers(v, el);
    const header = [...flyout.querySelectorAll<HTMLElement>('[role^="menuitem"]')];
    expect(header.length).toBeGreaterThan(0);
    // closed again: one menu at a time
    document.body.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    for (const m of [...el.querySelectorAll('[role="menu"]')]) m.remove();
    const menu = await legendMenu(el);
    expect(menu.className).toBe(flyout.className);
    const rows = [...menu.querySelectorAll<HTMLElement>('[role^="menuitem"]')];
    expect(rows[0].className).toBe(header[0].className);
    expect(rows.map((r) => r.getAttribute('role'))).toEqual(
      header.map((r) => r.getAttribute('role')).slice(0, rows.length),
    );
  });
});

const stack = (v: SeqViewer) => (v.spec.history as UiHistory).stack;
const hiddenOf = (v: SeqViewer) =>
  ['a', 'b', 'c'].map((id) => {
    let n = v.spec.root.value;
    for (const k of `${Pair}/${id}`.split('/').filter(Boolean)) n = n.children!.find((c) => c.id === k)!;
    return n.params?.hidden === true;
  });
const rowsOf = (menu: HTMLElement) => [...menu.querySelectorAll<HTMLElement>('[role^="menuitem"]')];
const labelsOf = (menu: HTMLElement) => rowsOf(menu).map((r) => r.textContent?.trim());
const row = (menu: HTMLElement, label: string) => rowsOf(menu).find((r) => r.textContent?.trim() === label)!;
const settled = async (v: SeqViewer, r?: unknown) => {
  await r;
  await v.settle();
  frames!.flush();
};
/** A menu closed: the next opens alone. */
const closeAll = () => {
  for (const m of document.querySelectorAll('[role="menu"]')) m.remove();
};

/** Both menus, the same tests: where each opens and how it is read. */
const menus = {
  legend: (el: HTMLElement) => legendMenu(el),
  header: (el: HTMLElement, v: SeqViewer) => headerLayers(v, el),
} as const;

describe.each(Object.keys(menus) as (keyof typeof menus)[])('the %s’s layers menu: Show all and Hide all', (which) => {
  const open = async (el: HTMLElement, v: SeqViewer) => {
    closeAll();
    return menus[which](el, v);
  };

  it('ends with a separator and the two entries, "Show all" disabled with its reason while every layer is shown', async () => {
    const { v, el } = await mounted();
    const menu = await open(el, v);
    expect(labelsOf(menu)).toEqual(['a', 'b', 'c', 'Show all', 'Hide all']);
    expect(rowsOf(menu).map((r) => r.hasAttribute('data-separated'))).toEqual([false, false, false, true, false]);
    const show = row(menu, 'Show all');
    expect(show.getAttribute('aria-disabled')).toBe('true');
    expect(show.getAttribute('aria-description')).toBe('Every layer is shown');
    expect(row(menu, 'Hide all').getAttribute('aria-disabled')).toBeNull();
    // a disabled one runs nothing
    show.click();
    await settled(v);
    expect(hiddenOf(v)).toEqual([false, false, false]);
  });

  it('Hide all hides every layer in one patch, one undo step; the row stays, dimmed; "Hide all" then says why not', async () => {
    const { v, el } = await mounted();
    const steps = stack(v).steps.value.length;
    const menu = await open(el, v);
    row(menu, 'Hide all').click();
    await settled(v);
    expect(hiddenOf(v)).toEqual([true, true, true]);
    expect(stack(v).steps.value.length).toBe(steps + 1);
    // the rule a single hide follows: the track's row stays in the legend, dimmed
    const group = [...el.querySelectorAll<HTMLElement>('[role="group"]')].find(
      (g) => g.getAttribute('aria-label') === 'Pair',
    )!;
    expect(group.dataset.hidden).toBe('true');
    const again = await open(el, v);
    expect(row(again, 'Hide all').getAttribute('aria-description')).toBe('Every layer is hidden');
    expect(row(again, 'Show all').getAttribute('aria-disabled')).toBeNull();
    await settled(v, stack(v).undo());
    expect(hiddenOf(v)).toEqual([false, false, false]);
    // undone: the step is on the redo side
    expect(stack(v).steps.value.length).toBe(steps);
  });

  it('Show all shows the hidden ones alone, one patch, one undo step; Hide all leaves the hidden as they are', async () => {
    const { v, el } = await mounted();
    await settled(v, v.execute('layer.hide', { address: `${Pair}/b` }));
    const steps = stack(v).steps.value.length;
    const menu = await open(el, v);
    expect(row(menu, 'Show all').getAttribute('aria-disabled')).toBeNull();
    row(menu, 'Show all').click();
    await settled(v);
    expect(hiddenOf(v)).toEqual([false, false, false]);
    expect(stack(v).steps.value.length).toBe(steps + 1);
    await settled(v, stack(v).undo());
    expect(hiddenOf(v)).toEqual([false, true, false]);
    // Hide all from the same state: a and c, b as it was; one step
    const next = await open(el, v);
    row(next, 'Hide all').click();
    await settled(v);
    expect(hiddenOf(v)).toEqual([true, true, true]);
    expect(stack(v).steps.value.length).toBe(steps + 1);
    await settled(v, stack(v).undo());
    expect(hiddenOf(v)).toEqual([false, true, false]);
  });

  it('the arrows reach every entry, the disabled one too (its reason read); Escape closes', async () => {
    const { v, el } = await mounted();
    const menu = await open(el, v);
    const rows = rowsOf(menu);
    const key = (k: string) =>
      (document.activeElement as HTMLElement).dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
    rows[0].focus();
    const reached = [document.activeElement];
    for (let i = 1; i < rows.length; i++) {
      key('ArrowDown');
      reached.push(document.activeElement);
    }
    expect(reached).toEqual(rows);
    expect(row(menu, 'Show all').getAttribute('aria-description')).toBe('Every layer is shown');
    // the header's flyout closes alone on the first Escape (its row focused), the menu on the next
    for (let i = 0; i < 2 && document.querySelector('[role="menu"]'); i++) key('Escape');
    expect(document.querySelectorAll('[role="menu"]')).toHaveLength(0);
  });
});
