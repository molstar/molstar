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
  port.show(
    menu,
    { x: 0, y: 0 },
    {
      preview() {},
      cancelPreview() {},
      closeMenu() {},
      submenuRestMs: 0,
      runItem: async (it: HeaderMenuItem) => {
        if (it.command) await v.execute(it.command, it.args);
      },
    } as never,
  );
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
    expect(rows.map((r) => r.getAttribute('role'))).toEqual(header.map((r) => r.getAttribute('role')).slice(0, rows.length));
  });
});
