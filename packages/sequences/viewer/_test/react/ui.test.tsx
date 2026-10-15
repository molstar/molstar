/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The React UI (UX3): component tests with a real context and a canvas stub. Parts read only their states (render
// counts on unrelated changes stay zero), slots replace parts, no `title` attributes, the tooltip moves without
// re-rendering, the default `Description` shows plain text, and keyboard focus is shown for Tab, not a click. Menus and
// focus are menus.test.tsx, the provider and `SeqViewer` viewer.test.tsx, the parts' content parts.test.tsx.

import { act } from 'react';
import { afterEach, describe, expect, it } from 'vitest';
import { userEvent } from 'vitest/browser';
import { SeqViewerConfig } from '../../src/index.js';
import type { RowProps } from '../../src/react/index.js';
import * as Parts from '../../src/react/index.js';
// The skin, for the legend's grid (its layout is CSS).
import '../../src/skin.css';
import { at, counted, keys, lines, mount, openedContext, pointer, settle, stubSurface, xOf } from './fixtures.js';

const cleanups: (() => Promise<void> | void)[] = [];
afterEach(async () => {
  for (const c of cleanups.splice(0)) await c();
});

async function opened(ui: Parameters<typeof mount>[1] = {}, strict = false) {
  const t = await openedContext();
  const m = await mount(t.ctx, { surface: stubSurface, ...ui }, strict);
  cleanups.push(m.unmount, () => t.ctx.dispose());
  return { ...t, ...m };
}

function RowToggle({ displayRow }: RowProps) {
  return (
    <button type="button" aria-label={`Row ${displayRow}`}>
      3D
    </button>
  );
}
const allParts = () =>
  counted({
    ViewTabs: Parts.ViewTabs,
    SectionHeader: Parts.SectionHeader,
    TrackActions: Parts.TrackActions,
    HiddenChips: Parts.HiddenChips,
    RowControl: RowToggle,
    Tooltip: Parts.Tooltip,
    Legend: Parts.Legend,
    SnapshotPlayer: Parts.SnapshotPlayer,
    Description: Parts.Description,
    Export: Parts.Export,
    Issues: Parts.IssueBadge,
    Busy: Parts.Busy,
    EmptyState: Parts.EmptyState,
    ErrorState: Parts.ErrorState,
  });
const changed = (counts: Record<string, number>) => Object.keys(counts).filter((k) => counts[k] > 0);

describe('the React UI (U1)', () => {
  it('parts read only their states: unrelated changes render nothing', async () => {
    const c = allParts();
    const t = await opened({ components: c.components });
    const { ctx, container } = t;
    const move = pointer(ctx);
    const res = at(ctx, 'chains.main', 't.res');
    c.reset();
    // A pan, a selection, a highlight, a pointer moving on the axis without a label change: nothing renders.
    await act(async () => {
      ctx.viewports.get('chains.main')?.panBy(40, 0, 'api');
      ctx.selection.set(ctx.refs.toLoci('a', { start: 3, stop: 9 }), 'api');
      ctx.highlight.set('host', ctx.refs.toLoci('dom'));
    });
    await settle();
    expect(changed(c.counts)).toEqual([]);
    // Hovering positions: the tooltip only (its label changes per position).
    await act(async () => {
      for (let i = 20; i < 30; i++) move('move', 'chains.main', xOf(ctx, 'chains.main', i), res.y);
    });
    expect(changed(c.counts)).toEqual(['Tooltip']);
    c.reset();
    // Busy, an issue, a snapshot: only the parts that show them.
    const busy = ctx.tasks.busy as unknown as { set(v: boolean): void };
    await act(async () => busy.set(true));
    expect(container.querySelector('.seqv-busy')).not.toBeNull();
    await act(async () => busy.set(false));
    expect(changed(c.counts)).toEqual(['Busy']);
    c.reset();
    await act(async () => ctx.report('seqv.test', 'a reported problem'));
    expect(changed(c.counts)).toEqual(['Issues']);
    c.reset();
    await act(async () => ctx.snapshots.apply('s2'));
    await settle();
    expect(changed(c.counts).sort()).toEqual(['SnapshotPlayer']);
    c.reset();
    // Snapshot 1 hides a track and has its own description: the chips (in the header's row), the description and the
    // legend follow.
    await act(async () => ctx.snapshots.apply('s1'));
    await settle();
    expect(changed(c.counts).sort()).toEqual([
      'Description',
      'HiddenChips',
      'Legend',
      'SectionHeader',
      'SnapshotPlayer',
    ]);
    // A view switch: the tabs, and the new view's parts; the player and the description stay.
    c.reset();
    await act(async () => ctx.actions.run('view.switch', 'msa'));
    await settle();
    for (const p of ['SnapshotPlayer', 'Description', 'Busy', 'Issues', 'EmptyState']) expect(c.counts[p]).toBe(0);
    expect(c.counts.ViewTabs).toBeGreaterThan(0);
    // The row controls: one beside each of the alignment's rows in view.
    expect(container.querySelectorAll('.seqv-row-control').length).toBeGreaterThan(0);
  });

  it('slots replace parts, and `none` leaves one out', async () => {
    const Tabs = () => <div data-testid="host-tabs">host tabs</div>;
    const t = await opened({ components: { ViewTabs: Tabs, Legend: 'none', Description: 'none' } });
    expect(t.container.querySelector('[data-testid=host-tabs]')).not.toBeNull();
    expect(t.container.querySelector('.seqv-tabs')).toBeNull();
    expect(t.container.querySelector('.seqv-legend')).toBeNull();
    expect(t.container.querySelector('.seqv-description')).toBeNull();
    expect(t.container.querySelector('.seqv-player')).not.toBeNull();
  });

  it("`clean` chrome leaves the controls out (the canvas parts stay); a slot's `default` brings one back", async () => {
    const t = await openedContext();
    t.ctx.config.set(SeqViewerConfig.UI.Chrome, 'clean');
    t.ctx.report('seqv.test', 'a reported problem');
    const m = await mount(t.ctx, { surface: stubSurface, components: { Legend: 'default' } });
    cleanups.push(m.unmount, () => t.ctx.dispose());
    for (const c of ['.seqv-tabs', '.seqv-section-header', '.seqv-badge'])
      expect(m.container.querySelector(c)).toBeNull();
    expect([m.container.querySelector('.seqv-legend'), m.container.querySelector('.seqv-player')]).not.toContain(null);
    const h = at(t.ctx, 'chains.main', 't.dom');
    await act(async () => pointer(t.ctx)('move', 'chains.main', h.header, h.y));
    expect(m.container.querySelector('.seqv-track-actions')).toBeNull();
  });

  it('no `title` attributes: hints, menus, tooltips, the issue list', async () => {
    const t = await opened({ components: { RowControl: RowToggle } });
    const { ctx, container } = t;
    const move = pointer(ctx);
    await act(async () => ctx.report('seqv.test', 'a reported problem'));
    await act(async () => (container.querySelector('.seqv-badge') as HTMLButtonElement).click());
    const h = at(ctx, 'chains.main', 't.dom');
    await act(async () => move('move', 'chains.main', h.header, h.y));
    const actions = container.querySelector('.seqv-track-actions button') as HTMLButtonElement;
    expect(actions.getAttribute('aria-label')).toBe('Actions for Domains');
    await act(async () => actions.click());
    expect([...container.querySelectorAll('[role=menuitem]')].map((m) => m.textContent)).toEqual([
      'Show only what this track covers',
      'Focus',
      'Colors ▸',
      'Hide track',
    ]);
    await act(async () => move('move', 'chains.main', xOf(ctx, 'chains.main', 12), at(ctx, 'chains.main', 't.res').y));
    expect(container.querySelector('.seqv-tooltip')).not.toBeNull();
    await act(async () => ctx.actions.run('view.switch', 'msa'));
    await settle();
    expect(container.querySelectorAll('button').length).toBeGreaterThan(10);
    expect(container.querySelectorAll('[title]').length).toBe(0);
  });

  it("Shift+F10 on the keyboard cursor's track opens its menu, the first item focused", async () => {
    const { ctx, container } = await opened();
    const key = keys(ctx);
    await act(async () => ctx.input.focus(ctx.layout.sectionIndex('chains.main'), true));
    for (let i = 0; i < 8 && ctx.ui.keyCursor.value?.track !== 't.dom'; i++)
      await act(async () => key('chains.main', 'ArrowDown'));
    expect(ctx.ui.keyCursor.value?.track).toBe('t.dom');
    expect(container.querySelector('[role=menu]')).toBeNull();
    const handled: boolean[] = [];
    await act(async () => void handled.push(key('chains.main', 'F10', true)));
    expect(handled).toEqual([true]);
    expect(document.activeElement?.textContent).toBe('Show only what this track covers');
  });

  it("a track's menu hides the track; its chip shows it again", async () => {
    const t = await opened();
    const { ctx, container } = t;
    const h = at(ctx, 'chains.main', 't.dom');
    await act(async () => pointer(ctx)('move', 'chains.main', h.header, h.y));
    await act(async () => (container.querySelector('.seqv-track-actions button') as HTMLButtonElement).click());
    const hide = [...container.querySelectorAll('[role=menuitem]')].find((m) => m.textContent === 'Hide track');
    await act(async () => (hide as HTMLButtonElement).click());
    await settle();
    const chip = container.querySelector('.seqv-chips button') as HTMLButtonElement;
    expect(chip.textContent).toBe('Show Domains');
    await act(async () => chip.click());
    await settle();
    expect(container.querySelector('.seqv-chips')).toBeNull();
    expect(ctx.history.canUndo.value).toBe(true);
  });

  it('the tooltip moves by its position properties without re-rendering', async () => {
    const c = counted({ Tooltip: Parts.Tooltip });
    const t = await opened({ components: c.components });
    const { ctx, container } = t;
    const move = pointer(ctx);
    const y = at(ctx, 'chains.main', 't.res').y;
    const x = xOf(ctx, 'chains.main', 12);
    await act(async () => move('move', 'chains.main', x, y));
    const tip = container.querySelector('.seqv-tooltip') as HTMLElement;
    expect(tip.textContent).toContain('13');
    const before = tip.style.getPropertyValue('--seqv-x');
    c.reset();
    // Within one cell: the label holds, the tooltip follows.
    await act(async () => {
      for (let d = 1; d <= 3; d++) move('move', 'chains.main', x + d, y + d);
    });
    expect(c.counts.Tooltip).toBe(0);
    expect(tip.style.getPropertyValue('--seqv-x')).not.toBe(before);
    expect(container.querySelector('.seqv-tooltip')).toBe(tip);
  });

  it("the default Description shows a plaintext description as paragraphs of text, a snapshot's over the document's, and a markdown one not at all", async () => {
    const t = await opened({}, true);
    const { ctx, container } = t;
    const paragraphs = () => [...container.querySelectorAll('.seqv-description p')].map((p) => p.textContent);
    expect(paragraphs()).toEqual(['The kinase and its domains, as given.', 'A second paragraph.']);
    await act(async () => ctx.snapshots.apply('s1'));
    await settle();
    expect(paragraphs()).toEqual(['Only the first domain is hidden.']);
    await t.unmount();
    cleanups.splice(cleanups.indexOf(t.unmount), 1);
    const md = await openedContext();
    const doc = md.ctx.document.shown.value;
    await md.ctx.setDocument({
      ...(doc as NonNullable<typeof doc>),
      metadata: { description: '# Title\n\nSome **bold** words.', description_format: 'markdown' },
      snapshots: [],
    });
    const m = await mount(md.ctx, { surface: stubSurface });
    cleanups.push(m.unmount, () => md.ctx.dispose());
    expect(m.container.querySelector('.seqv-description')).toBeNull();
  });

  it('keyboard focus shows on a surface for Tab, never for a click or the keys after it', async () => {
    const t = await opened({ components: { ViewTabs: 'none', SnapshotPlayer: 'none', Description: 'none' } });
    const { container } = t;
    const before = document.createElement('button');
    before.textContent = 'before';
    container.prepend(before);
    const surface = () => container.querySelector('.seqv-surface') as HTMLElement;
    const inner = () => surface().firstElementChild as HTMLElement;
    // Every focusable before the surface, then Tab onto it.
    const focusables = [...container.querySelectorAll<HTMLElement>('button, [tabindex]')];
    focusables[focusables.indexOf(inner()) - 1].focus();
    await userEvent.tab();
    expect(document.activeElement).toBe(inner());
    expect(surface().hasAttribute('data-focus-visible')).toBe(true);
    await userEvent.click(before);
    expect(surface().hasAttribute('data-focus-visible')).toBe(false);
    await userEvent.click(inner());
    expect(document.activeElement).toBe(inner());
    await userEvent.keyboard('{ArrowRight}{+}');
    expect(surface().hasAttribute('data-focus-visible')).toBe(false);
    expect(inner().getAttribute('aria-roledescription')).toBe('sequence view');
    expect(inner().getAttribute('aria-label')).toBe('Chains');
  });
  it("the track's menu sits on the hovered line, a rows track's row on every hover, without re-rendering", async () => {
    const c = counted({ TrackActions: Parts.TrackActions });
    const { ctx, container } = await opened({ components: c.components });
    const box = () => container.querySelector('.seqv-track-actions') as HTMLElement | null;
    const move = pointer(ctx);
    const dom = lines(ctx, 'chains.main', 't.dom')[0];
    await act(async () => move('move', 'chains.main', dom.header, dom.y));
    expect([box()?.style.top, box()?.style.height]).toEqual([`${dom.top}px`, `${dom.height}px`]);
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex('chains.main')));
    await act(async () => ctx.views.setCurrent('msa'));
    await settle();
    const section = 'msa.main';
    const rows = lines(ctx, section, 't.rows');
    expect(rows.length).toBeGreaterThan(4);
    await act(async () => move('move', section, rows[3].header, rows[3].y));
    const first = box();
    c.reset();
    // Row 3, then row 1, then row 4 of the same track: the same box on each hovered row, no render.
    for (const k of [3, 1, 4]) {
      await act(async () => move('move', section, rows[k].header, rows[k].y));
      expect([k, box()?.style.top, box()?.style.height]).toEqual([k, `${rows[k].top}px`, `${rows[k].height}px`]);
    }
    expect([box(), c.counts.TrackActions]).toEqual([first, 0]);
    // The pointer leaves the surface for the menu (one dispatch): it stays on row 4, also when the view moves under it.
    await act(async () => {
      ctx.input.leave(ctx.layout.sectionIndex(section));
      box()?.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    });
    await settle();
    await act(async () => ctx.viewports.get(section)?.panBy(10, 0, 'api'));
    expect([box(), box()?.style.top]).toEqual([first, `${rows[4].top}px`]);
    await act(async () => box()?.dispatchEvent(new PointerEvent('pointerout', { bubbles: true })));
    await settle();
    expect(box()).toBeNull();
  });

  it("the legend on the canvas's grid: its toggle at the labels' inset, its items from the data left", async () => {
    const { ctx, container } = await opened();
    const surface = container.querySelector('.seqv-surface') as HTMLElement;
    const left = (sel: string) => (container.querySelector(sel) as HTMLElement).getBoundingClientRect().left;
    const dataLeft = () => ctx.layout.get('chains.main')?.current.dataLeft;
    expect(container.querySelectorAll('.seqv-legend-item').length).toBeGreaterThan(0);
    for (const headerWidth of [150, 220]) {
      await act(async () => {
        ctx.setRenderProps({ headerWidth });
        ctx.layout.get('chains.main'); // As the canvas reads it on its next frame (the layout builds on read).
      });
      await settle();
      const x0 = surface.getBoundingClientRect().left;
      expect([headerWidth, left('.seqv-legend-item') - x0]).toEqual([headerWidth, dataLeft()]);
      expect(left('.seqv-legend summary') - x0).toBe(10);
    }
  });
});
