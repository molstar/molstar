/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { colorHex } from '@molstar/seqx-core/theme/index';
import {
  type LegendChannels,
  LegendChannelsModes,
  type LegendSelection,
  channelsRow,
  highlightColorOf,
  legendRows,
} from './legends.js';

export type { LegendChannels };

/** The channels the row names before "+n": it is one line, so its height holds whatever the channels or the text size. */
const MaxChannelsShown = 2;
import type { SeqViewer } from './viewer.js';

// The `Legend` part (_spec/README.md "Legends"): a part the host places (beside a track header, in a panel, in a Mol*
// host's legend control), never drawn on the canvas. One tab stop, a list with a roving focus: arrows move between
// entries, Enter or Space selects a value's items (as a click), Cmd/Ctrl+Enter toggles them, focus hovers them; the
// hovered item's entry is outlined; the active highlight channels are rows with a clear. Thin over the DOM-free rows
// (`legendRows`), which the headless twin reads. Everything it adds is taken away by the function it returns.

export interface LegendPartOptions extends LegendSelection {
  /**
   * `full` (the default): a group a legend, its entries labelled. `compact`: one line a legend, its title then its
   * swatches unlabelled (each named for the accessible name and on hover), a ramp a short bar with its domain's ends,
   * steps joined cells, a "+n" button past `maxItems` (it opens the full list in place); the same entries, focus and gestures.
   */
  readonly form?: 'full' | 'compact';
  /** `vertical` (the default): the legends stacked. `horizontal`: in a row, wrapping (a band under the canvas). */
  readonly orient?: 'vertical' | 'horizontal';
  /**
   * The highlight channels' row (_spec/README.md "Legends"): `reserve` (the default) holds its height whenever the part
   * is mounted, empty with no channel, so a highlight moves nothing below; `collapse` has the row only while a channel
   * is active.
   */
  readonly channels?: LegendChannels;
}

let mounts = 0;

/** An open legend's height: its list scrolls past it. */
const OpenMaxHeight = '14em';

/** Mounts the legends of a viewer into `element` (`ToolbarMount`'s shape with options). */
export function mountLegend(element: HTMLElement, viewer: SeqViewer, options: LegendPartOptions = {}): () => void {
  const channels = options.channels ?? 'reserve';
  if (!LegendChannelsModes.includes(channels))
    throw new RangeError(`legend part: channels must be one of ${LegendChannelsModes.join(', ')}, not "${channels}"`);
  const compact = options.form === 'compact';
  const horizontal = options.orient === 'horizontal';
  if (options.form !== undefined && !['full', 'compact'].includes(options.form))
    throw new RangeError(`legend part: form must be full or compact, not "${options.form}"`);
  if (options.orient !== undefined && !['vertical', 'horizontal'].includes(options.orient))
    throw new RangeError(`legend part: orient must be vertical or horizontal, not "${options.orient}"`);
  const doc = element.ownerDocument;
  // the part's root holds the legends' listbox and, outside it (a listbox holds options, not buttons), the channels row
  const root = doc.createElement('div');
  root.className = 'seqx-legend';
  root.style.cssText =
    'font: 12px var(--msui-font-sans, sans-serif); color: var(--msui-text, inherit); display: flex; flex-direction: column; align-items: flex-start; gap: 12px; padding: 12px 0 0;';
  const host = doc.createElement('div');
  host.tabIndex = 0;
  host.setAttribute('role', 'listbox');
  host.setAttribute('aria-label', 'Legend');
  host.style.cssText = horizontal
    ? 'display: flex; flex-flow: row wrap; align-items: flex-start; gap: 8px 24px;'
    : 'display: flex; flex-direction: column; align-items: flex-start; gap: 12px;';
  host.dataset.form = compact ? 'compact' : 'full';
  host.dataset.orient = horizontal ? 'horizontal' : 'vertical';
  const channelsEl = doc.createElement('div');
  root.append(host);
  element.append(root);
  let focus = 0;
  // the focused entry's name, said to a screen reader through `aria-activedescendant` and, in the compact form (where
  // the labels are not text), shown beside the groups while the part has the focus
  const idBase = `seqx-legend-${++mounts}`;
  const nameEl = doc.createElement('div');
  nameEl.setAttribute('aria-hidden', 'true');
  nameEl.dataset.focusName = '';
  nameEl.style.cssText = 'font-size: 11px; min-height: 1.4em; color: var(--msui-text-muted, gray);';
  if (compact) root.append(nameEl);
  const syncFocus = () => {
    const it = items[focus];
    if (it) host.setAttribute('aria-activedescendant', it.id);
    else host.removeAttribute('aria-activedescendant');
    nameEl.textContent = compact && it && doc.activeElement === host ? (it.getAttribute('aria-label') ?? '') : '';
  };
  let items: HTMLElement[] = [];
  /** What each entry's element hovers (its layer and value); the elements are made again on every draw. */
  const targets = new WeakMap<HTMLElement, { address: string; entry: string }>();

  /** What the highlight rows show: the channels, their styles and colors (a channel's place changes none of it). */
  const highlightKey = () =>
    viewer.highlight.names.value
      .map((n) => `${n}|${viewer.highlight.get(n)?.style ?? 'tint'}|${highlightColorOf(viewer, n)}`)
      .join('\n');
  let shownKey = '';

  const draw = () => {
    // a Clear button that has the focus keeps it across the redraw (its channel's new button, else the part itself)
    const active = doc.activeElement;
    const focusedClear = active instanceof HTMLElement && root.contains(active) ? active.dataset.clear : undefined;
    const focusedMore = active instanceof HTMLElement && host.contains(active) ? active.dataset.legendMore : undefined;
    shownKey = highlightKey();
    host.replaceChildren();
    items = [];
    const rows = legendRows(viewer, options);
    for (const row of rows) {
      const group = doc.createElement('div');
      group.setAttribute('role', 'group');
      group.setAttribute('aria-label', row.title);
      // an open legend wraps and is bounded (a legend over hundreds of values scrolls, never the whole page)
      const nowrap = compact && !horizontal && !row.open;
      group.style.cssText = `display: flex; flex-wrap: ${nowrap ? 'nowrap' : 'wrap'}; max-width: 100%;${nowrap ? ' overflow: hidden;' : ''}${row.open ? ` max-height: ${OpenMaxHeight}; overflow-y: auto;` : ''} align-items: center; gap: ${compact ? '0 2px' : '4px 8px'};`;
      group.dataset.open = String(row.open);
      const title = doc.createElement('span');
      title.textContent = row.title + (row.note ? ` (${row.note})` : '');
      // a group's title: muted, 600, 11 px (_spec/README.md "Legends")
      title.style.cssText = `color: var(--msui-text-muted, gray); font-weight: 600; font-size: 11px;${compact ? ' margin-right: 6px;' : ''}`;
      group.append(title);
      const largest = Math.max(1e-9, ...row.entries.map((x) => x.size ?? 0));
      for (const e of row.entries) {
        const item = doc.createElement('span');
        item.setAttribute('role', 'option');
        item.id = `${idBase}-${items.length}`;
        item.setAttribute('aria-label', e.name);
        item.dataset.legend = row.key;
        // an entry with an action alone is a target (its items lit, selected): a swatch alone keeps the default cursor
        if (e.action && e.value !== undefined) item.dataset.value = e.value;
        // the compact form: the label on hover and in the name, never as text (a ramp's ends aside)
        if (compact) item.title = e.name;
        item.style.cssText = `display: inline-flex; align-items: center; gap: 4px; cursor: ${e.action ? 'pointer' : 'default'}; border-radius: 3px; padding: 0 2px;${e.hovered ? ' outline: 2px solid var(--msui-accent, currentColor);' : ''}`;
        const swatch = doc.createElement('span');
        // the swatches match the canvas: never forced (forced colors keep data colors)
        // a size key's mark is a circle, its size relative to the row's largest (14 px across, in a 14 px cell): the
        // channel's range is in the layer's own units (a factor, a radius), so the key shows the ratio
        const dot = e.size === undefined ? 0 : Math.max(3, Math.round((14 * e.size) / largest));
        swatch.style.cssText = `forced-color-adjust: none; display: inline-block; width: ${dot || (e.ramp ? (compact ? 40 : 64) : 10)}px; height: ${dot || 10}px; ${dot ? `margin: ${(14 - dot) / 2}px 2px; border-radius: 50%; background: var(--msui-text-muted, gray)` : `background: ${e.ramp ? `linear-gradient(to right, ${e.ramp.map((c) => colorHex(c)).join(', ')})` : e.color !== undefined ? colorHex(e.color) : 'transparent'}`};`;
        if (!compact) item.append(swatch, doc.createTextNode(e.label));
        else if (e.ends) item.append(doc.createTextNode(e.ends[0]), swatch, doc.createTextNode(e.ends[1]));
        else item.append(swatch);
        if (e.action && e.value !== undefined) {
          const value = e.value;
          item.addEventListener(
            'click',
            (ev) =>
              void viewer.execute('legend.select', {
                address: row.layer,
                entry: value,
                toggle: ev.metaKey || ev.ctrlKey,
              }),
          );
          // hovered by the part's own listeners (below): an entry is replaced when the rows redraw under the pointer, and a
          // removed node never fires its `pointerleave`
          targets.set(item, { address: row.layer, entry: value });
        }
        // a step's cells join: the group's gap taken back
        if (compact && e.step) item.style.cssText += 'padding: 0; margin-right: -2px;';
        items.push(item);
        group.append(item);
      }
      if (row.beyond) {
        const rest = doc.createElement('span');
        rest.dataset.beyond = String(row.beyond);
        rest.textContent = `${row.beyond} more not shown`;
        rest.style.cssText = 'color: var(--msui-text-muted, gray); font-size: 11px;';
        group.append(rest);
      }
      if (row.more || row.open) {
        // "+n" opens the legend's full list in place, "Show fewer" closes it (the session's state, per legend)
        const more = doc.createElement('button');
        more.type = 'button';
        more.dataset.legendMore = row.id;
        more.setAttribute('aria-expanded', String(row.open));
        const label = row.open ? 'Show fewer' : `Show ${row.more} more`;
        // the accessible name starts with the visible text (WCAG 2.5.3)
        more.setAttribute(
          'aria-label',
          `${row.open ? 'Show fewer' : `+${row.more}, show ${row.more} more`}: ${row.title}`,
        );
        more.title = label;
        more.textContent = row.open ? 'Show fewer' : `+${row.more}`;
        more.style.cssText =
          'font: inherit; color: inherit; background: none; border: 0; padding: 0 2px; cursor: pointer; text-decoration: underline dotted;';
        more.addEventListener('click', () => viewer.legends.setOpen(row.id, !row.open));
        group.append(more);
      }
      host.append(group);
    }
    // the active highlight channels: their color (as the canvas draws it, by the channel's style) and a clear
    const names = viewer.highlight.names.value;
    const slot = channelsRow(channels, names.length);
    channelsEl.replaceChildren();
    channelsEl.dataset.channels = channels;
    // one line whatever the channels or the text size: the first `MaxChannelsShown`, the rest a "+n" chip; under `reserve`
    // with none, a hidden copy of a channel's row holds the line's real height
    channelsEl.style.cssText =
      'display: flex; flex-wrap: nowrap; align-items: center; gap: 12px; overflow: hidden; white-space: nowrap; max-width: 100%;';
    const channelItem = (h: string, hidden: boolean) => {
      const style = viewer.highlight.get(h)?.style ?? 'tint';
      const row = doc.createElement('span');
      row.style.cssText = 'display: inline-flex; align-items: center; gap: 4px;';
      if (hidden) {
        row.setAttribute('aria-hidden', 'true');
        row.style.visibility = 'hidden';
      } else row.dataset.highlight = h;
      const color = colorHex(highlightColorOf(viewer, h));
      const swatch = doc.createElement('span');
      swatch.setAttribute('aria-hidden', 'true');
      swatch.dataset.style = style;
      swatch.dataset.color = color;
      swatch.style.cssText = `forced-color-adjust: none; display: inline-block; width: 14px; height: 10px; box-sizing: border-box; ${
        style === 'outline'
          ? `border: 1px solid ${color};`
          : style === 'band'
            ? `border-bottom: 3px solid ${color};`
            : `background: ${color}; opacity: 0.5;`
      }`;
      row.append(swatch, doc.createTextNode(`Highlight: ${h} `));
      const clear = doc.createElement('button');
      clear.type = 'button';
      clear.textContent = 'Clear';
      if (hidden) clear.tabIndex = -1;
      else {
        clear.setAttribute('aria-label', `Clear the ${h} highlight`);
        clear.dataset.clear = h;
        clear.addEventListener('click', () => void viewer.execute('highlight.clear', { channel: h }));
      }
      row.append(clear);
      return row;
    };
    if (slot.present) {
      if (names.length) {
        for (const h of names.slice(0, MaxChannelsShown)) channelsEl.append(channelItem(h, false));
        if (names.length > MaxChannelsShown) {
          const more = doc.createElement('span');
          more.dataset.more = String(names.length - MaxChannelsShown);
          more.textContent = `+${names.length - MaxChannelsShown}`;
          more.title = names.slice(MaxChannelsShown).join(', ');
          channelsEl.append(more);
        }
      } else channelsEl.append(channelItem('legend', true));
      if (channelsEl.parentNode !== root) root.append(channelsEl);
    } else channelsEl.remove();
    if (focusedMore !== undefined)
      (
        [...host.querySelectorAll<HTMLElement>('button[data-legend-more]')].find(
          (b) => b.dataset.legendMore === focusedMore,
        ) ?? host
      ).focus();
    if (focusedClear !== undefined) {
      const again = [...root.querySelectorAll<HTMLElement>('button[data-clear]')].find(
        (b) => b.dataset.clear === focusedClear,
      );
      (again ?? host).focus();
    }
    focus = Math.min(focus, Math.max(0, items.length - 1));
    items.forEach((it, i) => it.setAttribute('aria-selected', String(i === focus)));
    syncFocus();
  };

  const entryAt = (i: number) => items[i];
  const onKey = (ev: KeyboardEvent) => {
    // a "+n" button has its own keys (Enter, Space: a click)
    if (ev.target instanceof HTMLButtonElement) return;
    if (!items.length) return;
    const it = entryAt(focus);
    if (ev.key === 'ArrowRight' || ev.key === 'ArrowDown') focus = Math.min(items.length - 1, focus + 1);
    else if (ev.key === 'ArrowLeft' || ev.key === 'ArrowUp') focus = Math.max(0, focus - 1);
    else if ((ev.key === 'Enter' || ev.key === ' ') && it?.dataset.value !== undefined) {
      const row = legendRows(viewer, options).find((r) => r.key === it.dataset.legend);
      void viewer.execute('legend.select', {
        address: row?.layer,
        entry: it.dataset.value,
        toggle: ev.key === 'Enter' && (ev.metaKey || ev.ctrlKey),
      });
    } else return;
    ev.preventDefault();
    items.forEach((x, i) => x.setAttribute('aria-selected', String(i === focus)));
    syncFocus();
    // focus on an entry hovers its items, as the pointer does
    const f = entryAt(focus);
    const row = f && legendRows(viewer, options).find((r) => r.key === f.dataset.legend);
    if (f?.dataset.value !== undefined && row)
      void viewer.execute('legend.hover', { address: row.layer, entry: f.dataset.value });
  };
  // the pointer: one pair of listeners on the part, so a redraw under the pointer cannot lose a leave. Over an entry it
  // hovers the entry's items, over anything else in the part (a title, a gap) or out of the part it clears the hover.
  host.addEventListener('pointerover', (ev) => {
    const el = ev.target instanceof Element ? ev.target.closest<HTMLElement>('[role="option"]') : null;
    const t = el && targets.get(el);
    void viewer.execute('legend.hover', t ? { address: t.address, entry: t.entry } : {});
  });
  host.addEventListener('pointerleave', () => void viewer.execute('legend.hover', {}));
  host.addEventListener('keydown', onKey);
  host.addEventListener('focus', syncFocus);
  host.addEventListener('blur', () => {
    nameEl.textContent = '';
    void viewer.execute('legend.hover', {});
  });
  const unlisten = [
    viewer.legends.legends.listen(draw),
    viewer.legends.hover.listen(draw),
    viewer.legends.cvdSafe.listen(draw),
    viewer.legends.open.listen(draw),
    // what is drawn moved (a pan, a tile arriving): an entry's action with it
    viewer.legends.drawn.listen(draw),
    // a channel's place changing is no change of what the rows show; a theme's tokens change their colors
    viewer.highlight.version.listen(() => highlightKey() !== shownKey && draw()),
    viewer.tokensVersion.listen(() => highlightKey() !== shownKey && draw()),
  ];
  draw();
  return () => {
    for (const u of unlisten) u();
    root.remove();
  };
}

/** The default mount (`SeqViewerSpec.ui.legends`): every legend, after whatever is in the target. */
export const legendsMount = (viewer: SeqViewer, target: HTMLElement): (() => void) => mountLegend(target, viewer);
