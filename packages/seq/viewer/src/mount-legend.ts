/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { colorHex } from '@molstar/seqx-core/theme/index';
import { type LegendChannels, LegendChannelsModes, channelsRow, highlightColorOf, legendRows } from './legends.js';

export type { LegendChannels };

/** The channels the row names before "+n": it is one line, so its height holds whatever the channels or the text size. */
const MaxChannelsShown = 2;
import type { SeqViewer } from './viewer.js';

// The `Legend` part (_spec/README.md "Legends"): a part the host places (beside a track header, in a panel, in a Mol*
// host's legend control), never drawn on the canvas. One tab stop, a list with a roving focus: arrows move between
// entries, Enter or Space selects a value's items (as a click), Cmd/Ctrl+Enter toggles them, focus hovers them; the
// hovered item's entry is outlined; the active highlight channels are rows with a clear. Thin over the DOM-free rows
// (`legendRows`), which the headless twin reads. Everything it adds is taken away by the function it returns.

export interface LegendPartOptions {
  /** A track's or a scale's address: only that one's legend (every legend when absent). */
  readonly of?: string;
  /** Entries a legend lists before "+n more" (12). */
  readonly maxItems?: number;
  /**
   * The highlight channels' row (_spec/README.md "Legends"): `reserve` (the default) holds its height whenever the part
   * is mounted, empty with no channel, so a highlight moves nothing below; `collapse` has the row only while a channel
   * is active.
   */
  readonly channels?: LegendChannels;
}

/** Mounts the legends of a viewer into `element` (`ToolbarMount`'s shape with options). */
export function mountLegend(element: HTMLElement, viewer: SeqViewer, options: LegendPartOptions = {}): () => void {
  const channels = options.channels ?? 'reserve';
  if (!LegendChannelsModes.includes(channels))
    throw new RangeError(`legend part: channels must be one of ${LegendChannelsModes.join(', ')}, not "${channels}"`);
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
  host.style.cssText = 'display: flex; flex-direction: column; align-items: flex-start; gap: 12px;';
  const channelsEl = doc.createElement('div');
  root.append(host);
  element.append(root);
  let focus = 0;
  let items: HTMLElement[] = [];

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
    shownKey = highlightKey();
    host.replaceChildren();
    items = [];
    const rows = legendRows(viewer, options);
    for (const row of rows) {
      const group = doc.createElement('div');
      group.setAttribute('role', 'group');
      group.setAttribute('aria-label', row.title);
      group.style.cssText = 'display: flex; flex-wrap: wrap; align-items: center; gap: 4px 8px;';
      const title = doc.createElement('span');
      title.textContent = row.title + (row.note ? ` (${row.note})` : '');
      // a group's title: muted, 600, 11 px (_spec/README.md "Legends")
      title.style.cssText = 'color: var(--msui-text-muted, gray); font-weight: 600; font-size: 11px;';
      group.append(title);
      for (const e of row.entries) {
        const item = doc.createElement('span');
        item.setAttribute('role', 'option');
        item.setAttribute('aria-label', e.name);
        item.dataset.legend = row.key;
        if (e.value !== undefined) item.dataset.value = e.value;
        item.style.cssText = `display: inline-flex; align-items: center; gap: 4px; cursor: ${e.value !== undefined ? 'pointer' : 'default'}; border-radius: 3px; padding: 0 2px;${e.hovered ? ' outline: 2px solid var(--msui-accent, currentColor);' : ''}`;
        const swatch = doc.createElement('span');
        // the swatches match the canvas: never forced (forced colors keep data colors)
        swatch.style.cssText = `forced-color-adjust: none; display: inline-block; width: ${e.ramp ? 64 : 10}px; height: 10px; background: ${e.ramp ? `linear-gradient(to right, ${e.ramp.map((c) => colorHex(c)).join(', ')})` : e.color !== undefined ? colorHex(e.color) : 'transparent'};`;
        item.append(swatch, doc.createTextNode(e.label));
        if (e.value !== undefined) {
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
          item.addEventListener(
            'pointerenter',
            () => void viewer.execute('legend.hover', { address: row.layer, entry: value }),
          );
          item.addEventListener('pointerleave', () => void viewer.execute('legend.hover', {}));
        }
        items.push(item);
        group.append(item);
      }
      if (row.more) {
        const more = doc.createElement('span');
        more.textContent = `+${row.more} more`;
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
    if (focusedClear !== undefined) {
      const again = [...root.querySelectorAll<HTMLElement>('button[data-clear]')].find(
        (b) => b.dataset.clear === focusedClear,
      );
      (again ?? host).focus();
    }
    focus = Math.min(focus, Math.max(0, items.length - 1));
    items.forEach((it, i) => it.setAttribute('aria-selected', String(i === focus)));
    if (items[focus]) host.setAttribute('aria-activedescendant', '');
  };

  const entryAt = (i: number) => items[i];
  const onKey = (ev: KeyboardEvent) => {
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
    // focus on an entry hovers its items, as the pointer does
    const f = entryAt(focus);
    const row = f && legendRows(viewer, options).find((r) => r.key === f.dataset.legend);
    if (f?.dataset.value !== undefined && row)
      void viewer.execute('legend.hover', { address: row.layer, entry: f.dataset.value });
  };
  host.addEventListener('keydown', onKey);
  host.addEventListener('blur', () => void viewer.execute('legend.hover', {}));
  const unlisten = [
    viewer.legends.legends.listen(draw),
    viewer.legends.hover.listen(draw),
    viewer.legends.cvdSafe.listen(draw),
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
