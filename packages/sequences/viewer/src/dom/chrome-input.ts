/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The chrome's input through the public streams (seqspec:interaction#streams): one delegated listener per mounted UI,
// at its root, reports a click, a context-menu request, a hover and a leave on a part (its root tagged
// `data-seqv-part`, as the default UIs' parts and a host's own are; a default part's class is the fallback), with the
// part's section, track and key. The surfaces are the engine's (`DomInput`). Edge-rate: nothing is looked up while a stream is unobserved.

import type { SeqViewerContext } from '../context.js';
import { Modifier } from '../input.js';
import type { SeqChromeKind, SeqChromeStream, SeqPart } from '../interaction.js';

/** The default UIs' parts by their root's class: the fallback where a root carries no `data-seqv-part`. */
const KINDS: readonly (readonly [string, SeqChromeKind])[] = [
  ['seqv-tabs', 'viewTab'],
  ['seqv-track-actions', 'trackAction'],
  ['seqv-switches', 'sectionSwitch'],
  ['seqv-chips', 'hiddenChip'],
  ['seqv-section-header', 'sectionHeader'],
  ['seqv-row-control', 'rowControl'],
  ['seqv-legend', 'legend'],
  ['seqv-player', 'snapshotPlayer'],
  ['seqv-search', 'search'],
  ['seqv-export', 'export'],
  ['seqv-issues', 'issues'],
  ['seqv-hover-card', 'hoverCard'],
  ['seqv-tooltip', 'tooltip'],
  ['seqv-description', 'description'],
];
const SELECTOR = ['[data-seqv-part]', '.seqv-surface', ...KINDS.map(([c]) => `.${c}`)].join(', ');

const modifiers = (e: MouseEvent) =>
  (e.shiftKey ? Modifier.shift : 0) |
  (e.ctrlKey ? Modifier.ctrl : 0) |
  (e.altKey ? Modifier.alt : 0) |
  (e.metaKey ? Modifier.meta : 0);

/** The part a node is in, and its root; null: none, or a surface (the engine's). */
export function chromePartOf(node: EventTarget | null): { part: SeqPart; root: Element } | null {
  const el = node instanceof Element ? node : null;
  const root = el?.closest(SELECTOR);
  if (!el || !root || root.classList.contains('seqv-surface')) return null;
  const tagged = root.getAttribute('data-seqv-part');
  const kind = (tagged ?? KINDS.find(([c]) => root.classList.contains(c))?.[1]) as SeqChromeKind | undefined;
  if (!kind || kind === ('surface' as string)) return null;
  const section = el.closest('[data-seqv-section]')?.getAttribute('data-seqv-section') ?? null;
  const track = el.closest('[data-seqv-track]')?.getAttribute('data-seqv-track') ?? null;
  const keyed = el.closest('[data-seqv-key]');
  const key = keyed && root.contains(keyed) ? keyed.getAttribute('data-seqv-key') : null;
  return { part: { kind, section, track, group: null, key }, root };
}

/**
 * Listens at `root` (a mounted UI's `.seqv-root`, or a host's box of parts) and reports its chrome's input to
 * `ctx.interaction`; returns the removal.
 */
export function listenChrome(root: HTMLElement, ctx: SeqViewerContext): () => void {
  const it = ctx.interaction;
  let over: Element | null = null;
  let overKey: string | null = null;
  let last: SeqPart | null = null;
  const report = (stream: SeqChromeStream, e: MouseEvent, count?: number) => {
    const at = chromePartOf(e.target);
    if (!at) return null;
    const r = at.root.getBoundingClientRect();
    const place = { x: e.clientX - r.left, y: e.clientY - r.top, clientX: e.clientX, clientY: e.clientY };
    it.chrome(stream, at.part, e, place, { count, buttons: e.buttons, button: e.button, modifiers: modifiers(e) });
    return at;
  };
  const click = (e: MouseEvent) => {
    if (it.observed('click')) report('click', e, Math.max(1, e.detail));
  };
  const menu = (e: MouseEvent) => {
    if (it.observed('menu')) report('menu', e);
  };
  const leave = (e: MouseEvent) => {
    if (last && it.observed('leave')) {
      const place = { x: Number.NaN, y: Number.NaN, clientX: e.clientX, clientY: e.clientY };
      it.chrome('leave', last, e, place);
    }
    over = null;
    overKey = null;
    last = null;
  };
  const hover = (e: PointerEvent) => {
    if (!it.observed('hover') && !it.observed('leave')) return;
    const at = chromePartOf(e.target);
    if (!at) {
      if (last) leave(e);
      return;
    }
    if (at.root === over && at.part.key === overKey) return;
    over = at.root;
    overKey = at.part.key;
    last = at.part;
    if (it.observed('hover')) report('hover', e);
  };
  const out = (e: PointerEvent) => {
    if (last) leave(e);
  };
  root.addEventListener('click', click);
  root.addEventListener('contextmenu', menu);
  root.addEventListener('pointerover', hover);
  root.addEventListener('pointerleave', out);
  return () => {
    root.removeEventListener('click', click);
    root.removeEventListener('contextmenu', menu);
    root.removeEventListener('pointerover', hover);
    root.removeEventListener('pointerleave', out);
  };
}
