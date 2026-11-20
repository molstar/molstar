/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { formatChord, isMac } from '@molstar/ui/commands/keys';
import type { KeyChord } from '@molstar/ui/commands/types';
import { type GestureBinding, parseChord } from './bindings.js';
import { changeKindOf } from './commands.js';
import type { SeqViewer } from './viewer.js';

// What the view answers to (_spec/README.md "Input", The keys panel), as data: generated from the viewer's commands with
// their keys (`view.keys`, else a command's shortcut, as the input reads them) and from the binding table's gestures,
// never hand-written; grouped (move, zoom, select, edit, view), each entry a command or a method once with every key
// and gesture that runs it. The keys panel draws it; a host draws its own from the same list.

export type KeyGroupId = 'move' | 'zoom' | 'select' | 'edit' | 'view';

/** One thing the view does, with what runs it: its keys as the platform writes them and its gestures. */
export interface KeyEntry {
  /** The command's id, or the method's name (`pan`, `drag-select`). */
  readonly action: string;
  /** What it does: the command's title. */
  readonly name: string;
  readonly keys: readonly string[];
  readonly gestures: readonly string[];
}

export interface KeyGroup {
  readonly id: KeyGroupId;
  readonly title: string;
  readonly entries: readonly KeyEntry[];
}

const GroupTitles: Readonly<Record<KeyGroupId, string>> = {
  move: 'Move',
  zoom: 'Zoom',
  select: 'Select',
  edit: 'Edit',
  view: 'View',
};

/** The model methods by name (a binding's action that is no command). */
const MethodNames: Readonly<Record<string, string>> = {
  pan: 'Pan',
  zoom: 'Zoom',
  'drag-select': 'Select a range',
  'drag-add': 'Add a range',
  'drag-subtract': 'Remove a range',
  'drag-extend': 'Extend the range',
};

const GestureWords: Readonly<Record<string, string>> = {
  click: 'Click',
  'double-click': 'Double-click',
  drag: 'Drag',
  wheel: 'Wheel',
  pinch: 'Pinch',
  'back-button': 'Back button',
  'forward-button': 'Forward button',
};

const TargetWords: Readonly<Record<string, string>> = {
  any: '',
  cell: 'on data',
  item: 'on an item',
  row: 'on a row',
  column: 'on a column',
  track: 'on a track',
  empty: 'on empty space',
  axis: 'on an axis',
  break: 'on a break',
  'header-slot': 'on a header’s action',
  'header-state': 'on a header’s state',
  'header-label': 'on a track’s label',
  'section-header': 'on a section’s header',
  window: 'on a window',
  'tree-node': 'on a tree node',
  'tree-toggle': 'on the tree’s toggle',
  'summary-row': 'on a folded clade',
  'column-title': 'on a column title',
  'columns-hidden': 'on the hidden columns',
};

/** The group of a command or a method, by what it changes. */
function groupOf(action: string): KeyGroupId {
  if (action === 'pan' || /^(view\.pan-|cursor\.(prev|next|up|down)$|window\.pan-)/.test(action)) return 'move';
  if (action === 'zoom' || /^(view\.(zoom-|fit)|rows\.zoom-|window\.(grow|shrink)|navigator\.zoom)/.test(action))
    return 'zoom';
  if (/^drag-|^selection\.(?!bookmark|copy)|^cursor\.|^clades\.select/.test(action) && action !== 'cursor.focus')
    return 'select';
  if (changeKindOf(action) === 'patch' || changeKindOf(action) === 'preview') return 'edit';
  return 'view';
}

/** A gesture as text: `Shift+Drag on data`, `⌘Wheel`. */
function gestureText(b: GestureBinding, mac: boolean): string {
  const word = GestureWords[b.gesture] ?? b.gesture;
  const mods = b.modifiers === '*' ? [] : (b.modifiers ?? []);
  const chord: KeyChord = {
    key: word,
    ...(mods.includes('shift') ? { shift: true } : {}),
    ...(mods.includes('mod') ? { mod: true } : {}),
    ...(mods.includes('alt') ? { alt: true } : {}),
  };
  const where = TargetWords[b.target] ?? `on ${b.target}`;
  return `${mods.includes('space') ? 'Space+' : ''}${formatChord(chord, mac)}${where ? ` ${where}` : ''}`;
}

/** The view's keys and gestures, grouped; empty groups left out. */
export function describeKeys(viewer: SeqViewer, mac = isMac()): KeyGroup[] {
  const entries = new Map<string, { action: string; name: string; keys: string[]; gestures: string[] }>();
  const entry = (action: string, name: string) => {
    let e = entries.get(action);
    if (!e) entries.set(action, (e = { action, name, keys: [], gestures: [] }));
    return e;
  };
  // the commands' keys, as the input matches them: the configured chords, else the shortcut
  const configured = viewer.config.view.keys;
  for (const c of viewer.commands) {
    const listed = configured[c.id];
    const chords = listed
      ? listed.map(parseChord).filter((k): k is KeyChord => k !== undefined)
      : c.shortcut
        ? [c.shortcut]
        : [];
    if (!chords.length) continue;
    const e = entry(c.id, c.title);
    for (const k of chords) {
      const text = formatChord(k, mac);
      if (!e.keys.includes(text)) e.keys.push(text);
    }
  }
  // the gestures: what each binding runs (the hover is the pointer itself: left out)
  for (const b of viewer.bindings.bindings) {
    if (!b.action || b.action === 'hover') continue;
    const name = MethodNames[b.action] ?? viewer.commands.find((c) => c.id === b.action)?.title ?? b.action;
    const e = entry(b.action, name);
    const text = gestureText(b, mac);
    if (!e.gestures.includes(text)) e.gestures.push(text);
  }
  const groups: KeyGroup[] = [];
  for (const id of ['move', 'zoom', 'select', 'edit', 'view'] as const) {
    const of = [...entries.values()].filter((e) => groupOf(e.action) === id);
    if (of.length) groups.push({ id, title: GroupTitles[id], entries: of });
  }
  return groups;
}
