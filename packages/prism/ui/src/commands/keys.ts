/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { KeyChord } from './types.js';

/** Is this an Apple platform (⌘ is the command key)? Elsewhere `mod` is Ctrl. */
export function isMac(): boolean {
  if (typeof navigator === 'undefined') return false;
  const platform =
    (navigator as Navigator & { userAgentData?: { platform?: string } }).userAgentData?.platform ?? navigator.platform;
  return /mac|iphone|ipad|ipod/i.test(platform ?? '');
}

const MacKeys: Record<string, string> = { mod: '⌘', alt: '⌥', shift: '⇧' };
const OtherKeys: Record<string, string> = { mod: 'Ctrl', alt: 'Alt', shift: 'Shift' };
const Named: Record<string, string> = {
  arrowup: '↑',
  arrowdown: '↓',
  arrowleft: '←',
  arrowright: '→',
  escape: 'Esc',
  enter: '↵',
  ' ': 'Space',
};

/** The key caps of a chord for this platform: `['⌘', 'K']` on macOS, `['Ctrl', 'K']` elsewhere. */
export function chordKeys(chord: KeyChord, mac = isMac()): string[] {
  const names = mac ? MacKeys : OtherKeys;
  const keys: string[] = [];
  // macOS order is ⌥⇧⌘; elsewhere Ctrl+Alt+Shift.
  if (mac) {
    if (chord.alt) keys.push(names.alt);
    if (chord.shift) keys.push(names.shift);
    if (chord.mod) keys.push(names.mod);
  } else {
    if (chord.mod) keys.push(names.mod);
    if (chord.alt) keys.push(names.alt);
    if (chord.shift) keys.push(names.shift);
  }
  const key = chord.key.toLowerCase();
  keys.push(Named[key] ?? (key.length === 1 ? key.toUpperCase() : chord.key));
  return keys;
}

/** A chord as text: `⌘K` on macOS, `Ctrl+K` elsewhere. */
export function formatChord(chord: KeyChord, mac = isMac()): string {
  return chordKeys(chord, mac).join(mac ? '' : '+');
}

/** A single character that is not a letter, a digit or a space: it implies its Shift (`?` is Shift+/ on most layouts). */
export function isSymbolKey(key: string): boolean {
  return /^[^\p{L}\p{N}\s]$/u.test(key);
}

/** A chord's identity: two chords with the same id are the same keys (a symbol's Shift is part of the symbol). */
export function chordId(chord: KeyChord): string {
  const key = chord.key.toLowerCase();
  const shift = !!chord.shift && !isSymbolKey(chord.key);
  return `${chord.mod ? 'mod+' : ''}${chord.alt ? 'alt+' : ''}${shift ? 'shift+' : ''}${key}`;
}

/** Keys that are only modifiers, or that a chord cannot be (Tab moves the focus). */
const NotChordKeys = new Set([
  'Shift',
  'Control',
  'Alt',
  'Meta',
  'OS',
  'CapsLock',
  'Dead',
  'Unidentified',
  'Process',
  'Tab',
]);

/**
 * The chord a keydown is, for recording a remap (commands.md §6): none for a modifier alone, Tab, or the platform's
 * other modifier (Ctrl on macOS, ⊞ elsewhere: commands never use it). ⌥T on macOS types "†": its letter is read from
 * `code`.
 */
export function chordFromEvent(
  e: { key: string; code?: string; altKey?: boolean; shiftKey?: boolean; metaKey?: boolean; ctrlKey?: boolean },
  mac = isMac(),
): KeyChord | undefined {
  if (NotChordKeys.has(e.key)) return undefined;
  if (mac ? e.ctrlKey : e.metaKey) return undefined;
  const mod = mac ? !!e.metaKey : !!e.ctrlKey;
  let key = e.key;
  const code = /^(?:Key([A-Z])|Digit([0-9]))$/.exec(e.code ?? '');
  if (code && (e.altKey || mod) && !/^[a-z0-9]$/i.test(key)) key = code[1] ?? code[2] ?? key;
  if (key.length === 1) key = key.toLowerCase();
  const shift = !!e.shiftKey && !isSymbolKey(key);
  return { key, ...(mod ? { mod: true } : {}), ...(e.altKey ? { alt: true } : {}), ...(shift ? { shift: true } : {}) };
}
