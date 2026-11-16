/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Chords (commands §6.5): written once with `mod`, shown per platform, compared by id, matched against key events,
// recorded from them. DOM-free: an event is anything with the fields read.

import { isApplePlatform } from '@molstar/core/util/browser';
import type { KeyChord } from './types.js';

/** The fields of a keydown the key layer reads. */
export interface KeyEventLike {
  key: string;
  altKey?: boolean;
  shiftKey?: boolean;
  metaKey?: boolean;
  ctrlKey?: boolean;
  /** Layout-independent key, e.g. "KeyT": read when `key` is not a plain letter (⌥T on macOS is "†"). */
  code?: string;
  /** Where the key went; a text-entry target keeps keys that type. */
  target?: unknown;
  /** A DOM event's composed path: across a shadow root `target` is its host, the path's first item the real one. */
  composedPath?(): EventTarget[];
  isComposing?: boolean;
  repeat?: boolean;
  defaultPrevented?: boolean;
}

/**
 * The platform the keys are for (commands §6.5): `mac` (⌘ is `mod`) or not (Ctrl is). Detected once from the browser
 * (core's `isApplePlatform`, false outside one); a host may set it (tests, a forced platform). Every default of `mac`
 * below reads it.
 */
export const KeyPlatform: { mac: boolean } = { mac: isApplePlatform() };

/** Is ⌘ the command key here? Elsewhere `mod` is Ctrl. */
export function isMac(): boolean {
  return KeyPlatform.mac;
}

const MacKeys = { mod: '⌘', alt: '⌥', shift: '⇧' };
const OtherKeys = { mod: 'Ctrl', alt: 'Alt', shift: 'Shift' };
const Named: Readonly<Record<string, string>> = {
  arrowup: '↑',
  arrowdown: '↓',
  arrowleft: '←',
  arrowright: '→',
  escape: 'Esc',
  enter: '↵',
  backspace: '⌫',
  delete: '⌦',
  ' ': 'Space',
};

/** The key caps of a chord for this platform: `['⌘', 'K']` on macOS, `['Ctrl', 'K']` elsewhere (a `kbd` each). */
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

/** No modifier at all (an arrow, a letter alone). */
export function isPlainChord(chord: KeyChord): boolean {
  return !chord.mod && !chord.alt && !chord.shift;
}

/** ⌘ or Ctrl without Alt never types a character (Alt does on macOS; Ctrl+Alt is AltGr). */
export function typesNothing(chord: KeyChord): boolean {
  return !!chord.mod && !chord.alt;
}

const PlainKey = /^[a-z0-9]$/;
const CodeKey = /^(?:Key([A-Z])|Digit([0-9]))$/;

/**
 * The letter or digit of the event's `code`, where the key itself is not one (⌥F on macOS is "ƒ"; a non-Latin
 * layout); `key` is the event's key in lower case.
 */
export function codeKey(e: KeyEventLike, key: string): string | undefined {
  if (!e.code || PlainKey.test(key)) return undefined;
  const m = CodeKey.exec(e.code);
  return m ? (m[1] ?? m[2] ?? '').toLowerCase() : undefined;
}

/**
 * Does the event press the chord? The modifiers exactly (the platform's other modifier never: commands do not use it);
 * the key by `key` first, then by `code` when the key is not a plain letter or digit.
 */
export function chordMatches(chord: KeyChord, e: KeyEventLike, mac: boolean): boolean {
  const mod = mac ? !!e.metaKey : !!e.ctrlKey;
  const other = mac ? !!e.ctrlKey : !!e.metaKey;
  if (!!chord.mod !== mod || other || !!chord.alt !== !!e.altKey) return false;
  // A symbol (`?` is Shift+/ on most layouts) already implies its modifiers; letters and digits match Shift strictly.
  if (!!chord.shift !== !!e.shiftKey && (chord.shift || !isSymbolKey(chord.key))) return false;
  const want = chord.key.toLowerCase();
  const key = e.key.toLowerCase();
  return key === want || codeKey(e, key) === want;
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
 * The chord a keydown is, for recording a remap (commands §6.3): none for a modifier alone, Tab, or the platform's other
 * modifier (Ctrl on macOS, ⊞ elsewhere). ⌥T on macOS types "†": its letter is read from `code`.
 */
export function chordFromEvent(e: KeyEventLike, mac = isMac()): KeyChord | undefined {
  if (NotChordKeys.has(e.key)) return undefined;
  if (mac ? e.ctrlKey : e.metaKey) return undefined;
  const mod = mac ? !!e.metaKey : !!e.ctrlKey;
  let key = e.key;
  const code = CodeKey.exec(e.code ?? '');
  if (code && (e.altKey || mod) && !/^[a-z0-9]$/i.test(key)) key = code[1] ?? code[2] ?? key;
  if (key.length === 1) key = key.toLowerCase();
  const shift = !!e.shiftKey && !isSymbolKey(key);
  return { key, ...(mod ? { mod: true } : {}), ...(e.altKey ? { alt: true } : {}), ...(shift ? { shift: true } : {}) };
}
