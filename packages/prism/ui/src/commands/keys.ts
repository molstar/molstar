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
