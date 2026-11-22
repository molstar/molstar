/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';

// The causes a genome runtime's read met (_spec/README.md "Data", Genome runtimes' internals): one issue a cause, with
// its count and its first ten places (EC-86), made by the format's own codes. Shared by the BED, GFF and VCF runtimes,
// not published (`./runtimes/internal/*`).

/** Places an issue keeps. */
export const MaxPlaces = 10;

/** The causes a read met, one issue each with its count and first places ("line 12"). */
export class Causes<Code extends string> {
  private readonly byKey = new Map<string, { code: Code; message: string; count: number; places: string[] }>();

  /** `issue` makes a cause's issue of the format's codes: its code, its message with count and places, its address. */
  constructor(private readonly issue: (code: Code, message: string, address: string) => Issue) {}

  add(code: Code, message: string, place?: string) {
    const key = `${code}|${message}`;
    let c = this.byKey.get(key);
    if (!c) {
      c = { code, message, count: 0, places: [] };
      this.byKey.set(key, c);
    }
    c.count++;
    if (place !== undefined && c.places.length < MaxPlaces) c.places.push(place);
  }

  /** Whether a cause with `code` was said. */
  has(code: Code): boolean {
    for (const c of this.byKey.values()) if (c.code === code) return true;
    return false;
  }

  get size() {
    return this.byKey.size;
  }

  issues(address: string): Issue[] {
    return [...this.byKey.values()].map((c) => {
      const where = c.places.length ? ` (${c.places.join(', ')}${c.count > c.places.length ? ', …' : ''})` : '';
      return this.issue(c.code, `${c.message}: ${c.count}${where}`, address);
    });
  }
}
