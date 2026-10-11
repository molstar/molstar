/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Up to `max` known values closest to `value` (prefix and containment first, then edit distance). */
export function suggest(value: string, known: readonly string[], max = 5): string[] {
  const scored = known.map((k) => ({ k, s: score(value, k) }));
  scored.sort((a, b) => a.s - b.s || a.k.localeCompare(b.k));
  return scored.slice(0, max).map((x) => x.k);
}

/**
 * The hint for an unknown id (errors say what would work): the known ones, closest to what was asked first, a few at
 * most.
 */
export function closestHint(what: string, id: string, known: readonly string[], max = 6): string {
  if (!known.length) return `There are no ${what} yet.`;
  const shown = suggest(id, known, max);
  const more = known.length - shown.length;
  return `${what[0]?.toUpperCase()}${what.slice(1)}, closest first: ${shown.join(', ')}${more > 0 ? ` (and ${more} more)` : ''}.`;
}

function score(a: string, b: string): number {
  const x = a.toLowerCase();
  const y = b.toLowerCase();
  const tail = (s: string) => s.split('.').pop() ?? s;
  if (y.includes(x) || tail(y).startsWith(tail(x))) return 0;
  return levenshtein(tail(x), tail(y));
}

function levenshtein(a: string, b: string): number {
  const row = Array.from({ length: b.length + 1 }, (_, i) => i);
  for (let i = 1; i <= a.length; i++) {
    let prev = row[0] as number;
    row[0] = i;
    for (let j = 1; j <= b.length; j++) {
      const tmp = row[j] as number;
      row[j] = Math.min((row[j] as number) + 1, (row[j - 1] as number) + 1, prev + (a[i - 1] === b[j - 1] ? 0 : 1));
      prev = tmp;
    }
  }
  return row[b.length] as number;
}
