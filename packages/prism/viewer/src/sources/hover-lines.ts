/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

const Entities: Record<string, string> = { amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ', '#39': "'" };

/**
 * Mol*'s loci labels (`plugin.behaviors.labels.highlight`: strings with `</br>` between lines and a `<small>` count)
 * as plain lines for the shell's hover label. Edge rate: once per hovered element, not per pointer move.
 */
export function hoverLines(labels: readonly unknown[]): readonly string[] {
  const out: string[] = [];
  for (const label of labels) {
    if (typeof label !== 'string') continue;
    const text = label
      .replace(/<\s*\/?\s*br\s*\/?\s*>/gi, '\n')
      .replace(/<[^>]*>/g, '')
      .replace(/&(#?\w+);/g, (m, name: string) => Entities[name] ?? m)
      .replace(/\s*\|\|\s*/g, ' ');
    for (const line of text.split('\n')) {
      const trimmed = line.trim();
      if (trimmed) out.push(trimmed);
    }
  }
  return out;
}
