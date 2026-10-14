/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The codon label (ST7, seqspec:viewer#labels): on a block a streamed layer reaches through a path with a codon step
// taken from the nucleotides (the block's elements each three of the path's system), the hovered element's three
// positions in reading order (`Paths.fromAt`: both parts of one split by a gap in the map, a partial codon's known
// ones), their letters (a lazy molecule's from its stream, asked when not held: `LabelContext.letters`), complemented by
// the registry's `complement` (`complement_rna`) letter map on the codon step's minus strand, and the positions in the
// system's numbering (the unit's own first). Label-rate: its strings are its only allocation.

import { letterMapTable, type StreamedData, StreamProp } from '@molstar/seq-core/format';
import { type ComposedPath, CoordinateSystemKind, Coords, Paths } from '@molstar/seq-core/model';
import { HoverKind } from './hit.js';
import type { LabelProvider } from './labels.js';

const at = new Int32Array(3);
const letters = new Uint8Array(3);
/** The `complement` entry's table, by entry. */
const tables = new WeakMap<object, Uint8Array | null>();

/** The path of a streamed layer of the hovered section on block b whose map steps three of its system an element. */
function codonPath(lc: Parameters<LabelProvider['label']>[1], b: number): ComposedPath | null {
  const r = lc.view.resolved;
  for (let k = 0; k < r.layers.count; k++) {
    if (r.tracks.section[r.layers.track[k]] !== lc.section) continue;
    const p = lc.stream(k)?.source.paths?.[b] ?? null;
    if (p && p.map.span === 3 && !Coords.sameSystem(p.map.from, p.from)) return p;
  }
  return null;
}

export const codonLabel: LabelProvider = {
  name: 'codon',
  priority: 90,
  applies: (h) => h.kind === HoverKind.seq && h.block >= 0,
  label(h, lc, out) {
    const p = codonPath(lc, h.block);
    if (p === null) return;
    // Its positions start on a unit or a molecule (a column or a row has no letters of its own here).
    const K = CoordinateSystemKind;
    if (p.from.kind !== K.unit && p.from.kind !== K.molecule) return;
    let known = 0;
    for (let k = 0; k < 3; k++) {
      at[k] = Paths.fromAt(p, h.index, k);
      if (at[k] >= 0) known++;
    }
    if (known === 0) return;
    const set = lc.set;
    const unit = p.from.kind === K.unit ? p.from.id : -1;
    const mol = set.molecules[unit >= 0 ? set.units.molecule[unit] : p.from.id];
    if (!mol) return;
    let held = 0;
    if (mol.lazy) {
      const data = mol.custom.props.get(StreamProp) as StreamedData | undefined;
      const s = data?.systems.get(set.strings.get(mol.key)) ?? -1;
      held = data && s >= 0 ? (lc.letters?.(data, s, at, 3, letters) ?? 0) : 0;
    } else
      for (let k = 0; k < 3; k++) {
        letters[k] = at[k] >= 0 ? mol.code[at[k]] : 0;
        if (letters[k] > 0) held++;
      }
    // The strand is the codon step's (a codon across a cyclic origin reads upward though its positions wrap).
    const strand = set.mappings[p.codon]?.codon?.strand ?? (at[2] >= 0 && at[2] < at[0] ? -1 : 1);
    const own = unit >= 0 ? set.units.numberings[unit].items[0] : undefined;
    const n = own ?? mol.numberings.items[0] ?? null;
    const pos = (x: number) => (n && x < n.length ? lc.format(n, x as never) : `${x + 1}`);
    const L = mol.length;
    /** Positions next to each other in reading order (across a cyclic origin too); else a gap (an intron). */
    const next = (a: number, b: number) => Math.abs(a - b) === 1 || (mol.cyclic && Math.abs(a - b) === L - 1);
    const name = strand < 0 ? (mol.alphabet === 'rna' ? 'complement_rna' : 'complement') : null;
    const entry = name ? lc.registry.svs.get('letterMaps', name) : undefined;
    let table = entry ? tables.get(entry) : null;
    if (entry && table === undefined) {
      table = letterMapTable(entry.map);
      tables.set(entry, table);
    }
    // The bases (an unknown one, past a clipped domain or not held: `·`) and the parts' positions, in reading order.
    let codon = '';
    let where = '';
    let first = -1;
    let last = -1;
    const part = () => `${pos(first)}${last !== first ? `–${pos(last)}` : ''}`;
    for (let k = 0; k < 3; k++) {
      const x = at[k];
      if (x >= 0 && last >= 0 && !next(x, last)) {
        codon += '|';
        where += `${part()} | `;
        first = -1;
      }
      const c = x < 0 || letters[k] === 0 ? 0 : table ? table[letters[k] & 127] : letters[k];
      codon += c > 0 ? String.fromCharCode(c).toUpperCase() : '·';
      if (x < 0) continue;
      if (first < 0) first = x;
      last = x;
    }
    where += part();
    const system = set.strings.get(unit >= 0 ? set.units.label[unit] : mol.label) || set.strings.get(mol.key);
    if (held === -1) out.line('codon', "failed: Retry in a track's … menu");
    else if (held < 0) out.line('codon', 'not served');
    else if (held < known) {
      out.line('codon', 'loading…');
      lc.waiting();
    } else out.line('codon', codon);
    out.line('codon at', `${system} ${where} (${strand < 0 ? '−' : '+'})`);
  },
};
