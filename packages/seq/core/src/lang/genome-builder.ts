/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder's genome helpers (_spec/README.md, Language builder, Genomes): an assembly's params from a list of
// sequences, a genome section (its ruler in bases and, given bands, an ideogram navigator: spelled, never inferred) and
// `colinearBlocks` (one sequence's blocks in another's order, from a map's runs read once). Their own module: a writer
// that never writes a genome never loads them.

import {
  type Address,
  ideogramNavigator,
  type Params,
  type RefValue,
  type SectionBuilder,
  type ViewBuilder,
} from './builder.js';

/** A sequence of an assembly as the helpers take it. */
export interface AssemblySequenceRow {
  readonly name: string;
  readonly length: number;
  readonly aliases?: readonly string[];
  readonly cyclic?: boolean;
}

/**
 * An assembly's params from its sequences (a chrom.sizes read by the writer): `name` its accession or name, the
 * sequences inline (columns written only when a row has them), any other param beside.
 */
export function assembly(
  name: string,
  sequences: readonly AssemblySequenceRow[],
  extra: Omit<Params<'assembly'>, 'name' | 'sequences'> = {},
): Params<'assembly'> {
  const table: Record<string, unknown> = { name: sequences.map((s) => s.name), length: sequences.map((s) => s.length) };
  if (sequences.some((s) => s.aliases)) table.aliases = sequences.map((s) => [...(s.aliases ?? [])]);
  if (sequences.some((s) => s.cyclic)) table.cyclic = sequences.map((s) => s.cyclic === true);
  return { name, sequences: table as Params<'assembly'>['sequences'], ...extra };
}

export interface GenomeOptions {
  /** The section's `include` (`['chr17']`, `{ min_length: 1e6 }`). */
  readonly include?: Params<'section'>['include'];
  /** `tracks` (the view's when absent) or `rings`. */
  readonly layout?: 'tracks' | 'rings';
  /** A bands annotation over the assembly: an ideogram navigator over the section (its id `ideogram`). */
  readonly bands?: RefValue<'ranges'>;
}

/**
 * A genome section of one assembly or several (a group each): its blocks the derived sequences, a ruler in bases
 * (track `ruler`), and with `bands` an ideogram navigator; more tracks are written in `body`.
 */
export function genome(
  view: ViewBuilder,
  id: string,
  of: RefValue<'assembly'> | readonly RefValue<'assembly'>[],
  options: GenomeOptions = {},
  body?: (section: SectionBuilder) => void,
): SectionBuilder {
  const params = {
    of: of as Params<'section'>['of'],
    ...(options.include !== undefined ? { include: options.include } : {}),
    ...(options.layout ? { layout: options.layout } : {}),
  };
  const section = view.section(id, params, (s) => {
    s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' });
    body?.(s);
  });
  if (options.bands !== undefined) view.navigator('ideogram', ideogramNavigator(section, options.bands));
  return section;
}

/** A map's runs as `colinearBlocks` reads them: a segments table (`from + i ↦ to + i`, or `to − i` on `-`). */
export interface Runs {
  readonly from: ArrayLike<number>;
  readonly to: ArrayLike<number>;
  readonly length: ArrayLike<number>;
  readonly strand?: ArrayLike<string | null> | { get(i: number): string | null };
}

/**
 * The blocks of `of` in `along`'s order, from a map's runs from `along` onto `of` (_spec/README.md, Language builder,
 * Genomes): a block a run of at least `minLength` (100 kb by default; shorter runs merged into the block before, else
 * after: their hull on `of`), in `along`'s order, `reversed` where its run is. The result is explicit data: written
 * as blocks (`section.block(id, params)`), exported as drawn. Why not a block's `{ along }`: a map read by region would
 * turn a block over as regions arrive; this reads the map once.
 */
export function colinearBlocks(
  runs: Runs,
  of: Address<'sequence'> | string,
  options: { readonly minLength?: number } = {},
): Params<'block'>[] {
  const min = options.minLength ?? 100_000;
  const strandAt = (i: number) => {
    const s = runs.strand;
    if (!s) return '+';
    return ('get' in s && typeof s.get === 'function' ? s.get(i) : (s as ArrayLike<string | null>)[i]) ?? '+';
  };
  const order = Array.from({ length: runs.from.length }, (_, i) => i).sort((a, b) => runs.from[a] - runs.from[b]);
  const blocks: { start: number; end: number; reversed: boolean; long: boolean }[] = [];
  for (const i of order) {
    const n = runs.length[i];
    if (!(n > 0)) continue;
    const reversed = strandAt(i) === '-';
    const start = reversed ? runs.to[i] - n + 1 : runs.to[i];
    const block = { start, end: start + n, reversed, long: n >= min };
    const last = blocks[blocks.length - 1];
    // a short run joins the block before it (a long one starts its own; a short one after a short one joins it too)
    if (last && (!block.long || !last.long)) {
      const keep = last.long || !block.long ? last.reversed : block.reversed;
      last.start = Math.min(last.start, block.start);
      last.end = Math.max(last.end, block.end);
      last.reversed = keep;
      last.long = last.long || block.long;
    } else blocks.push(block);
  }
  return blocks.map((b) => ({
    of,
    range: [b.start, b.end],
    ...(b.reversed ? { direction: 'reversed' as const } : {}),
  }));
}
