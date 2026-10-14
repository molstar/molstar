/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// BigBed: items (README, "BigBed"). A BED row is an item; a BED12 row with blocks is one row per block under one `id`
// (the viewer joins rows of an `id` into one feature: exons of a gene).

import type { StreamProvider } from "@molstar/seq-core/format";
import { openBbi } from "./bbi.ts";
import { spread, zoomRecords } from "./bigwig.ts";
import { ITEM_LEVELS, type ItemsRead, itemsProvider } from "./items.ts";
import { once, type RangeReader, text, view } from "./ranges.ts";

/** The BED columns after chrom, start and end, by position (a file without autoSql, or fields it does not name). */
const BED = [
  "name",
  "score",
  "strand",
  "thick_start",
  "thick_stop",
  "color",
  "block_count",
  "block_sizes",
  "block_starts",
];
const [NAME, BLOCKS] = [0, 6]; // field positions: name, and the three block fields of a BED12
/** The item levels counted from items (regions of 64 kb, cheap to read); coarser ones are the file's zoom levels. */
const FINE = 1;

/**
 * The names of the fields after the first three that an autoSql declaration lists: `string name;  "Name"`, an array
 * `int[blockCount] sizes;`, and an `enum(…)` or `set(…)`, whose values may run over lines (CB-03: skipped, they
 * shifted every later name onto the wrong column).
 */
export const sqlFields = (sql: string | null) =>
  sql === null
    ? []
    : [...sql.slice(sql.indexOf("(")).matchAll(/^\s*(?:\w+(?:\[\w*\])?|(?:enum|set)\([^)]*\))\s+(\w+)\s*;/gm)]
        .map((m) => m[1])
        .slice(3);

/**
 * A BigBed file as an `items` provider: one coordinate system per chromosome; the columns are the BED fields. The host
 * may name the field an item is labeled by (`label`, default the BED name) and the one that is its `category`, with
 * that field's values in order (`vocabulary`: a palette's order; ST4b).
 */
export function bigBedProvider(
  read: RangeReader,
  o: { maxItems?: number; label?: string; category?: string; vocabulary?: string[] } = {},
): StreamProvider {
  const open = once((signal: AbortSignal) => openBbi(read, signal));
  return itemsProvider(
    {
      vocabulary: o.category && o.vocabulary ? { category: o.vocabulary } : {},
      // The items' counts at the one fine level (bins of 256, a region of 64 kb, cheap to read), then the file's zoom
      // levels coarser than them, their summaries never reading an item (ST4b).
      async levels(signal) {
        const bbi = await open(signal);
        return [
          ...ITEM_LEVELS.slice(0, FINE),
          ...bbi.zooms.map((z) => z.reduction).filter((w) => w > ITEM_LEVELS[FINE - 1]),
        ];
      },
      async bins(system, start, stop, level, signal) {
        if (level < FINE) return null;
        const bbi = await open(signal);
        const chrom = bbi.chroms.get(system);
        const w = (await this.levels?.(signal))?.[level] ?? 0;
        const zoom = bbi.zooms.find((z) => z.reduction === w);
        if (!chrom || !zoom) throw new Error(`BigBed: no zoom level of ${w} on ${system}`);
        // Covered bases per bin over its width: the items' mean depth there, a density (a zoom record is not on the
        // level's grid: spread over the bins it meets by share, as BigWig's).
        const n = Math.ceil((stop - start) / w);
        const depth = new Float64Array(n);
        for (const b of await bbi.blocks(zoom.indexOffset, chrom, start, stop, signal))
          zoomRecords(b, bbi.le, chrom, (s0, e0, _valid, _lo, _hi, sum) =>
            spread(s0, e0, start, w, n, (j, share) => {
              depth[j] += (sum * share) / w;
            }),
          );
        return depth.map((d) => (d > 0 ? Math.max(1, Math.round(d)) : 0));
      },
      async systems(signal) {
        const bbi = await open(signal);
        if (bbi.kind !== "bigbed") throw new Error("not a BigBed file (this is a BigWig)");
        return [...bbi.chroms.values()].map((c) => ({ key: c.name, extent: c.size }));
      },
      async read(system, start, stop, signal): Promise<ItemsRead> {
        const bbi = await open(signal);
        const chrom = bbi.chroms.get(system);
        if (!chrom) throw new Error(`BigBed: no chromosome ${system}`);
        const named = sqlFields(bbi.autoSql);
        const names = Array.from(
          { length: Math.max(0, bbi.fieldCount - 3) },
          (_, k) => named[k] ?? BED[k] ?? `field${k + 4}`,
        );
        const blocked = bbi.definedFieldCount >= 12; // BED12: blocks are exons
        const labelAt = o.label ? names.indexOf(o.label) : NAME;
        const categoryAt = o.category ? names.indexOf(o.category) : -1;
        // Keyed by the file's names: a `__proto__` field is a column (null prototype), and a field named like a
        // column made here, or like an earlier field, is named by its position instead (one push per item, LA-39).
        const columns: Record<string, unknown[]> = Object.create(null);
        for (const c of ["start", "stop", "id", "label", ...(categoryAt >= 0 ? ["category"] : [])]) columns[c] = [];
        // The block fields are consumed.
        const kept = names.map((n, f) => {
          if (blocked && f >= BLOCKS && f < BLOCKS + 3) return null;
          const name = n in columns ? `field${f + 4}` : n;
          columns[name] = [];
          return name;
        });
        const [starts, stops]: number[][] = [[], []];
        for (const block of await bbi.blocks(bbi.dataIndexOffset, chrom, start, stop, signal))
          for (let at = 0, dv = view(block); at + 12 <= block.length; ) {
            // A last record with no NUL ends the block (GI-05: `indexOf` gave -1 and the loop read the block for ever).
            const nul = block.indexOf(0, at + 12);
            const end = nul < 0 ? block.length : nul;
            const [id, s, e] = [dv.getUint32(at, bbi.le), dv.getUint32(at + 4, bbi.le), dv.getUint32(at + 8, bbi.le)];
            const rest = text(block.subarray(at + 12, end)).split("\t");
            at = end + 1;
            if (id !== chrom.id || e <= start || s >= stop) continue;
            const key = `${rest[NAME] ?? ""}@${s}-${e}`;
            const sizes = blocked ? (rest[BLOCKS + 1] ?? "").split(",").filter(Boolean).map(Number) : [];
            const offsets = blocked ? (rest[BLOCKS + 2] ?? "").split(",").filter(Boolean).map(Number) : [];
            for (let k = 0; k < Math.max(1, sizes.length); k++) {
              columns.start.push(sizes.length ? s + offsets[k] : s);
              columns.stop.push(sizes.length ? s + offsets[k] + sizes[k] : e);
              columns.id.push(key);
              columns.label.push(rest[labelAt] ?? rest[NAME] ?? null);
              if (categoryAt >= 0) columns.category.push(rest[categoryAt] ?? null);
              kept.forEach((n, f) => {
                if (n) columns[n].push(rest[f]);
              });
            }
            starts.push(s);
            stops.push(e);
          }
        return { table: columns, start: starts, stop: stops };
      },
    },
    o.maxItems,
  );
}
