/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import * as plasmidGff from '../files/plasmid-gff.js';
import * as genbank from '../genomes/genbank.js';
import { nodeFixtureAssets } from '../node-assets.js';

// The plasmid from its GenBank record beside the GFF3 one (Features/Genome formats/GenBank): one file holds the letters
// and the features, the other two files, and the two views draw the same (their draw lists equal, the letters but their
// case: GenBank's ORIGIN is lowercase, as read), bla and the primer site one item across the origin in both.

function context() {
  return createSeqContext({
    registry: [
      { kind: 'format-runtime', runtime: SeqRuntime },
      { kind: 'format-runtime', runtime: RecordsRuntime },
      ...genbank.registry({ delayMs: 0 }),
    ],
    ports: { assets: () => nodeFixtureAssets() },
  });
}

const size = { width: 720, height: 900 };

describe('the plasmid from GenBank and from GFF3', () => {
  it('the two views’ draw lists are equal: the ring, its tracks by key, the detail, the labels', async () => {
    const ctx = context();
    const gff = await ctx.renderSvg(plasmidGff.build(), { size });
    const gb = await ctx.renderSvg(genbank.plasmid(), { size });
    expect(gff.issues).toEqual([]);
    expect(gb.issues.map((i) => i.code).sort()).toEqual([
      'seq.source.genbank.header',
      'seq.source.genbank.translation',
    ]);
    expect(gb.settled && gff.settled).toBe(true);
    const a = recordDrawLists(gff.lists);
    const b = recordDrawLists(gb.lists);
    expect(a.split('\n').length).toBeGreaterThan(100);
    expect(b).toBe(a);
    const count = (s: string) => gb.svg.split(`>${s}<`).length - 1;
    expect(count('bla (AmpR)')).toBe(2);
    ctx.retire();
  });

  it('zoomed to letters in the detail: the same marks, GenBank’s letters as read (lowercase)', async () => {
    const ctx = context();
    const viewport = { start: 2650, end: 2700 };
    const gff = await ctx.renderSvg(plasmidGff.build(), { size, viewport });
    const gb = await ctx.renderSvg(genbank.plasmid(), { size, viewport });
    const a = recordDrawLists(gff.lists).split('\n');
    const b = recordDrawLists(gb.lists).split('\n');
    expect(b.length).toBe(a.length);
    // every line equal but the letters' glyphs, which differ only by case (their colours and places the same)
    const differ = a.flatMap((l, k) => (l === b[k] ? [] : [[l, b[k]]]));
    expect(differ.length).toBe(1);
    expect(differ[0][0].startsWith('glyphs')).toBe(true);
    expect(differ[0][1].toUpperCase()).toBe(differ[0][0].toUpperCase());
    expect(gb.svg).toMatch(/>[acgt]</);
    ctx.retire();
  });

  it('a different file draws differently (the comparison sees a feature moved)', async () => {
    const ctx = context();
    const gff = await ctx.renderSvg(plasmidGff.build(), { size });
    const moved = await ctx.renderSvg(genbank.misread(), { size });
    expect(recordDrawLists(moved.lists)).not.toBe(recordDrawLists(gff.lists));
    ctx.retire();
  });
});
