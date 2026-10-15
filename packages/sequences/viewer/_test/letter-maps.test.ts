/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LM1 (seqspec:format#layer-types `letter_map`; seqspec:format#letter-maps): letters and rows drawn through a layer's
// map. The strongest check is equality: a mapped layer draws exactly what the same layer draws over the mapped letters
// stored (rects, texts, their colors), forward and reversed, linear and wrapped, for `letters`, `rows` and `of`; hits
// and hover
// labels name the letter shown and the stored one; the SVG export has the letters shown; the preset
// `reverseComplement` is a reversed piece with a complement.

import { createSeqViewSpec, translate } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import type { SeqViewerContext } from '../src/index.js';
import { renderSvg } from '../src/svg/index.js';
import { LayerFrame, RecordingBackend } from '../src/testing/index.js';
import { asText } from './conformance-fixtures.js';
import { drawn, entryOf, layerOf, record, rects, texts } from './draw-fixtures.js';
import { patch } from './graph-fixtures.js';
import { lit } from './overlay-fixtures.js';
import { all, parseXml } from './svg/xml.js';

const DNA = 'GAATTCGCCACCATGGCTTCCAAGGAACTGRYKN';
const COMPLEMENT: Record<string, string> = { A: 'T', T: 'A', C: 'G', G: 'C', R: 'Y', Y: 'R', K: 'M', N: 'N' };
const complement = (s: string) => [...s].map((c) => COMPLEMENT[c] ?? c).join('');
const reverse = (s: string) => [...s].reverse().join('');

/** `value` as DNA letters in one track "DNA", with the layer's `params`; the section's `direction`. */
function dna(value: string, params: Record<string, unknown> = {}, direction?: 1 | -1) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', label: 'Gene', alphabet: 'dna', value });
  const s = b
    .view({ ref: 'v', label: 'V' })
    .section({ ref: 'v.s', height: 'box', ...(direction ? { direction } : {}) });
  s.block({ sequence: a });
  s.track({ label: '' }).ruler();
  s.track({ label: 'DNA' }).letters(params as never);
  return b.getDocument();
}
/** A tile image's pixels as one number (FNV-1a). */
const hash = (c: unknown) => {
  let h = 2166136261;
  if (c instanceof Uint32Array) for (let i = 0; i < c.length; i++) h = Math.imul(h ^ c[i], 16777619);
  return h >>> 0;
};
/** What a track drew: its texts left to right and its rects, as numbers (a context's font ids differ). */
async function draw(doc: SeqViewSpecDoc, label = 'DNA', width = 1000) {
  const t = await drawn(doc, { width, height: 200 });
  const { rec } = record(t.ctx, label);
  const shown = texts(rec)
    .map((d) => [d.x, d.y, d.text, d.color] as const)
    .sort((p, q) => p[0] - q[0] || p[1] - q[1]);
  t.ctx.dispose();
  return { shown, rects: rects(rec), images: rec.images.content.map(hash) };
}

/** An alignment of `rows` in one `Rows` track. */
function aln(rows: string[], map: unknown = null) {
  const b = createSeqViewSpec();
  const a = b.data.alignment({ ref: 'aln', value: rows, flavor: 'aligned' } as never);
  const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ alignment: a });
  v.track({ label: 'Rows' }).rows({ letter_map: map } as never);
  return b.getDocument();
}

describe('letters through a map', () => {
  it("draws what the mapped letters draw when stored: the characters, their ink and their cells' colors", async () => {
    const mapped = await draw(dna(DNA, { letter_map: 'complement' }));
    expect(mapped.shown.map((p) => p[2]).join('')).toBe(complement(DNA));
    expect(mapped.rects.length).toBeGreaterThan(20);
    expect(mapped).toEqual(await draw(dna(complement(DNA))));
    // Not the stored letters' draw: the cells are colored by the letter shown (A is green, T blue).
    expect(mapped.rects).not.toEqual((await draw(dna(DNA))).rects);
  });

  it('an inline dict maps single letters, the others kept; a space hides a letter', async () => {
    const mapped = await draw(dna('GATTACA', { letter_map: { A: 'T', T: ' ' } }));
    expect(mapped.shown.map((p) => p[2]).join('')).toBe('GT  TCT'.replace(/ /g, ' '));
    expect((await draw(dna('GATTACA', { letter_map: null }))).shown.map((p) => p[2]).join('')).toBe('GATTACA');
  });

  it('a reversed piece with the complement draws the reverse complement stored', async () => {
    const mapped = await draw(dna(DNA, { letter_map: 'complement' }, -1));
    expect(mapped.shown.map((p) => p[2]).join('')).toBe(reverse(complement(DNA)));
    expect(mapped).toEqual(await draw(dna(reverse(complement(DNA)))));
  });

  it('in the text layout (wrapped), and zoomed out into cells', async () => {
    const mapped = await draw(asText(dna(DNA, { letter_map: 'complement' })), 'DNA', 400);
    const stored = await draw(asText(dna(complement(DNA))), 'DNA', 400);
    expect(mapped.shown.map((p) => p[2]).join('')).toContain(complement(DNA).slice(0, 10));
    expect(mapped).toEqual(stored);
    // 4,000 positions in 300 px: cells by their classes, from the mapped colors.
    const long = 'ACGTTGCAGGCT'.repeat(340);
    const a = await draw(dna(long, { letter_map: 'complement' }), 'DNA', 300);
    expect(a.rects.length).toBeGreaterThan(10);
    expect(a.rects).toEqual((await draw(dna(complement(long)), 'DNA', 300)).rects);
  });
});

describe('the preset', () => {
  it('reverseComplement draws the reverse complement of the section, its protein track untouched', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', label: 'Gene', alphabet: 'dna', value: 'ATGGCTTCCAAGGAACTGTAA' });
    const p = translate(b.data, { ref: 'p', label: 'Protein', from: a, cds: [[0, 21]] });
    const s = b.view({ ref: 'v', label: 'V', height: 'box' }).section({ ref: 'v.s' }).block({ sequence: a });
    s.track({ label: 'DNA' }).letters();
    s.track({ label: 'Protein' }).letters({ of: p });
    b.preset.reverseComplement(s);
    const t = await drawn(b.getDocument(), { width: 1000, height: 300 });
    const dnaLetters = texts(record(t.ctx, 'DNA').rec).sort((x, y) => x.x - y.x);
    expect(dnaLetters.map((d) => d.text).join('')).toBe('TTACAGTTCCTTGGAAGCCAT');
    // The protein's letters keep their order along the (reversed) axis: M at the right.
    const protein = texts(record(t.ctx, 'Protein').rec).sort((x, y) => x.x - y.x);
    expect(protein.at(-1)?.text.startsWith('M') || protein.at(-1)?.text === 'Met').toBe(true);
    t.ctx.dispose();
  });

  it("leaves a protein's own letters as stored (the map lists its alphabets), called before the tracks too", async () => {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: 'g', label: 'Gene', alphabet: 'dna', value: 'ATGGCTTCCAAG' });
    const p = b.data.sequence({ ref: 'p', label: 'Protein', value: 'MKTAYIAKQR' });
    const v = b.view({ ref: 'v', label: 'V', height: 'box' });
    const dnaSection = b.preset.reverseComplement(v.section({ ref: 'v.g' }).block({ sequence: g }));
    dnaSection.track({ label: 'DNA' }).letters();
    const proteinSection = b.preset.reverseComplement(v.section({ ref: 'v.p' }).block({ sequence: p }));
    proteinSection.track({ label: 'Protein' }).letters();
    const t = await drawn(b.getDocument(), { width: 1000, height: 300 });
    const read = (label: string) =>
      texts(record(t.ctx, label).rec)
        .sort((x, y) => x.x - y.x)
        .map((d) => d.text)
        .join('');
    expect(read('DNA')).toBe('CTTGGAAGCCAT');
    expect(read('Protein')).toBe(reverse('MKTAYIAKQR'));
    t.ctx.dispose();
  });
});

describe("another polymer's letters (`of`)", () => {
  function gene(map: unknown) {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: 'ATGGCTTCCAAGCACTGA' });
    const p = translate(b.data, { ref: 'p', from: g, cds: [[0, 18]] });
    const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: g });
    v.track({ label: 'Protein' }).letters({ of: p, letter_map: map } as never);
    return b.getDocument();
  }
  it('a residue the map changes is drawn as the one letter shown; the others keep their form', async () => {
    const wide = async (map: unknown) => {
      const t = await drawn(gene(map), { width: 1200, height: 200 });
      const out = texts(record(t.ctx, 'Protein').rec)
        .sort((x, y) => x.x - y.x)
        .map((d) => d.text);
      t.ctx.dispose();
      return out;
    };
    expect(await wide(null)).toEqual(['Met', 'Ala', 'Ser', 'Lys', 'His']);
    expect(await wide({ M: 'm', H: 'X' })).toEqual(['m', 'Ala', 'Ser', 'Lys', 'X']);
  });

  it('cells are colored by the residue shown: a mapped K draws the cells of a D stored', async () => {
    const cells = async (value: string, map: unknown) => {
      const b = createSeqViewSpec();
      const g = b.data.sequence({ ref: 'gene', alphabet: 'dna', value });
      const p = translate(b.data, { ref: 'p', from: g, cds: [[0, 18]] });
      const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: g });
      v.track({ label: 'Protein' }).letters({ of: p, letter_map: map } as never);
      return (await draw(b.getDocument(), 'Protein', 1200)).rects;
    };
    // AAG (Lys) mapped to D, against GAC (Asp) stored.
    const mapped = await cells('ATGGCTTCCAAGCACTGA', { K: 'D' });
    expect(mapped).toEqual(await cells('ATGGCTTCCGACCACTGA', null));
    expect(mapped).not.toEqual(await cells('ATGGCTTCCAAGCACTGA', null));
  });
});

describe('emphasis over mapped letters', () => {
  it("emphasizes the letters its track shows, through the track's letters layer's map", async () => {
    const doc = (map: string | null, value = DNA) => {
      const b = createSeqViewSpec();
      const a = b.data.sequence({ ref: 'a', label: 'Gene', alphabet: 'dna', value });
      const f = b.data.annotation({ ref: 'f', geometry: 'ranges', target: a, value: { start: [2], stop: [8] } });
      const s = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.s', height: 'box' }).block({ sequence: a });
      s.track({ label: 'DNA' }).letters({ ref: 'v.dna', letter_map: map }).emphasis({ data: f });
      return b.getDocument();
    };
    /** The emphasis layer's draw (the track's second layer): its letters left to right, and its rects. */
    const emphasized = (ctx: SeqViewerContext) => {
      const { k, section } = layerOf(ctx, 'DNA');
      const rec = new RecordingBackend();
      const entry = ctx.registry.get('renderers', 'emphasis');
      if (entry) new LayerFrame(ctx, entry, section, k + 1).draw(rec);
      return {
        shown: texts(rec)
          .sort((p, q) => p.x - q.x)
          .map((d) => [d.x, d.text, d.color]),
        rects: rects(rec),
      };
    };
    const once = async (d: SeqViewSpecDoc) => {
      const t = await drawn(d, { width: 1000, height: 200 });
      const out = emphasized(t.ctx);
      t.ctx.dispose();
      return out;
    };
    const t = await drawn(doc('complement'), { width: 1000, height: 200 });
    const mapped = emphasized(t.ctx);
    expect(mapped.shown.map((p) => p[1]).join('')).toBe(complement(DNA.slice(2, 8)));
    expect(mapped).toEqual(await once(doc(null, complement(DNA))));
    // A patch to the letters layer's map redraws the emphasis with it.
    for (const [map, letters] of [
      ['lower', DNA.slice(2, 8).toLowerCase()],
      [null, DNA.slice(2, 8)],
    ] as const) {
      expect((await t.ctx.applyPatch(patch([{ op: 'set', ref: 'v.dna', params: { letter_map: map } }]))).ok).toBe(true);
      expect(
        emphasized(t.ctx)
          .shown.map((p) => p[1])
          .join(''),
      ).toBe(letters);
    }
    t.ctx.dispose();
  });
});

describe('rows through a map', () => {
  it('rows draw what the mapped rows draw when stored: letters, gaps and cells', async () => {
    const rows = ['ACGT-ACGTN', 'ACGA-AC-TT', 'TCGTTAC-TA'];
    const mapped = await draw(aln(rows, 'complement'), 'Rows');
    expect(mapped.shown.length).toBeGreaterThan(10);
    expect(mapped).toEqual(await draw(aln(rows.map(complement)), 'Rows'));
    expect(mapped.images).not.toEqual((await draw(aln(rows), 'Rows')).images);
  });
});

describe('hover labels and the SVG export', () => {
  it('the label names the letter shown, and the stored one where it differs', async () => {
    const t = await lit(dna(DNA, { letter_map: 'complement' }), { width: 1000, height: 200 });
    const e = entryOf(t.ctx, 'DNA');
    const ax = t.ctx.layout.axis(t.section);
    const at = (i: number) => {
      const k = ax?.shownOf(0, i) ?? 0;
      return ax ? (ax.x(0, k as never) + ax.x(0, (k + 1) as never)) / 2 : 0;
    };
    t.hover(at(0), e.top + e.height / 2);
    const label = t.ctx.labels.compute(t.ctx.hover.live.current);
    expect(label?.title).toBe('Gene C1');
    expect(label?.lines.find((l) => l.key === 'stored')?.value).toBe('G');
    // N maps to N: nothing to say about the stored letter.
    t.hover(at(DNA.indexOf('N')), e.top + e.height / 2);
    const same = t.ctx.labels.compute(t.ctx.hover.live.current);
    expect([same?.title, same?.lines.some((l) => l.key === 'stored')]).toEqual([
      `Gene N${DNA.indexOf('N') + 1}`,
      false,
    ]);
    t.ctx.dispose();
  });

  it("a row's residue label says the same; a layer without a map shows the stored letter", async () => {
    for (const map of ['complement', null]) {
      const t = await lit(aln(['ACGT-ACGTN', 'ACGA-AC-TT'], map), { width: 1000, height: 300 });
      const e = entryOf(t.ctx, 'Rows');
      const ax = t.ctx.layout.axis(t.section);
      const k = ax?.shownOf(0, 0) ?? 0;
      t.hover(ax ? (ax.x(0, k as never) + ax.x(0, (k + 1) as never)) / 2 : 0, e.top + 6);
      const label = t.ctx.labels.compute(t.ctx.hover.live.current);
      expect(label?.lines.find((l) => l.key === 'residue')?.value).toBe(map ? 'T1' : 'A1');
      expect(label?.lines.find((l) => l.key === 'stored')?.value).toBe(map ? 'A' : undefined);
      t.ctx.dispose();
    }
  });

  it("hover reads the hovered layer's map: two strands in one section; a features track, its first letters layer", async () => {
    const doc = (tracks: readonly (string | null)[]) => {
      const b = createSeqViewSpec();
      const a = b.data.sequence({ ref: 'a', label: 'Gene', alphabet: 'dna', value: DNA });
      const f = b.data.annotation({ ref: 'f', geometry: 'ranges', target: a, value: { start: [20], stop: [24] } });
      const s = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.s', height: 'box' }).block({ sequence: a });
      tracks.forEach((map, k) => {
        s.track({ label: `L${k}` }).letters({ letter_map: map });
      });
      s.track({ label: 'Features' }).ranges({ data: f });
      return b.getDocument();
    };
    /** The title and `stored` line hovering position 0 of each track. */
    const hovered = async (tracks: readonly (string | null)[]) => {
      const t = await lit(doc(tracks), { width: 1000, height: 300 });
      const ax = t.ctx.layout.axis(t.section);
      const k = ax?.shownOf(0, 0) ?? 0;
      const x = ax ? (ax.x(0, k as never) + ax.x(0, (k + 1) as never)) / 2 : 0;
      const out = [...tracks.map((_, j) => `L${j}`), 'Features'].map((label) => {
        const e = entryOf(t.ctx, label);
        t.hover(x, e.top + e.height / 2);
        const l = t.ctx.labels.compute(t.ctx.hover.live.current);
        return `${l?.title} ${l?.lines.find((x) => x.key === 'stored')?.value ?? '-'}`;
      });
      t.ctx.dispose();
      return out;
    };
    // Both strands: each track names its own letter; the features track reads the first letters layer's (as stored).
    expect(await hovered([null, 'complement'])).toEqual(['Gene G1 -', 'Gene C1 G', 'Gene G1 -']);
    expect(await hovered(['complement', null])).toEqual(['Gene C1 G', 'Gene G1 -', 'Gene C1 G']);
  });

  it("a position an `of` layer shows no residue at reads the block's letters, not the `of` layer's map", async () => {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: 'g', label: 'Gene', alphabet: 'dna', value: 'GAATTCATGGCTTCCAAGTAA' });
    const p = translate(b.data, { ref: 'p', label: 'Protein', from: g, cds: [[6, 21]] });
    const s = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.s', height: 'box' }).block({ sequence: g });
    s.track({ label: 'DNA' }).letters();
    s.track({ label: 'Protein' }).letters({ of: p, letter_map: 'lower' });
    const t = await lit(b.getDocument(), { width: 1000, height: 300 });
    const ax = t.ctx.layout.axis(t.section);
    const x = (i: number) => {
      const k = ax?.shownOf(0, i) ?? 0;
      return ax ? (ax.x(0, k as never) + ax.x(0, (k + 1) as never)) / 2 : 0;
    };
    const e = entryOf(t.ctx, 'Protein');
    const at = (i: number) => {
      t.hover(x(i), e.top + e.height / 2);
      const l = t.ctx.labels.compute(t.ctx.hover.live.current);
      return `${l?.title} ${l?.lines.find((x) => x.key === 'stored')?.value ?? '-'}`;
    };
    // The UTR: the gene's own A, unmapped; the first codon: the protein's M through `lower`.
    expect(at(1)).toBe('Gene A2 -');
    expect(at(7)).toMatch(/^Protein m1 M$/);
    t.ctx.dispose();
  });

  it('the SVG export has the letters shown (reversed too)', async () => {
    for (const direction of [1, -1] as const) {
      const xml = await renderSvg(dna(DNA, { letter_map: 'complement' }, direction), { width: 1200 });
      const letters = all(parseXml(xml), 'text')
        .map((e) => [Number(e.attrs.x), e.text] as const)
        .sort((p, q) => p[0] - q[0])
        .map((p) => p[1])
        .filter((c) => /^[A-Za-z]$/.test(c))
        .join('');
      expect(letters).toBe(direction > 0 ? complement(DNA) : reverse(complement(DNA)));
    }
  });
});
