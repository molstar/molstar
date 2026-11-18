/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The emphasis passes (_spec/README.md "Render", Emphasis): the selection and the hover carried onto linked sections and
// drawn in the same look, lighter; the highlight channels as tints; the selection's tint over the data or under it (a
// config item); and an export projecting the selection it draws as the viewer does.

import { describe, expect, it } from 'vitest';
import { type Location, rangesLocation } from '@molstar/seqx-core/query';
import { projectOntoView } from '@molstar/seqx-core/view';
import { type DrawLists, type Frame, type PassName, recordDrawLists, renderSvg } from '../index.js';
import { p53GeneSpec, p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 240 };
const Protein = '/data/p53';
const Gene = '/data/tp53';

/** The rects of a pass as `[y, fill, stroke]`, from the recording. */
function rects(lists: DrawLists, pass: PassName): { y: number; line: string }[] {
  return recordDrawLists(lists, [pass])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((line) => ({ y: Number(line.split(' ')[2]), line }));
}

/** The rects of a pass within a section's band. */
function inSection(f: Frame, pass: PassName, section: number) {
  const s = f.layout.sections[section];
  return rects(f.lists, pass).filter((r) => r.y >= s.y && r.y < s.y + s.height);
}

describe('projectOntoView', () => {
  it('carries a residue onto its gene’s section (its codon), and keeps or drops the own part', async () => {
    const s = await settle(p53GeneSpec());
    const residue = rangesLocation(Protein, [174, 175]);
    const projected = projectOntoView(s.model, s.view, residue);
    expect(projected.parts.map((p) => p.system)).toEqual([Gene]);
    const codon = projected.parts[0].ranges;
    let n = 0;
    for (let i = 0; i < codon.length; i += 2) n += codon[i + 1] - codon[i];
    expect(n).toBe(3);
    // a hover's: closed and with its own part
    const own = projectOntoView(s.model, s.view, residue, { own: true, close: true });
    expect(own.parts.map((p) => p.system)).toEqual([Protein, Gene]);
    // not linked: the own part alone
    expect(projectOntoView(s.model, s.view, residue, { own: true, linked: false }).parts.map((p) => p.system)).toEqual([
      Protein,
    ]);
  });

  it('a location on a system no section shows reaches nothing without a path; an empty one is empty', async () => {
    const s = await settle(p53TracksSpec());
    expect(projectOntoView(s.model, s.view, rangesLocation('/data/none', [0, 1])).parts).toEqual([]);
    expect(projectOntoView(s.model, s.view, rangesLocation(Protein, [0, 1])).parts).toEqual([]);
  });
});

describe('the emphasis passes on linked sections', () => {
  async function geneFrame(emphasis: (s: Awaited<ReturnType<typeof settle>>) => object, config = {}) {
    const s = await settle(p53GeneSpec(), [], config);
    const input = { size: Size, viewport: { start: 170, end: 180 } };
    await s.frame(input);
    return { s, f: s.renderer.frame({ ...input, emphasis: emphasis(s) }) };
  }

  it('the selection projected onto the gene is the same look, lighter (half the opacity)', async () => {
    const { f } = await geneFrame((s) => {
      const selection = rangesLocation(Protein, [174, 176]);
      return { selection, projected: { selection: projectOntoView(s.model, s.view, selection) } };
    });
    expect(f.rebuilt).toBe(false);
    const protein = inSection(f, 'selection', 0);
    const gene = inSection(f, 'selection', 1);
    expect(protein.length).toBeGreaterThan(0);
    expect(gene.length).toBeGreaterThan(0);
    // the accent at 0.18 on its own section, at 0.09 projected; its outline opaque, then at half
    expect(protein[0].line).toContain('#8f4c062e');
    expect(protein.some((r) => r.line.includes('stroke #8f4c06ff'))).toBe(true);
    expect(gene[0].line).toContain('#8f4c0617');
    expect(gene.some((r) => r.line.includes('stroke #8f4c0680'))).toBe(true);
  });

  it('the hover projected onto the gene is drawn there, never on the hover’s own section', async () => {
    const { f } = await geneFrame((s) => ({
      hover: {
        layer: '/views/main/protein/letters/letters',
        section: 'protein',
        system: Protein,
        position: 175,
      },
      projected: {
        hover: projectOntoView(s.model, s.view, rangesLocation(Protein, [175, 176]), { own: true, close: true }),
      },
    }));
    const gene = inSection(f, 'hover', 1);
    expect(gene.length).toBeGreaterThan(0);
    // the hover's ink at half its opacity
    expect(gene.every((r) => r.line.includes('stroke #2d271e80'))).toBe(true);
    // the protein's hover is its layer's own outline alone, opaque
    const protein = inSection(f, 'hover', 0);
    expect(protein.length).toBe(1);
    expect(protein[0].line).toContain('stroke #2d271eff');
  });

  it('a hover from elsewhere (no layer: a structure’s) is projected onto every section', async () => {
    const { f } = await geneFrame((s) => ({
      projected: {
        hover: projectOntoView(s.model, s.view, rangesLocation(Protein, [175, 176]), { own: true, close: true }),
      },
    }));
    expect(inSection(f, 'hover', 0).length).toBeGreaterThan(0);
    expect(inSection(f, 'hover', 1).length).toBeGreaterThan(0);
  });

  it('highlight channels are tints of their colors in order, on every system their Location has', async () => {
    const { f } = await geneFrame((s) => {
      const own = rangesLocation(Protein, [172, 174]);
      const both: Location = projectOntoView(s.model, s.view, own, { own: true });
      return {
        highlights: [
          { location: both, color: 0xff000080 },
          { location: rangesLocation(Protein, [177, 178]), color: 0x00ff00ff },
        ],
      };
    });
    const protein = inSection(f, 'highlight', 0);
    // red at 0.5 · 0.3, green at 0.3, two tracks each
    expect(protein.filter((r) => r.line.includes('#ff000026')).length).toBe(2);
    expect(protein.filter((r) => r.line.includes('#00ff004d')).length).toBe(2);
    expect(protein.findIndex((r) => r.line.includes('#00ff004d'))).toBeGreaterThan(
      protein.findIndex((r) => r.line.includes('#ff000026')),
    );
    const gene = inSection(f, 'highlight', 1);
    expect(gene.length).toBeGreaterThan(0);
    expect(gene.every((r) => r.line.includes('#ff000026'))).toBe(true);
  });
});

describe('the selection’s tint: over the data (default) or under it', () => {
  it('under: the tint is in the background pass and the selection pass has the outlines alone', async () => {
    const s = await settle(p53TracksSpec(), [], { emphasis: { selectionPlace: 'under' } });
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    await s.frame(input);
    const a = rangesLocation(Protein, [3, 6]);
    const f = s.renderer.frame({ ...input, emphasis: { selection: a } });
    expect(f.rebuilt).toBe(false);
    // a changed selection draws the base again (its tint is there)
    expect(f.base).toBe(true);
    expect(rects(f.lists, 'background').filter((r) => r.line.includes('#8f4c062e')).length).toBe(3);
    expect(rects(f.lists, 'selection').every((r) => r.line.includes('stroke'))).toBe(true);
    // the same selection again: the emphasis alone
    expect(s.renderer.frame({ ...input, emphasis: { selection: a } }).base).toBe(false);
    // over (the default): the base never redraws for a selection
    const o = await settle(p53TracksSpec());
    await o.frame(input);
    const g = o.renderer.frame({ ...input, emphasis: { selection: a } });
    expect(g.base).toBe(false);
    expect(rects(g.lists, 'background').some((r) => r.line.includes('#8f4c062e'))).toBe(false);
  });
});

describe('an export draws the selection projected, as the viewer does', () => {
  it('a selection on the protein, given as an address, shows on the gene’s section too', async () => {
    const r = await renderSvg(p53GeneSpec(), {
      size: Size,
      viewport: { start: 170, end: 180 },
      selection: '/data/p53/(range :from 175 :to 176)',
    });
    expect(r.selection?.parts[0].ranges).toEqual(new Int32Array([174, 176]));
    const lines = rects(r.lists, 'selection');
    // the gene's tint, lighter
    expect(lines.some((l) => l.line.includes('#8f4c0617'))).toBe(true);
    expect(lines.some((l) => l.line.includes('#8f4c062e'))).toBe(true);
  });
});
