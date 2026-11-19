/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The emphasis passes (_spec/README.md "Render", Emphasis): the selection and the hover carried onto linked sections and
// drawn in the same look, lighter; the highlight channels as tints; the selection's tint over the data or under it (a
// config item); and an export projecting the selection it draws as the viewer does.

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { elements, type Location, rangesLocation } from '@molstar/seqx-core/query';
import { projectOntoView } from '@molstar/seqx-core/view';
import { type DrawLists, type Frame, type PassName, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
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

/** The rects of a pass within a section's band (an outline outside a mark starts a few px above it). */
function inSection(f: Frame, pass: PassName, section: number) {
  const s = f.layout.sections[section];
  return rects(f.lists, pass).filter((r) => r.y >= s.y - 4 && r.y < s.y + s.height);
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
    // the accent at 0.18 on its own section, at 0.09 projected; its thin outline at 0.6, then at half that
    expect(protein[0].line).toContain('#8f4c062e');
    expect(protein.some((r) => r.line.includes('stroke #8f4c0699 1'))).toBe(true);
    expect(gene[0].line).toContain('#8f4c0617');
    expect(gene.some((r) => r.line.includes('stroke #8f4c064d 1'))).toBe(true);
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
        hoverSection: '/views/main/protein',
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
    // red at 0.5 · 0.3, green at 0.3, one band each over the section's two tracks
    expect(protein.filter((r) => r.line.includes('#ff000026')).length).toBe(1);
    expect(protein.filter((r) => r.line.includes('#00ff004d')).length).toBe(1);
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
    expect(rects(f.lists, 'background').filter((r) => r.line.includes('#8f4c062e')).length).toBe(1);
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

describe('the order of the tints: highlights below the selection', () => {
  const sel = rangesLocation(Protein, [3, 6]);
  const highlights = [{ location: rangesLocation(Protein, [4, 8]), color: 0x00ff00ff }];

  it('under: the highlights’ tints in the background pass, before the selection’s; a changed channel redraws the base', async () => {
    const s = await settle(p53TracksSpec(), [], { emphasis: { selectionPlace: 'under' } });
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    await s.frame(input);
    const f = s.renderer.frame({ ...input, emphasis: { selection: sel, highlights } });
    const bg = rects(f.lists, 'background').map((r) => r.line);
    const h = bg.findIndex((l) => l.includes('#00ff004d'));
    const t = bg.findIndex((l) => l.includes('#8f4c062e'));
    expect(h).toBeGreaterThan(0);
    expect(t).toBeGreaterThan(h);
    expect(rects(f.lists, 'highlight')).toEqual([]);
    const again = s.renderer.frame({
      ...input,
      emphasis: { selection: sel, highlights: [{ ...highlights[0], color: 0xff0000ff }] },
    });
    expect(again.base).toBe(true);
  });

  it('over: the highlights over the data, the selection’s pass above them', async () => {
    const s = await settle(p53TracksSpec());
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    await s.frame(input);
    const f = s.renderer.frame({ ...input, emphasis: { selection: sel, highlights } });
    expect(rects(f.lists, 'highlight').length).toBeGreaterThan(0);
    expect(rects(f.lists, 'background').some((r) => r.line.includes('#00ff004d'))).toBe(false);
    expect(f.base).toBe(false);
  });
});

describe('the composition of hover and selection', () => {
  it('a cell’s hover is concentric and outside it by the offset: it never touches the selection’s outline', async () => {
    const s = await settle(p53TracksSpec());
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    await s.frame(input);
    const f = s.renderer.frame({
      ...input,
      emphasis: {
        selection: rangesLocation(Protein, [3, 6]),
        hover: { layer: '/views/main/protein/letters/letters', section: 'protein', system: Protein, position: 4 },
      },
    });
    // the cell of 4 is [189, 212) over the letters' [8, 28); 1.5 px wide, 1 px out: its middle 1.75 px outside; a
    // squircle, the cell's radius (a third of 20) plus the 1.75 (concentric)
    expect(rects(f.lists, 'hover').map((r) => r.line)).toEqual([
      'rect 187.25 6.25 26.5 23.5 none stroke #2d271eff 1.5 squircle 8.42',
    ]);
    // the selection's outline is inside its band, half its width in (the band starts on the tracks' top, 8)
    expect(rects(f.lists, 'selection').some((r) => / 8\.5 [\d.]+ [\d.]+ none stroke /.test(r.line))).toBe(true);
  });

  it('a selected item of zero length is a hairline at its place (EC-127)', async () => {
    const b = createSeqSpec({ title: 'boundary' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
    const sites = b.data.ranges('sites', {
      target: p,
      rows: { start: new Int32Array([4]), end: new Int32Array([4]), id: ['ins'] },
    });
    b.views.view('main', (v) => v.section('s', { of: p }).track('sites').layer('r', { type: 'ranges', data: sites }));
    const s = await settle(b.build());
    const input = { size: Size };
    await s.frame(input);
    const f = s.renderer.frame({ ...input, emphasis: { elements: elements([{ owner: '/data/sites', key: 'ins' }]) } });
    const lines = rects(f.lists, 'selection').map((r) => r.line);
    expect(lines.length).toBe(1);
    expect(lines[0]).toMatch(/^rect [\d.]+ [\d.]+ 2 14 #8f4c06ff$/);
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
