/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The items the hover contract closes (_spec/README.md "Render", Emphasis, The hover contract), each on its story's
// view as the story composes it, headless on a fake clock: the density row's highlight (Files › UniProt's entry), a
// ligand's score box lighting its arcs (Complex › As a complex), the strands' density bands on a ring and a strip
// (Plasmid › Ring at scale); and a legend's hover both ways (Genome formats › VCF › TP53 ClinVar).

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { type Frame, type RangesLay, walkPass } from '@molstar/seqx-core/render';
import { ComplexViews, ExampleViews, FilesViews, PlasmidViews } from '@molstar/seqx-examples';
import { LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import { describerFor } from '@molstar/seqx-viewer/describe';
import { legendRows } from '@molstar/seqx-viewer/legends';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { HoverOwner, openExample } from './draw-checks.js';
import { passBoxes, probesOf } from './hover-probes.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const open = (view: Parameters<typeof openExample>[0], options: Partial<Parameters<typeof openExample>[1]> = {}) =>
  openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000, awaitOpens: true, ...options });

const card = (v: SeqViewer) => describerFor({ ...v, layout: () => v.frame?.layout }).describe(v.hover.current);

/** The hover pass's paths (a pairs layer's lit arcs). */
function hoverPaths(f: Frame): number {
  let n = 0;
  walkPass(f.lists.passes.hover, {
    rect() {},
    path() {
      n++;
    },
    glyphs() {},
    text() {},
    image() {},
    clip() {},
    unclip() {},
  });
  return n;
}

/** The pointer on a probe: the frame drawn, the hover's own section's mark and its card. */
function hoverAt(v: SeqViewer, x: number, y: number) {
  v.input.move(x, y);
  const f = v.draw() ?? v.frame!;
  // the hover's mark over the pointer's point (its outline, the view's position box), within 3 px
  const marked = passBoxes(f.lists.passes.hover).some(
    (b) => b.x1 >= x - 3 && b.x0 <= x + 3 && b.y1 >= y - 3 && b.y0 <= y + 3,
  );
  return { f, marked, card: card(v) };
}

describe('a variant allele’s card', () => {
  it('a VCF allele’s hover card naming the record, not the allele (Genome formats › VCF › TP53 ClinVar): `pos:ref>alt`, its category and its INFO fields', {
    tags: ['medium'],
  }, async () => {
    const v = (viewer = await open(ExampleViews.find((x) => x.id === 'vcf-clinvar')!));
    const probe = probesOf(v.frame!).find((p) => p.type === 'ranges' && v.hit(p.x, p.y)?.element !== undefined)!;
    v.input.move(probe.x, probe.y);
    v.draw();
    const key = String(v.hover.current.element!.key);
    expect(key).toMatch(/^\d+:[ACGTN]+>[ACGTN]+$/);
    const c = card(v)!;
    expect(c.title).toBe(key);
    const rows = Object.fromEntries(c.rows.map((r) => [r.label, r.value]));
    expect(rows.Category).toMatch(/^(SNV|MNV|insertion|deletion|indel)$/);
    // the record's ClinVar id is a row, no longer the title
    expect(rows.ID).toMatch(/^\d+$/);
    // its INFO fields: ClinVar's significance among them
    const info = c.rows.filter((r) => !['Track', 'Extent', 'Category', 'ID', 'Strand'].includes(r.label));
    expect(info.length).toBeGreaterThan(0);
    expect(c.rows.map((r) => r.label)).toContain('CLNSIG');
  });
});

describe('legends both ways', () => {
  it('a legend’s hover working one way only (Genome formats › VCF › TP53 ClinVar): the “deletion” entry lights its variants, as a variant outlines it', {
    tags: ['medium'],
  }, async () => {
    const v = (viewer = await open(ExampleViews.find((x) => x.id === 'vcf-clinvar')!));
    const row = legendRows(v).find((r) => r.entries.some((e) => e.value === 'deletion'))!;
    const deletion = row.entries.find((e) => e.value === 'deletion')!;
    // a variant of it hovered outlines the entry (the way there)
    const variants = v
      .frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers))
      .find((l) => l.key === row.layer)!;
    const lay = variants.lay as RangesLay;
    // a deletion's mark the pointer reaches (variants of one site share a pixel: the one hit there)
    let found = false;
    for (let k = 0; k < lay.hn && !found; k++) {
      if (lay.painted[lay.hitem[k]] !== 'deletion') continue;
      v.input.move((lay.hx0[k] + lay.hx1[k]) / 2, lay.laneY + (lay.hlane[k] + 0.5) * lay.laneHeight);
      v.draw();
      found = v.legends.hover.value?.value === 'deletion';
      v.input.leave();
    }
    expect(found).toBe(true);
    // the way back: the entry is a target, its hover lights its variants in view
    expect(deletion.action).toBe(true);
    await v.execute('legend.hover', { address: row.layer, entry: 'deletion' });
    const f = v.draw() ?? v.frame!;
    expect(passBoxes(f.lists.passes.highlight).length).toBeGreaterThan(0);
    expect(v.highlight.get('legend')?.elements.keys.length).toBeGreaterThan(0);
  });
});

describe('the hover contract closes', () => {
  it('the density row’s missing highlight (Files › UniProt’s entry, Regions): the row lights the stretch hovered', async () => {
    // a track capped low with a density overflow: the Regions' overlapping items past its lanes are one density row
    const v = (viewer = await open(FilesViews.find((x) => x.id === 'uniprot')!, {
      config: [[SeqCoreConfig.render, { track: { maxPx: 44, overflow: 'density' } }]],
    }));
    const regions = v
      .frame!.layout.sections[0].tracks.flatMap((t) => t.layers)
      .find((l) => l.key.includes('/region/'))!;
    const lay = regions.lay as RangesLay;
    expect(lay.overflow).toBeGreaterThan(0);
    const probes = probesOf(v.frame!, 8).filter((p) => p.layer === regions.key && p.mark === 'image');
    expect(probes.length).toBeGreaterThan(0);
    for (const p of probes) {
      const r = hoverAt(v, p.x, p.y);
      expect(v.hover.current.layer).toBe(regions.key);
      expect(r.marked, `${p.x},${p.y}`).toBe(true);
      // its card names what lies there: the items under the pointer, else the position
      expect(r.card?.rows.length).toBeGreaterThan(0);
      v.input.leave();
    }
  });

  it('a ligand’s score box lighting its residues but not its arcs (Complex › As a complex): the box lights the arcs its position lights, with a card', async () => {
    const v = (viewer = await open(ComplexViews.find((x) => x.id === 'units')!, {
      behaviors: [{ provider: LinksHoverProvider, params: undefined }],
    }));
    const scores = v
      .frame!.layout.sections.flatMap((s) => s.lays.flat())
      .filter((l) => l && l.type.name === 'ranges' && (l.lay as RangesLay).m > 0 && l.key.includes('/scores/'));
    expect(scores.length).toBeGreaterThan(0);
    let checked = 0;
    for (const l of scores) {
      const probe = probesOf(v.frame!).find((p) => p.layer === l!.key && p.mark === 'rect');
      if (!probe) continue;
      // the box: its item
      const box = hoverAt(v, probe.x, probe.y);
      const item = v.hover.current;
      expect(item.element, l!.key).toBeDefined();
      const arcs = hoverPaths(box.f);
      expect(box.marked).toBe(true);
      expect(box.card?.rows.length).toBeGreaterThan(0);
      // the same position, hit as a position (the empty band's): the same arcs light
      const at = { layer: item.layer!, system: item.system!, position: item.extent!.parts[0].ranges[0] };
      v.input.leave();
      v.hover.set(at, HoverOwner);
      const position = hoverPaths(v.draw() ?? v.frame!);
      v.hover.set(undefined, HoverOwner);
      expect(position, l!.key).toBeGreaterThan(0);
      expect(arcs, l!.key).toBe(position);
      checked++;
    }
    expect(checked).toBeGreaterThan(0);
  });

  it('the strands’ density bands showing no hover (Plasmid › Ring at scale): on the ring and on the detail, a mark and a card', async () => {
    const v = (viewer = await open(PlasmidViews.find((x) => x.id === 'scale')!));
    const bands = probesOf(v.frame!).filter((p) => p.type === 'ranges' && p.mark === 'image');
    const sections = new Set<string>();
    for (const p of bands) {
      if (v.hit(p.x, p.y)?.target === 'window') continue;
      const r = hoverAt(v, p.x, p.y);
      expect(r.marked, `${p.layer} ${p.x},${p.y}`).toBe(true);
      expect(r.card?.title, p.layer).toBeTruthy();
      sections.add(p.section);
      v.input.leave();
    }
    // both views: the ring and its linear detail
    expect([...sections].sort()).toEqual(['/views/main/detail', '/views/main/ring']);
  });
});
