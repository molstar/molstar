/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Hit as RenderHit } from '@molstar/seqx-core/render';
import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { hoverBehavior } from '../behaviors/hover.js';
import { bindingTable, DefaultBindings } from '../bindings.js';
import { viewerConfig } from '../config.js';
import { type Hit, HoverModel, hitTarget } from '../models/hover.js';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';
import { type FixtureOptions, Gene, GeneSection, Protein, ProteinSection, partsOf, proteinOnGene } from './fixture.js';

function setup(options: FixtureOptions = {}, config = viewerConfig()) {
  const spec = new SpecModel(proteinOnGene(options));
  const view = new ViewModel(spec, config);
  const hover = new HoverModel(view, config);
  return { spec, view, hover };
}

const residue = (position: number): Hit => ({ layer: `${ProteinSection}/letters/l`, system: Protein, position });
const nucleotide = (position: number): Hit => ({ layer: `${GeneSection}/letters/l`, system: Gene, position });

describe('HoverModel: ownership (EC-141)', () => {
  it('one owner at a time: another owner takes it; a stale owner’s clear does nothing', () => {
    const { hover } = setup();
    const pointer = {};
    const structure = {};
    expect(hover.set(residue(1), pointer)).toBe(true);
    expect(hover.set(residue(2), structure)).toBe(true);
    expect(hover.owner).toBe(structure);
    // a late leave of the pointer never clears the 3D hover
    expect(hover.clear(pointer)).toBe(false);
    expect(hover.current.position).toBe(2);
    expect(hover.clear(structure)).toBe(true);
    expect(hover.active).toBe(false);
    expect(hover.location()).toBeUndefined();
    expect(hover.clear(structure)).toBe(false);
  });

  it('set copies the hit into one reused record, moves the version, and tells nothing for the same hit', () => {
    const { hover } = setup();
    const owner = {};
    const record = hover.current;
    let told = 0;
    const off = hover.version.listen(() => told++);
    hover.set(residue(1), owner);
    hover.set(residue(2), owner);
    expect(hover.set(residue(2), owner)).toBe(false);
    expect(hover.current).toBe(record);
    expect(record).toMatchObject({ owner, system: Protein, position: 2, section: ProteinSection });
    expect(told).toBe(2);
    off();
  });

  it('an element is the same hover by its key, whatever extent object comes with it; no owner is refused', () => {
    const { hover } = setup();
    const element = { owner: '/data/domains', key: 'TAD' };
    expect(hover.set({ system: Protein, element, extent: rangesLocation(Protein, [0, 3]) }, 'x')).toBe(true);
    expect(hover.set({ system: Protein, element: { ...element }, extent: rangesLocation(Protein, [0, 3]) }, 'x')).toBe(
      false,
    );
    expect(hover.set(residue(1), undefined)).toBe(false);
    expect(hover.owner).toBe('x');
  });

  it('an element’s hover is its extent; a hit’s target is what it hit', () => {
    const { hover } = setup();
    const extent = rangesLocation(Protein, [0, 3]);
    const element = { owner: '/data/domains', key: 'TAD' };
    hover.set({ system: Protein, element, extent }, 'x');
    expect(hover.location()).toBe(extent);
    expect(hitTarget({ system: Protein, element })).toBe('item');
    expect(hitTarget(residue(1))).toBe('cell');
    expect(hitTarget({ system: Protein })).toBe('empty');
    expect(hitTarget({ system: Protein, target: 'axis' })).toBe('axis');
  });

  it('a cell of two axes keeps its row flat in the record: a new row is a new hover, the same is none', () => {
    const { hover } = setup();
    const cell = (row: number): Hit => ({ ...residue(4), ySystem: '/data/aln/(rows)', yPosition: row, target: 'cell' });
    let told = 0;
    const off = hover.version.listen(() => told++);
    expect(hover.set(cell(7), 'x')).toBe(true);
    expect(hover.current).toMatchObject({ position: 4, ySystem: '/data/aln/(rows)', yPosition: 7 });
    expect(hover.set(cell(7), 'x')).toBe(false);
    expect(hover.set(cell(8), 'x')).toBe(true);
    expect(hover.set(residue(4), 'x')).toBe(true);
    expect(hover.current.ySystem).toBeUndefined();
    expect(Number.isNaN(hover.current.yPosition)).toBe(true);
    expect(told).toBe(3);
    expect(hitTarget({ system: Protein, element: { owner: 'a', key: 'c' }, target: 'tree-node' })).toBe('tree-node');
    expect(hitTarget({ system: Protein, element: { owner: 'a', key: 'c' }, target: 'summary-row' })).toBe(
      'summary-row',
    );
    off();
  });
});

describe('HoverModel: projection onto linked sections (EC-30, EC-144)', () => {
  it('a residue lights its codon on the gene, every piece of one split by the intron, on both strands', () => {
    const plus = setup();
    plus.hover.set(residue(1), 'p');
    expect(partsOf(plus.hover.projected(ProteinSection))).toEqual([[Protein, [1, 2]]]);
    expect(partsOf(plus.hover.projected(GeneSection))).toEqual([[Gene, [5, 7, 12, 13]]]);
    const minus = setup({ strand: '-' });
    minus.hover.set(residue(1), 'p');
    expect(partsOf(minus.hover.projected(GeneSection))).toEqual([[Gene, [13, 16]]]);
    minus.hover.set(residue(2), 'p');
    expect(partsOf(minus.hover.projected(GeneSection))).toEqual([[Gene, [5, 7, 12, 13]]]);
  });

  it('a nucleotide lifts to its whole codon on the gene and to its residue on the protein, on both strands', () => {
    const plus = setup();
    plus.hover.set(nucleotide(12), 'p');
    expect(partsOf(plus.hover.projected(GeneSection))).toEqual([[Gene, [5, 7, 12, 13]]]);
    expect(partsOf(plus.hover.projected(ProteinSection))).toEqual([[Protein, [1, 2]]]);
    const minus = setup({ strand: '-' });
    minus.hover.set(nucleotide(12), 'p');
    expect(partsOf(minus.hover.projected(GeneSection))).toEqual([[Gene, [5, 7, 12, 13]]]);
    expect(partsOf(minus.hover.projected(ProteinSection))).toEqual([[Protein, [2, 3]]]);
    // an untranslated nucleotide lights itself and nothing on the protein
    minus.hover.set(nucleotide(0), 'p');
    expect(partsOf(minus.hover.projected(GeneSection))).toEqual([[Gene, [0, 1]]]);
    expect(partsOf(minus.hover.projected(ProteinSection))).toEqual([]);
  });

  it('is computed when asked, once per version; configurable', () => {
    const { hover } = setup();
    hover.set(residue(1), 'p');
    const a = hover.projected(GeneSection);
    expect(hover.projected(GeneSection)).toBe(a);
    hover.set(residue(2), 'p');
    expect(hover.projected(GeneSection)).not.toBe(a);
    expect(hover.projected('/views/main/none')).toEqual({ parts: [], rects: [] });

    const off = setup({}, viewerConfig({ hover: { projectLinked: false, closeCodons: false } }));
    off.hover.set(nucleotide(12), 'p');
    expect(partsOf(off.hover.projected(GeneSection))).toEqual([[Gene, [12, 13]]]);
    expect(partsOf(off.hover.projected(ProteinSection))).toEqual([]);
  });
});

describe('the hover behavior', () => {
  it('moves set the hover through the bindings; a move over nothing or a leave clears only its own', () => {
    const { hover } = setup();
    const pointer = hoverBehavior(hover, bindingTable());
    pointer.move(residue(1));
    expect(hover.owner).toBe(pointer.owner);
    pointer.move(undefined);
    expect(hover.active).toBe(false);
    pointer.move(residue(1));
    hover.set(residue(3), 'structure');
    pointer.leave();
    expect(hover.current.position).toBe(3);
    // unbound on a target: no hover there
    const none = hoverBehavior(
      hover,
      bindingTable(DefaultBindings, [[{ gesture: 'hover', modifiers: '*', target: 'any', action: null }]]),
    );
    none.move(residue(1));
    expect(hover.owner).toBe('structure');
  });

  it('holds with any modifier held', () => {
    const { hover } = setup();
    const pointer = hoverBehavior(hover, bindingTable());
    for (const m of [{}, { shift: true }, { mod: true }, { alt: true }, { space: true }]) {
      pointer.leave();
      pointer.move(residue(1), m);
      expect(hover.owner).toBe(pointer.owner);
    }
  });
});

describe('the render hit', () => {
  it('is a viewer hit as the layout gives it', () => {
    const r: RenderHit = { layer: 'l', section: 's', system: 'p', position: 3 };
    const h: Hit = r;
    expect(hitTarget(h)).toBe('cell');
  });
});
