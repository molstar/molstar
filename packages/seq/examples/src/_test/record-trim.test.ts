/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { trimEntry, trimVariation } from '../files/record-trim.js';

// What the record script keeps (files.md, UniProt's entry as it is): members and records removed, none renamed or moved.

const feature = (type: string, extra = {}) => ({
  type,
  location: { start: { value: 1, modifier: 'EXACT' }, end: { value: 5, modifier: 'EXACT' } },
  description: 'd',
  featureId: 'F1',
  evidences: [{ id: 'x' }],
  featureCrossReferences: [{ id: 'y' }],
  ...extra,
});

describe('trimEntry', () => {
  it('keeps the members the views read, drops the evidence and the variants’ features, renames nothing', () => {
    const entry = {
      entryType: 'reviewed',
      primaryAccession: 'P1',
      comments: [1],
      references: [2],
      sequence: { value: 'MEEP', length: 4 },
      features: [feature('Region'), feature('Natural variant'), feature('Helix', { ligand: { name: 'Zn' } })],
    };
    const out = trimEntry(entry as never);
    expect(Object.keys(out)).toEqual(['entryType', 'primaryAccession', 'sequence', 'features']);
    expect((out.features as { type: string }[]).map((f) => f.type)).toEqual(['Region', 'Helix']);
    expect((out.features as object[]).map((f) => Object.keys(f))).toEqual([
      ['type', 'location', 'description', 'featureId'],
      ['type', 'location', 'description', 'featureId'],
    ]);
    expect((out.features as { location: unknown }[])[0].location).toEqual(entry.features[0].location);
  });
});

const variant = (begin: string, extra = {}) => ({
  type: 'VARIANT',
  begin,
  end: begin,
  wildType: 'M',
  mutatedType: 'K',
  alternativeSequence: 'K',
  consequenceType: 'missense',
  sourceType: 'mixed',
  xrefs: [
    { name: 'UniProt', id: 'VAR_1', url: 'u' },
    { name: 'dbSNP', id: 'rs1' },
    { name: 'ClinVar', id: 'RCV1' },
    { name: 'ClinVar', id: 'RCV2' },
  ],
  clinicalSignificances: [{ type: 'Pathogenic', sources: ['ClinVar'], reviewStatus: 'x' }],
  ...extra,
});

describe('trimVariation', () => {
  it('keeps UniProt’s single-residue changes inside the sequence, with xrefs by name and significances by type', () => {
    const response = {
      accession: 'P1',
      sequence: 'MEEPQ',
      other: 1,
      features: [
        variant('2'),
        variant('9'),
        variant('0'),
        variant('2', { end: '3' }),
        variant('3', { sourceType: 'large_scale_study' }),
        variant('3', { wildType: undefined }),
        variant('4', { xrefs: [{ name: 'dbSNP', id: 'rs2' }] }),
        variant('5', { sourceType: 'uniprot' }),
      ],
    };
    const out = trimVariation(JSON.parse(JSON.stringify(response)));
    expect(Object.keys(out)).toEqual(['accession', 'sequence', 'features']);
    const kept = out.features as Record<string, unknown>[];
    expect(kept.map((f) => f.begin)).toEqual(['2', '5']);
    expect(kept[0]).toEqual({
      begin: '2',
      end: '2',
      xrefs: [{ name: 'UniProt' }, { name: 'ClinVar' }],
      consequenceType: 'missense',
      wildType: 'M',
      mutatedType: 'K',
      sourceType: 'mixed',
      clinicalSignificances: [{ type: 'Pathogenic' }],
    });
  });
});
