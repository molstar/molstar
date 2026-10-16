/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import {
  DpdbFixtures,
  mutationsTable,
  parseArtifacts,
  parseEntry,
  parseModels,
  parseMutations,
  residueValues,
  valuesTable,
} from '../data.js';

// Dynamic PDB's responses as the story reads them (`fixtures/`, the API's as they are): parsed, and the sequence's
// tables made from them. Read from disk, no network.

const fixture = (name: string) => JSON.parse(readFileSync(new URL(`../fixtures/${name}`, import.meta.url), 'utf8'));

describe('Dynamic PDB data', () => {
  it('reads the entry: its PDB id, method, construct and mutations', () => {
    const e = parseEntry(fixture('entry.json'));
    expect(e).toMatchObject({ id: 'dpdb_gijqn2no', pdb: '3K7F', method: 'X-ray crystallography', resolution: 1.95 });
    expect(e.title.endsWith("'")).toBe(false);
    expect(e.entities[0]).toMatchObject({
      construct: 'deltaTM-FAAH',
      organism: 'Rattus norvegicus',
      uniprot: 'P97612',
    });
    expect(e.entities[0]!.mutations.map((m) => m.text)).toEqual(['L192F', 'F194Y', 'A377T', 'S435N', 'I491V', 'V495M']);
    expect(() => parseEntry({ data: { type: 'models' } })).toThrow();
  });

  it('reads the models, the deposited first, with metrics and residue data', () => {
    const [deposited, qfit] = parseModels(fixture('models.json'));
    expect(deposited).toMatchObject({ title: 'Deposited model', modelType: 'Deposited' });
    expect(deposited!.metrics).toMatchObject({ r_work: 0.1545, r_free: 0.1883, clashscore: 1.99 });
    expect(deposited!.residueData).toHaveLength(1088);
    expect(deposited!.residueData[0]).toMatchObject({
      asym: 'A',
      seq: 26,
      auth: 32,
      uniprot: { accession: 'P97612', position: 32 },
    });
    expect(qfit).toMatchObject({ title: 'qFit model', modelType: 'Multiconformer', residueData: [] });
    expect(qfit!.altlocFraction).toBeCloseTo(0.4, 2);
  });

  it("reads a model's files and maps them to what the story can read", () => {
    const [sf, model] = parseArtifacts(fixture('artifacts-m002.json'));
    expect(sf).toMatchObject({ type: 'structure_factors' });
    expect(model).toMatchObject({ type: 'model', name: '3k7f_qFit_010.cif' });
    // The bucket has no CORS: the copy; RCSB's files as they are.
    expect(DpdbFixtures.file(model!)).toMatchObject({ format: 'bcif' });
    expect(DpdbFixtures.file(model!).url).toMatch(/fixtures\/3k7f-qfit\.bcif$/);
    const [, deposited] = parseArtifacts(fixture('artifacts-m001.json'));
    expect(DpdbFixtures.file(deposited!)).toEqual({ url: 'https://files.rcsb.org/download/3K7F.cif', format: 'mmcif' });
  });

  it('reads mutations as written, leaving out what does not read', () => {
    expect(parseMutations('L192F, x, F194Y;A3T')).toEqual([
      { text: 'L192F', from: 'L', position: 192, to: 'F' },
      { text: 'F194Y', from: 'F', position: 194, to: 'Y' },
      { text: 'A3T', from: 'A', position: 3, to: 'T' },
    ]);
    expect(parseMutations(undefined)).toEqual([]);
  });

  it("makes the tracks' tables in the units' auth numbering, for the residues the model has", () => {
    const [deposited] = parseModels(fixture('models.json'));
    const residues = deposited!.residueData;
    // A model with chain A only (author chain A is label chain A), residues 32–600.
    const chains = new Map<string, string>();
    for (let auth = 32; auth <= 600; auth++) chains.set(`A:${auth}`, 'A');
    const b = valuesTable(chains, residueValues(residues, 'bIso'));
    expect(new Set(b.target)).toEqual(new Set(['A']));
    expect(b.pos[0]).toBe('32');
    expect(b.value[0]).toBe(72.42);
    expect(b.target.length).toBe(b.value.length);
    const mutations = parseEntry(fixture('entry.json')).entities[0]!.mutations;
    const v = mutationsTable(chains, residues, mutations);
    expect(v.label).toEqual(['L192F', 'F194Y', 'A377T', 'S435N', 'I491V', 'V495M']);
    expect(v.beg).toEqual(['192', '194', '377', '435', '491', '495']);
    expect(v.alt).toEqual(['F', 'Y', 'T', 'N', 'V', 'M']);
  });
});
