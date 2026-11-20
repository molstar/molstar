/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The molecule forms on a resolved spec, and `elements` over what the removed `component` and `atom` forms read:
// positions by component (alternatives, a component over several positions, a class through the compounds, the
// non-standard ones and the modified-type items beside them), molecules by identifier and class, a complex's units by
// chain, operator and author chain, atoms by name, the atoms system.

import { type Node, printAddress } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { type Evaluation, evaluateQuery, partOf } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { moleculesSpec } from '../../lang/_test/molecules.js';
import { resolveSpec } from '../model.js';

function withPtm(): Node {
  const spec = moleculesSpec();
  // a database's modification known by name only, on the Fc chain: an annotation, found beside the sequence's own
  const b = createSeqSpec();
  const fc = '/data/fc';
  b.data.ranges('ptm', {
    target: fc as never,
    rows: { start: [8], end: [9], category: ['glycosylation'], label: ['N-linked'] },
  });
  b.data.ranges('other', { target: fc as never, rows: { start: [10], end: [11], category: ['domain'] } });
  const data = spec.children!.find((c) => c.kind === 'data')!;
  const extra = b.build().children!.find((c) => c.kind === 'data')!.children!;
  return {
    ...spec,
    children: spec.children!.map((c) => (c === data ? { ...data, children: [...data.children!, ...extra] } : c)),
  };
}

const spec = withPtm();
const model = resolveSpec(spec);
const run = (address: string) => evaluateQuery(spec, model, address, { limit: 100 });
const loc = (r: Evaluation, system: string) =>
  r.result?.type === 'location' ? Array.from(partOf(r.result.location, system)?.ranges ?? []) : r.issues;
const nodes = (r: Evaluation) =>
  r.result?.type === 'nodes' ? r.result.nodes.map((m) => printAddress(m.address)) : r.issues;
const codes = (r: Evaluation) => r.issues.map((i) => i.code);

describe('(elements …) over components (the removed `component` form’s cases)', () => {
  it('positions by component: one or several codes, alternatives and a component over several positions', () => {
    expect(loc(run('/data/fc/(elements :component MSE)'), '/data/fc')).toEqual([0, 1]);
    expect(loc(run('/data/fc/(elements :component [CRO MSE])'), '/data/fc')).toEqual([0, 1, 3, 6]);
    // a standard position by its letter's standard component; an alternative counts at its position
    expect(loc(run('/data/trna/(elements :component U)'), '/data/trna')).toEqual([5, 8, 11, 13]);
    expect(loc(run('/data/fc/(elements :component LYS)'), '/data/fc')).toEqual([1, 2, 7, 8, 15, 16]);
    expect(loc(run('/data/glycan/(elements :label_comp_id NAG)'), '/data/glycan')).toEqual([0, 2]);
  });

  it('the non-standard positions (standard), and with the modified-type items of the annotations (modified)', () => {
    expect(loc(run('/data/fc/(elements :modified true)'), '/data/fc')).toEqual([0, 1, 3, 6, 8, 9]);
    expect(loc(run('/data/fc/(elements :standard false)'), '/data/fc')).toEqual([0, 1, 3, 6]);
    // a unit's are its sequence's, shifted and clipped to it
    expect(loc(run('/data/c/A/(elements :standard false)'), '/data/c/A')).toEqual([0, 1, 3, 6]);
    expect(loc(run('/data/fc/(except :of (elements :standard false))'), '/data/fc')).toEqual([1, 3, 6, 22]);
    expect(loc(run('/data/fc/(elements :component CRO :standard true)'), '/data/fc')).toEqual([]);
  });

  it('a class through the compound nodes, a word of the class; a code with none is no match', () => {
    expect(loc(run('/data/fc/(elements :component-class peptide)'), '/data/fc')).toEqual([3, 6]);
    expect(loc(run('/data/fc/(elements :component-class saccharide)'), '/data/fc')).toEqual([]);
  });

  it('refuses a field no property has, and a value of the wrong type', () => {
    expect(codes(run('/data/fc/(elements :class peptide)'))).toEqual(['seq.query.form.param']);
    expect(codes(run('/data/fc/(elements :standard 1)'))).toEqual(['seq.query.form.type']);
  });
});

describe('(molecule …), (chain …)', () => {
  it('molecules by identifier, matched exactly, and by class; under data or a complex', () => {
    expect(nodes(run('/data/(molecule :identifier "ccd:ATP")'))).toEqual(['/data/atp']);
    expect(nodes(run('/data/(molecule :identifier "ccd:atp")'))).toEqual([]);
    expect(nodes(run('/data/(molecule :class branched)'))).toEqual(['/data/glycan']);
    expect(nodes(run('/data/c/(molecule :class ligand)'))).toEqual(['/data/atp']);
  });

  it('a complex’s units by chain, operator and author chain, each by its id; the unit then a prefix', () => {
    expect(nodes(run('/data/c/(chain :id A)'))).toEqual(['/data/c/A']);
    expect(nodes(run('/data/c/(chain :auth H :operator 1_555)'))).toEqual(['/data/c/A']);
    const none = run('/data/c/(chain :id A :operator 2_555)');
    expect([nodes(none), codes(none)]).toEqual([[], ['seq.query.node.missing']]);
    expect(loc(run('/data/c/(chain :id A)/(at :label 240 :numbering auth)'), '/data/c/A')).toEqual([10, 11]);
    expect(codes(run('/data/c/(chain :operator 1_555)'))).toEqual(['seq.query.form.param']);
  });

  it('an entity is a source’s: a spec complex’s units have none, `(molecule :entity 1)` names nothing', () => {
    expect(nodes(run('/data/c/(molecule :entity 1)'))).toEqual([]);
    expect(codes(run('/data/c/(molecule)'))).toEqual(['seq.query.form.param']);
  });
});

describe('(atoms), (elements :atom …)', () => {
  it('a molecule’s or a unit’s atoms system, a prefix the location forms read', () => {
    expect(nodes(run('/data/atp/(atoms)'))).toEqual(['/data/atp/(atoms)']);
    expect(loc(run('/data/atp/(atoms)/(range :from 1 :to 2)'), '/data/atp/(atoms)')).toEqual([0, 2]);
    expect(nodes(run('/data/c/L/(atoms)'))).toEqual(['/data/c/L/(atoms)']);
    expect(codes(run('/data/fc/(atoms)'))).toEqual(['seq.query.node.missing']);
  });

  it('atoms by name, on the system’s atoms; on a polymer, the residues whose standard component has them', () => {
    expect(loc(run('/data/atp/(elements :atom O1G)'), '/data/atp/(atoms)')).toEqual([1, 2]);
    expect(loc(run('/data/c/L/(elements :atom [PG N1])'), '/data/c/L/(atoms)')).toEqual([0, 1, 2, 3]);
    // residues, as a structure viewer reads them: a CA on every amino acid, an SD on the MET under the MSE
    const ca = run('/data/fc/(elements :atom CA)');
    expect([loc(ca, '/data/fc'), codes(ca)]).toEqual([[0, 22], []]);
    expect(loc(run('/data/fc/(elements :atom SD)'), '/data/fc')).toEqual([0, 1]);
    expect(loc(run('/data/fc/(elements :element S)'), '/data/fc')).toEqual([0, 1]);
    expect(loc(run('/data/fc/(elements :atom O1G)'), '/data/fc')).toEqual([]);
    expect(loc(run('/data/trna/(elements :atom "O2\'" :component G)'), '/data/trna')).toEqual([0, 1, 2, 4, 14, 15]);
    // an atom's index is a held atom's alone
    expect(codes(run('/data/fc/(elements :atom-index 0)'))).toEqual(['seq.query.property.absent']);
  });
});

describe('the host’s identifier resolver and component names', () => {
  const host = resolveSpec(spec, {
    equivalents: (id) => (id === 'chembl:CHEMBL14249' ? ['ccd:ATP'] : []),
    componentNames: (code) => (code === 'NAG' ? { name: 'N-acetyl-D-glucosamine', class: 'D-saccharide' } : undefined),
  });
  const hosted = (address: string) => evaluateQuery(spec, host, address, { limit: 100 });

  it('a molecule by an identifier the host holds the same compound as one of its own', () => {
    expect(nodes(hosted('/data/(molecule :identifier "chembl:CHEMBL14249")'))).toEqual(['/data/atp']);
    expect(nodes(run('/data/(molecule :identifier "chembl:CHEMBL14249")'))).toEqual([]);
  });

  it('a class through the host’s names where no compound node describes the code', () => {
    expect(loc(hosted('/data/glycan/(elements :component-class saccharide)'), '/data/glycan')).toEqual([0, 2]);
    expect(loc(run('/data/glycan/(elements :component-class saccharide)'), '/data/glycan')).toEqual([]);
    expect(host.compound('NAG')?.name).toBe('N-acetyl-D-glucosamine');
    // a compound node wins over the host's names
    expect(host.compound('CRO')?.name).toBe('chromophore');
  });
});
