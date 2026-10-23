/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { EntitySourceColorThemeProvider } from '../color/entity-source.js';

// entity 1 is a zinc ion, entity 2 the only polymer: `sequence.sequences` has one entry, for entity 2
const Cif = `data_TEST
loop_
_entity.id
_entity.type
1 non-polymer
2 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
2 'polypeptide(L)'
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
2 1 GLY n
2 2 ALA n
loop_
_entity_src_gen.entity_id
_entity_src_gen.pdbx_src_id
_entity_src_gen.pdbx_gene_src_scientific_name
_entity_src_gen.pdbx_beg_seq_num
_entity_src_gen.pdbx_end_seq_num
_entity_src_gen.pdbx_gene_src_gene
_entity_src_gen.plasmid_name
2 1 'Homo sapiens' 1 2 TP53 ?
1 2 'Homo sapiens' . . ZN1 ?
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM   1 N  N  GLY A 2 1 0.0 0.0 0.0 1 A 1
ATOM   2 C  CA GLY A 2 1 1.4 0.0 0.0 1 A 1
ATOM   3 N  N  ALA A 2 2 3.0 0.5 0.0 2 A 1
ATOM   4 C  CA ALA A 2 2 4.2 0.9 0.0 2 A 1
HETATM 5 ZN ZN ZN  B 1 . 8.0 0.0 0.0 101 A 1
`;

describe('entity-source color theme (perf-audit-2 09/B2)', () => {
  it('colors a polymer whose entity follows a non-polymer one, and gives the ligand the default color', async () => {
    const parsed = await CIF.parseText(Cif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    expect(model.sequence.sequences.map((s) => s.entityId)).toEqual(['2']);
    const structure = Structure.ofModel(model);

    const params = EntitySourceColorThemeProvider.getParams({ structure });
    const theme = EntitySourceColorThemeProvider.factory({ structure }, PD.getDefaultValues(params));
    const colorOf = (atomId: number) => {
      const unit = structure.units.find((u) => u.elements.includes(atomId as any))!;
      return theme.color(StructureElement.Location.create(structure, unit, atomId as any), false);
    };
    const ligand = colorOf(4);
    expect(colorOf(0)).toBe(colorOf(3));
    expect(colorOf(0)).not.toBe(ligand);
    expect(Number.isNaN(ligand)).toBe(false);
  });
});
