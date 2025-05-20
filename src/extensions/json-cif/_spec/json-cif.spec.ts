/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { molfileToJSONCif } from '../utils';
import { CifFile } from '../../../mol-io/reader/cif';
import { trajectoryFromMmCIF } from '../../../mol-model-formats/structure/mmcif';
import { Task } from '../../../mol-task';
import { JSONCifLigandGraph } from '../ligand-graph';
import { parseJSONCif } from '../parser';
import { JSONCifDataBlock } from '../model';

describe('json-cif', () => {
    it('roundtrips', async () => {
        const { structure, jsoncif: file } = await molfileToJSONCif(MolString);

        expect(file.dataBlocks.length).toBe(1);
        expect(file.dataBlocks[0].categoryNames.length).toBe(2);
        expect(file.dataBlocks[0].categoryNames[0]).toBe('atom_site');
        expect(file.dataBlocks[0].categories['atom_site'].rows.length).toBe(structure.elementCount);

        const parsed = parseJSONCif(file);
        const parsedModel = await parseCifModel(parsed);
        expect(parsedModel.atomicHierarchy.atoms._rowCount).toBe(structure.elementCount);
    });

    it('ligand graph', async () => {
        const { structure, jsoncif: file } = await molfileToJSONCif(MolString);

        // remove atom
        let graph = new JSONCifLigandGraph(file.dataBlocks[0]);
        graph.removeAtom(graph.atoms[0]);
        let data = graph.getData().block;
        expect(data.categories.atom_site.rows.length).toBe(structure.elementCount - 1);

        // modify atom
        graph = new JSONCifLigandGraph(file.dataBlocks[0]);
        expect(file.dataBlocks[0].categories.atom_site.rows[0].type_symbol !== 'N').toBe(true);
        graph.modifyAtom(1, { type_symbol: 'N' });
        data = graph.getData().block;
        expect(data.categories.atom_site.rows[0].type_symbol).toBe('N');

        // add atom and bond
        graph = new JSONCifLigandGraph(file.dataBlocks[0]);
        const newAtom = graph.addAtom({ type_symbol: 'C', Cartn_x: 0, Cartn_y: 0, Cartn_z: 0 });
        graph.addOrUpdateBond(graph.atoms[0], newAtom, { value_order: 'sing', type_id: 'covale' });
        data = graph.getData().block;
        expect(data.categories.atom_site.rows.length).toBe(structure.elementCount + 1);
        expect(data.categories.molstar_bond_site.rows.length).toBe(file.dataBlocks[0].categories.molstar_bond_site.rows.length + 1);

        // remove bond
        graph.removeBond(graph.atoms[0], newAtom);
        data = graph.getData().block;
        expect(data.categories.atom_site.rows.length).toBe(structure.elementCount + 1);
        expect(data.categories.molstar_bond_site.rows.length).toBe(file.dataBlocks[0].categories.molstar_bond_site.rows.length);
    });

    it('ligand graph traversal', () => {
        const data: JSONCifDataBlock = {
            header: 'test',
            categoryNames: ['atom_site', 'molstar_bond_site'],
            categories: {
                atom_site: {
                    name: 'atom_site',
                    fieldNames: ['id', 'type_symbol'],
                    rows: [
                        { id: 1, type_symbol: 'C', Cartn_x: 0, Cartn_y: 0, Cartn_z: 0 },
                        { id: 2, type_symbol: 'C', Cartn_x: 1, Cartn_y: 0, Cartn_z: 0 },
                        { id: 3, type_symbol: 'C', Cartn_x: 2, Cartn_y: 0, Cartn_z: 0 },
                        { id: 4, type_symbol: 'C', Cartn_x: 2, Cartn_y: 0, Cartn_z: 0 },
                    ],
                },
                molstar_bond_site: {
                    name: 'molstar_bond_site',
                    fieldNames: ['atom_id_1', 'atom_id_2'],
                    rows: [
                        { atom_id_1: 1, atom_id_2: 4 },
                        { atom_id_1: 1, atom_id_2: 2 },
                        { atom_id_1: 2, atom_id_2: 3 },
                    ],
                },
            },
        };

        const graph = new JSONCifLigandGraph(data);

        const bfs = graph.traverse(1, 'bfs', [] as number[], (a, s) => s.push(a.row.id!));
        expect(bfs).toEqual([1, 4, 2, 3]);
        const dfs = graph.traverse(1, 'dfs', [] as number[], (a, s) => s.push(a.row.id!));
        expect(dfs).toEqual([1, 2, 3, 4]);
    });
});

async function parseCifModel(file: CifFile) {
    const models = await trajectoryFromMmCIF(file.blocks[0], file).run();
    const model = await Task.resolveInContext(models.getFrameAtIndex(0));
    return model;
}

const MolString = `2244
  -OEChem-04072009073D

 21 21  0     0  0  0  0  0  0999 V2000
    1.2333    0.5540    0.7792 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6952   -2.7148   -0.7502 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7958   -2.1843    0.8685 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7813    0.8105   -1.4821 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0857    0.6088    0.4403 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7927   -0.5515    0.1244 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7288    1.8464    0.4133 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1426   -0.4741   -0.2184 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0787    1.9238    0.0706 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7855    0.7636   -0.2453 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1409   -1.8536    0.1477 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1094    0.6715   -0.3113 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5305    0.5996    0.1635 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1851    2.7545    0.6593 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7247   -1.3605   -0.4564 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5797    2.8872    0.0506 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8374    0.8238   -0.5090 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7290    1.4184    0.8593 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2045    0.6969   -0.6924 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7105   -0.3659    0.6426 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2555   -3.5916   -0.7337 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  1 12  1  0  0  0  0
  2 11  1  0  0  0  0
  2 21  1  0  0  0  0
  3 11  2  0  0  0  0
  4 12  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  2  0  0  0  0
  6  8  2  0  0  0  0
  6 11  1  0  0  0  0
  7  9  1  0  0  0  0
  7 14  1  0  0  0  0
  8 10  1  0  0  0  0
  8 15  1  0  0  0  0
  9 10  2  0  0  0  0
  9 16  1  0  0  0  0
 10 17  1  0  0  0  0
 12 13  1  0  0  0  0
 13 18  1  0  0  0  0
 13 19  1  0  0  0  0
 13 20  1  0  0  0  0
M  END`;