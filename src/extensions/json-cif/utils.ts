/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseMol } from '../../mol-io/reader/mol/parser';
import { trajectoryFromMol } from '../../mol-model-formats/structure/mol';
import { Structure, to_mmCIF } from '../../mol-model/structure';
import { Task } from '../../mol-task';
import { JSONCifEncoder } from './encoder';

export async function molfileToJSONCif(molfile: string) {
    const parsed = await parseMol(molfile).run();
    if (parsed.isError) throw new Error(parsed.message);
    const models = await trajectoryFromMol(parsed.result).run();
    const model = await Task.resolveInContext(models.getFrameAtIndex(0));
    const structure = Structure.ofModel(model);
    const encoder = new JSONCifEncoder('Mol*', { formatJSON: true });

    to_mmCIF('mol', structure, false, {
        encoder,
        includedCategoryNames: new Set(['atom_site']),
        extensions: {
            molstar_bond_site: true,
        }
    });

    return {
        structure,
        molfile: parsed.result,
        jsoncif: encoder.getFile()
    };
}
