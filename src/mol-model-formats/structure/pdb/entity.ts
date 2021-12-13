/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Tokens } from '../../../mol-io/reader/common/text/tokenizer';
import { EntityCompound } from '../common/entity';

const Spec = {
    'MOL_ID': '',
    'MOLECULE': '',
    'CHAIN': '',
    'FRAGMENT': '',
    'SYNONYM': '',
    'EC': '',
    'ENGINEERED': '',
    'MUTATION': '',
    'OTHER_DETAILS': ''
};
type Spec = keyof typeof Spec

export function parseCmpnd(lines: Tokens, lineStart: number, lineEnd: number) {
    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);

    let currentSpec: Spec | undefined;
    let currentCompound: EntityCompound = { chains: [], description: '' };
    const compounds: EntityCompound[] = [];

    for (let i = lineStart; i < lineEnd; i++) {
        const line = getLine(i);
        // COLUMNS       DATA TYPE       FIELD         DEFINITION
        // ----------------------------------------------------------------------------------
        //  1 -  6       Record name     "COMPND"
        //  8 - 10       Continuation    continuation  Allows concatenation of multiple records.
        // 11 - 80       Specification   compound      Description of the molecular components.
        //               list

        const cmpnd = line.substr(10, 70).trim();
        const cmpndSpecEnd = cmpnd.indexOf(':');
        const cmpndSpec = cmpnd.substring(0, cmpndSpecEnd);

        let value: string;

        if (cmpndSpec in Spec) {
            currentSpec = cmpndSpec as Spec;
            value = cmpnd.substring(cmpndSpecEnd + 2);
        } else {
            value = cmpnd;
        }
        value = value.replace(/;$/, '');

        if (currentSpec === 'MOL_ID') {
            currentCompound = {
                chains: [],
                description: ''
            };
            compounds.push(currentCompound);
        } else if (currentSpec === 'MOLECULE') {
            if (currentCompound.description) currentCompound.description += ' ';
            currentCompound.description += value;
        } else if (currentSpec === 'CHAIN') {
            Array.prototype.push.apply(currentCompound.chains, value.split(/\s*,\s*/));
        }
    }

    // Define a seprate entity for each chain
    // --------------------------------------
    //
    // This is a workaround for how sequences are currently determined for PDB files.
    //
    // The current approach infers the "observed sequence" from the atomic hierarchy.
    // However, for example for PDB ID 3HHR, this approach fails, since chains B and C
    // belong to the same entity but contain different observed sequence, which causes display
    // errors in the sequence viewer (since the sequences are determined "per entity").
    //
    // A better approach could be to parse SEQRES categories and use it to construct
    // entity_poly_seq category. However, this would require constructing label_seq_id (with gaps)
    // from RES ID pdb column (auth_seq_id), which isn't a trivial exercise.
    //
    // (properly formatted) mmCIF structures do not exhibit this issue.
    const singletons: EntityCompound[] = [];
    for (const comp of compounds) {
        for (const chain of comp.chains) {
            singletons.push({
                description: comp.description,
                chains: [chain]
            });
        }
    }
    return singletons;
}

export function parseHetnam(lines: Tokens, lineStart: number, lineEnd: number) {
    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);

    const hetnams = new Map<string, string>();

    for (let i = lineStart; i < lineEnd; i++) {
        const line = getLine(i);
        // COLUMNS       DATA  TYPE    FIELD           DEFINITION
        // ----------------------------------------------------------------------------
        //  1 -  6       Record name   "HETNAM"
        //  9 - 10       Continuation  continuation    Allows concatenation of multiple records.
        // 12 - 14       LString(3)    hetID           Het identifier, right-justified.
        // 16 - 70       String        text            Chemical name.

        const het = line.substr(11, 3).trim();
        const name = line.substr(15).trim();

        if (hetnams.has(het)) {
            hetnams.set(het, `${hetnams.get(het)!} ${name}`);
        } else {
            hetnams.set(het, name);
        }
    }

    return hetnams;
}