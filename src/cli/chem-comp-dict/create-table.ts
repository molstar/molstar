#!/usr/bin/env node
/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import * as util from 'util';
import * as path from 'path';
import * as fs from 'fs';
require('util.promisify').shim();
const writeFile = util.promisify(fs.writeFile);

import { Database, Table, DatabaseCollection } from '../../mol-data/db';
import { CCD_Schema } from '../../mol-io/reader/cif/schema/ccd';
import { SetUtils } from '../../mol-util/set';
import { DefaultMap } from '../../mol-util/map';
import { mmCIF_chemCompBond_schema } from '../../mol-io/reader/cif/schema/mmcif-extras';
import { ccd_chemCompAtom_schema } from '../../mol-io/reader/cif/schema/ccd-extras';
import { DefaultDataOptions, ensureDataAvailable, getEncodedCif, readCCD, readPVCD } from './util';

type CCB = Table<CCD_Schema['chem_comp_bond']>
type CCA = Table<CCD_Schema['chem_comp_atom']>

function ccbKey(compId: string, atomId1: string, atomId2: string) {
    return atomId1 < atomId2 ? `${compId}:${atomId1}-${atomId2}` : `${compId}:${atomId2}-${atomId1}`;
}

function ccaKey(compId: string, atomId: string) {
    return `${compId}:${atomId}`;
}

function addChemCompBondToSet(set: Set<string>, ccb: CCB) {
    for (let i = 0, il = ccb._rowCount; i < il; ++i) {
        set.add(ccbKey(ccb.comp_id.value(i), ccb.atom_id_1.value(i), ccb.atom_id_2.value(i)));
    }
    return set;
}

function addChemCompAtomToSet(set: Set<string>, cca: CCA) {
    for (let i = 0, il = cca._rowCount; i < il; ++i) {
        set.add(ccaKey(cca.comp_id.value(i), cca.atom_id.value(i)));
    }
    return set;
}

function checkAddingBondsFromPVCD(pvcd: DatabaseCollection<CCD_Schema>) {
    const ccbSetByParent = DefaultMap<string, Set<string>>(() => new Set());

    for (const k in pvcd) {
        const { chem_comp, chem_comp_bond } = pvcd[k];
        if (chem_comp_bond._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length === 0) {
                const set = ccbSetByParent.getDefault(chem_comp.id.value(0));
                addChemCompBondToSet(set, chem_comp_bond);
            } else {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    const parentId = parentIds[i];
                    const set = ccbSetByParent.getDefault(parentId);
                    addChemCompBondToSet(set, chem_comp_bond);
                }
            }
        }
    }

    for (const k in pvcd) {
        const { chem_comp, chem_comp_atom, chem_comp_bond } = pvcd[k];
        if (chem_comp_bond._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length > 0) {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    const entryBonds = addChemCompBondToSet(new Set<string>(), chem_comp_bond);
                    const entryAtoms = addChemCompAtomToSet(new Set<string>(), chem_comp_atom);
                    const extraBonds = SetUtils.difference(ccbSetByParent.get(parentIds[i])!, entryBonds);
                    extraBonds.forEach(bk => {
                        const [a1, a2] = bk.split('|');
                        if (entryAtoms.has(a1) && entryAtoms.has(a2)) {
                            console.error(`Adding all PVCD bonds would wrongly add bond ${bk} for ${k}`);
                        }
                    });
                }
            }
        }
    }
}

function checkAddingAtomsFromPVCD(pvcd: DatabaseCollection<CCD_Schema>) {
    const ccaSetByParent = DefaultMap<string, Set<string>>(() => new Set());

    for (const k in pvcd) {
        const { chem_comp, chem_comp_atom } = pvcd[k];
        if (chem_comp_atom._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length === 0) {
                const set = ccaSetByParent.getDefault(chem_comp.id.value(0));
                addChemCompAtomToSet(set, chem_comp_atom);
            } else {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    const parentId = parentIds[i];
                    const set = ccaSetByParent.getDefault(parentId);
                    addChemCompAtomToSet(set, chem_comp_atom);
                }
            }
        }
    }
}

async function createBonds(
    ccd: DatabaseCollection<CCD_Schema>,
    pvcd: DatabaseCollection<CCD_Schema>,
    atomsRequested: boolean
) {
    const ccbSet = new Set<string>();

    const comp_id: string[] = [];
    const atom_id_1: string[] = [];
    const atom_id_2: string[] = [];
    const value_order: typeof mmCIF_chemCompBond_schema['value_order']['T'][] = [];
    const pdbx_aromatic_flag: typeof mmCIF_chemCompBond_schema['pdbx_aromatic_flag']['T'][] = [];
    const pdbx_stereo_config: typeof mmCIF_chemCompBond_schema['pdbx_stereo_config']['T'][] = [];
    const molstar_protonation_variant: string[] = [];

    function addBonds(compId: string, ccb: CCB, protonationVariant: boolean) {
        for (let i = 0, il = ccb._rowCount; i < il; ++i) {
            const atomId1 = ccb.atom_id_1.value(i);
            const atomId2 = ccb.atom_id_2.value(i);
            const k = ccbKey(compId, atomId1, atomId2);
            if (!ccbSet.has(k)) {
                atom_id_1.push(atomId1);
                atom_id_2.push(atomId2);
                comp_id.push(compId);
                value_order.push(ccb.value_order.value(i));
                pdbx_aromatic_flag.push(ccb.pdbx_aromatic_flag.value(i));
                pdbx_stereo_config.push(ccb.pdbx_stereo_config.value(i));
                molstar_protonation_variant.push(protonationVariant ? 'Y' : 'N');
                ccbSet.add(k);
            }
        }
    }

    // check adding bonds from PVCD
    checkAddingBondsFromPVCD(pvcd);

    // add bonds from PVCD
    for (const k in pvcd) {
        const { chem_comp, chem_comp_bond } = pvcd[k];
        if (chem_comp_bond._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length === 0) {
                addBonds(chem_comp.id.value(0), chem_comp_bond, false);
            } else {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    addBonds(parentIds[i], chem_comp_bond, true);
                }
            }
        }
    }

    // add bonds from CCD
    for (const k in ccd) {
        const { chem_comp, chem_comp_bond } = ccd[k];
        if (chem_comp_bond._rowCount) {
            addBonds(chem_comp.id.value(0), chem_comp_bond, false);
        }
    }

    const bondTable = Table.ofArrays(mmCIF_chemCompBond_schema, {
        comp_id, atom_id_1, atom_id_2, value_order,
        pdbx_aromatic_flag, pdbx_stereo_config, molstar_protonation_variant
    });

    const bondDatabase = Database.ofTables(
        CCB_TABLE_NAME,
        { chem_comp_bond: mmCIF_chemCompBond_schema },
        { chem_comp_bond: bondTable }
    );

    return { bonds: bondDatabase, atoms: atomsRequested ? createAtoms(ccd, pvcd) : void 0 };
}

function createAtoms(ccd: DatabaseCollection<CCD_Schema>, pvcd: DatabaseCollection<CCD_Schema>) {
    const ccaSet = new Set<string>();

    const comp_id: string[] = [];
    const atom_id: string[] = [];
    const charge: number[] = [];
    const pdbx_stereo_config: typeof CCD_Schema.chem_comp_atom['pdbx_stereo_config']['T'][] = [];

    function addAtoms(compId: string, cca: CCA) {
        for (let i = 0, il = cca._rowCount; i < il; ++i) {
            const atomId = cca.atom_id.value(i);
            const k = ccaKey(compId, atomId);
            if (!ccaSet.has(k)) {
                atom_id.push(atomId);
                comp_id.push(compId);
                charge.push(cca.charge.value(i));
                pdbx_stereo_config.push(cca.pdbx_stereo_config.value(i));
                ccaSet.add(k);
            }
        }
    }

    // check adding atoms from PVCD
    checkAddingAtomsFromPVCD(pvcd);

    // add atoms from PVCD
    for (const k in pvcd) {
        const { chem_comp, chem_comp_atom } = pvcd[k];
        if (chem_comp_atom._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length === 0) {
                addAtoms(chem_comp.id.value(0), chem_comp_atom);
            } else {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    addAtoms(parentIds[i], chem_comp_atom);
                }
            }
        }
    }

    // add atoms from CCD
    for (const k in ccd) {
        const { chem_comp, chem_comp_atom } = ccd[k];
        if (chem_comp_atom._rowCount) {
            addAtoms(chem_comp.id.value(0), chem_comp_atom);
        }
    }

    const atomTable = Table.ofArrays(ccd_chemCompAtom_schema, {
        comp_id, atom_id, charge, pdbx_stereo_config
    });

    return Database.ofTables(
        CCA_TABLE_NAME,
        { chem_comp_atom: ccd_chemCompAtom_schema },
        { chem_comp_atom: atomTable }
    );
}

async function run(out: string, binary = false, options = DefaultDataOptions, ccaOut?: string) {
    await ensureDataAvailable(options);
    const ccd = await readCCD();
    const pvcd = await readPVCD();

    const { bonds, atoms } = await createBonds(ccd, pvcd, !!ccaOut);

    const ccbCif = getEncodedCif(CCB_TABLE_NAME, bonds, binary);
    if (!fs.existsSync(path.dirname(out))) {
        fs.mkdirSync(path.dirname(out));
    }
    writeFile(out, ccbCif);

    if (!!ccaOut) {
        const ccaCif = getEncodedCif(CCA_TABLE_NAME, atoms, binary);
        if (!fs.existsSync(path.dirname(ccaOut))) {
            fs.mkdirSync(path.dirname(ccaOut));
        }
        writeFile(ccaOut, ccaCif);
    }
}

const CCB_TABLE_NAME = 'CHEM_COMP_BONDS';
const CCA_TABLE_NAME = 'CHEM_COMP_ATOMS';

const parser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Create a cif file with one big table of all chem_comp_bond entries from the CCD and PVCD.'
});
parser.add_argument('out', {
    help: 'Generated file output path.'
});
parser.add_argument('--forceDownload', '-f', {
    action: 'store_true',
    help: 'Force download of CCD and PVCD.'
});
parser.add_argument('--binary', '-b', {
    action: 'store_true',
    help: 'Output as BinaryCIF.'
});
parser.add_argument('--ccaOut', '-a', {
    help: 'Optional generated file output path for chem_comp_atom data.',
    required: false
});
parser.add_argument('--ccdUrl', '-c', {
    help: 'Fetch the CCD from a custom URL. This forces download of the CCD.',
    required: false
});
parser.add_argument('--pvcdUrl', '-p', {
    help: 'Fetch the PVCD from a custom URL. This forces download of the PVCD.',
    required: false
});
interface Args {
    out: string,
    forceDownload?: boolean,
    binary?: boolean,
    ccaOut?: string,
    ccdUrl?: string,
    pvcdUrl?: string
}
const args: Args = parser.parse_args();

run(args.out, args.binary, { forceDownload: args.forceDownload, ccdUrl: args.ccdUrl, pvcdUrl: args.pvcdUrl }, args.ccaOut);
