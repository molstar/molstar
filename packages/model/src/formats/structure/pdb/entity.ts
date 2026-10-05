/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ryan DiRisio <rjdiris@gmail.com>
 */

import { CifCategory, CifField } from '@molstar/io/reader/cif';
import type { Tokens } from '@molstar/io/reader/common/text/tokenizer';
import type { EntityCompound } from '../common/entity.js';
import type { mmCIF_Schema } from '@molstar/io/reader/cif/schema/mmcif';

const Spec = {
  MOL_ID: '',
  MOLECULE: '',
  CHAIN: '',
  FRAGMENT: '',
  SYNONYM: '',
  EC: '',
  ENGINEERED: '',
  MUTATION: '',
  OTHER_DETAILS: '',
};
type Spec = keyof typeof Spec;

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

    const cmpnd = line.substring(10, 80).trim();
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
        description: '',
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
        chains: [chain],
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

    const het = line.substring(11, 14).trim();
    const name = line.substring(15).trim(); // support any length

    if (hetnams.has(het)) {
      hetnams.set(het, `${hetnams.get(het)!} ${name}`);
    } else {
      hetnams.set(het, name);
    }
  }

  return hetnams;
}

export function parseSeqres(lines: Tokens, lineStart: number, lineEnd: number) {
  const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);

  const seqresMap = new Map<string, string[]>();

  for (let i = lineStart; i < lineEnd; i++) {
    const line = getLine(i);
    // COLUMNS        DATA TYPE      FIELD        DEFINITION
    // -------------------------------------------------------------------------------------
    // 1 -  6        Record name    "SEQRES"
    // 8 - 10        Integer        serNum       Serial number of the SEQRES record for  the
    //                                         current  chain. Starts at 1 and increments
    //                                         by one  each line. Reset to 1 for each chain.
    // 12             Character      chainID      Chain identifier. This may be any single
    //                                         legal  character, including a blank which is
    //                                         is  used if there is only one chain.
    // 14 - 17        Integer        numRes       Number of residues in the chain.
    //                                         This  value is repeated on every record.
    // 20 - 22        Residue name   resName      Residue name.
    // 24 - 26        Residue name   resName      Residue name.
    // 28 - 30        Residue name   resName      Residue name.
    // 32 - 34        Residue name   resName      Residue name.
    // 36 - 38        Residue name   resName      Residue name.
    // 40 - 42        Residue name   resName      Residue name.
    // 44 - 46        Residue name   resName      Residue name.
    // 48 - 50        Residue name   resName      Residue name.
    // 52 - 54        Residue name   resName      Residue name.
    // 56 - 58        Residue name   resName      Residue name.
    // 60 - 62        Residue name   resName      Residue name.
    // 64 - 66        Residue name   resName      Residue name.
    // 68 - 70        Residue name   resName      Residue name.
    const chainId = line.substring(11, 12);
    const residues = line.substring(19).trim().split(/\s+/); // support any number
    if (!seqresMap.has(chainId)) {
      seqresMap.set(chainId, []);
    }
    seqresMap.get(chainId)!.push(...residues);
  }

  return seqresMap;
}

export function getEntityPolySeq(
  seqresMap: Map<string, string[]>,
  getEntityIdForChain: (chainId: string) => string | undefined,
): CifCategory | undefined {
  const epsEntityIds: string[] = [];
  const epsNums: number[] = [];
  const epsMonIds: string[] = [];
  const epsHeteros: mmCIF_Schema['entity_poly_seq']['hetero']['T'][] = [];
  const processedEntities = new Set<string>();

  for (const [chainId, residues] of seqresMap) {
    const entityId = getEntityIdForChain(chainId);
    if (!entityId || processedEntities.has(entityId)) continue;
    processedEntities.add(entityId);

    for (let j = 0; j < residues.length; j++) {
      epsEntityIds.push(entityId);
      epsNums.push(j + 1);
      epsMonIds.push(residues[j]);
      epsHeteros.push('no');
    }
  }

  if (epsEntityIds.length > 0) {
    const entity_poly_seq: CifCategory.SomeFields<mmCIF_Schema['entity_poly_seq']> = {
      entity_id: CifField.ofStrings(epsEntityIds),
      num: CifField.ofNumbers(epsNums),
      mon_id: CifField.ofStrings(epsMonIds),
      hetero: CifField.ofStrings(epsHeteros),
    };
    return CifCategory.ofFields('entity_poly_seq', entity_poly_seq);
  }

  return undefined;
}

/**
 * Build pdbx_unobs_or_zero_occ_residues by comparing SEQRES with observed ATOM records
 * Collect observed (label_asym_id, label_seq_id) pairs per model, and auth_asym_id -> label_asym_id mapping.
 * Only include atoms belonging to the polymer entity for each SEQRES chain.
 * Non-polymer HETATMs (ions, ligands, water) share auth_asym_id in PDB format
 * and their sequential label_seq_id values can collide with unobserved SEQRES positions.
 */
export function getPdbxUnobsOrZeroOccResidues(
  seqresMap: Map<string, string[]>,
  getEntityIdForChain: (chainId: string) => string | undefined,
  atom_site: CifCategory.SomeFields<mmCIF_Schema['atom_site']>,
  modelCount: number,
) {
  const polymerEntityIds = new Map<string, string>();
  for (const chainId of seqresMap.keys()) {
    const entityId = getEntityIdForChain(chainId);
    if (entityId) polymerEntityIds.set(chainId, entityId);
  }

  const observedResidues = new Set<string>();
  const authToLabelAsym = new Map<string, string>();
  const rowCount = atom_site.label_asym_id!.rowCount;
  for (let i = 0; i < rowCount; ++i) {
    const authAsym = atom_site.auth_asym_id!.str(i);
    const entityId = atom_site.label_entity_id!.str(i);

    // Skip non-polymer atoms: their label_seq_id can collide with SEQRES positions
    const polymerEntityId = polymerEntityIds.get(authAsym);
    if (polymerEntityId && entityId !== polymerEntityId) continue;

    const labelAsym = atom_site.label_asym_id!.str(i);
    const labelSeq = atom_site.label_seq_id!.int(i);
    const modelN = atom_site.pdbx_PDB_model_num!.int(i);
    observedResidues.add(`${modelN}|${labelAsym}|${labelSeq}`);
    if (!authToLabelAsym.has(authAsym)) {
      authToLabelAsym.set(authAsym, labelAsym);
    }
  }

  const unobsIds: number[] = [];
  const unobsModelNums: number[] = [];
  const unobsLabelAsymIds: string[] = [];
  const unobsLabelCompIds: string[] = [];
  const unobsLabelSeqIds: number[] = [];
  const unobsAuthAsymIds: string[] = [];
  const unobsAuthCompIds: string[] = [];
  const unobsAuthSeqIds: number[] = [];
  const unobsPolymerFlags: mmCIF_Schema['pdbx_unobs_or_zero_occ_residues']['polymer_flag']['T'][] = [];
  const unobsOccupancyFlags: mmCIF_Schema['pdbx_unobs_or_zero_occ_residues']['occupancy_flag']['T'][] = [];
  let unobsId = 0;
  modelCount = Math.max(modelCount, 1);

  for (let m = 1; m <= modelCount; ++m) {
    for (const [chainId, residues] of seqresMap) {
      const labelAsymId = authToLabelAsym.get(chainId);
      if (!labelAsymId) continue;

      for (let j = 0; j < residues.length; j++) {
        const seqId = j + 1; // 1-based label_seq_id
        if (!observedResidues.has(`${m}|${labelAsymId}|${seqId}`)) {
          unobsId++;
          unobsIds.push(unobsId);
          unobsModelNums.push(m);
          unobsLabelAsymIds.push(labelAsymId);
          unobsLabelCompIds.push(residues[j]);
          unobsLabelSeqIds.push(seqId);
          unobsAuthAsymIds.push(chainId);
          unobsAuthCompIds.push(residues[j]);
          unobsAuthSeqIds.push(seqId);
          unobsPolymerFlags.push('y');
          unobsOccupancyFlags.push(1);
        }
      }
    }
  }

  if (unobsIds.length > 0) {
    const pdbx_unobs: CifCategory.SomeFields<mmCIF_Schema['pdbx_unobs_or_zero_occ_residues']> = {
      id: CifField.ofNumbers(unobsIds),
      PDB_model_num: CifField.ofNumbers(unobsModelNums),
      polymer_flag: CifField.ofStrings(unobsPolymerFlags),
      occupancy_flag: CifField.ofNumbers(unobsOccupancyFlags),
      label_asym_id: CifField.ofStrings(unobsLabelAsymIds),
      label_comp_id: CifField.ofStrings(unobsLabelCompIds),
      label_seq_id: CifField.ofNumbers(unobsLabelSeqIds),
      auth_asym_id: CifField.ofStrings(unobsAuthAsymIds),
      auth_comp_id: CifField.ofStrings(unobsAuthCompIds),
      auth_seq_id: CifField.ofNumbers(unobsAuthSeqIds),
    };
    return CifCategory.ofFields('pdbx_unobs_or_zero_occ_residues', pdbx_unobs);
  }

  return undefined;
}
