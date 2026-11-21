/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The services of the showcase (examples.md, Showcases: an application composing services): UniProt, PDBe (SIFTS and its
// entry files) and AlphaFold DB, each by the location its service answers at, the hosts a live read may reach, the
// licence its answers carry, and the recorded file that answers the same location in the fixture mode. A spec names the
// service's own location (never a fixture's): the mode is the application's port, so the spec reads the live service
// and the recorded file alike (a reader without Seq* can take the spec and the service).

/** A service: what an application composes and a policy lets through. */
export interface Service {
  readonly id: 'uniprot' | 'pdbe' | 'alphafold';
  readonly name: string;
  /** The hosts a live read may reach (the asset policy's `hosts`). */
  readonly hosts: readonly string[];
  readonly licence: string;
  readonly attribution: string;
}

export const UniProt: Service = {
  id: 'uniprot',
  name: 'UniProt',
  hosts: ['rest.uniprot.org'],
  licence: 'CC-BY-4.0',
  attribution: 'UniProt Consortium, UniProtKB',
};
export const Pdbe: Service = {
  id: 'pdbe',
  name: 'PDBe (SIFTS)',
  hosts: ['www.ebi.ac.uk'],
  licence: 'LicenseRef-EMBL-EBI-terms',
  attribution: 'PDBe, SIFTS',
};
export const AlphaFoldDb: Service = {
  id: 'alphafold',
  name: 'AlphaFold DB',
  hosts: ['alphafold.ebi.ac.uk'],
  licence: 'CC-BY-4.0',
  attribution: 'AlphaFold DB',
};

export const Services: readonly Service[] = [UniProt, Pdbe, AlphaFoldDb];

/** The hosts of every service: what the live mode's policy allows. */
export const ServiceHosts: readonly string[] = [...new Set(Services.flatMap((s) => s.hosts))];

/** The locations the services answer at. */
export const serviceUrl = {
  /** UniProtKB's REST entry. */
  uniprot: (accession: string) => `https://rest.uniprot.org/uniprotkb/${accession}.json`,
  /** SIFTS's mappings of a PDB entry's chains to UniProt (the entry's id in lower case). */
  sifts: (entry: string) => `https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/${entry.toLowerCase()}`,
  /** PDBe's mmCIF of an entry (the id in lower case). */
  entry: (entry: string) => `https://www.ebi.ac.uk/pdbe/entry-files/download/${entry.toLowerCase()}.cif`,
  /** AlphaFold DB's per-residue confidence (model version 6). */
  confidence: (accession: string) => `https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-confidence_v6.json`,
} as const;

/**
 * What the recorded fixtures answer, by location (relative to the examples package): the entries the showcase composes.
 * PDBe's mmCIF of 2OCJ is answered by the examples' committed file of the same entry (RCSB's copy, the same model).
 */
export const Recorded: ReadonlyMap<string, string> = new Map<string, string>([
  [serviceUrl.uniprot('P04637'), 'fixtures/p53/P04637.json'],
  [serviceUrl.uniprot('P24941'), 'fixtures/complex/P24941.json'],
  [serviceUrl.uniprot('P20248'), 'fixtures/complex/P20248.json'],
  [serviceUrl.sifts('2OCJ'), 'fixtures/services/sifts-2ocj.json'],
  [serviceUrl.sifts('1FIN'), 'fixtures/services/sifts-1fin.json'],
  [serviceUrl.entry('2OCJ'), 'fixtures/structures/2OCJ.cif.gz'],
  [serviceUrl.confidence('P04637'), 'fixtures/services/AF-P04637-F1-confidence_v6.json'],
  [serviceUrl.confidence('P24941'), 'fixtures/complex/AF-P24941-F1-confidence_v6.json'],
  [serviceUrl.confidence('P20248'), 'fixtures/complex/AF-P20248-F1-confidence_v6.json'],
]);

/** The service a location belongs to, by its host. */
export function serviceOf(url: string): Service | undefined {
  let host: string;
  try {
    host = new URL(url).host;
  } catch {
    return undefined;
  }
  return Services.find((s) => s.hosts.includes(host));
}
