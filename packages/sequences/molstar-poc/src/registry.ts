/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seq-molstar-poc`'s format half (seqspec:libraries#the-mol-adapter; seqspec:sources#data-by-reference): the
// structure formats, the `pdbe` and `afdb` services, the `molstar:structure` complex provider and the structure
// providers, as one registry entry.

import {
  type FormatEntry,
  type FormatResult,
  issueOf,
  type ServiceEntry,
  type SvsRegistryEntry,
} from '@molstar/seq-core/format';
import { StructureProvider } from './complex.js';
import { SiftsCoverageProvider, SiftsDepthProvider } from './coverage.js';
import { PaeFormat, PaeProvider } from './pae.js';
import { PlddtProvider } from './plddt.js';
import { SecondaryStructureProvider } from './secondary.js';
import { categoryTable, readStructure, type StructureFormat } from './source.js';

/** A structure format: the parsed structure (`object`), or with `category_name` a CIF category as a table. */
function structureFormat(name: StructureFormat, extensions: string[]): FormatEntry {
  const cif = name === 'mmcif' || name === 'bcif';
  return {
    name,
    extensions,
    gives: cif ? ['object', 'table'] : ['object'],
    reads: ['want', 'categoryName'],
    async read(bytes, request, ctx): Promise<FormatResult> {
      const source = await readStructure(bytes, name);
      if (request.categoryName == null && request.want !== 'table') return { kind: 'object', value: source };
      const table = source.frame && request.categoryName ? categoryTable(source.frame, request.categoryName) : null;
      if (!table) {
        ctx.report(
          issueOf(
            ctx.registries,
            'svs.source.category',
            '',
            `the ${name} file has no category "${request.categoryName}"`,
          ),
        );
        return { kind: 'table', columns: {}, length: 0 };
      }
      return { kind: 'table', ...table, text: true };
    },
  };
}

export const MolstarFormats: readonly FormatEntry[] = [
  structureFormat('mmcif', ['cif', 'mmcif']),
  structureFormat('bcif', ['bcif']),
  structureFormat('pdb', ['pdb', 'ent']),
  structureFormat('sdf', ['sdf', 'mol']),
  structureFormat('mol2', ['mol2']),
  PaeFormat,
];

/**
 * PDBe's BinaryCIF carries the updated mmCIF's SIFTS (`structure`); `sifts` is the updated mmCIF itself. AlphaFold DB's
 * models by entry id (`AF-P04637-F1`, the current version), or with a version (`AF-P04637-F1-model_v4`); `pae` the
 * model's predicted aligned error (pae.ts).
 */
export const MolstarServices: readonly ServiceEntry[] = [
  {
    name: 'pdbe',
    resources: ['structure', 'sifts'],
    resolve: (id, resource) => {
      const pdb = id.toLowerCase();
      if (!/^[0-9][a-z0-9]{3}$/.test(pdb)) return null;
      return resource === 'sifts'
        ? { uri: `https://www.ebi.ac.uk/pdbe/entry-files/download/${pdb}_updated.cif`, format: 'mmcif' }
        : { uri: `https://www.ebi.ac.uk/pdbe/entry-files/download/${pdb}.bcif`, format: 'bcif' };
    },
  },
  {
    name: 'afdb',
    resources: ['structure', 'pae'],
    resolve: (id, resource) => {
      if (!/^AF-[A-Z0-9-]+-F\d+(-model_v\d+)?$/.test(id)) return null;
      const file = id.includes('-model_v') ? id : `${id}-model_v6`;
      if (resource === 'pae')
        return {
          uri: `https://alphafold.ebi.ac.uk/files/${file.replace('-model_', '-predicted_aligned_error_')}.json`,
          format: 'afdb_pae',
        };
      return { uri: `https://alphafold.ebi.ac.uk/files/${file}.bcif`, format: 'bcif' };
    },
  },
];

/** The format half: register it beside `KernelRegistry` and `MvsRegistry` (`register(reg, [...])`). */
export const MolstarSvsEntry: SvsRegistryEntry = {
  formats: MolstarFormats,
  services: MolstarServices,
  complexProviders: [StructureProvider],
  providers: [PlddtProvider, PaeProvider, SecondaryStructureProvider, SiftsCoverageProvider, SiftsDepthProvider],
  issueCodes: [
    { code: 'molstar.unit_name.collision', severity: 'error' },
    { code: 'molstar.ins_code.invalid', severity: 'warning' },
  ],
};
