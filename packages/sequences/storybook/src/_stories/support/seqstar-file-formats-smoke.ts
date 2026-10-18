/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Seq*/File formats` stories' entries for the smoke test (seqstar-smoke.ts): every format a story offers, each a
// case, in the recorded mode only (the live mode fetches the public services; the smoke test never does). The Mol*
// adapter (`@molstar/prism-viewer/sequence/molstar`) is imported by a name held in a variable, as seqstar-molstar-smoke.ts
// does, so the smoke test loads Prism's adapter only when it runs these cases. No React.

import type { SvsRegistryEntry } from '@molstar/seq-core/format';
import ins1 from '../../../../../prism/viewer/src/sequence/molstar/_test/fixtures/structures/ins1.cif?raw';
import {
  ALIGNMENT_FORMATS,
  alignmentCase,
  archiveCase,
  DOCUMENT_FORMATS,
  documentCase,
  type FormatCase,
  genomeCase,
  pdbeAfdbDoc,
  pdbeAfdbFiles,
  recordedLoaders,
  STRUCTURE_FORMATS,
  sequencesCase,
  structuresCase,
  TABLE_FORMATS,
  tablesCase,
  tilesCase,
  uniprotDoc,
  uniprotFiles,
} from './seqstar-file-formats.js';
import type { AllowedWarning, SmokeCase, SmokeEnv, SmokeStory } from './seqstar-smoke.js';

const MOLSTAR = '@molstar/prism-viewer/sequence/molstar';
const molstar = async () => [
  { svs: ((await import(/* @vite-ignore */ MOLSTAR)) as { MolstarSvsEntry: SvsRegistryEntry }).MolstarSvsEntry },
];
const text = (b: Uint8Array) => new TextDecoder().decode(b);
// The structures are the Mol* adapter's fixtures (Prism's), the p53 data the demo's (this package's).
const fixture = (env: SmokeEnv, path: string) =>
  env.readFile(
    path.startsWith('p53/')
      ? `packages/sequences/storybook/_test/fixtures/${path}`
      : `packages/prism/viewer/src/sequence/molstar/_test/fixtures/${path}`,
  );
const of = (label: string, c: FormatCase, allow?: AllowedWarning[]): SmokeCase => ({
  label,
  doc: c.doc,
  loaders: c.loaders,
  allow,
});

export const FILE_FORMATS_SMOKE: SmokeStory[] = [
  { key: 'file-formats.stories.tsx#Sequences', cases: () => [of('FASTA, a value, residues', sequencesCase())] },
  {
    key: 'file-formats.stories.tsx#Alignments',
    cases: () => ALIGNMENT_FORMATS.map((f) => of(f, alignmentCase(f))),
  },
  { key: 'file-formats.stories.tsx#Tables', cases: () => TABLE_FORMATS.map((f) => of(f, tablesCase(f))) },
  {
    key: 'file-formats.stories.tsx#Structures',
    cases: async (env) => {
      const registry = await molstar();
      const files = { af: fixture(env, 'structures/AF-P04637-F1-model_v6.bcif'), ins1 };
      return STRUCTURE_FORMATS.map((f) => ({ ...of(f, structuresCase(f, files)), registry, fixtureStructures: false }));
    },
  },
  {
    key: 'file-formats.stories.tsx#UniProt',
    cases: (env) => {
      const files = uniprotFiles(
        text(fixture(env, 'p53/P04637.json')),
        text(fixture(env, 'p53/variation-P04637.json')),
      );
      return [{ label: 'P04637, recorded', doc: uniprotDoc(), loaders: recordedLoaders(files) }];
    },
  },
  {
    key: 'file-formats.stories.tsx#PdbeAfdb',
    cases: async (env) => {
      const files = pdbeAfdbFiles({
        pdbe: fixture(env, 'structures/2ocj.bcif'),
        afdb: fixture(env, 'structures/AF-P04637-F1-model_v6.bcif'),
      });
      return [
        {
          label: '2OCJ and AF-P04637-F1, recorded',
          doc: pdbeAfdbDoc(),
          loaders: recordedLoaders(files),
          registry: await molstar(),
          fixtureStructures: false,
        },
      ];
    },
  },
  {
    key: 'file-formats.stories.tsx#GenomeFiles',
    cases: () => [{ ...of('2bit, BigWig, BigBed, tabix VCF', genomeCase()), streams: true }],
  },
  { key: 'file-formats.stories.tsx#Tiles', cases: () => [{ ...of('seqstar-tiles/1', tilesCase()), streams: true }] },
  {
    key: 'file-formats.stories.tsx#Documents',
    cases: async (env) => [
      ...DOCUMENT_FORMATS.map((f) => of(`.${f}`, documentCase(f))),
      of('.svsx', (await archiveCase(env.readFile('packages/sequences/core/_test/format/fixtures/kinase.svsx'))).c),
    ],
  },
];
