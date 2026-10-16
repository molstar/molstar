/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Seq*/Mol* › p53` story's document for the smoke test (seqstar-smoke.ts): the demo's document over the committed
// fixtures, read through the test's `readFile`, with the Mol* adapter's format half registered (Prism's,
// `@molstar/prism-viewer/sequence/molstar`). No React. The adapter is imported by a name held in a variable, so the
// smoke test loads it only when it runs this case.

import type { SvsRegistryEntry } from '@molstar/seq-core/format';
import { P53_STRUCTURES, p53DemoDoc, p53Loaders } from './seqstar-p53-demo.js';
import type { SmokeCase, SmokeEnv } from './seqstar-smoke.js';

const MOLSTAR = '@molstar/prism-viewer/sequence/molstar';
interface Molstar {
  readonly MolstarSvsEntry: SvsRegistryEntry;
  structureUrl(ref: { source: { service: string; id: string } }): { url: string } | null;
}
const text = (b: Uint8Array) => new TextDecoder().decode(b);
const FILES: Record<string, string> = { 'AF-P04637-F1': 'AF-P04637-F1-model_v6' };

export async function p53MolstarCase(env: SmokeEnv): Promise<SmokeCase> {
  const { MolstarSvsEntry, structureUrl } = (await import(/* @vite-ignore */ MOLSTAR)) as Molstar;
  const at = (p: string) => env.readFile(`packages/sequences/core/_test/fixtures/${p}`);
  const sequence = (JSON.parse(text(at('p53.json'))) as { sequence: string }).sequence;
  const structures = new Map<string, Uint8Array>();
  for (const s of P53_STRUCTURES) {
    const url = structureUrl({ source: s.source })?.url;
    if (url)
      structures.set(
        url,
        env.readFile(
          `packages/prism/viewer/src/sequence/molstar/_test/fixtures/structures/${FILES[s.source.id] ?? s.source.id}.bcif`,
        ),
      );
  }
  return {
    label: 'p53, its features and five structures',
    doc: p53DemoDoc(sequence),
    loaders: p53Loaders({
      structures,
      fasta: `>sp|P04637|P53_HUMAN\n${sequence}\n`,
      features: text(env.readFile('packages/sequences/storybook/_test/fixtures/p53/P04637.json')),
      variants: text(env.readFile('packages/sequences/storybook/_test/fixtures/p53/variation-P04637.json')),
      ptms: text(env.readFile('packages/sequences/storybook/_test/fixtures/p53/P04637-ptms.tsv')),
    }),
    registry: [{ svs: MolstarSvsEntry }],
    fixtureStructures: false,
  };
}
