/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { readStructure, structureOf } from '@molstar/seq-molstar-poc';
import { manualScheduler } from '@molstar/seq-viewer';
import { fixtureMeasurer } from '@molstar/seq-viewer/testing';
import { fixedInk, observableState } from '@molstar/prism-ui';
import { createSequenceViewers, type HostStructure } from '../model.js';

// Headless host for the sequence model: Seq*'s molstar-poc recorded structures (the p53 demo's PDBe BinaryCIF files),
// read from disk; no network, no new fixture.

const FIXTURES = new URL('../../../../../sequences/molstar-poc/_test/fixtures/structures/', import.meta.url);

/** A recorded structure as the host would give it: its parse and the structure its 3D view draws. */
export async function hostStructure(id: '1tup' | '2ocj' | '1c26' = '1tup'): Promise<HostStructure> {
  const bytes = new Uint8Array(readFileSync(new URL(`${id}.bcif`, FIXTURES)));
  const source = await readStructure(bytes, 'bcif');
  return { key: id, title: id.toUpperCase(), source, structure: await structureOf(source, 0, null) };
}

export function headlessViewers(background = '#ffffff') {
  const scheduler = manualScheduler();
  const ink = observableState(fixedInk(background).value);
  const structure = observableState<HostStructure | undefined>(undefined);
  const viewers = createSequenceViewers({
    ink: ink.asReadonly,
    structure: structure.asReadonly,
    scheduler,
    text: fixtureMeasurer(),
  });
  return { viewers, viewer: viewers.get(), scheduler, ink, structure };
}
