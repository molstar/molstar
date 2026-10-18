/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import path from 'node:path';
import { manualScheduler } from '@molstar/seq-viewer';
import { fixtureMeasurer } from '@molstar/seq-viewer/testing';
import { fixedInk, observableState } from '@molstar/prism-ui';
import type { SequenceBehavior } from '../behaviors.js';
import { createSequenceViewers, type HostActions, type HostStructure } from '../model.js';
import { readStructure, structureOf } from '../molstar/index.js';

// Headless host for the sequence model: the Mol* adapter's recorded structures (the p53 PDBe BinaryCIF files),
// read from disk; no network, no new fixture.

const FIXTURES = path.resolve(import.meta.dirname, '../molstar/_test/fixtures/structures');

/** A recorded structure as the host would give it: its parse and the structure its 3D view draws. */
export async function hostStructure(id: '1tup' | '2ocj' | '1c26' = '1tup'): Promise<HostStructure> {
  const bytes = new Uint8Array(readFileSync(path.join(FIXTURES, `${id}.bcif`)));
  const source = await readStructure(bytes, 'bcif');
  return { key: id, title: id.toUpperCase(), source, structure: await structureOf(source, 0, null) };
}

/** The test behavior's track: the B-factors (read as the confidence) as values, a menu and a legend block. */
export const ValuesTrack = 'test.sequence.values-track';

/**
 * A behavior as an extension's would be (sequence.md §3): on every document, a computed node and a values track after the
 * viewer's tracks; its menu colors by `plddt-confidence`, its legend one block. `patched` lists each patch's outcome.
 */
export function valuesBehavior(patched: boolean[] = []): SequenceBehavior {
  return {
    id: 'test.sequence.values',
    title: 'Values',
    attach(ctx) {
      const sub = ctx.shown.subscribe((shown) => {
        if (!shown) return;
        void ctx
          .patch(shown, [
            {
              op: 'insert',
              parent: shown.data,
              node: {
                kind: 'computed',
                ref: 'test.values',
                params: { provider: 'molstar:plddt', options: { of: shown.complex } },
              },
            },
            {
              op: 'insert',
              parent: shown.section,
              after: shown.after,
              node: {
                kind: 'track',
                ref: ValuesTrack,
                params: { label: '' },
                children: [
                  {
                    kind: 'layer',
                    ref: 'test.values.layer',
                    params: { type: 'values', data: 'test.values', encoding: 'line', domain: [0, 100] },
                  },
                ],
              },
            },
          ])
          .then((ok) => patched.push(ok));
      });
      ctx.offer({
        menus: {
          [ValuesTrack]: [
            { id: 'color-values', label: 'Color by values', icon: 'paintbrush', theme: 'plddt-confidence' },
          ],
        },
        legend: [{ title: 'Values', items: [{ label: 'High', color: '#3f6fc4' }] }],
      });
      return () => sub.unsubscribe();
    },
  };
}

export function headlessViewers(
  o: { background?: string; actions?: HostActions; behaviors?: readonly SequenceBehavior[] } = {},
) {
  const background = o.background ?? '#ffffff';
  const scheduler = manualScheduler();
  const ink = observableState(fixedInk(background).value);
  const structure = observableState<HostStructure | undefined>(undefined);
  const viewers = createSequenceViewers({
    ink: ink.asReadonly,
    structure: structure.asReadonly,
    scheduler,
    text: fixtureMeasurer(),
    ...(o.actions ? { actions: o.actions } : {}),
    ...(o.behaviors ? { behaviors: o.behaviors } : {}),
  });
  return { viewers, viewer: viewers.get(), scheduler, ink, structure };
}
