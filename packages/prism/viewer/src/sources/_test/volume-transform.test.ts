/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { StateObjectCell } from '@molstar/core/state';
import { StateTransform } from '@molstar/core/state';
import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../_test/harness.js';
import { mapFile } from '../../api/_test/volume-fixtures.js';
import { isStreamedForStructure } from '../volume-transform.js';

// Which volumes move by themselves (`volume.transform`, the move handles): not density streamed for a structure, which
// hangs under it in the tree and moves with it.

/** A volume cell placed under `parent` (as streaming places its density under the structure) */
const under = (parent: string) => ({ transform: { parent } }) as unknown as StateObjectCell;

describe('isStreamedForStructure', () => {
  it('is a volume under a structure, not one loaded on its own', async () => {
    const v = await headlessViewer({ files: [cbs(), mapFile()] });
    const state = v.plugin.state.data;
    const structure = v.plugin.managers.structure.hierarchy.current.structures[0].cell.transform.ref;
    expect(isStreamedForStructure(state, under(structure))).toBe(true);
    expect(isStreamedForStructure(state, under(StateTransform.RootRef))).toBe(false);
    const own = v.plugin.managers.volume.hierarchy.current.volumes[0].cell;
    expect(isStreamedForStructure(state, own)).toBe(false);
  });

  it('leaves volume.transform moving a map loaded on its own', async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    // the map loaded on its own moves
    expect((await v.api.volume.transform({ translation: [1, 0, 0] })).ok).toBe(true);
  });
});
