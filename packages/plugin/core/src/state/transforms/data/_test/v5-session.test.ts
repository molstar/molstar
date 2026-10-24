/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it, vi } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultFormats } from '@molstar/plugin/default-registry';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';

// A session saved by Mol* 5.13.1 (`make-v5-session.cjs`): its `Download`s and `ReadFile` have no `asResponse`, so they
// load as Response Data (`.v6/data-response`), read by v5's own parsers from the session's assets.
const session = readFileSync(new URL('./v5-session.molx', import.meta.url));

describe('a 5.x session', () => {
  it('loads from its assets: text and gzipped downloads, a PDB and an opened file', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch');
    const plugin = new PluginContext({ behaviors: [], registry: [DefaultFormats] });
    await plugin.init();
    await plugin.managers.snapshot.open(new File([session], 'v5.molx'));

    const loaded: { from: string; data: string; atoms: number }[] = [];
    const cells = plugin.state.data.cells;
    cells.forEach((cell) => {
      if (!SO.Molecule.Trajectory.is(cell.obj)) return;
      // the data cell the trajectory was read from, past its parse (CIF) and inflate (DeflateData) cells
      let data = cells.get(cell.transform.parent)!;
      while (!(data.obj instanceof SO.Data.Response || data.obj instanceof SO.Data.String)) {
        data = cells.get(data.transform.parent)!;
      }
      const atoms = cell.obj.data.representative.atomicHierarchy.atoms._rowCount;
      loaded.push({ from: data.transform.transformer.id, data: data.obj.type.name, atoms });
    });
    expect(loaded).toEqual([
      { from: 'ms-plugin.download', data: 'Response Data', atoms: 2 },
      { from: 'ms-plugin.download', data: 'Response Data', atoms: 1 },
      // a gzip served without Content-Encoding: the download inflates it, DeflateData passes the text through
      { from: 'ms-plugin.defalate-data', data: 'String Data', atoms: 2 },
      { from: 'ms-plugin.read-file', data: 'Response Data', atoms: 2 },
    ]);
    expect(fetchSpy).not.toHaveBeenCalled();
    plugin.dispose();
  });
});
