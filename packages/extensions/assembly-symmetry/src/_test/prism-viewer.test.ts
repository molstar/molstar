/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import type { ApiModuleKey } from '@molstar/prism-viewer/api/area';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import type { AssemblySymmetryInfo, SymmetryInfo } from '../operations.js';
import { AssemblySymmetryExtension } from '../prism-viewer.js';
import { AssemblySymmetryData } from '../prop.js';

// Assembly symmetry in the Prism viewer on a headless viewer (parity.md §6): 1GRM (gramicidin, a dimer in assembly 1)
// and 1CBS (a monomer). RCSB PDB is never reached: `fetch` is stubbed with recorded-shape answers.

const fixture = (name: string) =>
  new File([readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url))], name);

/** What data.rcsb.org answers for 1GRM (shape of `rcsb_struct_symmetry`; values made up but plausible). */
const Gramicidin = [
  {
    kind: 'Global Symmetry',
    oligomeric_state: 'Homo 2-mer',
    stoichiometry: ['A2'],
    symbol: 'C2',
    type: 'Cyclic',
    clusters: [
      {
        avg_rmsd: 0.4,
        members: [
          { asym_id: 'A', pdbx_struct_oper_list_ids: ['1'] },
          { asym_id: 'B', pdbx_struct_oper_list_ids: ['1'] },
        ],
      },
    ],
    rotation_axes: [{ order: 2, start: [0, 0, -16], end: [0, 0, 16] }],
  },
  {
    kind: 'Local Symmetry',
    oligomeric_state: 'Monomer',
    stoichiometry: ['A'],
    symbol: 'C1',
    type: 'Asymmetric',
    clusters: [{ members: [{ asym_id: 'A' }] }],
  },
];

const Monomer = [
  {
    kind: 'Global Symmetry',
    oligomeric_state: 'Monomer',
    stoichiometry: ['A'],
    symbol: 'C1',
    type: 'Asymmetric',
    clusters: [{ members: [{ asym_id: 'A' }] }],
  },
];

let answer: (entry: string) => Response;
const fetchMock = vi.fn(async (_url: string, init?: RequestInit) => {
  const entry = JSON.parse(String(init?.body)).variables.entry_id as string;
  return answer(entry);
});

const json = (symmetry: unknown) =>
  new Response(JSON.stringify({ data: { assembly: { rcsb_struct_symmetry: symmetry } } }), { status: 200 });

beforeEach(() => {
  answer = (entry) => json(entry === '1GRM' ? Gramicidin : Monomer);
  fetchMock.mockClear();
  vi.stubGlobal('fetch', fetchMock);
});

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
  vi.unstubAllGlobals();
});

async function viewer(files: string[], modules?: readonly ApiModuleKey[]) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), AssemblySymmetryExtension] },
    ...(modules ? { modules } : {}),
  });
  viewers.push(v);
  for (const f of files) {
    const r = await v.api.load.files({ files: [fixture(f)] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

const axesCells = (v: PrismViewer) =>
  [...v.plugin.state.data.cells.values()].filter((c) =>
    c.transform.tags?.includes(AssemblySymmetryData.Tag.Representation),
  );

const themes = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures.flatMap((s) =>
    s.components.flatMap((c) =>
      c.representations.map((r) => (r.cell.transform.params as { colorTheme: { name: string } }).colorTheme.name),
    ),
  );

describe('the assembly-symmetry extension in the Prism viewer', () => {
  it('lists the symmetries RCSB PDB has, asking it once', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const r = await v.api['assembly-symmetry']!.get({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    const info = (r.ok ? r.value : undefined) as AssemblySymmetryInfo;
    expect(info.symmetries.map((s) => [s.index, s.symbol, s.kind])).toEqual([
      [1, 'C2', 'Global Symmetry'],
      [2, 'C1', 'Local Symmetry'],
    ]);
    expect(info.shown).toBeUndefined();
    expect(r.ok && r.did).toMatch(/C2 \(cyclic\), Homo 2-mer, stoichiometry A2/);
    // GraphQL to the configured server, with the entry and its assembly.
    expect(fetchMock).toHaveBeenCalledTimes(1);
    const [url, init] = fetchMock.mock.calls[0];
    expect(url).toBe(AssemblySymmetryData.DefaultServerUrl);
    expect(JSON.parse(String(init?.body)).variables).toEqual({ assembly_id: '1', entry_id: '1GRM' });
    await v.api['assembly-symmetry']!.get({});
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });

  it('shows the axes of the first symmetry that is not C1, as one step that Undo removes', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const steps = v.managers.history.steps.value.length;
    const r = await v.api['assembly-symmetry']!.show({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(r.ok && r.did).toBe("Showed the axes of 1GRM's Global Symmetry: C2 (cyclic).");
    expect((r.ok ? r.value : undefined) as SymmetryInfo).toMatchObject({ index: 1, symbol: 'C2', axes: 1 });
    const cells = axesCells(v);
    expect(cells.length).toBe(1);
    expect(cells[0].obj?.label).toBe('Global Symmetry');
    expect(cells[0].obj?.description).toBe('Cyclic (C2)');
    expect(v.managers.history.steps.value.length).toBe(steps + 1);
    // The chosen symmetry is in the structure's custom-properties node.
    const props = v.plugin.managers.structure.hierarchy.current.structures[0].properties!.cell.transform.params as {
      autoAttach: string[];
      properties: Record<string, { symmetryIndex: number }>;
    };
    expect(props.autoAttach).toContain('molstar_struct_symmetry');
    expect(props.properties.molstar_struct_symmetry.symmetryIndex).toBe(0);

    const again = await v.api['assembly-symmetry']!.show({ symmetry: 1 });
    expect(again.ok && again.unchanged).toBe(true);
    const got = await v.api['assembly-symmetry']!.get({});
    expect(got.ok && (got.value as AssemblySymmetryInfo).shown).toBe(1);

    await v.api.history.undo();
    expect(axesCells(v).length).toBe(0);
  });

  it('hides the axes, and says so when there are none', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const a = v.api['assembly-symmetry']!;
    const none = await a.hide({});
    expect(none.ok && none.unchanged).toBe(true);
    await a.show({});
    const r = await a.hide({});
    expect(r.ok && r.did).toBe('Hid the symmetry axes of 1GRM.');
    expect(axesCells(v).length).toBe(0);
  });

  it('colors the structure by symmetry cluster, as one step', async () => {
    const v = await viewer(['1grm_updated.cif']);
    const steps = v.managers.history.steps.value.length;
    const r = await v.api['assembly-symmetry']!.color({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(themes(v).length).toBeGreaterThan(0);
    expect(new Set(themes(v))).toEqual(new Set([AssemblySymmetryData.Tag.Cluster]));
    expect(v.managers.history.steps.value.length).toBe(steps + 1);
    // Coloring adds no axes.
    expect(axesCells(v).length).toBe(0);
  });

  it('refuses C1, numbers it does not have, and structures RCSB PDB does not annotate', async () => {
    const v = await viewer(['1grm_updated.cif', '1cbs_updated.cif', 'test.gro']);
    const a = v.api['assembly-symmetry']!;
    const c1 = await a.show({ symmetry: 2, structure: '1GRM' });
    expect(!c1.ok && c1.code).toBe('invalid');
    const beyond = await a.show({ symmetry: 5, structure: '1GRM' });
    expect(!beyond.ok && [beyond.code, beyond.hint]).toEqual([
      'invalid',
      'Give 1 to 2, or 0 for the first that is not C1.',
    ]);
    const monomer = await a.show({ structure: '1CBS' });
    expect(!monomer.ok && [monomer.code, monomer.message]).toEqual([
      'blocked',
      '1CBS has no symmetry but C1: there are no axes.',
    ]);
    const calls = fetchMock.mock.calls.length;
    const gro = v.plugin.managers.structure.hierarchy.current.structures[2].cell.transform.ref;
    const local = await a.show({ structure: gro });
    expect(!local.ok && local.code).toBe('blocked');
    expect(fetchMock.mock.calls.length).toBe(calls);
    expect(axesCells(v).length).toBe(0);
  });

  it('says when the service fails, and changes nothing', async () => {
    answer = () => new Response('oops', { status: 500 });
    const v = await viewer(['1grm_updated.cif']);
    const steps = v.managers.history.steps.value.length;
    const r = await v.api['assembly-symmetry']!.show({});
    expect(!r.ok && r.code).toBe('failed');
    expect(!r.ok && r.message).toMatch(/^Could not get the symmetry of 1GRM from RCSB PDB/);
    expect(axesCells(v).length).toBe(0);
    expect(v.managers.history.steps.value.length).toBe(steps);
  });

  it('needs the network, asks the agent before reaching it, and brings its commands and skill', async () => {
    const v = await viewer(['1grm_updated.cif'], ['files']);
    const r = await v.api['assembly-symmetry']!.show({});
    expect(!r.ok && r.code).toBe('unavailable');
    expect(fetchMock).not.toHaveBeenCalled();
    for (const id of ['get', 'show', 'color']) {
      const def = v.api.describe(`assembly-symmetry.${id}`)!.definition;
      expect(def.approval).toBe('network');
      expect(def.addresses?.({} as never)).toEqual([]);
    }
    expect(v.api.describe('assembly-symmetry.get')!.definition.read).toBe(true);
    expect(v.api.describe('assembly-symmetry.hide')!.definition.approval).toBeUndefined();
    for (const id of ['show', 'hide', 'color', 'get']) {
      expect(v.managers.commands.get(`assembly-symmetry.${id}`)).toBeDefined();
    }
    expect(AssemblySymmetryExtension.skills?.[0].tools.map((t) => t.operation)).toEqual([
      'assembly-symmetry.get',
      'assembly-symmetry.show',
      'assembly-symmetry.hide',
      'assembly-symmetry.color',
    ]);
  });
});
