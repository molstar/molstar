/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { toJSON } from '@molstar/seq-core/format';
import { Asset } from '@molstar/core/util/assets';
import { ReadFile } from '@molstar/plugin/state/transforms/data/fetch';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, until } from '../../_test/harness.js';
import { LoadersFeature } from '../../features/loaders/index.js';
import type { PrismViewer } from '../../viewer.js';

// The sequence's parse cache (sequence.md §2): `dataByUrl` indexes the plugin's data cells by address with their
// parses; a loaded SeqViewSpec reads a file the plugin holds from it (no refusal, no second read) and shares a structure
// the plugin parsed (no second parse). Headless, over 1CBS and 1CRN.

const seqOf = (v: PrismViewer) => v.views.sequence!.viewers.value!.get();
const trajectoryOf = (v: PrismViewer, title: string) =>
  v.plugin.managers.structure.hierarchy.current.structures.find(
    (s) => s.cell.obj?.data.model.entryId?.toLowerCase() === title,
  )?.model?.trajectory?.cell.obj?.data;

/** A plain-text file read into the plugin as a data cell (as a person's dropped file the 3D view cannot show). */
async function readText(v: PrismViewer, name: string, text: string) {
  await v.plugin
    .build()
    .toRoot()
    .apply(ReadFile, { file: Asset.File(new File([text], name)), isBinary: false })
    .commit();
}

describe('dataByUrl', () => {
  it("indexes the plugin's data cells and their parses by address, rebuilt after a change", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const index = v.views.sequence!.source.dataByUrl;
    const entry = index.entry('1cbs.cif');
    expect(entry?.data).toBeTruthy();
    expect(entry?.parse).toBeTruthy();
    expect(new TextDecoder().decode(index.bytes('1cbs.cif')).startsWith('data_1CBS')).toBe(true);
    const source = index.structure('1cbs.cif');
    expect(source?.format).toBe('mmcif');
    expect(source?.trajectory).toBe(trajectoryOf(v, '1cbs'));
    // Its CIF block, for the categories read beside the structure.
    expect(source?.frame?.header).toBe('1CBS');
    // One source per trajectory: its structures are made once whoever asks.
    expect(index.structure('1cbs.cif')).toBe(source);
    expect(index.entry('1crn.cif')).toBeUndefined();
    expectOk(await v.api.load.files({ files: [crambin()] }));
    expect(index.structure('1crn.cif')?.trajectory).toBe(trajectoryOf(v, '1crn'));
    expectOk(await v.api.session.clear());
    expect(index.entry('1cbs.cif')).toBeUndefined();
  });

  it("a loaded SeqViewSpec reads a file the plugin holds, and shares a second structure's parse", async () => {
    const v = await headlessViewer({ extensions: [LoadersFeature], files: [cbs(), crambin()] });
    await readText(v, 'p.fasta', '>p\nMKVLATGHEK\n');
    const b = createSeqViewSpec({ title: 'Two structures' });
    b.data.complex({ ref: 'cbs', provider: 'molstar:structure', source: { host_file: '1cbs.cif', format: 'mmcif' } });
    b.data.complex({ ref: 'crn', provider: 'molstar:structure', source: { host_file: '1crn.cif', format: 'mmcif' } });
    b.data.sequence({ ref: 'p', source: { host_file: 'p.fasta', format: 'fasta' } } as never);
    const section = b.view({ ref: 'svs-view', label: 'Two' }).section({ ref: 'svs-main', style: 'compact' });
    section.block({ complex: 'crn' } as never);
    section.track({ ref: 'svs-letters', label: 'Residues' }).letters({});
    expectOk(await v.api.load.svs!({ data: toJSON(b.getDocument()) }));
    await until(() => seqOf(v).structure.value?.complex === 'cbs');
    const seq = seqOf(v).seq;
    await until(
      () =>
        !!seq.graph.assembled.value?.set.molecules.some(
          (m) => m && seq.graph.assembled.value?.set.strings.get(m.key) === 'p',
        ),
    );
    const set = seq.graph.assembled.value!.set;
    const keys = (xs: readonly ({ key: number } | null | undefined)[]) =>
      xs.flatMap((x) => (x ? [set.strings.get(x.key)] : []));
    expect(keys(set.complexes)).toEqual(['cbs', 'crn']);
    expect(set.molecules.find((m) => m && set.strings.get(m.key) === 'p')?.length).toBe(10);
    expect(seq.issues.current.value.filter((i) => i.code === 'svs.source.load')).toEqual([]);
  });

  it('a file the plugin does not hold is refused, and says to load it first', async () => {
    const v = await headlessViewer({ extensions: [LoadersFeature], files: [cbs()] });
    const b = createSeqViewSpec({ title: 'Missing' });
    b.data.complex({ ref: 'cbs', provider: 'molstar:structure', source: { host_file: '1cbs.cif', format: 'mmcif' } });
    b.data.sequence({ ref: 'p', source: { host_file: 'missing.fasta', format: 'fasta' } } as never);
    const section = b.view({ ref: 'svs-view', label: 'Missing' }).section({ ref: 'svs-main', style: 'compact' });
    section.track({ ref: 'svs-letters', label: 'Residues' }).letters({});
    expectOk(await v.api.load.svs!({ data: toJSON(b.getDocument()) }));
    const seq = () => seqOf(v).seq;
    await until(() => seq().issues.current.value.some((i) => i.code === 'svs.source.load'));
    const issue = seq().issues.current.value.find((i) => i.code === 'svs.source.load');
    expect(issue?.message).toContain('load missing.fasta in the viewer first');
  });
});
