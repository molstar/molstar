/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { toJSON, zipFiles } from '@molstar/seq-core/format';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { cbs, crambin, exampleFile, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { LoadersFeature } from '../index.js';

// `load.svs` (sequence.md §6, David 2026-10-15): a SeqViewSpec replaces the followed structure's document until the
// person focuses or loads another structure; its first `molstar:structure` complex maps onto the matching loaded
// structure (loaded when none matches, in the same step); it is saved with the sequence's state. Headless, over 1CBS.

/** A small SeqViewSpec over one structure: a section of its own with a letters track (the fixture). */
function svsOver(source: Record<string, unknown>, title = 'Retinoic acid binding') {
  const b = createSeqViewSpec({ title });
  b.data.complex({ ref: 'cbs', provider: 'molstar:structure', source: source as never });
  const section = b.view({ ref: 'svs-view', label: title }).section({ ref: 'svs-main', style: 'compact' });
  section.track({ ref: 'svs-letters', label: 'Residues' }).letters({});
  return toJSON(b.getDocument());
}
const overFile = (title?: string) => svsOver({ host_file: '1cbs.cif', format: 'mmcif' }, title);

afterEach(() => {
  vi.unstubAllGlobals();
});

const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const structures = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures;
// Loaded by then (a structure shown, or load.svs).
const seqOf = (v: PrismViewer) => v.views.sequence!.viewers.value!.get();
/** The SVS is shown: the sequence's document is it (its section), over the structure the link maps. */
const loaded = (v: PrismViewer) => !!v.views.sequence!.viewers.value;
const showsSvs = (v: PrismViewer) =>
  loaded(v) &&
  seqOf(v).structure.value?.complex === 'cbs' &&
  !!seqOf(v).seq.document.base.value?.root.children?.some((n) => n.kind === 'view' && n.ref === 'svs-view');
const showsBuilt = (v: PrismViewer, title: string) =>
  loaded(v) && seqOf(v).structure.value?.title === title && seqOf(v).structure.value?.complex === 'structure';

async function viewer(files: File[] = [cbs()]) {
  const v = await headlessViewer({ extensions: [LoadersFeature], files });
  v.plugin.state.setSnapshotParams({ camera: false });
  return v;
}

describe('load.svs', () => {
  it('shows a SeqViewSpec over the loaded structure it names, as one step, and Undo takes it back', async () => {
    const v = await viewer();
    await until(() => showsBuilt(v, '1CBS'));
    const before = steps(v);
    const r = await v.api.load.svs!({ data: overFile() });
    expectOk(r);
    expect(r.did).toBe('The sequence shows Retinoic acid binding over 1CBS.');
    expect(structures(v)).toHaveLength(1);
    expect(steps(v)).toBe(before + 1);
    await until(() => showsSvs(v));
    // The source reads the 3D view's parse (the host file the loaders' `parsed` answers).
    const doc = seqOf(v).seq.document.base.value!;
    expect(JSON.stringify(doc)).toContain('"host_file":"structure:');
    const got = await v.api.sequence!.get();
    expectOk(got);
    expect(got.value).toMatchObject({ structure: '1CBS', document: 'Retinoic acid binding' });

    expectOk(await v.api.history.undo());
    await until(() => showsBuilt(v, '1CBS'));
    expect(seqOf(v).override.value).toBeUndefined();
  });

  it("applies the person's style and letter colors over the document's", async () => {
    const v = await viewer();
    expectOk(await v.api.sequence!.colors({ scheme: 'taylor' }));
    expectOk(await v.api.load.svs!({ data: overFile() }));
    await until(() => showsSvs(v));
    const letters = JSON.stringify(seqOf(v).seq.document.base.value);
    expect(letters).toContain('"name":"taylor"');
    // The person's Tracks over the document's own `compact`; then Text.
    expect(seqOf(v).seq.ui.style.value.get('svs-main')).toBe('prism_tracks');
    expectOk(await v.api.sequence!.style({ style: 'text' }));
    expect(seqOf(v).seq.ui.style.value.get('svs-main')).toBe('compact');
  });

  it('ends when the person focuses another structure', async () => {
    const v = await viewer([cbs(), crambin()]);
    expectOk(await v.api.load.svs!({ data: overFile() }));
    await until(() => showsSvs(v));
    expectOk(await v.api.focus.set({ target: 'polymer', structure: '1CRN' }));
    await until(() => showsBuilt(v, '1CRN'));
    expect(seqOf(v).override.value).toBeUndefined();
  });

  it('ends when another structure loads; Undo of that load brings it back', async () => {
    const v = await viewer();
    expectOk(await v.api.load.svs!({ data: overFile() }));
    await until(() => showsSvs(v));
    expectOk(await v.api.load.files({ files: [crambin()] }));
    await until(() => seqOf(v).override.value === undefined && showsBuilt(v, '1CBS'));

    expectOk(await v.api.history.undo());
    expect(structures(v)).toHaveLength(1);
    await until(() => showsSvs(v));
    expect(seqOf(v).override.value?.document.metadata?.title).toBe('Retinoic acid binding');
  });

  it('is kept by a session (version 4 of the sequence state) and restored with it', async () => {
    const v = await viewer();
    expectOk(await v.api.load.svs!({ data: overFile() }));
    await until(() => showsSvs(v));
    const saved = v.plugin.state.getSnapshot({ camera: false }).custom?.['prism-viewer.sequence'] as {
      version: number;
      override?: { structure: string };
    };
    expect(saved.version).toBe(4);
    expect(saved.override?.structure).toBe(structures(v)[0].cell.transform.ref);
    const molx = new Uint8Array(await (await v.plugin.managers.snapshot.serialize({ type: 'molx' })).arrayBuffer());

    const w = await viewer([]);
    await w.plugin.managers.snapshot.open(new File([molx], 'state.molx'));
    await until(() => structures(w).length === 1);
    await until(() => showsSvs(w));
  });

  it('reads version 1 of the sequence state as version 2 without a document', async () => {
    const v = await viewer();
    expectOk(await v.api.load.svs!({ data: overFile() }));
    await until(() => showsSvs(v));
    const snapshot = v.plugin.state.getSnapshot({ camera: false });
    const { override: _, ...rest } = snapshot.custom!['prism-viewer.sequence'] as Record<string, unknown>;
    await v.plugin.state.setSnapshot({
      ...snapshot,
      custom: { ...snapshot.custom, 'prism-viewer.sequence': { ...rest, version: 1, colors: 'zappo' } },
    });
    expect(seqOf(v).theme.value.letters).toBe('zappo');
    await until(() => showsBuilt(v, '1CBS'));
  });

  it('loads the structure the document names when none matches, in the same step; the next one matches it', async () => {
    const stub = vi.fn(async (input: string | URL | Request) => {
      const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
      if (url.includes('files.test/1cbs.cif')) return new Response(await cbs().arrayBuffer(), { status: 200 });
      if (url.includes('files.test/doc.svsj'))
        return new Response(svsOver({ uri: '1cbs.cif', format: 'mmcif' }), { status: 200 });
      return new Response('not found', { status: 404 });
    });
    vi.stubGlobal('fetch', stub);
    const v = await viewer([]);
    const r = await v.api.load.svs!({ url: 'https://files.test/doc.svsj' });
    expectOk(r);
    expect(r.value?.loaded).toEqual(['1CBS']);
    expect(steps(v)).toBe(1);
    await until(() => showsSvs(v));
    // The relative URI resolved against the document's address.
    expect(stub.mock.calls.map((c) => String(c[0]))).toContain('https://files.test/1cbs.cif');

    expectOk(
      await v.api.load.svs!({ data: svsOver({ uri: 'https://files.test/1cbs.cif', format: 'mmcif' }, 'Again') }),
    );
    expect(structures(v)).toHaveLength(1);
    await until(() => seqOf(v).override.value?.document.metadata?.title === 'Again');
  });

  it('matches an AlphaFold DB entry it names by accession: loaded by load.entry, or a file of its model', async () => {
    const stub = vi.fn(async (input: string | URL | Request) => {
      const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
      if (url.includes('alphafold.ebi.ac.uk/api/prediction/Q8W3K0'))
        return new Response(JSON.stringify([{ bcifUrl: 'https://af.test/AF-Q8W3K0-F1.bcif' }]), { status: 200 });
      if (url.includes('af.test/'))
        return new Response(await exampleFile('1cbs_full.bcif').arrayBuffer(), { status: 200 });
      return new Response('not found', { status: 404 });
    });
    vi.stubGlobal('fetch', stub);
    const v = await viewer([]);
    expectOk(await v.api.load.entry!({ source: 'afdb', id: 'Q8W3K0' }));
    const requests = stub.mock.calls.length;
    // The download's label is load.entry's identity: the structure is that one, not loaded again.
    const r = await v.api.load.svs!({ data: svsOver({ service: 'afdb', id: 'Q8W3K0' }) });
    expectOk(r);
    expect(structures(v)).toHaveLength(1);
    expect(stub.mock.calls.length).toBe(requests);
    await until(() => showsSvs(v));

    // A model's file, whose entry id is AF-<accession>-F1.
    const text = (await cbs().text()).replace(/1CBS/g, 'AF-Q8W3K0-F1');
    const w = await viewer([new File([text], 'AF-Q8W3K0-F1.cif')]);
    expectOk(await w.api.load.svs!({ data: svsOver({ service: 'afdb', id: 'Q8W3K0' }) }));
    expect(structures(w)).toHaveLength(1);
    expect(stub.mock.calls.length).toBe(requests);
  });

  it('opens .svsj files with the structure they name, dropped together', async () => {
    const v = await viewer([]);
    const doc = new File([overFile()], 'doc.svsj');
    const r = await v.api.load.files({ files: [cbs(), doc] });
    expectOk(r);
    expect(r.did).toBe('Opened 1cbs.cif (1 structure); the sequence shows doc.svsj.');
    await until(() => showsSvs(v));
  });

  it('reads an archive (.svsx): its document, and the structure it holds, opened as a file', async () => {
    const enc = new TextEncoder();
    /** An archive with the document at `index.svsj` and 1CBS at `data/1cbs.cif` (a relative URI in it). */
    const archive = async (name = 'cbs.svsx') =>
      new File(
        [
          (await zipFiles({
            'index.svsj': enc.encode(svsOver({ uri: 'data/1cbs.cif', format: 'mmcif' })),
            'data/1cbs.cif': new Uint8Array(await cbs().arrayBuffer()),
          })) as Uint8Array<ArrayBuffer>,
        ],
        name,
      );
    const v = await viewer([]);
    const r = await v.api.load.svs!({ file: await archive() });
    expectOk(r);
    expect(r.value?.loaded).toEqual(['1CBS']);
    expect(steps(v)).toBe(1);
    await until(() => showsSvs(v));
    // The override names the archive its relative URIs resolve in (kept in memory, out of the saved state).
    expect(seqOf(v).override.value?.archive).toMatch(/^svsx-/);

    // Dropped with nothing else: the same.
    const w = await viewer([]);
    const dropped = await w.api.load.files({ files: [await archive('dropped.svsx')] });
    expectOk(dropped);
    await until(() => showsSvs(w));

    // A zip that is not an archive of ours.
    const junk = new File([(await zipFiles({ 'a.txt': enc.encode('x') })) as Uint8Array<ArrayBuffer>], 'junk.svsx');
    expect(await w.api.load.svs!({ file: junk })).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('fails on neither or both, on what is not a SeqViewSpec, on no structure, and on a file not loaded', async () => {
    const v = await viewer();
    expect(await v.api.load.svs!({})).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.svs!({ url: 'https://a.test/x.svsj', data: '{}' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await v.api.load.svs!({ data: '{"not": "svs"}' })).toMatchObject({ ok: false, code: 'invalid' });
    const plain = toJSON(createSeqViewSpec({ title: 'none' }).getDocument());
    expect(await v.api.load.svs!({ data: plain })).toMatchObject({ ok: false, code: 'invalid' });
    const elsewhere = svsOver({ host_file: 'other.cif', format: 'mmcif' });
    expect(await v.api.load.svs!({ data: elsewhere })).toMatchObject({ ok: false, code: 'not-found' });
    expect(seqOf(v).override.value).toBeUndefined();
  });
});
