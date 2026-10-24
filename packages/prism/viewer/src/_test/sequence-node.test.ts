/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '@molstar/model/model/structure';
import type { DocNode } from '@molstar/seq-core/format';
import { describe, expect, it, vi } from 'vitest';
import { type SceneNode, sequenceNodeId } from '../model/scene.js';
import type { SequenceViewer } from '../sequence/model.js';
import type { SequenceTrackProvider } from '../sequence/tracks.js';
import type { ViewerExtension } from '../spec.js';
import type { PrismViewer } from '../viewer.js';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { cbs, crambin, expectOk, headlessViewer, until } from './harness.js';

// Each structure's sequence (sequence.md §3, §7.1, §8.1): track providers put in and taken out (an extension's
// added, the viewer's own removed), the chain filter as the section's blocks by one patch (the document and the link
// kept), and the structure tree's Sequence node: its eye the band, its children the tracks, its menu and inspector.
// Headless.

/** An extension's track, out by default: what Add track offers. */
const Extra: SequenceTrackProvider = {
  id: 'test.sequence.extra',
  label: 'Extra',
  order: 40,
  track: 'test.extra',
  shown: false,
  build: (ctx) => ({
    track: {
      kind: 'track',
      ref: 'test.extra',
      params: { label: '' },
      children: [{ kind: 'layer', ref: 'test.extra.letters', params: { type: 'letters', color: { ...ctx.letters } } }],
    },
  }),
};
const ExtraExtension: ViewerExtension = { id: 'test', title: 'Test', sequence: { tracks: [Extra] } };

async function opened(files = [cbs()]) {
  const v = await headlessViewer({ files, extensions: [ExtraExtension] });
  const viewer = (await v.views.sequence!.load()).viewers.get();
  await until(() => viewer.tracks.value.length === 3);
  return { v, viewer };
}
const refs = (viewer: SequenceViewer) => viewer.tracks.value.map((t) => t.ref);
const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
/** The section's blocks in the document shown: the units it names. */
function blocks(viewer: SequenceViewer): string[] {
  const out: string[] = [];
  const walk = (n: DocNode) => {
    if (n.kind === 'block') out.push(String(n.params?.unit));
    for (const c of n.children ?? []) walk(c);
  };
  const doc = viewer.seq.document.shown.value;
  if (doc) walk(doc.root);
  return out;
}
const structureNodes = (v: PrismViewer) => v.views.scene.roots.value.filter((n) => n.kind === 'structure');
const sequenceOf = (n: SceneNode) => n.children.value.find((c) => c.kind === 'sequence');

describe("each structure's sequence", () => {
  it("adds an extension's track and removes the viewer's own, each a step Undo takes back", async () => {
    const { v, viewer } = await opened();
    const shown = viewer.structure.value!;
    expect(viewer.addableTracks({ ref: shown.ref, structure: shown.structure }).map((a) => a.label)).toEqual(['Extra']);
    const before = steps(v);
    expectOk(await v.api.sequence!.addTrack!({ track: 'Extra' }));
    expect(steps(v)).toBe(before + 1);
    expect(refs(viewer)).toEqual(['ruler', 'letters', 'ss-track', 'test.extra']);
    // Put in by a patch: the same document shown, its link kept.
    expect(viewer.structure.value).toBe(shown);
    expect(await v.api.sequence!.addTrack!({ track: 'Extra' })).toMatchObject({ ok: true, unchanged: true });

    expectOk(await v.api.sequence!.removeTrack!({ track: 'Ruler' }));
    expect(refs(viewer)).toEqual(['letters', 'ss-track', 'test.extra']);
    expect(viewer.edits.value[shown.ref]?.tracks).toEqual({ ruler: { in: false }, 'test.extra': { in: true } });
    const missing = await v.api.sequence!.addTrack!({ track: 'nothing' });
    expect(missing).toMatchObject({ ok: false, code: 'not-found' });
    expect((missing as { hint?: string }).hint).toContain('Ruler');

    expectOk(await v.api.history.undo());
    await until(() => refs(viewer).length === 4);
    expect(refs(viewer)).toEqual(['ruler', 'letters', 'ss-track', 'test.extra']);
    expectOk(await v.api.history.undo());
    await until(() => refs(viewer).length === 3);
    expect(viewer.edits.value).toEqual({});
  });

  it('filters the chains by blocks in one patch, the document kept; every copy apart; every chain again', async () => {
    const { v, viewer } = await opened();
    const shown = viewer.structure.value!;
    expect(blocks(viewer)).toEqual([]);
    expectOk(await v.api.sequence!.chains!({ polymer: true }));
    expect(blocks(viewer)).toEqual(['structure/A']);
    expect(viewer.structure.value).toBe(shown);
    expect(await v.api.sequence!.chains!({ polymer: true })).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.sequence!.chains!({ chains: 'A', each: true }));
    expect(blocks(viewer)).toEqual(['structure/A']);
    const section = (n: DocNode): DocNode | undefined =>
      n.kind === 'section' ? n : (n.children ?? []).reduce<DocNode | undefined>((f, c) => f ?? section(c), undefined);
    expect(section(viewer.seq.document.shown.value!.root)?.params?.copies).toBe('each');
    expect(await v.api.sequence!.chains!({ chains: 'Z' })).toMatchObject({ ok: false, code: 'not-found' });
    expectOk(await v.api.sequence!.chains!({}));
    expect(blocks(viewer)).toEqual([]);
    expect(viewer.edits.value).toEqual({});
  });

  it('puts a Sequence node under each structure: its eye the band, its tracks its children, its menu Add track', async () => {
    const { v, viewer } = await opened([cbs(), crambin()]);
    await until(() => structureNodes(v).every((n) => !!sequenceOf(n)));
    const [first, second] = structureNodes(v);
    const seqA = sequenceOf(first)!;
    const seqB = sequenceOf(second)!;
    expect(seqA.id).toBe(sequenceNodeId(first.id));
    // The band shows the first, by the rule (Auto); the other is not shown.
    expect(seqA.isHidden.value).toBe(false);
    expect(seqA.description.value).toContain('auto');
    expect(seqB.isHidden.value).toBe(true);
    expect(seqA.children.value.map((c) => c.label.value)).toEqual(['Ruler', 'Residues', 'Secondary structure']);

    // The second's eye: the band follows it (chosen, not Auto).
    seqB.toggleVisibility();
    await until(() => viewer.structure.value?.title === '1CRN' && viewer.tracks.value.length === 3);
    await until(() => !sequenceOf(structureNodes(v)[1])!.isHidden.value);
    expect(sequenceOf(structureNodes(v)[0])!.isHidden.value).toBe(true);
    expect(sequenceOf(structureNodes(v)[1])!.description.value).not.toContain('auto');

    // A track's eye and Remove, through the API.
    const ruler = sequenceOf(structureNodes(v)[1])!.children.value[0];
    ruler.toggleVisibility();
    await until(() => viewer.tracks.value[0].hidden);
    expect(ruler.removable).toBe(true);
    ruler.remove();
    await until(() => viewer.tracks.value.length === 2);

    // Add track lists the extension's; Follow automatically while one is chosen.
    const menu = sequenceOf(structureNodes(v)[1])!.menu!;
    const add = menu.find((m) => m.id === 'add-track')!;
    expect(add.items!.map((i) => i.label)).toEqual(['Ruler', 'Extra']);
    add.items!.find((i) => i.label === 'Extra')!.run!();
    await until(() => refs(viewer).includes('test.extra'));
    menu.find((m) => m.id === 'auto')!.run!();
    await until(() => v.views.sequence!.source.followed.value === undefined);

    // Its inspector: the chain filter; a chain left out is a filter.
    const node = sequenceOf(structureNodes(v)[1])!;
    const inspection = v.views.scene.inspect(node);
    expect(inspection.slot).toBe('inspector:sequence');
    expect(inspection.sequence).toBeDefined();
  });

  it('takes edits right after a load, its document still being made: they are in it once it is', async () => {
    const v = await headlessViewer({ files: [cbs()], extensions: [ExtraExtension] });
    expectOk(await v.api.sequence!.chains!({ polymer: true }));
    expectOk(await v.api.sequence!.hideTrack!({ track: 'Ruler' }));
    expectOk(await v.api.sequence!.addTrack!({ track: 'Extra' }));
    const viewer = (await v.views.sequence!.load()).viewers.get();
    await until(() => viewer.tracks.value.length === 4 && blocks(viewer).length === 1);
    expect(viewer.tracks.value.filter((t) => t.hidden).map((t) => t.ref)).toEqual(['ruler']);
    expect(blocks(viewer)).toEqual(['structure/A']);
  });

  it('puts the node right under Contents; `sequence.tree: false` leaves it out, the band and the API kept', async () => {
    const { v } = await opened();
    await until(() => !!sequenceOf(structureNodes(v)[0]));
    expect(
      structureNodes(v)[0]
        .children.value.map((c) => c.kind)
        .slice(0, 2),
    ).toEqual(['contents', 'sequence']);

    const spec = DefaultPrismViewerSpec();
    const w = await headlessViewer({
      files: [cbs()],
      spec: { sequence: spec.sequence ? { ...spec.sequence, tree: false } : false },
    });
    const viewer = (await w.views.sequence!.load()).viewers.get();
    await until(() => viewer.tracks.value.length === 3);
    expect(structureNodes(w)[0].children.value.some((c) => c.kind === 'sequence')).toBe(false);
    expectOk(await w.api.sequence!.removeTrack!({ track: 'Ruler' }));
  });

  it("the inspector's chains: presets, an entity left out, none; the node's style and letter colors", async () => {
    const { v, viewer } = await opened();
    await until(() => !!sequenceOf(structureNodes(v)[0]));
    const chains = v.views.scene.inspect(sequenceOf(structureNodes(v)[0])!).sequence!;
    expect(chains.entities.value.map((e) => e.kind)).toEqual(['polymer', 'ligand']);
    expect(chains.preset.value).toBe('all');
    const ligand = chains.entities.value.find((e) => e.kind === 'ligand')!;
    chains.toggleEntity(ligand.id);
    await until(() => chains.preset.value === 'polymer');
    expect(blocks(viewer)).toEqual(['structure/A']);
    chains.apply('none');
    await until(() => chains.preset.value === 'none');
    expect(chains.shown.value.size).toBe(0);
    chains.apply('all');
    await until(() => chains.preset.value === 'all');
    expect(blocks(viewer)).toEqual([]);

    // Style and letter colors from the node's menu: the band's, through the API.
    const menu = sequenceOf(structureNodes(v)[0])!.menu!;
    menu.find((m) => m.id === 'style')!.items!.find((i) => i.id === 'text')!.run!();
    await until(() => viewer.style.value === 'compact');
    menu.find((m) => m.id === 'colors')!.items!.find((i) => i.id === 'colors:zappo')!.run!();
    await until(() => viewer.theme.value.letters === 'zappo');
    await until(() => sequenceOf(structureNodes(v)[0])!.children.value[1].description.value === 'Letters · Zappo');
  });

  it('a refused edit builds the document again and leaves the API free; None empties the band', async () => {
    // A provider whose track's ref is not a ref (`@`): Seq* refuses its insert.
    const Broken: SequenceTrackProvider = {
      ...Extra,
      id: 'test.sequence.broken',
      label: 'Broken',
      track: 'test.broken@1',
      build: () => ({
        track: {
          kind: 'track',
          ref: 'test.broken@1',
          params: { label: '' },
          children: [{ kind: 'layer', ref: 'test.broken.values', params: { type: 'values', data: 'test.nothing' } }],
        },
      }),
    };
    const v = await headlessViewer({
      files: [cbs()],
      extensions: [{ id: 'test', title: 'Test', sequence: { tracks: [Broken] } }],
    });
    const viewer = (await v.views.sequence!.load()).viewers.get();
    await until(() => viewer.tracks.value.length === 3);
    const warned = vi.spyOn(console, 'warn').mockImplementation(() => {});
    await v.api.sequence!.addTrack!({ track: 'Broken' });
    expect(warned.mock.calls.some((c) => String(c[0]).includes('refused'))).toBe(true);
    // The queue is free: the next calls run.
    expectOk(await v.api.sequence!.style({ style: 'text' }));
    expectOk(await v.api.sequence!.chains!({ chains: 'none' }));
    await until(() => viewer.empty.value);
    expectOk(await v.api.sequence!.chains!({}));
    await until(() => !viewer.empty.value);
  });

  it('removing the focused structure the band follows: the tree rebuilt once the removal is done, no error', async () => {
    const { v, viewer } = await opened([cbs(), crambin()]);
    const second = v.plugin.managers.structure.hierarchy.current.structures[1];
    const ref = second.cell.transform.ref;
    const errors: unknown[] = [];
    const caught = (e: unknown) => errors.push(e);
    process.on('uncaughtException', caught);
    try {
      v.plugin.managers.structure.focus.setFromLoci(Structure.toStructureElementLoci(second.cell.obj!.data));
      await until(() => viewer.structure.value?.ref === ref);
      expectOk(await v.api.state!.remove({ node: ref }));
      await until(() => structureNodes(v).length === 1 && viewer.structure.value?.title === '1CBS');
      await new Promise((r) => setTimeout(r, 50));
      expect(errors).toEqual([]);
      expect(!!sequenceOf(structureNodes(v)[0])).toBe(true);
    } finally {
      process.off('uncaughtException', caught);
    }
  });

  it('edits a structure not shown: its own tracks, kept for when it is shown', async () => {
    const { v, viewer } = await opened([cbs(), crambin()]);
    await until(() => structureNodes(v).every((n) => !!sequenceOf(n)));
    expectOk(await v.api.sequence!.removeTrack!({ track: 'Ruler', structure: '1crn' }));
    expectOk(await v.api.sequence!.addTrack!({ track: 'Extra', structure: '1crn' }));
    expect(refs(viewer)).toEqual(['ruler', 'letters', 'ss-track']);
    await until(() => sequenceOf(structureNodes(v)[1])!.children.value.length === 3);
    expect(sequenceOf(structureNodes(v)[1])!.children.value.map((c) => c.label.value)).toEqual([
      'Residues',
      'Secondary structure',
      'Extra',
    ]);
    expectOk(await v.api.sequence!.follow({ structure: '1crn' }));
    await until(() => viewer.structure.value?.title === '1CRN' && viewer.tracks.value.length === 3);
    expect(refs(viewer)).toEqual(['letters', 'ss-track', 'test.extra']);
  });

  it("follows the 3D view's secondary structure when the person sets it otherwise (the structure's properties)", async () => {
    const { v, viewer } = await opened();
    const state = v.plugin.state.data;
    const props = [...state.cells.values()].find((c) =>
      c.transform.transformer.id.endsWith('custom-structure-properties'),
    );
    expect(props).toBeDefined();
    const ss = () => JSON.stringify(viewer.seq.document.shown.value).match(/"assignment":"([^"\\]|\\.)*"/)?.[0] ?? '';
    expect(ss()).toContain('auto');
    const params = props!.transform.params as { properties: Record<string, unknown> };
    await v.plugin.runTask(
      state.updateTree(
        state
          .build()
          .to(props!.transform.ref)
          .update({
            ...params,
            properties: {
              ...params.properties,
              molstar_computed_secondary_structure: { type: { name: 'model', params: {} } },
            },
          }),
      ),
    );
    await until(() => ss().includes('model'));
  });
});
