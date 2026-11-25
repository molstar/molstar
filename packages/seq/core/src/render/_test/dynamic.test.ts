/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInLayerTypes, BuiltInLayouts } from '../built-ins.js';
import { DataStore, SeqRuntime, type StoreSettings } from '@molstar/seqx-core/data';
import { itemAddress, nodeAddress, OverlayBuilder, resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveView, resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import {
  attachSources,
  bundledTextMetrics,
  type Frame,
  type FrameInput,
  LightTokens,
  NoAssets,
  OverlaySources,
  type RenderConfig,
  Renderer,
  renderConfig,
} from '../index.js';
import { renderSvg } from '../render-svg.js';
import { Genes, transcriptsSpec } from './transcripts.js';

// Dynamic nodes drawn: transcripts arriving with the gene's regions, one track each from a template; their domains
// from nested sources only when a template needs them, through each transcript's codon map; the overlay's version in
// the draw lists' key while the static tiers stay.

const size = { width: 640, height: 240 };

async function live(
  spec = transcriptsSpec({ domains: true }),
  settings: Partial<StoreSettings> = {},
  config: Partial<RenderConfig> = {},
) {
  const model = resolveSpec(spec, { runtimes: [SeqRuntime] });
  const view = resolveViews(model)[0];
  const store = new DataStore({ settings });
  await attachSources(model, store);
  const overlays = new OverlayBuilder({ runtimes: [SeqRuntime] });
  const sources = new OverlaySources(store, NoAssets);
  const renderer = new Renderer({
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    config: renderConfig(config),
  });
  const frame = async (input: Partial<FrameInput> = {}, rounds = 6) => {
    let f = renderer.frame({ size, ...input, overlay: overlays.of(model, store) });
    for (let i = 0; i < rounds; i++) {
      await sources.sync(overlays.of(model, store), view);
      sources.hold(f.layout.shown);
      store.want(f.wants);
      await store.whenSettled();
      f = renderer.frame({ size, ...input, overlay: overlays.of(model, store) });
    }
    sources.hold(f.layout.shown);
    return f;
  };
  return { model, view, store, overlays, sources, renderer, frame };
}

describe('dynamic nodes through templates', () => {
  it('a template expands into a track per transcript, in key order, as the regions bring them', async () => {
    const r = await renderSvg(transcriptsSpec(), { size, runtimes: [SeqRuntime] });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    expect(r.svg).toContain('Transcript T1');
    expect(r.svg).toContain('Transcript T2');
    expect(r.svg).not.toContain('Domains');
  });

  it('a nested source is not requested until a template needs it; then its domains show through the codon map', async () => {
    const without = await live(transcriptsSpec());
    await without.frame();
    expect(without.store.nodes.size).toBe(4);
    expect(without.store.coverage(nodeAddress(Genes, 'T1-domains'))).toBeUndefined();

    const withDomains = await live();
    const f = await withDomains.frame();
    const nested = nodeAddress(Genes, 'T1-domains');
    expect(withDomains.store.coverage(nested)?.systems[0].name).toBe('T1');
    const tracks = f.layout.sections[0].tracks.map((t) => t.track.label);
    expect(tracks).toEqual(['letters', 'Transcript T1', 'Transcript T2', 'Domains T1-domains', 'Domains T2-domains']);
    // planning through the dynamic map: the gene's view carried back onto T1's residues, all of them in view
    const wants = f.wants.filter((w) => w.source === nested && w.priority === 'view');
    expect(wants.map((w) => [w.system, w.start, w.end])).toEqual([['T1', 0, 240]]);
    // a domain's marks are pieces of its exons on the gene
    const domains = f.layout.sections[0].tracks[3].layers[0];
    expect(domains.layer.path?.from).toBe(itemAddress(nodeAddress(Genes, 'T1'), 'T1'));
    expect(f.layout.issues).toEqual([]);
  });

  it('the overlay’s version is part of the key; the model and the resolved view stay', async () => {
    const spec = transcriptsSpec();
    const { model, view, store, overlays, renderer, frame } = await live(spec);
    // the static tiers before anything arrived
    const hash = model.hash;
    const tracks = view.sections[0].tracks;
    const printed = JSON.stringify(view, (_, v) => (typeof v === 'function' ? undefined : v));
    expect(store.nodes.size).toBe(0);
    await frame();
    expect(store.nodes.size).toBe(4);
    // after the arrivals: the same model and resolved view, unchanged, and a fresh resolve gives the same hash
    expect(model.hash).toBe(hash);
    expect(resolveSpec(spec, { runtimes: [SeqRuntime] }).hash).toBe(hash);
    expect(resolveView(model, view.address)).toBe(view);
    expect(view.sections[0].tracks).toBe(tracks);
    expect(view.sections[0].tracks.map((t) => t.template !== undefined)).toEqual([false, true]);
    expect(JSON.stringify(view, (_, v) => (typeof v === 'function' ? undefined : v))).toBe(printed);
    const o1 = overlays.of(model, store);
    const a = renderer.frame({ size, overlay: o1 });
    const same = renderer.frame({ size, overlay: o1 });
    expect(same.rebuilt).toBe(false);
    // the same nodes under a later version: only the key moves
    const o2 = Object.assign(Object.create(o1), { version: o1.version + 1 });
    const b = renderer.frame({ size, overlay: o2 });
    expect(b.rebuilt).toBe(true);
    expect(b.key).not.toBe(a.key);
    expect(resolveViews(model)[0]).toBe(view);
    expect(renderer.options.model).toBe(model);
  });

  it('a nested source that sends its nodes by region is asked for its regions in view', async () => {
    const { store, frame } = await live(transcriptsSpec({ domains: true, nestedByRegion: true }));
    const f = await frame();
    const nested = nodeAddress(Genes, 'T1-domains');
    expect(store.getNode(nested, 'domains')).toBeDefined();
    expect(f.wants.filter((w) => w.source === nested).map((w) => [w.system, w.start, w.end])).toContainEqual([
      'T1',
      0,
      240,
    ]);
    expect(f.layout.sections[0].tracks.map((t) => t.track.label)).toContain('Domains T1-domains');
  });

  it('a template instance holds its node only while in view: out of it, the node goes with its region', async () => {
    const T1 = nodeAddress(Genes, 'T1');
    const item = itemAddress(T1, 'T1');
    const labels = (f: Frame) => f.layout.sections[0].tracks.map((t) => t.track.label);
    // a byte bound nothing unwanted fits under: a region out of view goes at once
    const run = await live(transcriptsSpec({ domains: true }), { maxBytes: 1 });
    await run.frame();
    expect(run.store.getNode(Genes, 'T1')).toBeDefined();
    // on T1's region (the first): its nodes are shown and held, its domains through its item system
    const on = await run.frame({ viewport: { start: 600, end: 700 } });
    expect(on.layout.shown.map((n) => n.address)).toContain(T1);
    expect(labels(on)).toContain('Domains T1-domains');
    expect(run.overlays.of(run.model, run.store).system(item)?.length).toBe(240);
    // T1's region leaves the view and is evicted: its instances hold nothing, T1 and its tracks go; T2's stay
    const away = await run.frame({ viewport: { start: 1400, end: 1500 } });
    expect(run.store.get('/data/genes|tp53|1|0')).toBeUndefined();
    expect(away.layout.shown.map((n) => n.address)).not.toContain(T1);
    expect(run.store.getNode(Genes, 'T1')).toBeUndefined();
    expect(run.overlays.of(run.model, run.store).system(item)).toBeUndefined();
    expect(labels(away)).toEqual(['letters', 'Transcript T2', 'Domains T2-domains']);
    // back on it: the region and its nodes come again
    const back = await run.frame({ viewport: { start: 600, end: 700 } });
    expect(labels(back)).toContain('Transcript T1');

    // without holdShown, the node goes with its region
    const off = await live(transcriptsSpec({ domains: true }), { maxBytes: 1 }, { nodes: { holdShown: false } });
    await off.frame();
    const f = await off.frame({ viewport: { start: 1400, end: 1500 } });
    expect(f.layout.shown).toEqual([]);
    expect(off.store.getNode(Genes, 'T1')).toBeUndefined();
  });

  it('a template makes at most maxTemplateTracks tracks, and says so past it', async () => {
    const r = await renderSvg(transcriptsSpec(), { size, runtimes: [SeqRuntime], config: { maxTemplateTracks: 1 } });
    expect(r.issues.map((i) => [i.code, i.address])).toEqual([
      ['seq.render.template.limit', '/views/main/gene/transcripts'],
    ]);
    expect(r.svg).toContain('Transcript T1<');
    expect(r.svg).not.toContain('Transcript T2<');
  });

  it('past the nesting limit a nested source is an issue, never opened (EC-179)', async () => {
    const r = await renderSvg(transcriptsSpec({ domains: true }), {
      size,
      runtimes: [SeqRuntime],
      resolve: { maxNestedDepth: 0 },
    });
    expect(r.issues.map((i) => i.code)).toContain('seq.resolve.node.depth');
    expect(r.svg).not.toContain('Domains T1');
  });
});
