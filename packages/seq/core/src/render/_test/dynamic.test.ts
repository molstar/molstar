/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import { itemAddress, nodeAddress, OverlayBuilder, resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import {
  attachSources,
  bundledTextMetrics,
  LightTokens,
  NoAssets,
  OverlaySources,
  Renderer,
  renderSvg,
} from '../index.js';
import { Genes, transcriptsSpec } from './transcripts.js';

// Dynamic nodes drawn: transcripts arriving with the gene's regions, one track each from a template; their domains
// from nested sources only when a template needs them, through each transcript's codon map; the overlay's version in
// the draw lists' key while the static tiers stay.

const size = { width: 640, height: 240 };

async function live(spec = transcriptsSpec({ domains: true })) {
  const model = resolveSpec(spec, { runtimes: [SeqRuntime] });
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store);
  const overlays = new OverlayBuilder({ runtimes: [SeqRuntime] });
  const sources = new OverlaySources(store, NoAssets);
  const renderer = new Renderer({ model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() });
  const frame = async (rounds = 6) => {
    let f = renderer.frame({ size, overlay: overlays.of(model, store) });
    for (let i = 0; i < rounds; i++) {
      await sources.sync(overlays.of(model, store), view);
      store.want(f.wants);
      await store.whenSettled();
      f = renderer.frame({ size, overlay: overlays.of(model, store) });
    }
    return f;
  };
  return { model, view, store, overlays, renderer, frame };
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
    const { model, view, store, overlays, renderer, frame } = await live(transcriptsSpec());
    await frame();
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
