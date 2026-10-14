/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqViewerConfig } from "@molstar/seq-viewer";
import { useSeqState, useSeqViewer } from "@molstar/seq-viewer/react";
import { streamsOf } from "@molstar/seq-viewer/testing";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { useEffect, useMemo } from "react";
import {
  dynamicDoc,
  dynamicKnobs,
  dynamicLoaders,
  items,
  keyedDoc,
  keyedLoaders,
  PEAKS,
} from "./support/seqstar-docs-streamed.js";
import f from "./support/seqstar-format.module.css";
import { PrismPage } from "./support/seqstar-prism.js";
import { ViewerCard } from "./support/seqstar-viewer.js";
import { TextAction } from "./ui/index.js";

/**
 * Dynamic nodes: a streamed source's manifest and patches add data nodes to the set,
 * never to the document. A fixture provider over one 2 Mb chromosome serves a signal; its manifest brings the genes, and
 * each region of data (at base resolution, zoomed in past 10 positions a pixel) brings a peak annotation. Three tracks
 * name peaks A, B and C and stay hidden until their node arrives; nodes no track shows leave with their regions once
 * the cache is full. Latency and the cache are knobs.
 */
const meta = { title: "Seq*/Streaming/Streamed data", parameters: { layout: "fullscreen" } } satisfies Meta;
export default meta;
type Story = StoryObj<{ latency: number; cacheKb: number }>;

/** The buttons, the cache knob applied, and what the set holds of the stream. */
function Panel({ cacheKb }: { cacheKb: number }) {
  const { context: ctx } = useSeqViewer();
  const a = useSeqState(ctx.graph.assembled);
  useEffect(() => ctx.config.set(SeqViewerConfig.Stream.CacheBytes, cacheKb * 1024), [ctx, cacheKb]);
  const held: string[] = [];
  for (const x of a?.set.annotations ?? []) if (x) held.push(a?.set.strings.get(x.key) ?? "");
  const dynamic = held.filter((k) => k.startsWith("sig/"));
  const go = (at: number | null) => {
    const vp = ctx.viewports.get("dyn.main");
    if (at === null) vp?.fit();
    else vp?.zoomTo(0, at - 4000, at + 4000);
  };
  return (
    <div className={f.stack}>
      <div className={f.row}>
        <TextAction onClick={() => go(null)}>Whole chromosome</TextAction>
        {Object.entries(PEAKS).map(([k, at]) => (
          <TextAction key={k} onClick={() => go(at)}>{`Zoom to peak ${k}`}</TextAction>
        ))}
      </div>
      <div className={f.code}>
        Held from the stream: {dynamic.length} node{dynamic.length === 1 ? "" : "s"} ({dynamic.slice(0, 8).join(", ")}
        {dynamic.length > 8 ? ", …" : ""})
      </div>
    </div>
  );
}

function Card({ latency, cacheKb }: { latency: number; cacheKb: number }) {
  dynamicKnobs.latency = latency;
  const d = useMemo(dynamicDoc, []);
  const options = useMemo(() => ({ loaders: dynamicLoaders() }), []);
  return (
    <ViewerCard doc={d} options={options} label="A 2 Mb chromosome; peaks arrive as you zoom in and pan">
      <Panel cacheKb={cacheKb} />
    </ViewerCard>
  );
}

export const FirstPart: Story = {
  name: "Dynamic nodes: manifest, patch and eviction",
  args: { latency: 300, cacheKb: 256 },
  argTypes: {
    latency: { control: { type: "range", min: 0, max: 2000, step: 50 }, description: "ms before each region's answer" },
    cacheKb: { control: { type: "range", min: 64, max: 4096, step: 64 }, description: "Stream.CacheBytes, KiB" },
  },
  render: (args) => (
    <PrismPage
      title="Dynamic nodes: from a manifest, with a patch, and out again"
      caption={
        <>
          The genes come with the source's manifest, before any data. Zoom in past about 10 positions a pixel (or press
          a button): each region of data brings a peak, and a track that names one (Peak A, B, C) appears when it
          arrives. Pan on at that zoom and the regions behind you leave under the cache's bound, taking the peaks no
          track shows with them; the shown ones stay, and one that returns keeps its id. The document never changes.
        </>
      }
    >
      <Card latency={args.latency} cacheKb={args.cacheKb} />
    </PrismPage>
  ),
};

// --- Track templates over a keyed source, and sources nested two deep ---------------------------------------

/** What arrived: the manifests asked for, the instances shown, and a button that adds a gene as a patch would. */
function TemplatePanel({ asked }: { asked: string[] }) {
  const { context: ctx } = useSeqViewer();
  const a = useSeqState(ctx.graph.assembled);
  const issues = useSeqState(ctx.issues.current);
  const rv = a?.view?.resolved;
  let instances = 0;
  for (let k = 0; rv && k < rv.layers.count; k++) if (rv.layers.node[k] !== null) instances++;
  const deep = issues.filter((i) => i.code === "svs.stream.depth").length;
  const add = () => {
    const n = streamsOf(ctx).nodes;
    const k = n.entries.size + 1;
    n.add("genes", [items(`X${k}`, `Gene X${k} (a patch's)`, 50_000 + ((k * 137_000) % 1_800_000), 60_000)]);
  };
  return (
    <div className={f.stack}>
      <div className={f.row}>
        <TextAction onClick={add}>Add a gene, as a patch would</TextAction>
      </div>
      <div className={f.code}>
        Instance tracks: {instances}. Manifests asked for: {asked.length} ({asked.slice(-6).join(", ")}
        {asked.length > 6 ? ", …" : ""}). Sources deeper than max_depth: {deep}.
      </div>
    </div>
  );
}

function TemplateCard({ latency, maxDepth }: { latency: number; maxDepth: number }) {
  dynamicKnobs.latency = latency;
  // Remounted per max_depth (its key), so the document, the log and the loaders are made once each.
  const d = useMemo(() => keyedDoc(maxDepth), [maxDepth]);
  const asked = useMemo<string[]>(() => [], []);
  const options = useMemo(() => ({ loaders: keyedLoaders(asked) }), [asked]);
  return (
    <ViewerCard doc={d} options={options} label="Genes, their transcripts and their domains, by template">
      <TemplatePanel asked={asked} />
    </ViewerCard>
  );
}

export const SecondPart: StoryObj<{ latency: number; maxDepth: number }> = {
  name: "Dynamic nodes: templates and nested sources",
  args: { latency: 400, maxDepth: 4 },
  argTypes: {
    latency: {
      control: { type: "range", min: 0, max: 2000, step: 50 },
      description: "ms before a nested source's manifest",
    },
    maxDepth: { control: { type: "range", min: 1, max: 4, step: 1 }, description: "the genes source's max_depth" },
  },
  render: (args) => (
    <PrismPage
      title="Dynamic nodes: rows templated from a keyed source, sources nested two deep"
      caption={
        <>
          One track template per level: a track for every gene the genes source's manifest keys, for every transcript of
          the source each gene nests, and for every domain of the source each transcript nests. A nested source is asked
          for only because a template reaches its nodes; at max_depth 2 the domains' sources are too deep, reported and
          never asked. A gene added as a patch would add it gets its track at once, in key order.
        </>
      }
    >
      <TemplateCard key={args.maxDepth} latency={args.latency} maxDepth={args.maxDepth} />
    </PrismPage>
  ),
};
