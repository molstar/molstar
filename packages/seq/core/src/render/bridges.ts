/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bridge chunk's entry (_spec/README.md "Render", Bridges): each placed bridge laid out after its sides, its band
// between them, its `pairs` layers as ribbons (`pairs-ribbons.ts`), drawn over the page and hit before the sections.
// Its key is its sides' layouts (a pan of one side relays that side and the bridge only) and its layers' preps. Loaded
// when a resolved view has a bridge (`bridge-load.ts`): a view without one loads none of it.

import type { Issue } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import type { ResolvedBridge, ResolvedLayer } from '@molstar/seqx-core/view';
import type { ThemeContext } from './channels.js';
import type { RenderConfig } from './config.js';
import type { PassName, PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { Box, Emphasis, Hit, LayoutContext } from './layer-type.js';
import { issueOf } from '@molstar/core/repr/tree';
import { BridgePlacement } from './layout/bridge.js';
import { type BridgeLayoutContext, RibbonsLay, RibbonsType } from './pairs-ribbons.js';
import type { PairsPrep } from './pairs-lay.js';
import type {
  BridgeLayout,
  BridgesLayout,
  LayerLayout,
  LayoutInput,
  LayoutState,
  SectionLayout,
} from './tracks-types.js';
import { escapeXml } from './svg.js';

/** What laying out a view's bridges reads of its layout pass. */
export interface BridgesInput {
  readonly input: LayoutInput;
  readonly state: LayoutState;
  readonly config: RenderConfig;
  readonly theme: ThemeContext;
  /** The view's sections laid out. */
  readonly sections: readonly SectionLayout[];
  /** The layout's wants: a bridge's reads join them. */
  readonly wants: Want[];
  report(i: Issue): void;
}

/** What a bridge's layer was laid out for: its sides' layouts and its prep (its key). */
interface Kept {
  a: SectionLayout;
  b: SectionLayout;
  prep: PairsPrep;
  lay: RibbonsLay;
  /** What it said, said again when it is drawn again. */
  said: Issue[];
}

const KeptOf = new WeakMap<LayoutState, Map<string, Kept>>();

/** The view's placed bridges laid out (_spec/README.md "Render", Bridges), in spec order; hidden ones left out. */
export function layoutBridges(x: BridgesInput): BridgesLayout | undefined {
  const { input, state } = x;
  const view = input.view;
  let kept = KeptOf.get(state);
  if (!kept) KeptOf.set(state, (kept = new Map()));
  const used = new Set<string>();
  const bySection = new Map(x.sections.map((s) => [s.section.address, s]));
  const out: Bridge[] = [];
  const all = view.bridges ?? [];
  for (let i = 0; i < all.length; i++) {
    const br = all[i];
    if (!br.placed || br.hidden) continue;
    const a = bySection.get(br.a);
    const b = bySection.get(br.b);
    if (!a || !b) continue;
    const placement = new BridgePlacement().place(a, b, br.above === 'a', x.input.dpr ?? 1);
    const id = br.address.slice(br.address.lastIndexOf('/') + 1);
    const lays: LayerLayout[] = [];
    for (const layer of br.layers) {
      if (!layer.placed || layer.params.hidden === true || !layer.address) continue;
      if (!input.types?.pairs) {
        x.report(
          issueOf(RenderIssues, 'seq.render.layer.type', 'No type of the layer’s fallbacks is drawn here', {
            address: layer.address,
            hint: RenderIssues.codes['seq.render.layer.type'].hint,
          }),
        );
        continue;
      }
      used.add(layer.address);
      lays.push(layBridgeLayer(x, kept, i, id, layer, placement, a, b));
    }
    out.push(new Bridge(i, br, id, placement, lays));
  }
  for (const k of [...kept.keys()]) if (!used.has(k)) kept.delete(k);
  return out.length ? new Bridges(out) : undefined;
}

/** A view's bridges: their draw, hit and export groups (`BridgesLayout`). */
class Bridges implements BridgesLayout {
  constructor(readonly list: readonly Bridge[]) {}
  draw(pass: PassName, out: PassWriter, emphasis: Emphasis | undefined): void {
    for (let i = 0; i < this.list.length; i++) this.list[i].draw(pass, out, emphasis);
  }
  hit(at: Float64Array): Hit | null | undefined {
    const i = this.band(at[0], at[1]);
    return i < 0 ? undefined : (this.list[i].hit(at) ?? null);
  }
  at(x: number, y: number): number {
    const i = this.band(x, y);
    return i < 0 ? -1 : this.list[i].index;
  }
  private band(x: number, y: number): number {
    for (let i = 0; i < this.list.length; i++) {
      const b = this.list[i].box;
      if (y >= b.y && y < b.y + b.height && x >= b.x && x < b.x + b.width) return i;
    }
    return -1;
  }
  svg(document: string): string {
    const groups = this.list
      .map(
        (b) =>
          `<g data-bridge="${escapeXml(b.address)}" data-a="${escapeXml(b.a)}" data-b="${escapeXml(b.b)}" data-y="${b.box.y}" data-height="${b.box.height}" data-count="${b.count}"/>`,
      )
      .join('');
    // before the passes' groups, after the sections'
    const at = document.indexOf('<g data-pass=');
    return at < 0 ? document : document.slice(0, at) + groups + document.slice(at);
  }
}

/** One layer of a bridge: laid out again when its sides' layouts or its prep changed, else drawn again as kept. */
function layBridgeLayer(
  x: BridgesInput,
  kept: Map<string, Kept>,
  index: number,
  id: string,
  layer: ResolvedLayer,
  placement: BridgePlacement,
  a: SectionLayout,
  b: SectionLayout,
): LayerLayout {
  const { input, state } = x;
  const key = layer.address!;
  let prep = state.preps.get(layer) as PairsPrep | undefined;
  // the composition's pairs type prepares it (its data, its channels): the chunk shares its prepared form
  if (prep === undefined) state.preps.set(layer, (prep = input.types!.pairs.prepare(layer, input.model) as PairsPrep));
  const k = kept.get(key);
  // a source read by region reads the store each frame: never kept by its sides alone
  if (k && k.a === a && k.b === b && k.prep === prep && !prep.region) {
    for (const i of k.said) x.report(i);
    k.lay.band = placement.box;
    return { key, layer, type: RibbonsType, lay: k.lay };
  }
  const said: Issue[] = [];
  const ctx: BridgeLayoutContext = {
    box: placement.box,
    axis: a.axis!,
    section: placement.a.section,
    key,
    dpr: x.input.dpr ?? 1,
    config: x.config,
    tokens: input.tokens,
    text: input.text,
    settled: input.settled ?? false,
    theme: x.theme,
    want: (w) => x.wants.push(w),
    report: (i) => {
      said.push(i);
      x.report(i);
    },
    bridge: { placement, index, id },
  } satisfies LayoutContext & { bridge: unknown };
  const lay = RibbonsType.layout(prep, ctx, input.data, input.overlay!, k?.lay);
  kept.set(key, { a, b, prep, lay, said });
  return { key, layer, type: RibbonsType, lay };
}

/** A bridge laid out: its band, its layers, its draw and its hit (_spec/README.md "Render", Bridges). */
class Bridge implements BridgeLayout {
  readonly address: string;
  readonly a: string;
  readonly b: string;
  readonly label: string;
  constructor(
    readonly index: number,
    bridge: ResolvedBridge,
    readonly id: string,
    readonly placement: BridgePlacement,
    readonly lays: readonly LayerLayout[],
  ) {
    this.address = bridge.address;
    this.a = bridge.a;
    this.b = bridge.b;
    this.label = bridge.label;
  }

  get box(): Box {
    return this.placement.box;
  }

  /** The ribbons drawn (its ARIA description's count). */
  get count(): number {
    let n = 0;
    for (const l of this.lays) n += (l.lay as RibbonsLay).n;
    return n;
  }

  draw(pass: PassName, out: PassWriter, emphasis: Emphasis | undefined): void {
    for (let i = 0; i < this.lays.length; i++) this.lays[i].type.draw(this.lays[i].lay, pass, out, emphasis);
  }

  hit(at: Float64Array): Hit | undefined {
    for (let i = this.lays.length - 1; i >= 0; i--) {
      const h = this.lays[i].type.hit(this.lays[i].lay, at);
      if (h) return h;
    }
    return undefined;
  }
}
