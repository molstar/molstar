/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lucideBox } from '@molstar/ui/components/icon/_generated/lucide/box';
import type { ComponentExpressionT, HexColorT } from '@molstar/mvs-builder/tree/mvs/param-types';
import { location } from '@molstar/seqx-core/query';
import {
  type ColorRun,
  colorRuns,
  type StructureBinding,
  type StructureFacts,
  type StructureTarget,
  structureFacts,
  type UnitFacts,
} from '@molstar/seqx-core/structure';
import type { BehaviorProvider, SeqViewerBehavior } from '@molstar/seqx-viewer';
import type { ViewerCommand } from '@molstar/seqx-viewer/viewer';
import type { HeaderActionEntry, HeaderFacts } from '@molstar/seqx-viewer/header-entries';
import type { StructureSide } from '@molstar/seqx-viewer/link-structure';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { FakeStructureViewer, Selectors } from './fake-structure/fake-viewer.js';
import { type ColorNode, structureState } from './fake-structure/mvs-state.js';

// The host's adapter for the fake structure viewer (structure-syncing, The fake host's adapter): the one place that
// knows both vocabularies. Seq*'s words are Locations on the units of a complex (`structureFacts` says what a unit is);
// the fake's are MolViewSpec selectors. It translates on the facts once and looks up per message, as a real host's
// adapter will, and gives `linkStructure` its `side`. A run on a polymer unit is `{ label_asym_id, instance_id,
// beg_label_seq_id, end_label_seq_id }` (`label_seq_id` is the unit's position on its molecule plus one); a molecule
// unit is `{ label_asym_id, instance_id }`; inward the reverse. Not placed: a branched unit's residue by author number,
// `label_entity_id` and class selectors, MolQL (no dataset of the fake has a glycan, and the fake emits none).

/** What a host's adapter is: the translation in both directions and the side a link takes. */
export interface HostAdapter {
  readonly side: StructureSide;
  readonly facts: StructureFacts;
  selectorsOfTarget(target: StructureTarget): Selectors;
  targetOfSelectors(selectors: Selectors): StructureTarget | null;
}

/** The selector fields the adapter places. */
const Placed = new Set(['label_asym_id', 'instance_id', 'label_seq_id', 'beg_label_seq_id', 'end_label_seq_id']);

/** A unit's start on its molecule, the base `label_seq_id` counts from; a ligand has one residue. */
const baseOf = (u: UnitFacts) => (u.kind === 'polymer' ? u.start : 0);

/** The adapter of `fake` over a complex's facts. */
export function hostAdapter(fake: FakeStructureViewer, facts: StructureFacts): HostAdapter {
  const byAddress = new Map(facts.units.map((u) => [u.address, u]));

  const selectorsOfTarget = (target: StructureTarget): ComponentExpressionT[] => {
    const out: ComponentExpressionT[] = [];
    for (const p of target.parts) {
      const u = byAddress.get(p.system);
      if (!u || u.chain === undefined) continue;
      const unit = { label_asym_id: u.chain, ...(u.operator !== undefined ? { instance_id: u.operator } : {}) };
      if (u.kind !== 'polymer') {
        out.push(unit);
        continue;
      }
      for (let i = 0; i + 1 < p.ranges.length; i += 2) {
        const beg = baseOf(u) + p.ranges[i] + 1;
        const end = baseOf(u) + p.ranges[i + 1];
        out.push(
          beg === end ? { ...unit, label_seq_id: beg } : { ...unit, beg_label_seq_id: beg, end_label_seq_id: end },
        );
      }
    }
    return out;
  };

  const targetOfSelectors = (selectors: Selectors): StructureTarget | null => {
    const parts: { system: string; ranges: number[] }[] = [];
    // a field this adapter does not place would widen the target silently: the whole is not placed
    if (selectors.some((s) => Object.keys(s).some((k) => !Placed.has(k)))) return null;
    for (const s of selectors) {
      for (const u of facts.units) {
        if (s.label_asym_id !== undefined && s.label_asym_id !== u.chain) continue;
        if (s.instance_id !== undefined && u.operator !== undefined && s.instance_id !== u.operator) continue;
        const length = u.length ?? 1;
        if (u.kind !== 'polymer') parts.push({ system: u.address, ranges: [0, length] });
        else {
          const lo = s.label_seq_id ?? s.beg_label_seq_id;
          const hi = s.label_seq_id ?? s.end_label_seq_id;
          const from = Math.max(0, (lo === undefined ? 1 : lo) - 1 - baseOf(u));
          const to = Math.min(length, (hi === undefined ? baseOf(u) + length : hi) - baseOf(u));
          if (to > from) parts.push({ system: u.address, ranges: [from, to] });
        }
      }
    }
    return parts.length ? location({ parts }) : null;
  };

  const side: StructureSide = {
    highlight: (t) => fake.setHighlight(t ? selectorsOfTarget(t) : null),
    select: (t) => fake.setSelection(t ? selectorsOfTarget(t) : null),
    focus: (t, { durationMs }) => fake.focusOn(selectorsOfTarget(t), durationMs),
    onHover: (fn) => fake.onHover((s, at) => fn(targetOfSelectors(s), at)),
    onSelect: (fn) => fake.onSelect((s) => fn(targetOfSelectors(s))),
    onFocus: (fn) =>
      fake.onFocus((s) => {
        const t = targetOfSelectors(s);
        if (t) fn(t);
      }),
    signal: fake.signal,
  };
  return { side, facts, selectorsOfTarget, targetOfSelectors };
}

/** The model a viewer's translation reads: its static model with the store's current overlay. */
export function viewerModel(viewer: SeqViewer) {
  const model = viewer.spec.model.value;
  return { model, overlay: viewer.store.currentOverlay(model) };
}

/** The facts of a complex of the viewer's model, or undefined when the address is no complex. */
export function viewerFacts(viewer: SeqViewer, complex: string): StructureFacts | undefined {
  return structureFacts(viewerModel(viewer), complex);
}

const hex = (rgba: number): HexColorT => `#${((rgba >>> 8) & 0xffffff).toString(16).padStart(6, '0')}`;

/** A track's color runs as the `color` nodes of the fake's MolViewSpec state, with the light or the dark tokens. */
export function colorNodesOfRuns(
  adapter: HostAdapter,
  runs: readonly ColorRun[],
  scheme: 'light' | 'dark',
): ColorNode[] {
  return runs.map((r) => ({
    selector: adapter.selectorsOfTarget(r.target),
    color: hex(scheme === 'dark' ? r.dark : r.color),
  }));
}

/** The track each fake is colored by (a snapshot of it). */
const ColoredBy = new WeakMap<FakeStructureViewer, string>();
/** The layer each viewer's structure is colored by: what the header's pressed icon shows (one at a time). */
const SourceOf = new WeakMap<SeqViewer, string>();

/**
 * "Color by this track": the layer's colors as runs on the units (`colorRuns`, the drawn track's: the viewer's renderer
 * compiles them), written as the fake's MolViewSpec state and taken by it. A snapshot: it does not follow the track.
 * Run again on the same track it restores the base colors. Returns the message of what happened.
 */
export function colorByLayer(
  viewer: SeqViewer,
  fake: FakeStructureViewer,
  adapter: HostAdapter,
  layerAddress: string,
  scheme: 'light' | 'dark',
  complexes?: readonly StructureBinding[],
): string {
  const renderer = viewer.renderer;
  const resolved = viewer.view.resolved.value;
  if (!renderer || !resolved) return 'The view is not drawn yet';
  const layer = resolved.sections
    .flatMap((s) => s.tracks.flatMap((t) => t.layers))
    .find((l) => l.address === layerAddress);
  if (!layer) return `No layer ${layerAddress}`;
  // the same track again: back to the base state
  if (ColoredBy.get(fake) === layerAddress) {
    ColoredBy.delete(fake);
    SourceOf.delete(viewer);
    fake.setState(structureState(fake.entry, []));
    return `Restored the base colors`;
  }
  const r = colorRuns(viewerModel(viewer), viewer.store.store, layer, complexes, { colors: renderer.layerColors() });
  const issues = fake.setState(structureState(fake.entry, colorNodesOfRuns(adapter, r.runs, scheme), layerAddress));
  if (issues) return `The structure refused the state: ${issues[0]}`;
  ColoredBy.set(fake, layerAddress);
  SourceOf.set(viewer, layerAddress);
  return `Colored ${r.runs.length} runs${r.complete ? '' : ' (some of the track’s data is still pending)'}`;
}

const ColorByCommand = 'fake-host.color-by';

/** The host's behavior that provides the command "Color the structure by this track" (what its header entry runs). */
export function colorByBehavior(run: (layer: string) => string): SeqViewerBehavior {
  const provider: BehaviorProvider<undefined, { commands: readonly ViewerCommand[] }> = {
    name: 'fake-host.color-by',
    create: () => ({
      commands: [
        {
          id: ColorByCommand,
          title: 'Color the structure by this track',
          category: 'Structure',
          when: (ctx: { layer?: string }) => (ctx?.layer ? true : 'No layer to color by'),
          run: (_args: unknown, ctx: { layer: string }) => ({ message: run(ctx.layer) }),
        },
      ],
    }),
  };
  return { provider, params: undefined };
}

const colorLayer = (f: HeaderFacts) => f.layers.find((l) => !l.address.includes(':') && l.schemes.length > 0);
/** Whether the structure is colored by this header's track. */
const colorSource = (f: HeaderFacts) => !!f.viewer && SourceOf.get(f.viewer) === colorLayer(f)?.address;

/** The host's header entry: on a track with a letters, values or ranges layer, "Color the structure by this track". */
export const ColorByEntry: HeaderActionEntry = {
  kind: 'header-action',
  entry: {
    name: 'fake-host:color-by',
    on: 'track',
    icon: 'box',
    iconMarkup: lucideBox,
    label: (f) =>
      colorSource(f)
        ? 'The structure is colored by this track (click to restore)'
        : 'Color the structure by this track',
    command: ColorByCommand,
    order: 0,
    when: (f) => f.layers.some((l) => !l.address.includes(':') && l.schemes.length > 0),
    args: (f) => ({ layer: colorLayer(f)!.address }),
    pressed: (f) => colorSource(f),
  },
};
