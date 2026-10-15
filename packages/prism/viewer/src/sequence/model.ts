/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasInk, type IconName, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { Structure } from '@molstar/model/model/structure';
import type { FormatResult, Loaders } from '@molstar/seq-core/format';
import { MolstarSvsEntry, type StructureSource } from '@molstar/seq-molstar-poc';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  type Scheduler,
  type SeqViewerContext,
  rafScheduler,
  SeqViewerOptionalBehaviors,
  SeqViewerSpec,
  type TextMeasurer,
} from '@molstar/seq-viewer';
import {
  ComplexRef,
  DataRef,
  LettersRef,
  SecondaryTrack,
  DefaultSequenceTheme,
  hostFileOf,
  SectionRef,
  type SequenceTheme,
  sequenceDocument,
} from './document.js';
import type { SequenceBehavior, SequenceOffer, SequenceShown, SequenceTrackMenuItem } from './behaviors.js';
import { legendOf, type SequenceLegend, withTracks } from './legend.js';
import { PrismRegistry } from './palette.js';
import { DefaultTracksRail, TextPreset, type TracksRail, TracksPreset, tracksPreset } from './presets.js';
import { type SequenceFonts, skinFor } from './skin.js';

// The sequence model (prism-sequence.md §2.1): owns the Seq* contexts, one per id ('default' only for now), and keeps
// each in step with what the viewer follows: the canvas ink (its skin), the structure the 3D view shows (its
// document), and the person's choices (the style, the color theme). Plain arguments, no DOM and no plugin: testable
// headless.

export type SequenceStyle = 'compact' | 'tracks';
export type { SequenceTheme };

/** A structure the host (the 3D view) shows and has parsed. */
export interface HostStructure {
  /** Tells one structure from another: a new key rebuilds the view. */
  readonly key: string;
  /** What the person calls it ("1TQN"). */
  readonly title: string;
  /** The host's parse (`StructureSource.of(trajectory, format)`): shared with the 3D view, never read again. */
  readonly source: StructureSource;
  /** The structure the 3D view draws: what the link maps. */
  readonly structure: Structure;
  /** The ligands the 3D view draws, by chemical component id; absent: every ligand. The sequence lists the same ones. */
  readonly ligands?: readonly string[];
}

/** What a viewer shows, once its document is in: the link's input. */
export interface SequenceStructure {
  readonly title: string;
  readonly structure: Structure;
  /** The Seq* complex made from it. */
  readonly complex: string;
}

/**
 * How the 3D view colors the polymer: a color theme's name (`secondary-structure`, an extension's `plddt-confidence`),
 * `default` for the viewer's default coloring.
 */
export type StructureColoring = string;

/** What the sequence view may ask of the 3D view (the host implements it over its plugin, `sources/`). */
export interface HostActions {
  colorStructure(coloring: StructureColoring): void | Promise<void>;
  /** How the 3D view's polymer is colored now (read from the plugin's state, so a change made elsewhere shows); undefined: no one theme. */
  coloring?(): StructureColoring | undefined;
}

/** An action a track offers in its context menu. */
export interface TrackAction {
  readonly id: string;
  readonly label: string;
  readonly icon?: IconName;
  run(): void;
}

/** The coloring a track's icon in the gutter stands for (the track's first action), and whether the 3D view has it now. */
export interface TrackColoring extends TrackAction {
  readonly coloring: StructureColoring;
  /** The 3D polymer is colored this way now (read from the host when asked). */
  active(): boolean;
}

/**
 * What the viewer's own tracks offer, by ref: the model's decision, as it set the tracks and their colors (design §6);
 * behaviors offer their tracks' (`SequenceOffer.menus`). The 3D coloring each track stands for: its own track colors the
 * structure the same way.
 */
const TrackMenus: Readonly<Record<string, readonly SequenceTrackMenuItem[]>> = {
  [LettersRef]: [{ id: 'color-default', label: 'Apply default coloring', icon: 'paintbrush', theme: 'default' }],
  [SecondaryTrack]: [
    { id: 'color-secondary', label: 'Color by secondary structure', icon: 'paintbrush', theme: 'secondary-structure' },
  ],
};

export interface SequenceViewer {
  readonly id: string;
  readonly seq: SeqViewerContext;
  readonly style: ReadonlyState<SequenceStyle>;
  /** UI state only (`ctx.ui.setStyle`): the document does not change. */
  setStyle(style: SequenceStyle): void;
  /** What it shows, once the document is in; undefined: nothing (no structure, or the document is being made). */
  readonly structure: ReadonlyState<SequenceStructure | undefined>;
  readonly theme: ReadonlyState<SequenceTheme>;
  /** Written into the document the model builds. */
  setTheme(theme: SequenceTheme): void;
  /** What the theme's colors mean, for the ink's theme. */
  readonly legend: ReadonlyState<SequenceLegend>;
  /** The tracks version's rail: which navigator, how tall (the model registers the preset again). */
  readonly rail: ReadonlyState<TracksRail>;
  setRail(rail: TracksRail): void;
  /**
   * The context menu of a track (by its ref): its actions on the 3D view; none: the track has no menu.
   */
  trackActions(track: string): readonly TrackAction[];
  /** The coloring icon of a track (by its ref); none: the track has none, or the host offers no actions. */
  trackColoring(track: string): TrackColoring | undefined;
  /** The host's fonts (the control reads its `--msp-font-*` once): part of the skin, not of the ink. */
  setFonts(fonts: SequenceFonts): void;
  /** What the document shows, for behaviors: a new value per document. */
  readonly shown: ReadonlyState<SequenceShown | undefined>;
  /** Attaches a behavior (an extension's, sequence.md §3); returns its removal. An id taken throws. */
  addBehavior(behavior: SequenceBehavior): () => void;
}

export interface SequenceViewers {
  /** The viewers by id; the first slice creates only 'default'. */
  get(id?: string): SequenceViewer;
  readonly ids: ReadonlyState<readonly string[]>;
  /** Disposes every viewer's context. */
  dispose(): void;
}

export interface SequenceViewersDeps {
  /** The canvas ink: each viewer's skin follows it. */
  ink: ReadonlyState<CanvasInk>;
  /** The structure the 3D view shows (the structure source, §3); undefined: none. */
  structure: ReadonlyState<HostStructure | undefined>;
  /** Frames and time; default the window's (`requestAnimationFrame`): a headless host passes its own. */
  scheduler?: Scheduler;
  /** The text version's longest line, px: its lines wrap there, centered in the band (default 0: the whole width). */
  textWidth?: number;
  /** What the sequence view may ask of the 3D view: the tracks' menus run through it; none: no track has a menu. */
  actions?: HostActions;
  /** Canvas widths; default Seq*'s monospace estimate (a surface mounted on the context measures with its canvas). */
  text?: TextMeasurer;
  /** Behaviors attached to each viewer from the start (sequence.md §3); more come with `addBehavior`. */
  behaviors?: readonly SequenceBehavior[];
}

/** The letters' size until the host's tokens are read (`setFonts`): Prism's panel text size, `--msp-font-size-sm`. */
const DefaultFontPx = 12;

const windowFrames = (): Scheduler =>
  rafScheduler({
    requestAnimationFrame: (cb) => requestAnimationFrame(cb),
    cancelAnimationFrame: (id) => cancelAnimationFrame(id),
    setTimeout: (cb, ms) => setTimeout(cb, ms) as unknown as number,
    clearTimeout: (id) => clearTimeout(id),
    now: () => performance.now(),
  });

const presetOf = (style: SequenceStyle) => (style === 'tracks' ? TracksPreset : TextPreset);

class Viewer implements SequenceViewer {
  readonly seq: SeqViewerContext;
  private readonly styleState = observableState<SequenceStyle>('tracks');
  private readonly themeState = observableState<SequenceTheme>(DefaultSequenceTheme);
  private readonly structureState = observableState<SequenceStructure | undefined>(undefined);
  private readonly legendState = observableState<SequenceLegend>(legendOf(DefaultSequenceTheme.letters, false));
  readonly style = this.styleState.asReadonly;
  readonly theme = this.themeState.asReadonly;
  readonly structure = this.structureState.asReadonly;
  readonly legend = this.legendState.asReadonly;
  private readonly railState = observableState<TracksRail>(DefaultTracksRail);
  readonly rail = this.railState.asReadonly;
  private readonly shownState = observableState<SequenceShown | undefined>(undefined);
  readonly shown = this.shownState.asReadonly;
  /** The behaviors attached, by id: their dispose and their offer. */
  private readonly behaviors = new Map<string, { dispose: () => void; offer?: SequenceOffer }>();
  private dropPreset: (() => void) | undefined;
  private fonts: SequenceFonts | undefined;
  private builds = 0;
  private disposed = false;
  private readonly subs: { unsubscribe(): void }[] = [];

  constructor(
    readonly id: string,
    private readonly deps: SequenceViewersDeps,
  ) {
    const loaders: Loaders = {
      // The structure comes from the host (`parsed`); nothing is fetched.
      load: (r) => Promise.reject(new Error(`the sequence view loads no files (${r.uri ?? r.hostFile ?? 'a file'})`)),
      parsed: async (r): Promise<FormatResult | null> => {
        const host = deps.structure.value;
        return host && r.hostFile === hostFileOf(host.key) ? { kind: 'object', value: host.source } : null;
      },
    };
    const base = DefaultSeqViewerSpec();
    const behaviors = [...base.behaviors];
    // What another view lights comes into view (text scrolls its lines in, tracks pan at the zoom): Seq*'s optional
    // `revealHighlight`, following the link channels (`links:*`) and the 3D view's hover, which the link writes as an
    // external hover owned by `3d` (seqspec:interaction#highlight).
    behaviors.push(
      SeqViewerSpec.Behavior(SeqViewerOptionalBehaviors.RevealHighlight, { channels: ['links:*'], hovers: ['3d'] }),
    );
    this.seq = createSeqViewer(
      {
        ...base,
        behaviors,
        registry: [DefaultRegistry, { svs: MolstarSvsEntry }, { svs: PrismRegistry }],
        // The text version (Seq*'s `compact` preset leaves these to the host): whole blocks per line over the full width,
        // centered, the letters in the layer's colors and no larger than the panel text (set from the host's tokens:
        // `setFonts`); `textWidth` limits the line.
        render: {
          wrappedAlign: 'center',
          wrappedColors: true,
          unobservedGray: true,
          wrappedMaxWidth: deps.textWidth ?? 0,
          wrappedFontPx: DefaultFontPx,
          maxLetterPx: DefaultFontPx,
        },
      },
      { loaders, scheduler: deps.scheduler ?? windowFrames(), ...(deps.text ? { text: deps.text } : {}) },
    );
    this.dropPreset = this.seq.register({ presets: [tracksPreset(DefaultTracksRail, !!this.deps.actions)] });
    this.applySkin(deps.ink.value);
    this.subs.push(
      deps.ink.subscribe((ink) => {
        this.applySkin(ink);
        this.refreshLegend();
      }),
      deps.structure.subscribe(() => void this.build()),
    );
    for (const b of deps.behaviors ?? []) this.addBehavior(b);
  }

  addBehavior(behavior: SequenceBehavior): () => void {
    if (this.behaviors.has(behavior.id)) throw new Error(`A sequence behavior '${behavior.id}' is attached already`);
    const entry: { dispose: () => void; offer?: SequenceOffer } = { dispose: () => {} };
    this.behaviors.set(behavior.id, entry);
    let attached = true;
    entry.dispose = behavior.attach({
      shown: this.shown,
      patch: async (shown, ops) => {
        // Only onto the document it was made for: a later one has none of the nodes it refers to.
        if (!attached || this.disposed || shown !== this.shownState.value) return false;
        const r = await this.seq.applyPatch(
          { format: 'seq_view_spec_patch', ops: [...ops] },
          { layer: 'base', origin: `host:${behavior.id}` as never },
        );
        return r.ok;
      },
      offer: (offer) => {
        if (!attached) return;
        entry.offer = offer;
        this.refreshLegend();
      },
    });
    return () => {
      if (!attached) return;
      attached = false;
      this.behaviors.delete(behavior.id);
      entry.dispose();
      this.refreshLegend();
    };
  }

  /** A track's menu: the viewer's own, or a behavior's offer. */
  private menuOf(track: string): readonly SequenceTrackMenuItem[] {
    const own = TrackMenus[track];
    if (own) return own;
    for (const { offer } of this.behaviors.values()) {
      const menu = offer?.menus?.[track];
      if (menu) return menu;
    }
    return [];
  }

  /** What the colors mean: the letters' scheme, the secondary structure, and the behaviors' blocks (a pLDDT track's bands). */
  private refreshLegend() {
    const dark = this.deps.ink.value.theme === 'dark';
    const letters = legendOf(this.themeState.value.letters, dark);
    if (!this.deps.structure.value) {
      this.legendState.set(letters);
      return;
    }
    const offered = [...this.behaviors.values()].flatMap(({ offer }) => offer?.legend ?? []);
    this.legendState.set(withTracks(letters, offered));
  }

  private applySkin(ink: CanvasInk) {
    if (!this.disposed) this.seq.setSkin(skinFor(ink, this.fonts));
  }

  trackActions(track: string, o?: { unapplied?: boolean }): readonly TrackAction[] {
    const host = this.deps.actions;
    if (!host) return [];
    const now = o?.unapplied ? host.coloring?.() : undefined;
    return this.menuOf(track)
      .filter((a) => a.theme !== now)
      .map((a) => ({
        id: a.id,
        label: a.label,
        icon: a.icon,
        run: () => void host.colorStructure(a.theme),
      }));
  }

  trackColoring(track: string): TrackColoring | undefined {
    const host = this.deps.actions;
    const a = this.menuOf(track)[0];
    if (!host || !a) return undefined;
    const active = () => host.coloring?.() === a.theme;
    // Applying the coloring the 3D view has is nothing to do: the icon marks it, a click leaves it.
    return {
      id: a.id,
      label: a.label,
      icon: a.icon,
      coloring: a.theme,
      active,
      run: () => active() || void host.colorStructure(a.theme),
    };
  }

  setFonts(fonts: SequenceFonts) {
    this.fonts = fonts;
    if (fonts.size) this.seq.setRenderProps({ wrappedFontPx: fonts.size, maxLetterPx: fonts.size });
    this.applySkin(this.deps.ink.value);
  }

  setStyle(style: SequenceStyle) {
    // Only the preset (and with it the render props) changes: one document serves both styles, so the link survives.
    this.styleState.set(style);
    this.seq.ui.setStyle(SectionRef, presetOf(style));
  }

  setRail(rail: TracksRail) {
    const now = this.railState.value;
    if (now.navigator === rail.navigator && now.height === rail.height) return;
    this.railState.set(rail);
    this.dropPreset?.();
    this.dropPreset = this.seq.register({ presets: [tracksPreset(rail, !!this.deps.actions)] });
    // The section takes the preset again (the layout reads it when the style is set).
    this.seq.ui.setStyle(SectionRef, null);
    this.seq.ui.setStyle(SectionRef, presetOf(this.styleState.value));
  }

  setTheme(theme: SequenceTheme) {
    if (theme.letters === this.themeState.value.letters) return;
    this.themeState.set(theme);
    this.refreshLegend();
    void this.build();
  }

  /** The document for what the host shows now; the latest build wins, an earlier one is dropped when it ends. */
  private async build(): Promise<void> {
    const host = this.deps.structure.value;
    const build = ++this.builds;
    // A document is being made: what the behaviors patch for the last one is dropped meanwhile (`patch` checks it).
    if (this.shownState.value) this.shownState.set(undefined);
    this.refreshLegend();
    if (!host) {
      this.structureState.set(undefined);
      this.shownState.set(undefined);
      return;
    }
    const doc = sequenceDocument({
      title: host.title,
      key: host.key,
      format: host.source.format,
      theme: this.themeState.value,
      ligands: host.ligands,
    });
    // A new structure is not shown as the old one while its document is made.
    if (this.structureState.value && this.structureState.value.structure !== host.structure)
      this.structureState.set(undefined);
    await this.seq.setDocument(doc);
    if (this.disposed || build !== this.builds) return;
    this.seq.ui.setStyle(SectionRef, presetOf(this.styleState.value));
    this.structureState.set({ title: host.title, structure: host.structure, complex: ComplexRef });
    // The behaviors add theirs to this document (sequence.md §3).
    this.shownState.set({
      key: host.key,
      title: host.title,
      structure: host.structure,
      complex: ComplexRef,
      data: DataRef,
      section: SectionRef,
      after: SecondaryTrack,
    });
  }

  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    this.builds++;
    for (const { dispose } of [...this.behaviors.values()]) dispose();
    this.behaviors.clear();
    this.dropPreset?.();
    for (const s of this.subs) s.unsubscribe();
    this.subs.length = 0;
    this.seq.dispose();
  }
}

export function createSequenceViewers(deps: SequenceViewersDeps): SequenceViewers {
  const viewers = new Map<string, Viewer>();
  const ids = observableState<readonly string[]>([]);
  const make = (id: string) => {
    const v = new Viewer(id, deps);
    viewers.set(id, v);
    ids.set([...viewers.keys()]);
    return v;
  };
  make('default');
  return {
    get(id = 'default') {
      const v = viewers.get(id);
      if (!v) throw new Error(`no sequence viewer "${id}"`);
      return v;
    },
    ids: ids.asReadonly,
    dispose() {
      for (const v of viewers.values()) v.dispose();
      viewers.clear();
      ids.set([]);
    },
  };
}
