/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasInk, type IconName, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { Structure } from '@molstar/model/model/structure';
import { archiveLoaders, type FormatResult, type Loaders, type SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { StructureSource } from '@molstar/seq-molstar-poc';
import {
  createSeqViewer,
  type Scheduler,
  type GraphResult,
  type SeqViewerContext,
  rafScheduler,
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
  svsDocument,
} from './document.js';
import type { SequenceBehavior, SequenceOffer, SequenceShown, SequenceTrackMenuItem } from './behaviors.js';
import { legendOf, type SequenceLegend, withTracks } from './legend.js';
import { PrismSequenceSpec } from './spec.js';
import { DefaultTracksRail, TextPreset, type TracksRail, TracksPreset, tracksPreset } from './presets.js';
import { type SequenceFonts, skinFor, plainCellHex } from './skin.js';

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
  /** The structure's state ref (a loaded SeqViewSpec names the one it maps onto by it); absent: none. */
  readonly ref?: string;
  /** What the person calls it ("1TQN"). */
  readonly title: string;
  /** The host's parse (`StructureSource.of(trajectory, format)`): shared with the 3D view, never read again. */
  readonly source: StructureSource;
  /** The structure the 3D view draws: what the link maps. */
  readonly structure: Structure;
  /** The ligands the 3D view draws, by chemical component id; absent: every ligand. The sequence lists the same ones. */
  readonly ligands?: readonly string[];
}

/**
 * A SeqViewSpec document shown instead of the one the model builds (`load.svs`, sequence.md §6; David, 2026-10-15):
 * while the followed structure is the one it maps its first `molstar:structure` complex onto. Saved with the sequence's
 * state (§7.1); the viewer ends it when the person focuses or loads another structure.
 */
export interface SequenceOverride {
  readonly document: SeqViewSpecDoc;
  /** The loaded structure (its state ref) the document's first `molstar:structure` complex maps onto. */
  readonly structure: string;
  /** Where it came from (its address); absent: given as data. */
  readonly url?: string;
  /**
   * The `.svsx` it came from (`addArchive`'s id): its relative URIs resolve inside it. The files stay in memory for the
   * page's life, out of the saved state; after a reload the document's archived data does not load.
   */
  readonly archive?: string;
}

/** An `.svsx`'s files and its main file's directory (relative URIs resolve against it). */
export interface SequenceArchive {
  readonly files: ReadonlyMap<string, Uint8Array>;
  readonly base: string;
}

/** What a viewer shows, once its document is in: the link's input. */
export interface SequenceStructure {
  readonly title: string;
  readonly structure: Structure;
  /** The Seq* complex made from it. */
  readonly complex: string;
  /** The section the band shows: the built document's, or a loaded SeqViewSpec's first with tracks (sequence.md §6). */
  readonly section: string;
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
  /** The loaded SeqViewSpec shown instead of the built document; undefined: none. */
  readonly override: ReadonlyState<SequenceOverride | undefined>;
  /**
   * Shows a loaded SeqViewSpec over its structure (when that is the followed one), or the built document again
   * (undefined); resolves with what Seq* made of the document shown (an override not shown yet: ok).
   */
  setOverride(override: SequenceOverride | undefined): Promise<GraphResult>;
  /** Keeps an `.svsx`'s files for an override that names it (`archive`); returns its id. */
  addArchive(archive: SequenceArchive): string;
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
  // Each document is its own value (identity): the 3D link is made again for it, though the structure and the complex
  // ref are the same, as its Seq* set is not.
  private readonly structureState = observableState<SequenceStructure | undefined>(undefined, {
    compare: (a, b) => a === b,
  });
  private readonly legendState = observableState<SequenceLegend>(legendOf(DefaultSequenceTheme.letters, false));
  readonly style = this.styleState.asReadonly;
  readonly theme = this.themeState.asReadonly;
  readonly structure = this.structureState.asReadonly;
  readonly legend = this.legendState.asReadonly;
  private readonly railState = observableState<TracksRail>(DefaultTracksRail);
  readonly rail = this.railState.asReadonly;
  private readonly shownState = observableState<SequenceShown | undefined>(undefined);
  readonly shown = this.shownState.asReadonly;
  private readonly overrideState = observableState<SequenceOverride | undefined>(undefined, {
    compare: (a, b) => a === b || (!!a && !!b && a.document === b.document && a.structure === b.structure),
  });
  readonly override = this.overrideState.asReadonly;
  /** The sections of the document shown: the person's style applies to each. */
  private sections: readonly string[] = [SectionRef];
  /** The behaviors attached, by id: their dispose and their offer. */
  private readonly behaviors = new Map<string, { dispose: () => void; offer?: SequenceOffer }>();
  private dropPreset: (() => void) | undefined;
  private fonts: SequenceFonts | undefined;
  private builds = 0;
  private disposed = false;
  /** The `.svsx` files overrides name (`archive`), for the page's life. */
  private readonly archives = new Map<string, SequenceArchive>();
  private archiveIds = 0;
  private readonly subs: { unsubscribe(): void }[] = [];

  constructor(
    readonly id: string,
    private readonly deps: SequenceViewersDeps,
  ) {
    const refuse = (r: { uri: string | null; hostFile?: string | null }) =>
      Promise.reject(new Error(`the sequence view loads no files (${r.uri ?? r.hostFile ?? 'a file'})`));
    const loaders: Loaders = {
      // The structure comes from the host (`parsed`); nothing is fetched, but what an `.svsx` holds (its relative URIs).
      load: (r, signal) => {
        const id = this.overrideState.value?.archive;
        const archive = id ? this.archives.get(id) : undefined;
        if (!archive) return refuse(r);
        return archiveLoaders(archive.files, { load: (q) => refuse(q) }, archive.base).load(r, signal);
      },
      parsed: async (r): Promise<FormatResult | null> => {
        const host = deps.structure.value;
        return host && r.hostFile === hostFileOf(host.key) ? { kind: 'object', value: host.source } : null;
      },
    };
    this.seq = createSeqViewer(PrismSequenceSpec({ textWidth: deps.textWidth }), {
      loaders,
      scheduler: deps.scheduler ?? windowFrames(),
      ...(deps.text ? { text: deps.text } : {}),
    });
    this.dropPreset = this.seq.register({ presets: [tracksPreset(DefaultTracksRail, !!this.deps.actions)] });
    this.applySkin(deps.ink.value);
    this.subs.push(
      deps.ink.subscribe((ink) => {
        this.applySkin(ink);
        this.refreshLegend();
        // No colors: the letters' cells are the surface's color, which changed.
        if (this.themeState.value.letters === null) void this.build();
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
    for (const s of this.sections) this.seq.ui.setStyle(s, presetOf(style));
  }

  setRail(rail: TracksRail) {
    const now = this.railState.value;
    if (now.navigator === rail.navigator && now.height === rail.height) return;
    this.railState.set(rail);
    this.dropPreset?.();
    this.dropPreset = this.seq.register({ presets: [tracksPreset(rail, !!this.deps.actions)] });
    // The section takes the preset again (the layout reads it when the style is set).
    for (const s of this.sections) {
      this.seq.ui.setStyle(s, null);
      this.seq.ui.setStyle(s, presetOf(this.styleState.value));
    }
  }

  setTheme(theme: SequenceTheme) {
    if (theme.letters === this.themeState.value.letters) return;
    this.themeState.set(theme);
    // Text mode: the letters in the scheme's colors, or plain text.
    this.seq.setRenderProps({ wrappedColors: theme.letters !== null });
    this.refreshLegend();
    void this.build();
  }

  addArchive(archive: SequenceArchive): string {
    const id = `svsx-${++this.archiveIds}`;
    this.archives.set(id, archive);
    return id;
  }

  setOverride(override: SequenceOverride | undefined): Promise<GraphResult> {
    const now = this.overrideState.value;
    this.overrideState.set(override);
    // Only a change to what the followed structure shows rebuilds.
    const host = this.deps.structure.value;
    const applies = (o: SequenceOverride | undefined) => !!o && !!host && host.ref === o.structure;
    if (this.overrideState.value === now || (!applies(now) && !applies(override)))
      return Promise.resolve({ ok: true, issues: [] });
    return this.build();
  }

  /** The document for what the host shows now; the latest build wins, an earlier one is dropped when it ends. */
  private async build(): Promise<GraphResult> {
    const host = this.deps.structure.value;
    const build = ++this.builds;
    // A document is being made: what the behaviors patch for the last one is dropped meanwhile (`patch` checks it).
    if (this.shownState.value) this.shownState.set(undefined);
    this.refreshLegend();
    if (!host) {
      this.structureState.set(undefined);
      this.shownState.set(undefined);
      return { ok: true, issues: [] };
    }
    // A loaded SeqViewSpec over its structure (sequence.md §6), else the document built for the structure.
    const override = this.overrideState.value;
    const svs =
      override && host.ref === override.structure
        ? svsDocument(override.document, {
            key: host.key,
            format: host.source.format,
            theme: this.themeState.value,
            plain: plainCellHex(this.deps.ink.value),
          })
        : undefined;
    const doc =
      svs?.doc ??
      sequenceDocument({
        title: host.title,
        key: host.key,
        format: host.source.format,
        theme: this.themeState.value,
        plain: plainCellHex(this.deps.ink.value),
        ligands: host.ligands,
      });
    // A new structure is not shown as the old one while its document is made.
    if (this.structureState.value && this.structureState.value.structure !== host.structure)
      this.structureState.set(undefined);
    const result = await this.seq.setDocument(doc);
    if (this.disposed || build !== this.builds) return result;
    // The person's style over the document's sections (a loaded SeqViewSpec's presets too).
    this.sections = svs?.sections ?? [SectionRef];
    for (const s of this.sections) this.seq.ui.setStyle(s, presetOf(this.styleState.value));
    const complex = svs?.complex ?? ComplexRef;
    const section = svs?.section ?? SectionRef;
    this.structureState.set({ title: host.title, structure: host.structure, complex, section });
    // The behaviors add theirs to this document (sequence.md §3).
    this.shownState.set({
      key: host.key,
      title: host.title,
      structure: host.structure,
      complex,
      data: svs?.data ?? DataRef,
      section,
      after: svs?.after ?? SecondaryTrack,
    });
    return result;
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
