/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasInk, type IconName, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { deepEqual } from '@molstar/core/util';
import type { Structure } from '@molstar/model/model/structure';
import {
  archiveLoaders,
  type DocNode,
  type FormatResult,
  type Loaders,
  type Op,
  type SeqViewSpecDoc,
} from '@molstar/seq-core/format';
import { STRUCTURE_FORMATS, type StructureSource } from './molstar/source.js';
import {
  createSeqViewer,
  type Scheduler,
  type GraphResult,
  type SeqViewerContext,
  rafScheduler,
  type TextMeasurer,
} from '@molstar/seq-viewer';
import {
  blocksOf,
  type BuiltTrack,
  ComplexRef,
  DataRef,
  DefaultSequenceTheme,
  findRef,
  hostFileOf,
  insertTracks,
  lettersColor,
  SectionRef,
  type SequenceTheme,
  sequenceDocument,
  svsDocument,
} from './document.js';
import type { SequenceBehavior, SequenceOffer, SequenceShown, SequenceTrackMenuItem } from './behaviors.js';
import { filteredUnits, sequenceChains } from './chain-filter.js';
import { colorScheme, legendOf, type SequenceLegend, withTracks } from './legend.js';
import { PrismSequenceSpec } from './spec.js';
import { DefaultTracksRail, TextPreset, type TracksRail, TracksPreset, tracksPreset } from './presets.js';
import { type SequenceFonts, skinFor, plainCellHex } from './skin.js';
import {
  applying,
  type ChainFilter,
  DefaultSequenceTracks,
  defaultIn,
  isFiltered,
  type SequenceEdits,
  type SequenceTrackContext,
  type SequenceTrackProvider,
  type StructureSequenceEdits,
  trackIn,
  withStructureEdits,
  withTrackEdit,
} from './tracks.js';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/prism-ui/reactive/frame-state';

// The sequence model (sequence.md §1): owns the Seq* contexts, one per id ('default' only for now), and keeps
// each in step with what the viewer follows: the canvas ink (its skin), the structure the 3D view shows (its
// document), and the person's choices (the style, the color theme, each structure's tracks and chains). Plain
// arguments, no DOM and no plugin: testable headless.

export type { ChainFilter, SequenceEdits, StructureSequenceEdits };

export type SequenceStyle = 'compact' | 'tracks';
export type { SequenceTheme };

/** A structure the host (the 3D view) shows and has parsed. */
export interface HostStructure {
  /**
   * Tells one structure from another: a new key rebuilds the view. A frame of the same structure (the same atoms, new
   * coordinates: a trajectory step) keeps the key with a new `structure`: the document stays (sequence.md §2).
   */
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
  /**
   * The settings of the secondary structure the 3D view draws (`secondaryAssignmentOf`): told again when they change,
   * so the tracks' data that name them follow (a patch, not a new document).
   */
  readonly assignment?: string;
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

/** A track of a structure's sequence, for the structure tree (sequence.md §8.1), the commands and the gutter's hint. */
export interface SequenceTrack {
  readonly ref: string;
  /** Its name for people: the document's label, else its provider's or the behavior's name for it, else its ref. */
  readonly label: string;
  readonly hidden: boolean;
  /** A provider's track: Remove takes it out (a loaded SeqViewSpec's own tracks only hide). */
  readonly removable: boolean;
  /** What it draws, in a few words (its provider's `describe`); none: a document's own or a behavior's track. */
  readonly description?: string;
}

/** A track provider that applies to a structure and is out: what Add track offers. */
export interface AddableTrack {
  readonly id: string;
  readonly label: string;
  readonly track: string;
}

/** A data cell the plugin holds for an address, and its parse (state refs). */
export interface DataByUrlEntry {
  /** A Download's URL, or a read file's name. */
  readonly address: string;
  /** The data cell. */
  readonly data: string;
  /** The trajectory parsed from it; absent: not parsed (or not a structure). */
  readonly parse?: string;
}

/**
 * The parse cache (sequence.md §2): the plugin's data cells by address with their parses, so a document's file the
 * plugin holds is read from it, and a structure it parsed is shared, never parsed again.
 */
export interface DataByUrl {
  entry(address: string): DataByUrlEntry | undefined;
  /** The bytes the plugin holds for an address (a string cell encoded); undefined: none. */
  bytes(address: string): Uint8Array | undefined;
  /** The structure the plugin parsed from an address, as the adapter's source over its trajectory; undefined: none. */
  structure(address: string): StructureSource | undefined;
}

/** What a viewer shows, once its document is in: the link's input. */
export interface SequenceStructure {
  /** What its edits are kept under (the structure's state ref, else its key; `structureKey`). */
  readonly ref: string;
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

/** What a host structure's edits are kept under: its state ref, else (a host without refs: tests) its key. */
export const structureKey = (host: { readonly ref?: string; readonly key: string }) => host.ref ?? host.key;

/** The origin of the edits' patches on the document's base. */
const EditsOrigin = 'host:prism-viewer.sequence.edits';

/** Every track node under `node` (a section, a group, an inserted subtree), in document order. */
function tracksUnder(node: DocNode, out: DocNode[] = []): DocNode[] {
  if (node.kind === 'track') out.push(node);
  for (const c of node.children ?? []) tracksUnder(c, out);
  return out;
}
const hiddenIn = (n: DocNode) => n.params?.hidden === true;

export interface SequenceViewer {
  readonly id: string;
  readonly seq: SeqViewerContext;
  readonly style: ReadonlyState<SequenceStyle>;
  /** UI state only (`ctx.ui.setStyle`): the document does not change. */
  setStyle(style: SequenceStyle): void;
  /** What it shows, once the document is in; undefined: nothing (no structure, or the document is being made). */
  readonly structure: ReadonlyState<SequenceStructure | undefined>;
  /**
   * The Structure the 3D view draws now: `structure`'s, or a later frame of it (same atoms); the link maps onto it. A
   * frame step changes only this (sequence.md §2).
   */
  readonly frame: ReadonlyState<Structure | undefined>;
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
  /** The tracks of the section the band shows, in order, hidden ones included; empty while a document is made. */
  readonly tracks: ReadonlyState<readonly SequenceTrack[]>;
  /** The person's edits over each structure's document, by structure (sequence.md §7.1). */
  readonly edits: ReadonlyState<SequenceEdits>;
  /** Replaces the edits (a saved state applied); resolves once the document shown has them. */
  setEdits(edits: SequenceEdits): Promise<void>;
  /**
   * Hides or shows a track, by ref: of the section shown, or of another structure (by its ref, with its Structure);
   * resolves false when no such track is there, or it is already so (nothing changes).
   */
  setTrackHidden(track: string, hidden: boolean, of?: StructureOf): Promise<boolean>;
  /**
   * Hides or shows several tracks in one edit (one patch, one update of the document); resolves how many changed (the
   * unknown and the already so are left).
   */
  setTracksHidden(tracks: readonly string[], hidden: boolean, of?: StructureOf): Promise<number>;
  /** The track providers (sequence.md §3): the viewer's own, then the extensions' as they come. */
  readonly providers: ReadonlyState<readonly SequenceTrackProvider[]>;
  /** Adds a track provider (an extension's); returns its removal. An id taken throws. */
  addTrackProvider(provider: SequenceTrackProvider): () => void;
  /** A structure's tracks (sequence.md §8.1): the shown one's from its document, another's from the providers. */
  tracksOf(of: StructureOf): readonly SequenceTrack[];
  /** What Add track offers for a structure: the providers that apply and are out. */
  addableTracks(of: StructureOf): readonly AddableTrack[];
  /** Puts a provider's track in or takes it out (by provider id or track ref); false: no such provider, or already so. */
  setTrackIn(track: string, on: boolean, of?: StructureOf): Promise<boolean>;
  /** A structure's chain filter; undefined: none. */
  chainsOf(ref: string): ChainFilter | undefined;
  /** Sets a structure's chain filter (undefined: every chain); false: already so. */
  setChains(filter: ChainFilter | undefined, of?: StructureOf): Promise<boolean>;
  /** A letter scheme's colors as a CSS gradient of hard stops, for the ink's theme; undefined: none (no colors). */
  swatchOf(scheme: string | null): string | undefined;
  /** The band shows no chain (the filter leaves none): it has nothing to draw. */
  readonly empty: ReadonlyState<boolean>;
}

/** A structure named for an edit: its state ref and the Structure the 3D view draws (what providers apply to). */
export interface StructureOf {
  readonly ref: string;
  readonly structure: Structure;
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
  /**
   * The track providers from the start (sequence.md §3; the spec's `sequence.tracks`); more come with
   * `addTrackProvider`. Absent: the viewer's own (`DefaultSequenceTracks`).
   */
  tracks?: readonly SequenceTrackProvider[];
  /** The plugin's data cells and parses by address (sequence.md §2); none: a document's files are refused. */
  dataByUrl?: DataByUrl;
}

/** The letters' size until the host's tokens are read (`setFonts`): Prism's panel text size, `--msp-font-size-sm`. */

const windowFrames = (): Scheduler =>
  rafScheduler({
    requestAnimationFrame: (cb) => requestFrameCallback(cb),
    cancelAnimationFrame: (id) => cancelFrameCallback(id),
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
  private readonly tracksState = observableState<readonly SequenceTrack[]>([], {
    compare: (a, b) =>
      a.length === b.length &&
      a.every(
        (t, i) =>
          t.ref === b[i].ref &&
          t.label === b[i].label &&
          t.hidden === b[i].hidden &&
          t.removable === b[i].removable &&
          t.description === b[i].description,
      ),
  });
  readonly tracks = this.tracksState.asReadonly;
  private readonly editsState = observableState<SequenceEdits>({}, { compare: (a, b) => a === b });
  readonly edits = this.editsState.asReadonly;
  private readonly providersState = observableState<readonly SequenceTrackProvider[]>([]);
  readonly providers = this.providersState.asReadonly;
  /** The providers' `whenChanged` listeners, by id. */
  private readonly providerSubs = new Map<string, { unsubscribe(): void }>();
  /**
   * The providers' tracks in the document shown, by track ref, with their data nodes' refs (what Remove takes out with
   * them); made with each document, kept by `sync`.
   */
  private placed = new Map<string, readonly string[]>();
  /** The document shown is a loaded SeqViewSpec's: its sections keep their blocks, the base providers stay out. */
  private svsShown = false;
  /** Where a loaded SeqViewSpec's section takes the providers' tracks: after its own last one. */
  private svsAfter: string | undefined;
  /** The sections of the document shown: the person's style applies to each. */
  private sections: readonly string[] = [SectionRef];
  /** The behaviors attached, by id: their dispose and their offer. */
  private readonly behaviors = new Map<string, { dispose: () => void; title: string; offer?: SequenceOffer }>();
  private dropPreset: (() => void) | undefined;
  private fonts: SequenceFonts | undefined;
  private builds = 0;
  /** The host key the document shown (or being made) is for: a frame of it changes only `frame`. */
  private documentKey: string | undefined;
  /** The secondary structure's settings the tracks' data were last brought to (`HostStructure.assignment`). */
  private assignment: string | undefined;
  private readonly frameState = observableState<Structure | undefined>(undefined, { compare: (a, b) => a === b });
  readonly frame = this.frameState.asReadonly;
  private disposed = false;
  /** The `.svsx` files overrides name (`archive`), for the page's life. */
  private readonly archives = new Map<string, SequenceArchive>();
  private archiveIds = 0;
  private readonly subs: { unsubscribe(): void }[] = [];
  /** The edits' patches, one after another: each reads the document the one before it left. */
  private syncing: Promise<void> = Promise.resolve();

  constructor(
    readonly id: string,
    private readonly deps: SequenceViewersDeps,
  ) {
    const refuse = (r: { uri: string | null; hostFile?: string | null }) =>
      Promise.reject(
        new Error(`the sequence view fetches no files; load ${r.uri ?? r.hostFile ?? 'it'} in the viewer first`),
      );
    /** A request's address: a URI resolved against the loaded document's own, else the host file. */
    const addressOf = (r: { uri: string | null; hostFile?: string | null }) => {
      const base = this.overrideState.value?.url;
      if (r.uri === null) return r.hostFile ?? undefined;
      try {
        return base ? new URL(r.uri, base).href : r.uri;
      } catch {
        return r.uri;
      }
    };
    // What the plugin holds (sequence.md §2): its bytes, else refused; nothing is fetched.
    const held = (r: { uri: string | null; hostFile?: string | null }): Promise<Uint8Array> => {
      const address = addressOf(r);
      const bytes = address === undefined ? undefined : deps.dataByUrl?.bytes(address);
      return bytes ? Promise.resolve(bytes) : refuse(r);
    };
    const loaders: Loaders = {
      // The structure comes from the host (`parsed`); a file the plugin holds is read from it; an `.svsx`'s relative
      // URIs from the archive; nothing is fetched.
      load: (r, signal) => {
        const id = this.overrideState.value?.archive;
        const archive = id ? this.archives.get(id) : undefined;
        if (!archive) return held(r);
        return archiveLoaders(archive.files, { load: (q) => held(q) }, archive.base).load(r, signal);
      },
      parsed: async (r, format): Promise<FormatResult | null> => {
        const host = deps.structure.value;
        if (host && r.hostFile === hostFileOf(host.key)) return { kind: 'object', value: host.source };
        // Another structure the plugin parsed (a document's second complex, a file named by URL): shared, not read again.
        if (!(STRUCTURE_FORMATS as readonly string[]).includes(format)) return null;
        const address = addressOf(r);
        const source = address === undefined ? undefined : deps.dataByUrl?.structure(address);
        return source ? { kind: 'object', value: source } : null;
      },
    };
    this.seq = createSeqViewer(PrismSequenceSpec({ textWidth: deps.textWidth }), {
      loaders,
      scheduler: deps.scheduler ?? windowFrames(),
      ...(deps.text ? { text: deps.text } : {}),
    });
    this.dropPreset = this.seq.register({ presets: [tracksPreset(DefaultTracksRail, !!this.deps.actions)] });
    this.applySkin(deps.ink.value);
    // The providers before the first document: it is made with their tracks.
    for (const p of deps.tracks ?? DefaultSequenceTracks) this.addTrackProvider(p);
    this.subs.push(
      deps.ink.subscribe((ink) => {
        this.applySkin(ink);
        this.refreshLegend();
        // No colors: the letters' cells are the surface's color, which changed.
        if (this.themeState.value.letters === null) void this.build();
      }),
      deps.structure.subscribe((host) => {
        // A frame of the structure shown (sequence.md §2): no new document; the link follows the new Structure. Its
        // secondary structure's settings changed: the tracks' data follow (`syncNow`), not at every frame.
        if (host && host.key === this.documentKey) {
          this.frameState.set(host.structure);
          if (host.assignment !== this.assignment) {
            this.assignment = host.assignment;
            void this.sync();
          }
        } else {
          this.assignment = host?.assignment;
          void this.build();
        }
      }),
      // The tracks the band shows, for the tree and the commands: read from the document shown (behaviors' patches
      // included).
      this.seq.document.shown.subscribe(() => this.refreshTracks()),
      this.structure.subscribe(() => {
        this.refreshTracks();
        this.refreshEmpty();
      }),
      this.edits.subscribe(() => this.refreshEmpty()),
    );
    for (const b of deps.behaviors ?? []) this.addBehavior(b);
  }

  /** An id is one of a behavior or a provider (one namespace, sequence.md §3). */
  private taken(id: string): boolean {
    return this.behaviors.has(id) || this.providersState.value.some((p) => p.id === id);
  }

  addBehavior(behavior: SequenceBehavior): () => void {
    if (this.taken(behavior.id)) throw new Error(`A sequence behavior '${behavior.id}' is attached already`);
    const entry: { dispose: () => void; title: string; offer?: SequenceOffer } = {
      dispose: () => {},
      title: behavior.title,
    };
    this.behaviors.set(behavior.id, entry);
    let attached = true;
    entry.dispose = behavior.attach({
      shown: this.shown,
      patch: async (shown, ops) => {
        // Only onto the document it was made for: a later one has none of the nodes it refers to.
        if (!attached || this.disposed || shown !== this.shownState.value) return false;
        // The person's edits over the tracks it adds (sequence.md §7.1): hidden from the start, no flash.
        const r = await this.seq.applyPatch(
          { format: 'seq_view_spec_patch', ops: ops.map((op) => this.withEdits(op)) },
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

  addTrackProvider(provider: SequenceTrackProvider): () => void {
    if (this.taken(provider.id)) throw new Error(`A sequence track provider '${provider.id}' is there already`);
    this.providersState.set([...this.providersState.value, provider]);
    // Its answer changed (a property computed later): its track goes in or out of the document shown.
    if (provider.whenChanged)
      this.providerSubs.set(
        provider.id,
        provider.whenChanged.subscribe(() => this.resync()),
      );
    this.resync();
    let added = true;
    return () => {
      if (!added) return;
      added = false;
      this.providerSubs.get(provider.id)?.unsubscribe();
      this.providerSubs.delete(provider.id);
      this.providersState.set(this.providersState.value.filter((p) => p !== provider));
      this.resync();
    };
  }

  /** The document shown brought to the providers and the edits, once one is in. */
  private resync() {
    if (this.shownState.value) void this.sync();
  }

  private providerOf(track: string): SequenceTrackProvider | undefined {
    const list = this.providersState.value;
    const t = track.trim().toLowerCase();
    return list.find((p) => p.id === track || p.track === track) ?? list.find((p) => p.label.toLowerCase() === t);
  }

  /** A track's menu: its provider's, or a behavior's offer. */
  private menuOf(track: string): readonly SequenceTrackMenuItem[] {
    const own = this.providersState.value.find((p) => p.track === track)?.menu;
    if (own) return own;
    for (const { offer } of this.behaviors.values()) {
      const menu = offer?.menus?.[track];
      if (menu) return menu;
    }
    return [];
  }

  /** What a track draws, by its provider (the letters' scheme named as the menu names it). */
  private describe(p: SequenceTrackProvider | undefined): string | undefined {
    if (!p?.describe) return undefined;
    const name = this.themeState.value.letters;
    const letters = name === null ? 'No colors' : (colorScheme(name)?.label ?? name);
    return p.describe({ letters });
  }

  /** A track's `description`, where its provider has one (spread into a `SequenceTrack`). */
  private describedAs(track: string): { description?: string } {
    const d = this.describe(this.providersState.value.find((p) => p.track === track));
    return d ? { description: d } : {};
  }

  /** A track's name for people: the document's label, its provider's, a behavior's offer or title, else its ref. */
  private labelOf(node: DocNode): string {
    const ref = node.ref ?? '';
    const own = typeof node.params?.label === 'string' ? node.params.label.trim() : '';
    if (own) return own;
    const provided = this.providersState.value.find((p) => p.track === ref);
    if (provided) return provided.label;
    for (const [id, { offer, title }] of this.behaviors) {
      const offered = offer?.labels?.[ref];
      if (offered) return offered;
      if (ref.startsWith(id)) return title;
    }
    return ref;
  }

  /** The loaded SeqViewSpec shown hides this track itself (its own value, under the person's edits). */
  private sourceHidden(ref: string): boolean {
    const o = this.overrideState.value;
    const host = this.deps.structure.value;
    if (!o || !host || host.ref !== o.structure) return false;
    // The override's hidden tracks, by ref, once an override (a tree walk a track was tracks × nodes).
    if (this.hiddenOf?.document !== o.document) {
      const refs = new Set<string>();
      for (const t of tracksUnder(o.document.root)) if (t.ref && hiddenIn(t)) refs.add(t.ref);
      this.hiddenOf = { document: o.document, refs };
    }
    return this.hiddenOf.refs.has(ref);
  }
  private hiddenOf: { document: unknown; refs: ReadonlySet<string> } | undefined;

  /** The ref the shown structure's edits are under; undefined: none shown. */
  private shownRef(): string | undefined {
    const host = this.deps.structure.value;
    return host ? structureKey(host) : undefined;
  }

  /** The edits of the structure shown. */
  private shownEdits(): StructureSequenceEdits | undefined {
    const ref = this.shownRef();
    return ref === undefined ? undefined : this.editsState.value[ref];
  }

  /** The edits written into a document's tracks (a fresh copy the model owns: written in place). */
  private applyEdits(root: DocNode, edits = this.shownEdits()): void {
    const tracks = edits?.tracks;
    if (!tracks) return;
    for (const t of tracksUnder(root)) {
      const entry = t.ref ? tracks[t.ref] : undefined;
      if (entry?.hidden === undefined) continue;
      const { hidden: _, ...rest } = t.params ?? {};
      t.params = entry.hidden ? { ...rest, hidden: true } : rest;
    }
  }

  /** A behavior's insert with the person's edits over the tracks it adds. */
  private withEdits(op: Op): Op {
    if (op.op !== 'insert' || !this.shownEdits()?.tracks) return op;
    const node = structuredClone(op.node);
    this.applyEdits(node);
    return { ...op, node };
  }

  /** What providers build over: the structure shown and the document's complex. */
  private trackContext(structure: Structure, complex: string): SequenceTrackContext {
    return {
      structure,
      complex,
      letters: lettersColor(this.themeState.value, plainCellHex(this.deps.ink.value)),
    };
  }

  /** A provider's track and data, or none when it throws (one provider never takes the document down). */
  private built(p: SequenceTrackProvider, ctx: SequenceTrackContext): BuiltTrack | undefined {
    try {
      return p.build(ctx);
    } catch (e) {
      console.warn(`The sequence track '${p.id}' did not build:`, e);
      return undefined;
    }
  }

  /**
   * The providers whose track belongs in the document shown, in order: they apply, the person has them in (or they
   * are in by default), and a loaded SeqViewSpec has no node of their ref.
   */
  private wanted(structure: Structure, edits: StructureSequenceEdits | undefined, doc?: DocNode) {
    return applying(this.providersState.value, structure).filter(
      (p) => trackIn(p, edits, this.svsShown) && (!doc || this.placed.has(p.track) || !findRef(doc, p.track)),
    );
  }

  /**
   * The document shown brought to the providers and the edits (sequence.md §7.1): one base patch for what differs, the
   * providers' tracks put in or taken out, the chain filter's blocks, the tracks hidden or shown. Nothing while a
   * document is made (its build writes them).
   */
  private sync(): Promise<void> {
    this.syncing = this.syncing.then(() => this.syncNow()).catch(() => undefined);
    return this.syncing;
  }

  private async syncNow(): Promise<void> {
    const base = this.seq.document.base.value;
    const shown = this.shownState.value;
    const host = this.deps.structure.value;
    if (!base || !shown || !host || this.disposed) return;
    const edits = this.shownEdits();
    const ops: Op[] = [];
    const section = findRef(base.root, shown.section);
    if (!section) return;
    // The providers' tracks: out first (their data with them), then in, each after the one before it in order.
    const wanted = this.wanted(host.structure, edits, base.root);
    const keep = new Set(wanted.map((p) => p.track));
    const placed = new Map(this.placed);
    for (const [track, data] of this.placed) {
      if (keep.has(track)) continue;
      if (findRef(base.root, track)) ops.push({ op: 'remove', ref: track });
      for (const d of data) if (findRef(base.root, d)) ops.push({ op: 'remove', ref: d });
      placed.delete(track);
    }
    const ctx = this.trackContext(host.structure, shown.complex);
    const blocks = (section.children ?? []).filter((n) => n.kind === 'block' && n.ref).map((n) => n.ref!);
    let anchor: string = this.svsShown ? (this.svsAfter ?? 'first') : (blocks[blocks.length - 1] ?? 'first');
    for (const p of wanted) {
      if (placed.has(p.track)) {
        anchor = p.track;
        continue;
      }
      const b = this.built(p, ctx);
      if (!b) continue;
      for (const d of b.data ?? []) ops.push({ op: 'insert', parent: shown.data, after: 'last', node: d });
      const node = structuredClone(b.track);
      this.applyEdits(node, edits);
      ops.push({ op: 'insert', parent: shown.section, after: anchor, node });
      placed.set(
        p.track,
        (b.data ?? []).flatMap((d) => (d.ref ? [d.ref] : [])),
      );
      anchor = p.track;
    }
    // The tracks already in: their data as the providers build it now (an option can follow the 3D view, as the secondary
    // structure's assignment): a param that differs is set, which recomputes the node.
    // Those the document was built with (a base provider's) as those a sync put in; not the ones this sync puts in.
    for (const p of applying(this.providersState.value, host.structure)) {
      if (!trackIn(p, edits, this.svsShown) || placed.get(p.track) !== this.placed.get(p.track)) continue;
      // Only the viewer's own nodes: under a loaded SeqViewSpec, those a sync put in (its own are the document's).
      if (!findRef(base.root, p.track) || (this.svsShown && !this.placed.has(p.track))) continue;
      const b = this.built(p, ctx);
      for (const d of b?.data ?? []) {
        const n = d.ref ? findRef(base.root, d.ref) : undefined;
        if (!n || !d.params) continue;
        const changed = Object.keys(d.params).filter((k) => !deepEqual(n.params?.[k], d.params![k]));
        if (changed.length)
          ops.push({ op: 'set', ref: d.ref!, params: Object.fromEntries(changed.map((k) => [k, d.params![k]])) });
      }
    }
    // The chain filter (not under a loaded SeqViewSpec: its sections keep their own blocks).
    if (!this.svsShown) {
      const units = filteredUnits(edits?.chains, sequenceChains(host.structure));
      const want = units ? blocksOf(units) : [];
      const same = want.length === blocks.length && want.every((b, i) => b.ref === blocks[i]);
      if (!same) {
        for (const r of blocks) ops.push({ op: 'remove', ref: r });
        let after = 'first';
        for (const b of want) {
          ops.push({ op: 'insert', parent: shown.section, after, node: b });
          after = b.ref!;
        }
      }
      const each = !!edits?.chains?.each;
      if (each !== (section.params?.copies === 'each'))
        ops.push({ op: 'set', ref: shown.section, params: { copies: each ? 'each' : null } });
    }
    // The tracks hidden or shown. The override's hidden tracks read afresh for this sync (once, not once a track).
    this.hiddenOf = undefined;
    for (const t of tracksUnder(section)) {
      if (!t.ref || !findRef(base.root, t.ref)) continue;
      if (!keep.has(t.ref) && this.placed.has(t.ref)) continue;
      const entry = edits?.tracks?.[t.ref]?.hidden;
      const want = entry ?? this.sourceHidden(t.ref);
      if (want !== hiddenIn(t)) ops.push({ op: 'set', ref: t.ref, params: { hidden: want ? true : null } });
    }
    this.placed = placed;
    if (!ops.length) return;
    const r = await this.seq.applyPatch(
      { format: 'seq_view_spec_patch', ops },
      { layer: 'base', origin: EditsOrigin as never },
    );
    // Refused (a node the patch named went meanwhile): the document made again from the edits, once: inside this sync,
    // so the build does not queue a sync behind the one running (it would wait for itself).
    if (!r.ok && !this.disposed) {
      console.warn('The sequence edits were refused; the document is made again:', r.issues);
      await this.build({ sync: false });
    } else {
      // `shown` stays as it is: a new one would tell behaviors this is another document (they add their nodes again).
      this.refreshLegend();
      this.refreshTracks();
    }
  }

  /** The tracks of the section shown, from the document shown. */
  private refreshTracks() {
    const doc = this.seq.document.shown.value;
    const section = this.structureState.value?.section;
    const node = doc && section ? findRef(doc.root, section) : undefined;
    const list = node
      ? tracksUnder(node).flatMap((t) =>
          t.ref
            ? [
                {
                  ref: t.ref,
                  label: this.labelOf(t),
                  hidden: hiddenIn(t),
                  removable: this.placed.has(t.ref),
                  ...this.describedAs(t.ref),
                },
              ]
            : [],
        )
      : [];
    this.tracksState.set(list);
  }

  setEdits(edits: SequenceEdits): Promise<void> {
    this.editsState.set(edits);
    this.refreshLegend();
    return this.sync();
  }

  /** Whether `of` names the structure shown (none named: it does). */
  private isShown(of: StructureOf | undefined): boolean {
    return !of || of.ref === this.shownRef();
  }

  tracksOf(of: StructureOf): readonly SequenceTrack[] {
    if (this.isShown(of) && this.structureState.value) return this.tracksState.value;
    const edits = this.editsState.value[of.ref];
    return applying(this.providersState.value, of.structure)
      .filter((p) => trackIn(p, edits, false))
      .map((p) => ({
        ref: p.track,
        label: p.label,
        hidden: edits?.tracks?.[p.track]?.hidden ?? false,
        removable: true,
        ...this.describedAs(p.track),
      }));
  }

  addableTracks(of: StructureOf): readonly AddableTrack[] {
    const shown = this.isShown(of) && !!this.structureState.value;
    const edits = this.editsState.value[of.ref];
    const doc = shown ? this.seq.document.base.value?.root : undefined;
    return applying(this.providersState.value, of.structure)
      .filter((p) => (shown ? !this.placed.has(p.track) && !(doc && findRef(doc, p.track)) : !trackIn(p, edits, false)))
      .map((p) => ({ id: p.id, label: p.label, track: p.track }));
  }

  async setTrackIn(track: string, on: boolean, of?: StructureOf): Promise<boolean> {
    const p = this.providerOf(track);
    const host = this.deps.structure.value;
    const ref = of?.ref ?? this.shownRef();
    const structure = of?.structure ?? host?.structure;
    if (!p || ref === undefined || !structure) return false;
    if (on && p.applies && !p.applies(structure)) return false;
    const shown = this.isShown(of) && !!this.structureState.value;
    const edits = this.editsState.value[ref];
    const now = shown ? this.placed.has(p.track) : trackIn(p, edits, false);
    if (now === on) return false;
    if (shown && on && this.seq.document.base.value && findRef(this.seq.document.base.value.root, p.track))
      return false;
    const svs = shown && this.svsShown;
    // Only a change from the default is an edit; taken out, its eye goes with it.
    const entry = {
      ...(on === defaultIn(p, svs) ? {} : { in: on }),
      ...(on ? { hidden: edits?.tracks?.[p.track]?.hidden } : {}),
    };
    const clean = Object.fromEntries(Object.entries(entry).filter(([, v]) => v !== undefined));
    this.editsState.set(withStructureEdits(this.editsState.value, ref, withTrackEdit(edits, p.track, clean)));
    if (shown) await this.sync();
    return true;
  }

  async setTrackHidden(track: string, hidden: boolean, of?: StructureOf): Promise<boolean> {
    return (await this.setTracksHidden([track], hidden, of)) > 0;
  }

  async setTracksHidden(refs: readonly string[], hidden: boolean, of?: StructureOf): Promise<number> {
    const shown = this.isShown(of);
    const host = this.deps.structure.value;
    const target = of ?? (host ? { ref: structureKey(host), structure: host.structure } : undefined);
    if (!target) return 0;
    const ref = target.ref;
    // The structure shown, its document in: its tracks as shown; else (another, or one being made) its providers'.
    const list = shown && this.structureState.value ? this.tracksState.value : this.tracksOf(target);
    const now = new Map(list.map((t) => [t.ref, t.hidden]));
    let edits = this.editsState.value[ref];
    let changed = 0;
    for (const t of new Set(refs)) {
      if (!now.has(t) || now.get(t) === hidden) continue;
      changed++;
      // Only a change from the document's own value is an edit: showing a track the person hid drops its entry.
      const { hidden: _, ...rest } = edits?.tracks?.[t] ?? {};
      edits = withTrackEdit(edits, t, hidden === (shown && this.sourceHidden(t)) ? rest : { ...rest, hidden });
    }
    if (!changed) return 0;
    this.editsState.set(withStructureEdits(this.editsState.value, ref, edits ?? {}));
    if (shown) await this.sync();
    return changed;
  }

  swatchOf(scheme: string | null): string | undefined {
    if (scheme === null) return undefined;
    const legend = legendOf(scheme, this.deps.ink.value.theme === 'dark');
    const colors = legend.items.map((i) => i.color).slice(0, 8);
    if (!colors.length) {
      const ramp = legend.ramps[0];
      return ramp
        ? `linear-gradient(to right, ${ramp.stops.map((s) => `${s.color} ${s.t * 100}%`).join(', ')})`
        : undefined;
    }
    const w = 100 / colors.length;
    return `linear-gradient(to right, ${colors.map((c, i) => `${c} ${i * w}% ${(i + 1) * w}%`).join(', ')})`;
  }

  /** The structure shown leaves out every chain. */
  private readonly emptyState = observableState(false);
  readonly empty = this.emptyState.asReadonly;
  private refreshEmpty() {
    const host = this.deps.structure.value;
    const f = host && !this.svsShown ? this.editsState.value[structureKey(host)]?.chains : undefined;
    const units = host && f ? filteredUnits(f, sequenceChains(host.structure)) : undefined;
    this.emptyState.set(!!units && units.length === 0);
  }

  chainsOf(ref: string): ChainFilter | undefined {
    return this.editsState.value[ref]?.chains;
  }

  async setChains(filter: ChainFilter | undefined, of?: StructureOf): Promise<boolean> {
    const ref = of?.ref ?? this.shownRef();
    if (ref === undefined) return false;
    const next: ChainFilter | undefined = isFiltered(filter)
      ? {
          ...(filter.polymer ? { polymer: true } : {}),
          ...(filter.each ? { each: true } : {}),
          ...(filter.chains !== undefined ? { chains: [...filter.chains] } : {}),
        }
      : undefined;
    const edits = this.editsState.value[ref];
    if (JSON.stringify(edits?.chains ?? null) === JSON.stringify(next ?? null)) return false;
    const { chains: _, ...rest } = edits ?? {};
    this.editsState.set(withStructureEdits(this.editsState.value, ref, next ? { ...rest, chains: next } : rest));
    if (this.isShown(of)) await this.sync();
    return true;
  }

  /**
   * What the colors mean: the letters' scheme, the secondary structure, the providers' blocks for their tracks in, and
   * the behaviors' (a pLDDT track's bands).
   */
  private refreshLegend() {
    const dark = this.deps.ink.value.theme === 'dark';
    const letters = legendOf(this.themeState.value.letters, dark);
    if (!this.deps.structure.value) {
      this.legendState.set(letters);
      return;
    }
    const provided = this.providersState.value.filter((p) => this.placed.has(p.track)).flatMap((p) => p.legend ?? []);
    const offered = [...this.behaviors.values()].flatMap(({ offer }) => offer?.legend ?? []);
    this.legendState.set(withTracks(letters, [...provided, ...offered]));
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
    const from = this.styleState.value;
    this.styleState.set(style);
    for (const s of this.sections) this.seq.ui.setStyle(s, presetOf(style));
    // Back from Text, Tracks show the whole sequence: the text's scrolling and reveals moved the same viewport.
    if (from === 'compact' && style === 'tracks') for (const s of this.sections) this.seq.viewports.get(s)?.fit({});
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
  private async build(o: { readonly sync?: boolean } = {}): Promise<GraphResult> {
    const host = this.deps.structure.value;
    const build = ++this.builds;
    // A document is being made: what the behaviors patch for the last one is dropped meanwhile (`patch` checks it).
    if (this.shownState.value) this.shownState.set(undefined);
    this.refreshLegend();
    this.documentKey = host?.key;
    if (!host) {
      this.placed = new Map();
      this.structureState.set(undefined);
      this.frameState.set(undefined);
      this.shownState.set(undefined);
      return { ok: true, issues: [] };
    }
    const ref = structureKey(host);
    const edits = this.editsState.value[ref];
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
    const svsShown = !!svs;
    // The providers' tracks in (sequence.md §3): their data under the data node, their tracks in order.
    const ctx = this.trackContext(host.structure, svs?.complex ?? ComplexRef);
    const placed = new Map<string, readonly string[]>();
    const tracks: BuiltTrack[] = [];
    const before = this.svsShown;
    this.svsShown = svsShown;
    for (const p of applying(this.providersState.value, host.structure)) {
      if (!trackIn(p, edits, svsShown) || (svs && findRef(svs.doc.root, p.track))) continue;
      const b = this.built(p, ctx);
      if (!b) continue;
      tracks.push(b);
      placed.set(
        p.track,
        (b.data ?? []).flatMap((d) => (d.ref ? [d.ref] : [])),
      );
    }
    this.svsShown = before;
    let doc: SeqViewSpecDoc;
    if (svs) {
      const data = findRef(svs.doc.root, svs.data);
      const section = findRef(svs.doc.root, svs.section);
      if (data && section) insertTracks(data, section, tracks, svs.after);
      doc = svs.doc;
    } else {
      const units = filteredUnits(edits?.chains, sequenceChains(host.structure));
      doc = sequenceDocument({
        title: host.title,
        key: host.key,
        format: host.source.format,
        ligands: host.ligands,
        ...(units ? { units } : {}),
        ...(edits?.chains?.each ? { each: true } : {}),
        tracks,
      });
    }
    // The person's edits (sequence.md §7.1): written into every document built, so a rebuild keeps them.
    this.applyEdits(doc.root, edits);
    // A new structure is not shown as the old one while its document is made.
    if (this.structureState.value && this.structureState.value.structure !== host.structure)
      this.structureState.set(undefined);
    const result = await this.seq.setDocument(doc);
    if (this.disposed || build !== this.builds) return result;
    this.placed = placed;
    this.svsShown = svsShown;
    this.svsAfter = svs?.after;
    // The person's style over the document's sections (a loaded SeqViewSpec's presets too).
    this.sections = svs?.sections ?? [SectionRef];
    for (const s of this.sections) this.seq.ui.setStyle(s, presetOf(this.styleState.value));
    const complex = svs?.complex ?? ComplexRef;
    const section = svs?.section ?? SectionRef;
    // The frame shown now (frames stepped while the document was made), before the link reads it.
    const now = this.deps.structure.value;
    this.frameState.set(now?.key === host.key ? now.structure : host.structure);
    this.structureState.set({ ref, title: host.title, structure: host.structure, complex, section });
    // The behaviors add theirs to this document (sequence.md §3), after the section's last track.
    const shownSection = findRef(doc.root, section);
    const lastTrack = (shownSection?.children ?? []).filter((n) => n.kind === 'track').pop()?.ref;
    this.shownState.set({
      key: host.key,
      title: host.title,
      structure: host.structure,
      complex,
      data: svs?.data ?? DataRef,
      section,
      ...(lastTrack ? { after: lastTrack } : {}),
    });
    this.refreshLegend();
    // Edits changed while the document was made: brought to them now (not from inside a sync: `syncNow`'s refusal).
    if (o.sync !== false) await this.sync();
    return result;
  }

  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    this.builds++;
    for (const { dispose } of [...this.behaviors.values()]) dispose();
    this.behaviors.clear();
    for (const sub of this.providerSubs.values()) sub.unsubscribe();
    this.providerSubs.clear();
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
