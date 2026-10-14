/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasInk, observableState, type ReadonlyState } from '@molstar/prism-ui';
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
  type TextMeasurer,
} from '@molstar/seq-viewer';
import {
  ComplexRef,
  DefaultSequenceTheme,
  hostFileOf,
  SectionRef,
  type SequenceTheme,
  sequenceDocument,
} from './document.js';
import { legendOf, type SequenceLegend } from './legend.js';
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
}

/** What a viewer shows, once its document is in: the link's input. */
export interface SequenceStructure {
  readonly title: string;
  readonly structure: Structure;
  /** The Seq* complex made from it. */
  readonly complex: string;
}

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
  /** The host's fonts (the control reads its `--msp-font-*` once): part of the skin, not of the ink. */
  setFonts(fonts: SequenceFonts): void;
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
  /** Canvas widths; default Seq*'s monospace estimate (a surface mounted on the context measures with its canvas). */
  text?: TextMeasurer;
}

const windowFrames = (): Scheduler =>
  rafScheduler({
    requestAnimationFrame: (cb) => requestAnimationFrame(cb),
    cancelAnimationFrame: (id) => cancelAnimationFrame(id),
    setTimeout: (cb, ms) => setTimeout(cb, ms) as unknown as number,
    clearTimeout: (id) => clearTimeout(id),
    now: () => performance.now(),
  });

class Viewer implements SequenceViewer {
  readonly seq: SeqViewerContext;
  private readonly styleState = observableState<SequenceStyle>('compact');
  private readonly themeState = observableState<SequenceTheme>(DefaultSequenceTheme);
  private readonly structureState = observableState<SequenceStructure | undefined>(undefined);
  private readonly legendState = observableState<SequenceLegend>(legendOf(DefaultSequenceTheme.letters, false));
  readonly style = this.styleState.asReadonly;
  readonly theme = this.themeState.asReadonly;
  readonly structure = this.structureState.asReadonly;
  readonly legend = this.legendState.asReadonly;
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
    this.seq = createSeqViewer(
      { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: MolstarSvsEntry }] },
      { loaders, scheduler: deps.scheduler ?? windowFrames(), ...(deps.text ? { text: deps.text } : {}) },
    );
    this.applySkin(deps.ink.value);
    this.subs.push(
      deps.ink.subscribe((ink) => {
        this.applySkin(ink);
        this.legendState.set(legendOf(this.themeState.value.letters, ink.theme === 'dark'));
      }),
      deps.structure.subscribe(() => void this.build()),
    );
  }

  private applySkin(ink: CanvasInk) {
    if (!this.disposed) this.seq.setSkin(skinFor(ink, this.fonts));
  }

  setFonts(fonts: SequenceFonts) {
    this.fonts = fonts;
    this.applySkin(this.deps.ink.value);
  }

  setStyle(style: SequenceStyle) {
    this.styleState.set(style);
    this.seq.ui.setStyle(SectionRef, style);
  }

  setTheme(theme: SequenceTheme) {
    if (theme.letters === this.themeState.value.letters) return;
    this.themeState.set(theme);
    this.legendState.set(legendOf(theme.letters, this.deps.ink.value.theme === 'dark'));
    void this.build();
  }

  /** The document for what the host shows now; the latest build wins, an earlier one is dropped when it ends. */
  private async build(): Promise<void> {
    const host = this.deps.structure.value;
    const build = ++this.builds;
    if (!host) {
      this.structureState.set(undefined);
      return;
    }
    const doc = sequenceDocument({
      title: host.title,
      key: host.key,
      format: host.source.format,
      theme: this.themeState.value,
    });
    // A new structure is not shown as the old one while its document is made.
    if (this.structureState.value && this.structureState.value.structure !== host.structure)
      this.structureState.set(undefined);
    await this.seq.setDocument(doc);
    if (this.disposed || build !== this.builds) return;
    this.seq.ui.setStyle(SectionRef, this.styleState.value);
    this.structureState.set({ title: host.title, structure: host.structure, complex: ComplexRef });
  }

  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    this.builds++;
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
