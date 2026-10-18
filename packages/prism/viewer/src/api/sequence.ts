/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphResult } from '@molstar/seq-viewer';
import type { SequenceArchive, SequenceOverride, SequenceTrack, SequenceViewer } from '../sequence/model.js';
import type { TracksRail } from '../sequence/presets.js';
import type { LoadedSequence } from '../viewer.js';
import { LetterSchemes } from '../sequence/palette.js';
import type { SequenceSource } from '../sources/sequence.js';
import { type ApiAreaContext, ApiArea, ApiResult } from './area.js';
import { OptionalParam, StructureParam } from './params.js';
import { isFailure, labelOf, oneStructure } from './scope.js';

// `sequence` (api.md §2, sequence.md §6, §8): the one sequence viewer's style (Tracks or Text), the letters' colors,
// which structure it follows, the tracks' rail, its hidden tracks (`sequence.track.hide`, `sequence.track.show`). Its
// style, colors, rail and edits are the sequence's state, saved with snapshots (`prism-viewer.sequence`), so changing
// them is a step; following is view state. Present only when the spec has a sequence (viewer-spec.md §2.3).

const Schemes = LetterSchemes.map(
  (s) => [s ?? 'none', s === null ? 'None' : s.replace(/^prism:/, '')] as [string, string],
);

/** The tracks' rail navigators, as `sequence.rail` takes them (`TracksRail`). */
const Rails: [TracksRail['navigator'], string][] = [
  ['ruler', 'Ruler'],
  ['strip', 'Strip'],
  ['band', 'Band'],
  ['rail', 'Rail'],
  ['none', 'None'],
];

export class SequenceApi extends ApiArea {
  readonly id = 'sequence';

  constructor(
    ctx: ApiAreaContext,
    /** The sequence's runtime and viewers, loaded on first need (its own chunk, viewer-spec.md §2.3). */
    private readonly loaded: () => Promise<LoadedSequence>,
    private readonly source: SequenceSource,
  ) {
    super(ctx);
  }

  /** The viewer, once its runtime is in. */
  async #viewer(): Promise<SequenceViewer> {
    return (await this.loaded()).viewers.get();
  }

  /** The runtime's module: Seq*'s readers for the loaders (`load.svs`). */
  async runtime(): Promise<LoadedSequence['runtime']> {
    return (await this.loaded()).runtime;
  }

  readonly style = this.operation({
    id: 'sequence.style',
    title: 'Show the sequence as tracks or text',
    surface: {
      kind: 'agent-only',
      reason: "The sequence's menu and its commands (sequence.text, sequence.tracks, …) set it.",
    },
    params: {
      style: PD.Select<'tracks' | 'text'>('tracks', [
        ['tracks', 'Tracks'],
        ['text', 'Text'],
      ]),
    },
    run: async ({ style }) => {
      const viewer = await this.#viewer();
      const want = style === 'text' ? 'compact' : 'tracks';
      if (viewer.style.value === want) return ApiResult.already(`The sequence is already shown as ${style}.`);
      viewer.setStyle(want);
      return ApiResult.ok(`The sequence is shown as ${style}.`, { sequence: true });
    },
  });

  readonly colors = this.operation({
    id: 'sequence.colors',
    title: "Color the sequence's letters",
    surface: {
      kind: 'agent-only',
      reason: "The sequence's menu and its commands (sequence.text, sequence.tracks, …) set it.",
    },
    params: { scheme: PD.Select(Schemes[0][0], Schemes, { description: 'A letter color scheme' }) },
    run: async ({ scheme }) => {
      const viewer = await this.#viewer();
      const letters = scheme === 'none' ? null : scheme;
      if (viewer.theme.value.letters === letters)
        return ApiResult.already(`The letters are already colored by ${scheme}.`);
      viewer.setTheme({ letters });
      return ApiResult.ok(`The letters are colored by ${scheme}.`, { sequence: true });
    },
  });

  readonly rail = this.operation({
    id: 'sequence.rail',
    title: "Set the tracks' rail",
    surface: {
      kind: 'agent-only',
      reason: "The sequence's commands (sequence.rail.ruler, sequence.rail.strip, …) set it.",
    },
    // The rail is in the sequence's saved state (`savedOf`), so it is a step like the style and colors, and Undo
    // restores it.
    params: {
      navigator: OptionalParam(
        PD.Select<TracksRail['navigator']>('ruler', Rails, {
          description:
            "What the tracks' axis shows along the sequence: the ruler's own hairline, or a strip, a band or a rail row, or nothing; absent: kept",
        }),
      ),
      height: OptionalParam(
        PD.Numeric(
          30,
          { min: 10, max: 120, step: 1 },
          {
            description:
              "The rail row's height in px (not the ruler's or none, which have no row); absent: kept while the navigator stays, else the navigator's own",
          },
        ),
      ),
    },
    run: async ({ navigator, height }) => {
      if (navigator === undefined && height === undefined) {
        return ApiResult.fail('invalid', 'Change what?', 'It takes navigator and height.');
      }
      const viewer = await this.#viewer();
      const now = viewer.rail.value;
      const next = navigator ?? now.navigator;
      if (height !== undefined && (next === 'ruler' || next === 'none')) {
        return ApiResult.fail(
          'invalid',
          `The ${next === 'none' ? 'empty' : 'ruler'} rail has no row to give a height.`,
          'Use the strip, band or rail navigator.',
        );
      }
      // The height belongs to its navigator: a new navigator starts at its own unless one is given.
      const h = height ?? (next === now.navigator ? now.height : undefined);
      const label = Rails.find((r) => r[0] === next)?.[1].toLowerCase();
      const said = `the ${label}${h === undefined ? '' : `, ${h} px`}`;
      if (next === now.navigator && h === now.height) return ApiResult.already(`The tracks' rail is already ${said}.`);
      viewer.setRail(h === undefined ? { navigator: next } : { navigator: next, height: h });
      return ApiResult.ok(`The tracks' rail is ${said}.`, { sequence: true });
    },
  });

  readonly follow = this.operation({
    id: 'sequence.follow',
    title: 'Show the sequence of a structure',
    surface: {
      kind: 'agent-only',
      reason: "The sequence's menu and its commands (sequence.text, sequence.tracks, …) set it.",
    },
    params: {
      structure: StructureParam('A structure by label or state ref; empty: the one the person works with'),
    },
    history: false,
    run: async ({ structure }, ctx) => {
      if (!structure) {
        this.source.setFollowed(undefined);
        return ApiResult.ok('The sequence follows the structure in use.', { sequence: true });
      }
      const s = oneStructure(ctx, structure);
      if (isFailure(s)) return s;
      this.source.setFollowed(s.cell.transform.ref);
      return ApiResult.ok(`The sequence shows ${labelOf(ctx, s)}.`, { sequence: true });
    },
  });

  /** A track of the section shown by ref, else by name (case aside, one only); a failure says which there are. */
  #track(viewer: SequenceViewer, track: string): SequenceTrack | ApiResult {
    const all = viewer.tracks.value;
    const listed = all.map((t) => `${t.label} (${t.ref})`).join(', ');
    if (!all.length) return ApiResult.fail('unavailable', 'No sequence tracks are shown.');
    const byRef = all.find((t) => t.ref === track);
    if (byRef) return byRef;
    const byName = all.filter((t) => t.label.toLowerCase() === track.trim().toLowerCase());
    if (byName.length === 1) return byName[0];
    if (byName.length > 1)
      return ApiResult.fail('ambiguous', `Several tracks are called ${track}.`, `Name one by its ref: ${listed}.`);
    return ApiResult.fail('not-found', `No sequence track ${track}.`, `The tracks: ${listed}.`);
  }

  readonly hideTrack = this.operation({
    id: 'sequence.track.hide',
    title: 'Hide a sequence track',
    surface: {
      kind: 'agent-only',
      reason: "A track's menu (Hide track) and the Sequence menu's Tracks run its commands.",
    },
    // The edits are in the sequence's saved state (`savedOf`, version 3): a step, which Undo takes back (sequence.md §7).
    params: {
      track: PD.Text('', { description: "A track of the sequence's band, by ref or name (sequence.get lists them)" }),
    },
    run: async ({ track }) => {
      const viewer = await this.#viewer();
      const t = this.#track(viewer, track);
      if (ApiResult.isFailure(t)) return t;
      const { ref, label } = t as SequenceTrack;
      if (!(await viewer.setTrackHidden(ref, true))) return ApiResult.already(`The ${label} track is already hidden.`);
      return ApiResult.ok(`The ${label} track is hidden.`, { sequence: true });
    },
  });

  readonly showTrack = this.operation({
    id: 'sequence.track.show',
    title: 'Show a hidden sequence track',
    surface: {
      kind: 'agent-only',
      reason: "The Sequence menu's Tracks and the band's hidden-tracks hint run its commands.",
    },
    params: {
      track: OptionalParam(
        PD.Text('', { description: 'A hidden track, by ref or name; absent: every hidden track (show all)' }),
      ),
    },
    run: async ({ track }) => {
      const viewer = await this.#viewer();
      if (track === undefined || track === '') {
        const hidden = viewer.tracks.value.filter((t) => t.hidden);
        if (!hidden.length) return ApiResult.already('No sequence track is hidden.');
        await viewer.setTracksHidden(
          hidden.map((t) => t.ref),
          false,
        );
        const n = hidden.length;
        return ApiResult.ok(`${n === 1 ? 'The hidden track is' : `The ${n} hidden tracks are`} shown.`, {
          sequence: true,
        });
      }
      const t = this.#track(viewer, track);
      if (ApiResult.isFailure(t)) return t;
      const { ref, label } = t as SequenceTrack;
      if (!(await viewer.setTrackHidden(ref, false))) return ApiResult.already(`The ${label} track is already shown.`);
      return ApiResult.ok(`The ${label} track is shown.`, { sequence: true });
    },
  });

  /**
   * Shows a loaded SeqViewSpec over its structure (sequence.md §6): the loaders' `load.svs` calls it inside its run, so
   * it is part of that step. Not an operation of its own: the document comes from a load.
   */
  /** Shows a loaded SeqViewSpec over its structure; an `.svsx`'s files come with it (its relative URIs). */
  async showDocument(override: SequenceOverride, archive?: SequenceArchive): Promise<GraphResult> {
    const viewer = await this.#viewer();
    this.source.show(override.structure);
    return viewer.setOverride(archive ? { ...override, archive: viewer.addArchive(archive) } : override);
  }

  readonly get = this.operation({
    id: 'sequence.get',
    title: 'Describe the sequence view',
    params: {},
    read: true,
    run: async () => {
      const viewer = await this.#viewer();
      const shown = viewer.structure.value;
      const style = viewer.style.value === 'compact' ? 'text' : 'tracks';
      const override = shown ? viewer.override.value : undefined;
      const document = override
        ? (override.document.metadata?.title ?? override.url ?? 'a loaded document')
        : undefined;
      const value = {
        structure: shown?.title,
        style,
        colors: viewer.theme.value.letters ?? 'none',
        rail: viewer.rail.value,
        tracks: viewer.tracks.value.map((t) => ({ ref: t.ref, name: t.label, hidden: t.hidden })),
        ...(document ? { document } : {}),
      };
      const what = document ? `The SeqViewSpec ${document} over ${shown?.title}` : `The sequence of ${shown?.title}`;
      return ApiResult.ok(
        shown ? `${what}, as ${style}, colored by ${value.colors}.` : 'No sequence is shown.',
        {},
        { value },
      );
    },
  });
}
