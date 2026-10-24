/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphResult } from '@molstar/seq-viewer';
import type { Structure } from '@molstar/model/model/structure';
import { sequenceChains } from '../sequence/chain-filter.js';
import type {
  ChainFilter,
  SequenceArchive,
  SequenceOverride,
  SequenceTrack,
  SequenceViewer,
  StructureOf,
} from '../sequence/model.js';
import type { TracksRail } from '../sequence/presets.js';
import type { LoadedSequence } from '../viewer.js';
import { LetterSchemes } from '../sequence/palette.js';
import type { SequenceSource } from '../sources/sequence.js';
import { type ApiAreaContext, ApiArea, type ApiContext, ApiResult } from './area.js';
import { OptionalParam, StructureParam } from './params.js';
import { isFailure, labelOf, oneStructure } from './scope.js';

// `sequence` (api.md §2, sequence.md §6, §8): the one sequence viewer's style (Tracks or Text), the letters' colors,
// which structure it follows, the tracks' rail, and each structure's tracks (`sequence.track.hide`, `.show`, `.add`,
// `.remove`) and chains (`sequence.chains`). Its style, colors, rail and edits are the sequence's state, saved with
// snapshots (`prism-viewer.sequence`), so changing them is a step; following is view state. Present only when the spec
// has a sequence (viewer-spec.md §2.3).

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
      reason: "The eye of a structure's Sequence node in the structure tree and its Follow automatically set it.",
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

  /** The structure named (another than the one shown: its own tracks and chains); none named: the one shown. */
  #of(ctx: ApiContext, structure: string | undefined): StructureOf | undefined | ApiResult {
    if (!structure) return undefined;
    const s = oneStructure(ctx, structure);
    if (isFailure(s)) return s;
    const data = s.cell.obj?.data as Structure | undefined;
    if (!data) return ApiResult.fail('unavailable', `${labelOf(ctx, s)} is not built yet.`);
    return { ref: s.cell.transform.ref, structure: data };
  }

  /**
   * The structure an edit is for: the one named, else the one the sequence follows (its document may still be being
   * made: the edit is in it once it is); a failure when there is none.
   */
  #target(of: StructureOf | undefined): StructureOf | ApiResult {
    if (of) return of;
    const host = this.source.structure.value;
    if (!host?.ref) return ApiResult.fail('unavailable', 'No sequence is shown.', 'Name a structure.');
    return { ref: host.ref, structure: host.structure };
  }

  /**
   * A track of a structure by ref, else by name (case aside, one only); a failure says which there are. None named:
   * the one shown, or while its document is made, the one the sequence follows.
   */
  #track(viewer: SequenceViewer, track: string, of?: StructureOf): SequenceTrack | ApiResult {
    const host = this.source.structure.value;
    const pending =
      !of && !viewer.structure.value && host?.ref ? { ref: host.ref, structure: host.structure } : undefined;
    const target = of ?? pending;
    const all = target ? viewer.tracksOf(target) : viewer.tracks.value;
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
      structure: StructureParam('A structure by label or state ref; empty: the one the sequence shows'),
    },
    run: async ({ track, structure }, ctx) => {
      const viewer = await this.#viewer();
      const of = this.#of(ctx, structure);
      if (ApiResult.isFailure(of as ApiResult)) return of as ApiResult;
      const t = this.#track(viewer, track, of as StructureOf | undefined);
      if (ApiResult.isFailure(t)) return t;
      const { ref, label } = t as SequenceTrack;
      if (!(await viewer.setTrackHidden(ref, true, of as StructureOf | undefined)))
        return ApiResult.already(`The ${label} track is already hidden.`);
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
      structure: StructureParam('A structure by label or state ref; empty: the one the sequence shows'),
    },
    run: async ({ track, structure }, ctx) => {
      const viewer = await this.#viewer();
      const named = this.#of(ctx, structure);
      if (ApiResult.isFailure(named as ApiResult)) return named as ApiResult;
      const of = named as StructureOf | undefined;
      if (track === undefined || track === '') {
        const hidden = (of ? viewer.tracksOf(of) : viewer.tracks.value).filter((t) => t.hidden);
        if (!hidden.length) return ApiResult.already('No sequence track is hidden.');
        await viewer.setTracksHidden(
          hidden.map((t) => t.ref),
          false,
          of,
        );
        const n = hidden.length;
        return ApiResult.ok(`${n === 1 ? 'The hidden track is' : `The ${n} hidden tracks are`} shown.`, {
          sequence: true,
        });
      }
      const t = this.#track(viewer, track, of);
      if (ApiResult.isFailure(t)) return t;
      const { ref, label } = t as SequenceTrack;
      if (!(await viewer.setTrackHidden(ref, false, of)))
        return ApiResult.already(`The ${label} track is already shown.`);
      return ApiResult.ok(`The ${label} track is shown.`, { sequence: true });
    },
  });

  readonly addTrack = this.operation({
    id: 'sequence.track.add',
    title: 'Add a sequence track',
    surface: {
      kind: 'agent-only',
      reason: "The Sequence node's Add track in the structure tree and the sequence.track.add command run it.",
    },
    // The edits are in the sequence's saved state (`savedOf`, version 4): a step, which Undo takes back (sequence.md §7).
    params: {
      track: PD.Text('', {
        description: 'A track to add, by its provider id or name (sequence.get lists what can be added)',
      }),
      structure: StructureParam('A structure by label or state ref; empty: the one the sequence shows'),
    },
    run: async ({ track, structure }, ctx) => {
      const viewer = await this.#viewer();
      const named = this.#of(ctx, structure);
      if (ApiResult.isFailure(named as ApiResult)) return named as ApiResult;
      const of = this.#target(named as StructureOf | undefined);
      if (ApiResult.isFailure(of as ApiResult)) return of as ApiResult;
      const addable = viewer.addableTracks(of as StructureOf);
      const t = track.trim().toLowerCase();
      const p = addable.find((a) => a.id === track || a.track === track || a.label.toLowerCase() === t);
      if (!p) {
        const there = viewer.tracksOf(of as StructureOf).find((x) => x.ref === track || x.label.toLowerCase() === t);
        if (there) return ApiResult.already(`The ${there.label} track is already there.`);
        const listed = addable.map((a) => `${a.label} (${a.id})`).join(', ');
        return ApiResult.fail(
          'not-found',
          `No sequence track ${track} to add.`,
          listed ? `What can be added: ${listed}.` : 'Every track this structure has is there.',
        );
      }
      if (!(await viewer.setTrackIn(p.id, true, named as StructureOf | undefined)))
        return ApiResult.already(`The ${p.label} track is already there.`);
      return ApiResult.ok(`The ${p.label} track is added.`, { sequence: true });
    },
  });

  readonly removeTrack = this.operation({
    id: 'sequence.track.remove',
    title: 'Remove a sequence track',
    surface: {
      kind: 'agent-only',
      reason: "A track's Remove in the structure tree and the sequence.track.remove command run it.",
    },
    params: {
      track: PD.Text('', { description: "A track of the sequence's band, by ref or name (sequence.get lists them)" }),
      structure: StructureParam('A structure by label or state ref; empty: the one the sequence shows'),
    },
    run: async ({ track, structure }, ctx) => {
      const viewer = await this.#viewer();
      const named = this.#of(ctx, structure);
      if (ApiResult.isFailure(named as ApiResult)) return named as ApiResult;
      const of = named as StructureOf | undefined;
      const t = this.#track(viewer, track, of);
      if (ApiResult.isFailure(t)) return t;
      const { ref, label, removable } = t as SequenceTrack;
      if (!removable)
        return ApiResult.fail('invalid', `The ${label} track is the document's own.`, 'Hide it: sequence.track.hide.');
      if (!(await viewer.setTrackIn(ref, false, of))) return ApiResult.already(`The ${label} track is not there.`);
      return ApiResult.ok(`The ${label} track is removed.`, { sequence: true });
    },
  });

  readonly chains = this.operation({
    id: 'sequence.chains',
    title: 'Choose the chains the sequence shows',
    surface: {
      kind: 'agent-only',
      reason: "The Sequence node's Chains in the structure tree and the sequence.chains.* commands set it.",
    },
    // A structure's own (sequence.md §3, §7.1): a step, which Undo takes back. Nothing given: every chain.
    params: {
      structure: StructureParam('A structure by label or state ref; empty: the one the sequence shows'),
      chains: OptionalParam(
        PD.Text('', {
          description:
            'The chains to show, comma-separated, as sequence.get lists them (label_asym_id: the copies of an assembly are their chain once) or by author id; "none": no chain (to add them one by one); absent: every chain',
        }),
      ),
      polymer: OptionalParam(PD.Boolean(false, { description: 'Only the polymer chains (no ligands)' })),
      each: OptionalParam(
        PD.Boolean(false, {
          description: 'Chains of one molecule each a block of their own (merged by default: one block ×N)',
        }),
      ),
    },
    run: async ({ structure, chains, polymer, each }, ctx) => {
      const viewer = await this.#viewer();
      const named = this.#of(ctx, structure);
      if (ApiResult.isFailure(named as ApiResult)) return named as ApiResult;
      const of = this.#target(named as StructureOf | undefined);
      if (ApiResult.isFailure(of as ApiResult)) return of as ApiResult;
      const all = sequenceChains((of as StructureOf).structure);
      let kept: string[] | undefined;
      if (chains !== undefined && chains.trim().toLowerCase() === 'none') kept = [];
      else if (chains !== undefined && chains.trim() !== '') {
        // By name (`label_asym_id`), else by the author's id (which may name several): indexed once.
        const byName = new Map(all.map((c) => [c.name, c]));
        const byAuth = new Map<string, string[]>();
        for (const c of all) {
          const list = byAuth.get(c.auth);
          if (list) list.push(c.name);
          else byAuth.set(c.auth, [c.name]);
        }
        const want = new Set<string>();
        for (const name of String(chains)
          .split(',')
          .map((c: string) => c.trim())
          .filter(Boolean)) {
          const hit = byName.has(name) ? [name] : byAuth.get(name);
          if (!hit) {
            return ApiResult.fail(
              'not-found',
              `No chain ${name}.`,
              `The chains: ${all.map((c) => c.label).join(', ')}.`,
            );
          }
          for (const n of hit) want.add(n);
        }
        // In the structure's order.
        kept = all.filter((c) => want.has(c.name)).map((c) => c.name);
      }
      const filter: ChainFilter = {
        ...(polymer ? { polymer: true } : {}),
        ...(each ? { each: true } : {}),
        ...(kept ? { chains: kept } : {}),
      };
      const said = kept
        ? kept.length
          ? `chains ${kept.join(', ')}`
          : 'no chain'
        : polymer
          ? 'the polymer chains'
          : 'every chain';
      if (!(await viewer.setChains(filter, named as StructureOf | undefined)))
        return ApiResult.already(`The sequence already shows ${said}.`);
      return ApiResult.ok(`The sequence shows ${said}${each ? ', every copy apart' : ''}.`, { sequence: true });
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
      const of = shown ? { ref: shown.ref, structure: shown.structure } : undefined;
      const filter = shown ? viewer.chainsOf(shown.ref) : undefined;
      const value = {
        structure: shown?.title,
        ...(shown ? { ref: shown.ref } : {}),
        followed: this.source.followed.value ? 'chosen' : 'automatic',
        style,
        colors: viewer.theme.value.letters ?? 'none',
        rail: viewer.rail.value,
        tracks: viewer.tracks.value.map((t) => ({
          ref: t.ref,
          name: t.label,
          hidden: t.hidden,
          removable: t.removable,
        })),
        addable: of ? viewer.addableTracks(of).map((a) => ({ id: a.id, name: a.label })) : [],
        chains: of
          ? sequenceChains(of.structure).map((c) => ({
              name: c.name,
              label: c.label,
              entity: c.entityLabel,
              kind: c.kind,
            }))
          : [],
        ...(filter ? { filter } : {}),
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
