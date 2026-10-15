/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { SequenceViewer } from '../sequence/model.js';
import { LetterSchemes } from '../sequence/palette.js';
import type { SequenceSource } from '../sources/sequence.js';
import { type ApiAreaContext, ApiArea, ApiResult } from './area.js';
import { isFailure, labelOf, oneStructure } from './scope.js';

// `sequence` (operations.md §2, sequence.md §6): the one sequence viewer's style (Tracks or Text), the letters' colors,
// which structure it follows. Its style and colors are the sequence's state, saved with snapshots
// (`prism-viewer.sequence`), so changing them is a step; following is view state. Present only when the spec has a
// sequence (viewer-spec.md §2.3).

const Schemes = LetterSchemes.map(
  (s) => [s ?? 'none', s === null ? 'None' : s.replace(/^prism:/, '')] as [string, string],
);

export class SequenceApi extends ApiArea {
  readonly id = 'sequence';

  constructor(
    ctx: ApiAreaContext,
    private readonly viewer: SequenceViewer,
    private readonly source: SequenceSource,
  ) {
    super(ctx);
  }

  readonly style = this.operation({
    id: 'sequence.style',
    title: 'Show the sequence as tracks or text',
    params: {
      style: PD.Select<'tracks' | 'text'>('tracks', [
        ['tracks', 'Tracks'],
        ['text', 'Text'],
      ]),
    },
    run: async ({ style }) => {
      const want = style === 'text' ? 'compact' : 'tracks';
      if (this.viewer.style.value === want) return ApiResult.already(`The sequence is already shown as ${style}.`);
      this.viewer.setStyle(want);
      return ApiResult.ok(`The sequence is shown as ${style}.`, { sequence: true });
    },
  });

  readonly colors = this.operation({
    id: 'sequence.colors',
    title: "Color the sequence's letters",
    params: { scheme: PD.Select(Schemes[0][0], Schemes, { description: 'A letter color scheme' }) },
    run: async ({ scheme }) => {
      const letters = scheme === 'none' ? null : scheme;
      if (this.viewer.theme.value.letters === letters)
        return ApiResult.already(`The letters are already colored by ${scheme}.`);
      this.viewer.setTheme({ letters });
      return ApiResult.ok(`The letters are colored by ${scheme}.`, { sequence: true });
    },
  });

  readonly follow = this.operation({
    id: 'sequence.follow',
    title: 'Show the sequence of a structure',
    params: {
      structure: PD.Text('', { description: 'A structure by label or ref; empty: the one the person works with' }),
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

  readonly get = this.operation({
    id: 'sequence.get',
    title: 'Describe the sequence view',
    params: {},
    read: true,
    run: async () => {
      const shown = this.viewer.structure.value;
      const style = this.viewer.style.value === 'compact' ? 'text' : 'tracks';
      const value = { structure: shown?.title, style, colors: this.viewer.theme.value.letters ?? 'none' };
      return ApiResult.ok(
        shown ? `The sequence of ${shown.title}, as ${style}, colored by ${value.colors}.` : 'No sequence is shown.',
        {},
        { value },
      );
    },
  });
}
