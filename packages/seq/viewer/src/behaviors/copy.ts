/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `selection.copyLetters` (_spec/README.md "Behaviors", Copy): the selection's letters as a section shows them, in its
// reading order (its display order, a reversed section's included), complemented as drawn, only the positions shown (a
// collapsed stretch's are not), a run of consecutive display units a line; a section of several letters rows gives
// them as FASTA lines (`>label`, then the row's letters). With `stored`, the stored letters of the selection's span on
// the section's system, in position order, as stored. The section: the cursor's, else the first the selection projects
// onto with letters. Letters not held yet are loaded first; the text is written through the behavior's `write` (the
// platform's clipboard by default), the command's message says how many letters ("Copied 393 letters").

import { type ResolvedLayer, type ResolvedSection, projectOntoView } from '@molstar/seqx-core/view';
import type { Want } from '@molstar/seqx-core/data';
import { LettersType, shownLetters } from '@molstar/seqx-core/render';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { SelectContext } from './select.js';
import { parseChord } from '@molstar/seqx-viewer/bindings';
import type { ViewerConfig } from '@molstar/seqx-viewer/config';
import {
  axisTotal,
  elementAt,
  nearestOffsetOf,
  type SelectionModel,
  type SpecModel,
  type StoreModel,
  type ViewModel,
} from '@molstar/seqx-viewer/models';

/** Writes copied text where the platform keeps it (the clipboard): what `selection.copyLetters` hands its letters to. */
export type CopyWrite = (text: string) => void | Promise<void>;

/** The platform's clipboard, when there is one (a page); headless, nothing. */
export const clipboardWrite: CopyWrite = (text) => {
  const c = (globalThis as { navigator?: { clipboard?: { writeText?(t: string): Promise<void> } } }).navigator
    ?.clipboard;
  return c?.writeText ? c.writeText(text).catch(() => undefined) : undefined;
};

/** How many times letters not held are loaded, and how long each waits, before the copy gives what it has. */
const LoadRounds = 3;
const LoadMs = 10_000;

export function copyLettersCommand(
  deps: {
    readonly spec: SpecModel;
    readonly view: ViewModel;
    readonly store: StoreModel;
    readonly selection: SelectionModel;
    readonly config: ViewerConfig;
  },
  write: CopyWrite = clipboardWrite,
): CommandDefinition<SelectContext> {
  const { spec, view, store, selection } = deps;
  const shortcut = parseChord(deps.config.view.keys['selection.copy-letters']?.[0] ?? '');

  const lettersLayers = (s: ResolvedSection) => {
    const out: { label: string; layer: ResolvedLayer }[] = [];
    for (const t of s.tracks)
      for (const l of t.layers) if (l.placed && l.type === 'letters') out.push({ label: t.label, layer: l });
    return out;
  };

  /** The positions to copy on a section: runs in reading order, or the span in position order (`stored`). */
  const runsOf = (s: ResolvedSection, ranges: ArrayLike<number>, stored: boolean): number[][] => {
    if (stored) {
      const span: number[] = [];
      for (let p = ranges[0]; p < ranges[ranges.length - 1]; p++) span.push(p);
      return [span];
    }
    const total = axisTotal(s.axis);
    if (total === undefined) return [];
    // the display units the selection's hull spans, walked in display order: a unit is copied when its position is
    // selected; a unit not copied ends a run
    let lo = Infinity;
    let hi = -Infinity;
    for (let i = 0; i < ranges.length; i++) {
      const o = nearestOffsetOf(s.axis, ranges[i]);
      if (o === undefined) continue;
      lo = Math.min(lo, o);
      hi = Math.max(hi, o);
    }
    const runs: number[][] = [];
    let run: number[] | undefined;
    for (let u = Math.max(0, Math.floor(lo) - 1); u <= Math.min(total - 1, Math.ceil(hi)); u++) {
      const p = elementAt(s.axis, u);
      if (p !== undefined && inRanges(ranges, p)) {
        if (!run) runs.push((run = []));
        run.push(p);
      } else run = undefined;
    }
    return runs;
  };

  /** The text, or undefined when nothing selected shows letters; letters not held are loaded first. */
  const lettersOf = async (stored: boolean): Promise<string | undefined> => {
    const shown = selection.value.value.shown;
    const resolved = view.resolved.value;
    if (!resolved || !shown.parts.length) return undefined;
    const model = spec.model.value;
    const loc = projectOntoView(model, resolved, shown, { own: true });
    const cursor = selection.cursor.value?.section;
    const order = [...resolved.sections].sort((a, b) => (a.address === cursor ? -1 : b.address === cursor ? 1 : 0));
    for (const s of order) {
      const part = loc.parts.find((p) => p.system === s.axis.system);
      const rows = lettersLayers(s);
      if (!part?.ranges.length || !rows.length) continue;
      const runs = runsOf(s, part.ranges, stored);
      if (!runs.length) continue;
      const preps = rows.map((r) => LettersType.prepare(r.layer, model));
      for (let round = 0; ; round++) {
        const wants: Want[] = [];
        const texts = preps.map((prep) =>
          runs
            .map((positions) => {
              const r = shownLetters(prep, s.axis.dir, store.store, positions, stored);
              wants.push(...r.wants);
              return r.text;
            })
            .filter((t) => t.length)
            .join('\n'),
        );
        if (!wants.length || round >= LoadRounds) {
          if (rows.length === 1) return texts[0];
          return rows.map((r, i) => `>${r.label}\n${texts[i]}`).join('\n');
        }
        const w = store.store.wanter();
        w.want(wants);
        await w.whenSettled(LoadMs);
        w.dispose();
      }
    }
    return undefined;
  };

  return {
    id: 'selection.copy-letters',
    title: 'Copy the selection’s letters',
    category: 'Selection',
    ...(shortcut ? { shortcut: { ...shortcut }, scope: 'canvas' as const } : {}),
    when: () => (selection.value.value.shown.parts.length ? true : 'Nothing selected'),
    run: async (_, ctx): Promise<CommandResult> => {
      const text = await lettersOf(!!ctx?.stored);
      if (!text) return { message: '' };
      await write(text);
      // a count, not the text: a message is shown and announced (a chromosome's letters would be read out)
      let n = 0;
      for (const line of text.split('\n')) if (!line.startsWith('>')) n += line.length;
      return { message: `Copied ${n.toLocaleString('en-US')} letter${n === 1 ? '' : 's'}` };
    },
  };
}

/** Whether a position is in sorted `[start, end)` pairs. */
function inRanges(r: ArrayLike<number>, p: number): boolean {
  let lo = 0;
  let hi = r.length >> 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (r[2 * mid + 1] <= p) lo = mid + 1;
    else hi = mid;
  }
  return lo < r.length >> 1 && r[2 * lo] <= p;
}
