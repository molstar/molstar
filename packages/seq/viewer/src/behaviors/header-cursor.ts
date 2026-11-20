/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewLayout } from '@molstar/seqx-core/render';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { elementAt, type Cursor } from '@molstar/seqx-viewer/models';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The keyboard cursor on track headers (_spec/README.md "Behaviors", Headers: the cursor), the `headers` behavior's: it
// takes the `cursor` behavior's moves over (so a viewer without it carries none of this). `cursor.prev` at a section's
// first position goes onto its first track's header (`part: 'header'`, `track` its address; a section's line its
// section's address), ↑ ↓ step between the section's headers and on to the next section's, → goes back to the data
// where it left; Enter opens the header's menu; Space and the extends do nothing there.

/** A section's headers the cursor stands on, in display order: its line when drawn, then its tracks' label lines. */
function headersOf(viewer: SeqViewer, address: string): string[] {
  const s = viewer.frame?.layout.sections.find((x) => x.section.address === address);
  const out: string[] = [];
  if (!s) return out;
  if (s.header) out.push(address);
  for (const t of s.tracks)
    if (t.header?.hit?.target === 'header-label' && !out.includes(t.track.address)) out.push(t.track.address);
  return out;
}

/** ↑ ↓ on a header: the next header of its section, else the next section's first (last going up), else none. */
function headerStep(viewer: SeqViewer, c: Cursor, delta: -1 | 1): Cursor | undefined {
  const list = headersOf(viewer, c.section);
  const i = list.indexOf(c.track!);
  if (i >= 0 && list[i + delta] !== undefined) return { ...c, track: list[i + delta] };
  const sections = viewer.view.resolved.value?.sections ?? [];
  for (let k = sections.findIndex((x) => x.address === c.section) + delta; k >= 0 && k < sections.length; k += delta) {
    const next = headersOf(viewer, sections[k].address);
    const to = sections[k];
    const v = viewer.view.get(to.address);
    const p = next.length && v && !v.pending ? elementAt(to.axis, Math.ceil(v.start)) : undefined;
    if (p === undefined) continue;
    return {
      section: to.address,
      system: to.axis.system,
      position: p,
      part: 'header',
      track: next.at(delta > 0 ? 0 : -1),
    };
  }
  return undefined;
}

/**
 * The cursor's moves taken over, in the viewer's command list: `menu` opens the header's menu (Enter). Returns what the
 * frame calls: a layout that no longer has the cursor's track moves it to the next track's header, else the previous
 * one's, else the section's line.
 */
export function headerCursor(viewer: SeqViewer, menu: () => CommandResult): () => void {
  const selection = viewer.selection;
  const moved = (c: Cursor | undefined): CommandResult => {
    if (c) selection.setCursor({ ...c, keyboard: true });
    return { message: c ? 'Moved' : '' };
  };
  const over = (id: string, header: (c: Cursor) => CommandResult, after?: (before: Cursor | undefined) => void) => {
    const i = viewer.commands.findIndex((x) => x.id === id);
    if (i < 0) return;
    const was = viewer.commands[i] as CommandDefinition<unknown>;
    viewer.commands[i] = {
      ...was,
      when: (ctx) => (selection.cursor.value?.track !== undefined ? true : (was.when?.(ctx) ?? true)),
      run: (args, ctx, signal) => {
        const c = selection.cursor.value;
        if (c?.track !== undefined) return header(c);
        const r = was.run(args, ctx, signal);
        after?.(c);
        return r;
      },
    };
  };
  const nothing = () => ({ message: '' });
  // ← from the first position (a cell of one axis, the cursor shown) onto the section's first track's header
  over('cursor.prev', nothing, (before) => {
    const c = selection.cursor.value;
    if (!before?.keyboard || !c || before.ySystem !== undefined || c.position !== before.position) return;
    const first = headersOf(viewer, c.section).find((h) => h !== c.section);
    if (first !== undefined) selection.setCursor({ ...c, part: 'header', track: first, keyboard: true });
  });
  over('cursor.next', (c) => {
    const { part: _, track: __, ...cell } = c;
    return moved(cell);
  });
  over('cursor.up', (c) => moved(headerStep(viewer, c, -1)));
  over('cursor.down', (c) => moved(headerStep(viewer, c, 1)));
  for (const id of ['prev', 'next', 'up', 'down']) {
    over(`cursor.extend-${id}`, nothing);
    over(`cursor.add-${id}`, nothing);
  }
  over('cursor.toggle', nothing);
  over('cursor.focus', menu);

  /** The headers' order the last layout drew, by section: where the cursor goes when its track is hidden. */
  let order = new Map<string, string[]>();
  let orderFor: unknown;
  const relaid = (layout: ViewLayout | undefined) => {
    const c = selection.cursor.value;
    const s = c && viewer.view.resolved.value?.sections.find((x) => x.address === c.section);
    if (c?.track !== undefined && c.track !== c.section && s && !s.tracks.some((t) => t.address === c.track)) {
      const was = order.get(c.section) ?? [];
      const i = was.indexOf(c.track);
      const there = (a: string) => a === c.section || s.tracks.some((t) => t.address === a);
      const next = was.slice(i + 1).find(there) ?? was.slice(0, Math.max(0, i)).reverse().find(there) ?? c.section;
      selection.setCursor({ ...c, track: next });
    }
    if (!layout) return;
    orderFor = layout;
    order = new Map();
    for (const x of layout.sections) order.set(x.section.address, headersOf(viewer, x.section.address));
  };
  // a frame that kept its layout (a hover's): nothing read, nothing made (the work in a function of its own: what its
  // closures capture is allocated on entry)
  return () => {
    const layout = viewer.frame?.layout;
    if (layout !== orderFor) relaid(layout);
  };
}
