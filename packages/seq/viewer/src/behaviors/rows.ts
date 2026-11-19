/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `rows` behavior (_spec/README.md "Behaviors", Rows): an alignment's rows band as commands. Folds are session
// state (the band's `folds` in the view model: no undo step); a toggle is a hand's move, the bulk folds and a reveal are
// jumps (back restores the folds and the rows' window together); `clades.keep` writes the folds shown into the layer's
// `collapsed` (one patch), and a patch or an undo that changes `collapsed` resets the session's folds to it. Pins,
// hidden rows and the order are spec edits. The rows' scroll and zoom are the band's (`bandTop` in display rows,
// `rowScale` px a row); a fold keeps the row at the band's top there (the rows under the eye do not jump), and the
// view's history keeps a band by its top row's key (`BandRows`).

import type { Op } from '@molstar/core/repr/patch';
import type { Node, Value } from '@molstar/core/repr/tree';
import {
  allClades,
  cladesAtDepth,
  cladesHolding,
  cladesWithin,
  type RowsLay,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { EditOutcome, SpecModel, ViewModel } from '@molstar/seqx-viewer/models';
import type { EditOrigin } from '@molstar/seqx-viewer/ports';
import { CladeActions } from '@molstar/seqx-viewer/bindings';

/** What a rows command is given: the rows layer's track (the first rows band when absent) and what each needs. */
export interface RowsContext {
  readonly track?: string;
  /** `clades.toggle`: the clade by name. */
  readonly clade?: string;
  /** `clades.reveal`, `rows.pin`, `rows.hide`: data rows. */
  readonly rows?: readonly number[];
  /** `clades.collapseTo`: a depth, or a height by branch length. */
  readonly depth?: number;
  readonly height?: number;
  /** `rows.order`: `data`, `tree`, or a field of the rows table. */
  readonly order?: Value;
  readonly origin?: EditOrigin;
}

export const RowsCommandIds = {
  toggle: 'clades.toggle',
  collapseAll: 'clades.collapseAll',
  expandAll: 'clades.expandAll',
  collapseTo: 'clades.collapseTo',
  reveal: 'clades.reveal',
  keep: 'clades.keep',
  pin: 'rows.pin',
  hide: 'rows.hide',
  show: 'rows.show',
  order: 'rows.order',
  zoomIn: 'rows.zoom-in',
  zoomOut: 'rows.zoom-out',
} as const;

/** The actions a binding may name when the viewer has the behavior (a tree node's click). */
export const RowsActions: readonly string[] = CladeActions;

/**
 * The rows' zoom steps a key takes (×2) and its bounds, px a row (fixed: 1/64 shows 10,000 rows in about 160 px, a
 * band's whole family; 40 is twice a letter row, past which a row reads no better).
 */
const ZoomStep = 2;
const MinRowScale = 1 / 64;
const MaxRowScale = 40;

export interface RowsBehavior {
  readonly commands: readonly CommandDefinition<RowsContext>[];
  /** The wheel's vertical delta over a rows band scrolls its rows (px); false when there is none to scroll. */
  wheel(track: string, dy: number): boolean;
  /** The rows layer's last layout of a track (what the commands read). */
  lay(track?: string): RowsLay | undefined;
  dispose(): void;
}

/** A layer's track: its address less its last step. */
export function trackOfLayer(layer: string): string {
  return layer.slice(0, layer.lastIndexOf('/'));
}

export function rowsBehavior(spec: SpecModel, view: ViewModel, layout: () => ViewLayout | undefined): RowsBehavior {
  let disposed = false;
  /** The rows layers' laid out state by track, from the last frame. */
  const lays = (): Map<string, RowsLay> => {
    const out = new Map<string, RowsLay>();
    for (const s of layout()?.sections ?? [])
      for (const band of s.lays)
        for (const l of band ?? []) if (l.type.name === 'rows') out.set(trackOfLayer(l.key), l.lay as RowsLay);
    return out;
  };
  const lay = (track?: string): RowsLay | undefined => {
    const all = lays();
    return track !== undefined ? all.get(track) : all.values().next().value;
  };
  const trackOf = (ctx: RowsContext): string | undefined => {
    if (ctx.track !== undefined) return ctx.track;
    for (const t of lays().keys()) return t;
    return undefined;
  };
  const needs =
    (tree: boolean) =>
    (ctx: RowsContext): true | string => {
      if (disposed) return 'The viewer is retired';
      const l = lay(trackOf(ctx));
      if (!l) return 'No rows there';
      if (tree && !(l.clades && l.foldsApply)) return 'The rows are not in tree order';
      return true;
    };
  /** New folds, the row at the band's top kept there (a fold above it moves it up, never off). */
  const setFolds = (track: string, names: readonly string[]) => {
    const l = lay(track);
    const at = l?.keyAt(l.top);
    const top = l && at ? l.topOf(at.key, at.offset, names) : undefined;
    view.setBand(track, top === undefined ? { folds: names } : { folds: names, bandTop: top });
  };
  const namesOf = (l: RowsLay, ids: readonly number[]) => ids.map((c) => l.clades!.name[c]);
  const folded = (track: string, l: RowsLay): readonly string[] => view.band(track)?.folds ?? l.folds;
  /** A bulk fold: a jump (back restores the folds and the rows' window together). */
  const jumpFolds = (track: string, names: readonly string[], done: string): CommandResult => {
    const o = view.jump(() => {
      const before = view.band(track);
      setFolds(track, names);
      return view.band(track) !== before;
    });
    return { message: o === 'none' || o === 'same' ? 'Already so' : done };
  };
  const said =
    (done: string) =>
    (o: EditOutcome): CommandResult =>
      o.ok
        ? { message: o.changed ? done : 'Already so', undoable: o.changed }
        : { message: o.issues.find((i) => i.severity === 'error')?.message ?? 'Refused' };
  /** One patch on the rows layer's params. */
  const edit = (label: string, l: RowsLay, ops: Op[], ctx: RowsContext, done: string) => {
    const o = spec.patch(
      label,
      (_root, revision) => ({ base: revision, ops }),
      ctx.origin ? { origin: ctx.origin } : {},
    );
    return o instanceof Promise ? o.then(said(done)) : said(done)(o);
  };
  const at = (l: RowsLay) => l.key.split('/').filter(Boolean);
  const labelsOf = (l: RowsLay, rows: readonly number[]) => rows.map((r) => l.label(r));
  /** A list param of the rows layer as the spec has it now (a command acts on the current spec, not the last frame). */
  const listParam = (l: RowsLay, key: string): string[] => {
    let n: Node | undefined = spec.root.value;
    for (const id of at(l)) n = n?.children?.find((c) => c.id === id);
    const v = n?.params?.[key] ?? l.params?.[key];
    return Array.isArray(v) ? (v.filter((x) => typeof x === 'string') as string[]) : [];
  };
  /** The rows' zoom about the band's middle row: `rowScale` by a factor, `bandTop` keeping that row where it is. */
  const zoom = (track: string, l: RowsLay, factor: number): CommandResult => {
    const scale = Math.max(MinRowScale, Math.min(MaxRowScale, l.rowScale * factor));
    if (scale === l.rowScale) return { message: 'Already so' };
    const area = Math.max(0, l.box.height - l.pinnedHeight);
    // the middle row in display rows below the pinned band, at the band's middle before and after
    const mid = l.top + area / 2 / l.rowScale;
    const top = Math.max(0, Math.min(mid - area / 2 / scale, l.maxTop(scale, l.display.count)));
    view.setBand(track, { rowScale: scale, bandTop: top });
    return { message: factor > 1 ? 'Zoomed in' : 'Zoomed out' };
  };

  const commands: CommandDefinition<RowsContext>[] = [
    {
      id: RowsCommandIds.toggle,
      title: 'Fold or unfold a clade',
      category: 'Rows',
      when: (ctx) => {
        const ok = needs(true)(ctx);
        if (ok !== true) return ok;
        return ctx.clade !== undefined && lay(trackOf(ctx))!.clades!.name.includes(ctx.clade) ? true : 'No clade there';
      },
      run: (_, ctx): CommandResult => {
        const track = trackOf(ctx)!;
        const l = lay(track)!;
        const now = folded(track, l);
        const next = now.includes(ctx.clade!) ? now.filter((c) => c !== ctx.clade) : [...now, ctx.clade!];
        setFolds(track, next);
        return { message: next.length > now.length ? 'Folded' : 'Unfolded' };
      },
    },
    {
      id: RowsCommandIds.collapseAll,
      title: 'Collapse every clade',
      category: 'Rows',
      when: needs(true),
      run: (_, ctx) => {
        const track = trackOf(ctx)!;
        const l = lay(track)!;
        return jumpFolds(track, namesOf(l, allClades(l.clades!)), 'Collapsed');
      },
    },
    {
      id: RowsCommandIds.expandAll,
      title: 'Expand every clade',
      category: 'Rows',
      when: needs(true),
      run: (_, ctx) => jumpFolds(trackOf(ctx)!, [], 'Expanded'),
    },
    {
      id: RowsCommandIds.collapseTo,
      title: 'Collapse to a depth or a distance',
      category: 'Rows',
      when: (ctx) => {
        const ok = needs(true)(ctx);
        if (ok !== true) return ok;
        return ctx.depth !== undefined || ctx.height !== undefined ? true : 'No depth or height given';
      },
      run: (_, ctx) => {
        const track = trackOf(ctx)!;
        const l = lay(track)!;
        const ids =
          ctx.depth !== undefined ? cladesAtDepth(l.clades!, ctx.depth) : cladesWithin(l.clades!, ctx.height!);
        return jumpFolds(track, namesOf(l, ids), 'Collapsed');
      },
    },
    {
      id: RowsCommandIds.reveal,
      title: 'Open the clades holding rows',
      category: 'Rows',
      when: (ctx) => {
        const ok = needs(true)(ctx);
        if (ok !== true) return ok;
        return ctx.rows?.length ? true : 'No rows given';
      },
      run: (_, ctx) => {
        const track = trackOf(ctx)!;
        const l = lay(track)!;
        const open = new Set<string>();
        for (const r of ctx.rows!) for (const c of cladesHolding(l.clades!, r)) open.add(l.clades!.name[c]);
        const now = folded(track, l);
        const next = now.filter((c) => !open.has(c));
        if (next.length === now.length) return { message: 'Already so' };
        return jumpFolds(track, next, 'Revealed');
      },
    },
    {
      id: RowsCommandIds.keep,
      title: 'Keep the folds as the opening state',
      category: 'Rows',
      when: needs(true),
      run: (_, ctx) => {
        const track = trackOf(ctx)!;
        const l = lay(track)!;
        const folds = [...folded(track, l)];
        return edit('Keep the folds', l, [{ op: 'set', at: at(l), key: 'collapsed', value: folds }], ctx, 'Kept');
      },
    },
    {
      id: RowsCommandIds.pin,
      title: 'Pin or unpin rows',
      category: 'Rows',
      when: (ctx) => {
        const ok = needs(false)(ctx);
        if (ok !== true) return ok;
        return ctx.rows?.length ? true : 'No rows given';
      },
      run: (_, ctx) => {
        const l = lay(trackOf(ctx))!;
        const pinned = listParam(l, 'pinned');
        const labels = labelsOf(l, ctx.rows!);
        // a pinned row is unpinned (a toggle: its inverse is itself)
        const all = labels.every((x) => pinned.includes(x));
        const next = all
          ? pinned.filter((x) => !labels.includes(x))
          : [...pinned, ...labels.filter((x) => !pinned.includes(x))];
        return edit(
          'Pin rows',
          l,
          [{ op: 'set', at: at(l), key: 'pinned', value: next }],
          ctx,
          all ? 'Unpinned' : 'Pinned',
        );
      },
    },
    {
      id: RowsCommandIds.hide,
      title: 'Hide rows',
      category: 'Rows',
      when: (ctx) => {
        const ok = needs(false)(ctx);
        if (ok !== true) return ok;
        return ctx.rows?.length ? true : 'No rows given';
      },
      run: (_, ctx) => {
        const l = lay(trackOf(ctx))!;
        const hidden = listParam(l, 'hidden_rows');
        const next = [...hidden, ...labelsOf(l, ctx.rows!).filter((x) => !hidden.includes(x))];
        return edit('Hide rows', l, [{ op: 'set', at: at(l), key: 'hidden_rows', value: next }], ctx, 'Hidden');
      },
    },
    {
      id: RowsCommandIds.show,
      title: 'Show the hidden rows',
      category: 'Rows',
      when: (ctx) => {
        const ok = needs(false)(ctx);
        if (ok !== true) return ok;
        return listParam(lay(trackOf(ctx))!, 'hidden_rows').length ? true : 'No hidden rows';
      },
      run: (_, ctx) => {
        const l = lay(trackOf(ctx))!;
        return edit('Show the hidden rows', l, [{ op: 'unset', at: at(l), key: 'hidden_rows' }], ctx, 'Shown');
      },
    },
    {
      id: RowsCommandIds.order,
      title: 'Order the rows',
      category: 'Rows',
      when: (ctx) => {
        const ok = needs(false)(ctx);
        if (ok !== true) return ok;
        return ctx.order !== undefined ? true : 'No order given';
      },
      run: (_, ctx) => {
        const l = lay(trackOf(ctx))!;
        return edit('Order the rows', l, [{ op: 'set', at: at(l), key: 'order', value: ctx.order! }], ctx, 'Ordered');
      },
    },
    {
      id: RowsCommandIds.zoomIn,
      title: 'Zoom the rows in',
      category: 'Rows',
      when: needs(false),
      run: (_, ctx) => zoom(trackOf(ctx)!, lay(trackOf(ctx))!, ZoomStep),
    },
    {
      id: RowsCommandIds.zoomOut,
      title: 'Zoom the rows out',
      category: 'Rows',
      when: needs(false),
      run: (_, ctx) => zoom(trackOf(ctx)!, lay(trackOf(ctx))!, 1 / ZoomStep),
    },
  ];

  // a patch or an undo that changes a rows layer's `collapsed` resets the session's folds to it
  const seen = new Map<string, Value | undefined>();
  const scan = (root: Node, tell: boolean) => {
    const walk = (n: Node, address: string) => {
      for (const c of n.children ?? []) {
        if (c.id === undefined) continue;
        const a = `${address}/${c.id}`;
        if (c.kind === 'layer') {
          if (c.params?.type !== 'rows') continue;
          const now = c.params?.collapsed;
          if (tell && seen.has(a) && seen.get(a) !== now) resetFolds(trackOfLayer(a));
          seen.set(a, now);
        } else walk(c, a);
      }
    };
    walk(root, '');
  };
  const resetFolds = (track: string) => {
    if (view.band(track)?.folds !== undefined) view.setBand(track, { folds: undefined });
  };
  scan(spec.model.value.spec, false);
  const stop = spec.model.listen((model) => scan(model.spec, true));

  // the history keeps a band by its top row's key; restored under the step's folds and zoom, clamped
  view.setBandRows({
    keyAt(track, band) {
      const l = lay(track);
      if (!l) return undefined;
      const rs = band.rowScale ?? l.rowScale;
      return l.keyAt(Math.min(band.bandTop, l.maxTop(rs, l.display.count)));
    },
    topOf(track, band, key, offset) {
      // a band with no folds of its own shows its layer's `collapsed`
      const l = lay(track);
      return l?.topOf(key, offset, band.folds ?? l.collapsed, band.rowScale);
    },
    axisOf(system) {
      for (const l of lays().values()) if (l.rowsSystem === system) return l.rowsAxis();
      return undefined;
    },
  });

  return {
    commands,
    wheel(track, dy) {
      const l = lay(track);
      if (!l || disposed) return false;
      // the wheel's px as display rows at the band's zoom
      const top = Math.max(0, Math.min(l.maxTop(l.rowScale, l.display.count), l.top + dy / l.rowScale));
      if (top === l.top) return false;
      view.setBand(track, { bandTop: top });
      return true;
    },
    lay,
    dispose() {
      disposed = true;
      view.setBandRows(undefined);
      stop();
    },
  };
}
