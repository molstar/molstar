/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import {
  bind,
  Button,
  Choices,
  type DockId,
  Icon,
  Legend,
  observableState,
  onCleanup,
  type ShellModel,
  Show,
} from '@molstar/prism-ui';
import type { StructureView3D } from '@molstar/seq-molstar-poc/viewer';
import type { HoverLabel } from '@molstar/seq-viewer';
import { mountSection } from '@molstar/seq-viewer/canvas';
import { combineLatest, map, of, switchMap } from 'rxjs';
import { SectionRef } from './document.js';
import { linkSequence } from './link.js';
import type { SequenceStyle, SequenceViewer } from './model.js';
import { classes } from './_generated/sequence-control.module.css.js';

// The sequence control (prism-sequence.md §4): Seq*'s section surface on the canvas, the switch between compact text
// and tracks, the legend, the link with the 3D view while it is shown, and the Seq* hover label for the shell's one
// cursor-following label. Seq* draws the letters, hover, selection and focus itself and handles the pointer.

export interface SequenceControlProps {
  viewer: SequenceViewer;
  /**
   * The 3D view of a structure, for the link (hover and selection both ways, focus); none: the view stands alone.
   * Over the real plugin it is `pluginView3D(plugin, structure)` (sources/).
   */
  view3d?: (structure: Structure) => StructureView3D | undefined;
  /** The shell and its dock, for the height rules (§4): tracks follow their content, compact keeps the dock's size. */
  shell?: ShellModel;
  dock?: DockId;
  /** The Seq* hover label as lines for the shell's hover label (`[]`: none), once a frame and only when it changes. */
  onHover?: (lines: readonly string[]) => void;
  onCollapse?: () => void;
}

const Styles = [
  { id: 'compact', label: 'Text', title: 'Compact text: the letters in groups of 10' },
  { id: 'tracks', label: 'Tracks', title: 'Tracks: one zoomable track' },
] as const;

const NoLines: readonly string[] = [];

/** A Seq* hover label as the shell's lines: the title, then `key: value` per line. */
export function labelLines(label: HoverLabel | null): readonly string[] {
  if (!label) return NoLines;
  const lines = [label.title];
  for (const l of label.lines) lines.push(l.value ? `${l.key}: ${l.value}` : l.key);
  return lines;
}

const sameLines = (a: readonly string[], b: readonly string[]) =>
  a.length === b.length && a.every((x, i) => x === b[i]);

export function SequenceControl(props: SequenceControlProps) {
  const { viewer } = props;
  const seq = viewer.seq;

  const surface = (<div class={classes.surface} />) as HTMLElement;
  const legendOpen = observableState(false);
  const bar = (
    <div class={classes.bar}>
      <div class={classes.switch}>
        <Choices choices={Styles} active={viewer.style} onPick={(id) => viewer.setStyle(id as SequenceStyle)} />
      </div>
      <span class={classes.spacer} />
      <span data-legend-button="">
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="palette" />}
          title="What the colors mean"
          pressed={legendOpen.asReadonly}
          onClick={() => legendOpen.set(!legendOpen.value)}
        />
      </span>
      {props.onCollapse ? (
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="chevron-up" />}
          title="Hide the sequence"
          onClick={() => props.onCollapse?.()}
        />
      ) : null}
    </div>
  ) as HTMLElement;

  const section = mountSection(seq, SectionRef, surface);
  onCleanup(() => section.dispose());

  // The fonts: the host's tokens, read once at the start of a frame (the control is in the document by then).
  const fontFrame = requestAnimationFrame(() => {
    const css = getComputedStyle(surface);
    const mono = css.getPropertyValue('--msp-font-mono').trim();
    const ui = css.getPropertyValue('--msp-font-sans').trim();
    if (mono && ui) viewer.setFonts({ mono, ui });
  });
  onCleanup(() => cancelAnimationFrame(fontFrame));

  // The link while this control is rendered (link.ts): subscribed on mount, dropped on cleanup.
  if (props.view3d) onCleanup(linkSequence(viewer, props.view3d));

  // The hover label: the latest label of a frame, handed on when its lines differ from the last handed on.
  if (props.onHover) {
    const handOn = props.onHover;
    let pending: HoverLabel | null = null;
    let last = NoLines;
    let frame = 0;
    const flush = () => {
      frame = 0;
      const lines = labelLines(pending);
      if (sameLines(lines, last)) return;
      last = lines;
      handOn(lines);
    };
    const sub = seq.labels.current.subscribe((label) => {
      pending = label;
      if (!frame) frame = requestAnimationFrame(flush);
    });
    // Leaving the surface ends the label at once (Seq* clears its hover on the same event).
    const leave = () => {
      pending = null;
      if (!frame) frame = requestAnimationFrame(flush);
    };
    surface.addEventListener('pointerleave', leave);
    onCleanup(() => {
      sub.unsubscribe();
      surface.removeEventListener('pointerleave', leave);
      if (frame) cancelAnimationFrame(frame);
      if (last.length > 0) handOn(NoLines);
    });
  }

  // The height rules: tracks follow Seq*'s content height (the bar's, from its ResizeObserver entry, plus the section's);
  // compact keeps the size the dock had, the text scrolling inside it. Written only when the size changes.
  if (props.shell) {
    const { shell } = props;
    const dockId = props.dock ?? 'top';
    const compactSize = shell.dock(dockId).state.value.size;
    const barHeight = observableState(0);
    const ro = new ResizeObserver((entries) => barHeight.set(entries[entries.length - 1].contentRect.height));
    ro.observe(bar);
    onCleanup(() => {
      ro.disconnect();
      shell.setSize(dockId, compactSize);
    });
    const content = combineLatest([viewer.style, viewer.structure]).pipe(
      switchMap(([style]) => (style === 'tracks' ? (seq.layout.contentHeight(SectionRef) ?? of(0)) : of(0))),
    );
    bind(combineLatest([viewer.style, content, barHeight.asReadonly]), ([style, height, bar]) => {
      if (style === 'tracks' && height > 0) shell.setSize(dockId, height + bar);
      else shell.setSize(dockId, compactSize);
    });
  }

  // The legend: opened by its button, closed by Escape or a press elsewhere (listeners only while it is open).
  const outside = (e: PointerEvent) => {
    if (!(e.target as Element | null)?.closest?.(`.${classes.popup}, [data-legend-button]`)) legendOpen.set(false);
  };
  const escape = (e: KeyboardEvent) => {
    if (e.key === 'Escape') legendOpen.set(false);
  };
  const watching = legendOpen.asReadonly.subscribe((open) => {
    if (open) {
      document.addEventListener('pointerdown', outside, true);
      document.addEventListener('keydown', escape);
    } else {
      document.removeEventListener('pointerdown', outside, true);
      document.removeEventListener('keydown', escape);
    }
  });
  onCleanup(() => {
    watching.unsubscribe();
    document.removeEventListener('pointerdown', outside, true);
    document.removeEventListener('keydown', escape);
  });

  return (
    <div class={classes.root}>
      {bar}
      {surface}
      <Show when={legendOpen}>
        {() => (
          <div class={classes.popup} role="dialog" aria-label="Legend">
            <Legend
              title={viewer.legend.pipe(map((l) => l.title))}
              items={viewer.legend.pipe(map((l) => l.items))}
              ramps={viewer.legend.pipe(map((l) => l.ramps))}
            />
          </div>
        )}
      </Show>
    </div>
  );
}
