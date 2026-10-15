/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { activity, For } from '../../dom/control-flow.js';
import { playClose } from '../../dom/motion.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel, ShellSection, ShellViewContext } from '../../shell/shell.js';
import { classes } from './_generated/section-list.module.css.js';
import { ViewHost } from './view-host.js';

/**
 * The sections of a slot as a list of collapsible blocks (the Export panel). Each opens as the user last left it; until
 * then, open, or folded with `collapsed` (the Selection tab: its sections read as a list of what can be done). The
 * focused section (a command opened the panel at it) is expanded and scrolled into view; the others fold. A focus
 * `section/inner` opens that section and gives it `inner` as its own `focus` (the Image section's Animation).
 */
export function SectionList(props: {
  shell: ShellModel;
  slot: string;
  focus?: string;
  subject?: unknown;
  /** Inside a surface that already has a background (the inspector): no backgrounds of its own, smaller headers. */
  flat?: boolean;
  /** Sections the user has not opened yet start folded. */
  collapsed?: boolean;
}) {
  const { shell, subject } = props;
  const ctx: ShellViewContext = { shell, subject };
  const all = shell.sections(props.slot);
  const shown = (list: readonly ShellSection[]) => list.filter((s) => !s.when || s.when(ctx));
  const sections = all.pipe(map(shown));
  const focus = props.focus;
  const focuses = (id: string) => !!focus && (focus === id || focus.startsWith(`${id}/`));
  const inner = (id: string) => (focus?.startsWith(`${id}/`) ? focus.slice(id.length + 1) : undefined);
  return (
    <div class={classes.list} data-flat={props.flat ? '' : undefined}>
      <For each={sections} key={(s) => s.id}>
        {(section) => (
          <SectionBlock
            ctx={inner(section.id) ? { ...ctx, focus: inner(section.id) } : ctx}
            slot={props.slot}
            section={section}
            // Opened at a section: only that one open (the rest fold away); otherwise as the user left it, or the default.
            open={
              focus
                ? focuses(section.id)
                : (shell.sectionOpen(props.slot, section.id) ?? (!props.collapsed || !!section.open))
            }
            focused={focuses(section.id)}
          />
        )}
      </For>
    </div>
  );
}

function SectionBlock(props: {
  ctx: ShellViewContext;
  slot: string;
  section: ShellSection;
  open: boolean;
  focused: boolean;
}) {
  const { ctx, section } = props;
  const open = observableState(props.open);
  // What the header says (its chevron turns at once), ahead of the body while it closes.
  const expanded = observableState(props.open);
  // Opened by the user from here on: reveal it (as Show's `reveal`); the first render stays still.
  const toggled = observableState(false);
  let body: HTMLElement | undefined;
  let closing: Animation | undefined;
  const toggle = () => {
    toggled.set(true);
    closing?.cancel();
    closing = undefined;
    expanded.set(!expanded.value);
    // Closing plays the reveal back first (the body stays rendered, so it only hides after).
    const anim = !expanded.value && open.value && body ? playClose(body, false) : undefined;
    if (anim) {
      closing = anim;
      ctx.shell.setSectionOpen(props.slot, section.id, false);
      anim.finished.then(
        () => {
          if (closing !== anim) return;
          closing = undefined;
          open.set(false);
          anim.cancel();
        },
        () => {},
      );
      return;
    }
    open.set(expanded.value);
    ctx.shell.setSectionOpen(props.slot, section.id, open.value);
  };
  const block = (
    <section class={classes.section} data-section={section.id}>
      <button type="button" class={classes.header} aria-expanded={expanded} onClick={toggle}>
        <span class={classes.title}>{section.title}</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      <div
        class={classes.body}
        ref={(el) => {
          body = el;
        }}
        hidden={open.asReadonly.pipe(map((o) => !o))}
        data-reveal={toggled.asReadonly.pipe(map((t) => (t ? '' : undefined)))}
      >
        {/* A folded section stays mounted but quiet: its views follow no model until it opens (`activity`). The `hidden`
            above is bound outside it, or it would pause with it. */}
        {activity(open.asReadonly, () => (
          <ViewHost renderer={section.render} ctx={ctx} />
        ))}
      </div>
    </section>
  ) as HTMLElement;
  if (props.focused) {
    // Asked for (a command, the selection bar's Measure): as if the user opened it.
    ctx.shell.setSectionOpen(props.slot, section.id, true);
    queueMicrotask(() => block.scrollIntoView({ block: 'nearest' }));
  }
  return block;
}
