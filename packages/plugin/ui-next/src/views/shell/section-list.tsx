/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { For } from '../../dom/control-flow.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel, ShellSection, ShellViewContext } from '../../shell/shell.js';
import { classes } from './_generated/section-list.module.css.js';
import { ViewHost } from './view-host.js';

/**
 * The sections of a slot as a list of collapsible blocks (the Export panel). Each opens as the user last left it; until
 * then, open, or folded with `collapsed` (the Selection tab: its sections read as a list of what can be done). The
 * focused section (a command opened the panel at it) is expanded and scrolled into view; the others fold.
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
  return (
    <div class={classes.list} data-flat={props.flat ? '' : undefined}>
      <For each={sections} key={(s) => s.id}>
        {(section) => (
          <SectionBlock
            ctx={ctx}
            slot={props.slot}
            section={section}
            // Opened at a section: only that one open (the rest fold away); otherwise as the user left it, or the default.
            open={
              props.focus ? props.focus === section.id : (shell.sectionOpen(props.slot, section.id) ?? !props.collapsed)
            }
            focused={props.focus === section.id}
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
  const toggle = () => {
    open.set(!open.value);
    ctx.shell.setSectionOpen(props.slot, section.id, open.value);
  };
  const block = (
    <section class={classes.section} data-section={section.id}>
      <button type="button" class={classes.header} aria-expanded={open} onClick={toggle}>
        <span class={classes.title}>{section.title}</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      <div class={classes.body} hidden={open.asReadonly.pipe(map((o) => !o))}>
        <ViewHost renderer={section.render} ctx={ctx} />
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
