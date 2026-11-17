/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, State } from '@molstar/core/reactive';
import { Expandable, ExpandableList } from '../components/expandable/expandable.js';
import { effect, For } from '../dom/index.js';
import type { ShellModel, ShellSection, ShellViewContext } from './shell-model.js';
import { ViewHost } from './view-frame.js';

/**
 * The sections of a slot (shell §2: extensions add sections to a view, not views), as expandable panels in order: each
 * open as the person left it this session, else open, or folded with `collapsed` (a list read as what can be done). A
 * `show` at a section (`ctx.focus`) opens it and brings it into view. A section shows while its `when` holds.
 */
export function SectionList(props: {
  shell: ShellModel;
  slot: string;
  ctx: ShellViewContext;
  /** Sections not touched yet start folded (unless their own `open`). */
  collapsed?: boolean;
  /** `flat` inside a view that is already a well (an inspector). */
  variant?: 'well' | 'flat';
}) {
  const { shell, ctx } = props;
  const view = ctx.ref;
  const sections = derive(shell.sections(props.slot), (all) => all.filter((s) => !s.when || s.when(ctx)));
  return (
    <ExpandableList>
      <For each={sections} by={(s) => s.id}>
        {(section) => <SectionBlock shell={shell} ctx={ctx} view={view} section={section} props={props} />}
      </For>
    </ExpandableList>
  );
}

function SectionBlock(p: {
  shell: ShellModel;
  ctx: ShellViewContext;
  view: string;
  section: ShellSection;
  props: { collapsed?: boolean; variant?: 'well' | 'flat' };
}) {
  const { shell, ctx, view, section } = p;
  const open = new State(shell.sectionOpen(view, section.id) ?? (!p.props.collapsed || !!section.open));
  let el: HTMLElement | undefined;
  effect(ctx.focus, (f) => {
    if (f !== section.id) return;
    open.set(true);
    queueMicrotask(() => el?.scrollIntoView({ block: 'nearest' }));
  });
  const block = (
    <Expandable
      title={section.heading ?? section.title}
      open={open}
      variant={p.props.variant}
      onToggle={(o) => shell.setSectionOpen(view, section.id, o)}
    >
      {() => <ViewHost renderer={section.render} ctx={ctx} />}
    </Expandable>
  ) as HTMLElement;
  el = block;
  block.setAttribute('data-section', section.id);
  return block;
}
