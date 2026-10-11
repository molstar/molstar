/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { distinctUntilChanged, map } from 'rxjs';
import { Show } from '../../dom/control-flow.js';
import { onCleanup } from '../../dom/scope.js';
import type { ShellModel, ShellPage } from '../../shell/shell.js';
import { classes } from './_generated/page.module.css.js';
import { ViewHost } from './view-host.js';

/** The page open over the app (ShellPage), if any. */
export function PageHost(props: { shell: ShellModel }) {
  const { shell } = props;
  // Built again only for another page: opened again at another focus, the page stays and follows `shell.page` itself.
  const page = shell.page.pipe(
    map((open) => open?.page),
    distinctUntilChanged(),
  );
  return <Show when={page}>{(p) => <Page shell={shell} page={p} />}</Show>;
}

function Page(props: { shell: ShellModel; page: ShellPage }) {
  const { shell } = props;
  const open = { page: props.page, focus: shell.page.value?.focus };
  const ctx = { shell, ...(open.focus ? { focus: open.focus } : {}) };
  const page = (
    <section class={classes.page} role="region" aria-label={open.page.title} tabIndex={-1} data-msp-page="">
      <header class={classes.header}>
        <span class={classes.title}>
          {open.page.icon ? <Icon name={open.page.icon} /> : null}
          {open.page.title}
        </span>
        <div class={classes.headerContent}>
          {open.page.header ? <ViewHost renderer={open.page.header} ctx={ctx} /> : null}
        </div>
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="x" />}
          title="Close (Esc)"
          onClick={() => shell.closePage()}
        />
      </header>
      <div class={classes.body}>
        <ViewHost renderer={open.page.render} ctx={ctx} />
      </div>
    </section>
  ) as HTMLElement;
  // Esc closes it when it has the focus (or nothing has), after anything nearer has had it (an open list, a field's
  // own Escape). With the focus in the panel beside it, Esc is the panel's.
  const before = document.activeElement as HTMLElement | null;
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented) return;
    const active = document.activeElement;
    if (active && active !== document.body && !page.contains(active)) return;
    shell.closePage();
    e.preventDefault();
  };
  document.addEventListener('keydown', onKey);
  queueMicrotask(() => {
    if (!page.contains(document.activeElement)) page.focus({ preventScroll: true });
  });
  onCleanup(() => {
    document.removeEventListener('keydown', onKey);
    if (before?.isConnected) before.focus({ preventScroll: true });
  });
  return page;
}
