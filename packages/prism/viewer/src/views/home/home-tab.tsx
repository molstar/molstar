/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { formatChord } from '@molstar/prism-ui/commands/keys';
import { ActionFooter } from '@molstar/prism-ui';
import { Button } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ParamsForm } from '@molstar/prism-ui';
import { ScrollArea } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import type { AnnotationOfferModel, HomeAction, HomeModel, RemoteStatesModel, SessionModel } from '../../model/home.js';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { classes } from './_generated/home-tab.module.css.js';
import { KeptRecordsList } from './kept-records.js';

export interface HomeTabProps {
  home: HomeModel;
  /** More sections after the actions (the `home.import` slot). */
  more?: () => Child;
  /** Opens the command menu (the hint at the bottom). */
  onOpenCommandMenu?: () => void;
}

/**
 * The Home tab: the state actions applicable to the root as forms (download structure, open files, download density,
 * load trajectory, …; the first one open), remote states when on, and a pointer to the command menu. Files dropped
 * anywhere on the plugin open too (the shell's drop overlay).
 */
export function HomeTab(props: HomeTabProps) {
  const { home } = props;
  return (
    <div class={classes.home}>
      {/* What scrolls; the hint below stays in view. */}
      <ScrollArea class={classes.content} contentClass={classes.contentInner}>
        {home.session ? <Session session={home.session} /> : null}
        {home.annotationOffer ? <AnnotationOfferLine model={home.annotationOffer} /> : null}
        <section class={classes.group} aria-label="Actions">
          {/* The first action (Download Structure, open), then Open Files, then the rest (folded). */}
          <For each={home.actions.pipe(map((a) => a.slice(0, 1)))} key={(a) => a.id}>
            {(action) => <ActionBlock action={action} open />}
          </For>
          <OpenFiles home={home} />
          <For each={home.actions.pipe(map((a) => a.slice(1)))} key={(a) => a.id}>
            {(action) => <ActionBlock action={action} open={false} />}
          </For>
          {/* Extensions' loaders (Zenodo): among the actions, folded, as one of them. */}
          {props.more?.()}
        </section>
        {home.remoteStates ? <RemoteStatesList model={home.remoteStates} caption /> : null}
        {home.keptRecords ? <KeptRecordsList model={home.keptRecords} /> : null}
      </ScrollArea>
      <button type="button" class={classes.hint} onClick={() => props.onOpenCommandMenu?.()}>
        <kbd>{formatChord({ key: 'k', mod: true })}</kbd>
        <span>Search every command</span>
      </button>
    </div>
  );
}

/**
 * Open files: a block that folds, like the actions (parity.md §1), holding a drop area: drop files on it, or click it
 * for the file picker; they open as files dropped anywhere do. Its hint says what it reads (the formats on hover).
 */
function OpenFiles(props: { home: HomeModel }) {
  const { home } = props;
  const open = observableState(false);
  const over = observableState(false);
  const pick = (files: readonly File[]) => {
    if (files.length) home.openFiles(files);
  };
  const input = (
    <input
      type="file"
      multiple
      hidden
      onChange={(e) => {
        const el = e.currentTarget as HTMLInputElement;
        const files = [...(el.files ?? [])];
        el.value = '';
        pick(files);
      }}
    />
  ) as HTMLInputElement;
  // Only drags that carry files (the page's drop overlay, when it is on, takes them before this area does).
  const carriesFiles = (e: DragEvent) => !!e.dataTransfer?.types.includes('Files');
  return (
    <div class={classes.action} data-open={open.asReadonly.pipe(map((o) => (o ? '' : undefined)))}>
      <button type="button" class={classes.actionHeader} aria-expanded={open} onClick={() => open.set(!open.value)}>
        <span class={classes.actionTitle}>Open Files</span>
        <span class={classes.actionDescription}>Structures, maps, trajectories and sessions from your computer</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      <Show when={open} reveal>
        {() => (
          <div class={classes.actionBody}>
            <button
              type="button"
              class={classes.drop}
              data-over={over.asReadonly.pipe(map((o) => (o ? '' : undefined)))}
              title={home.formats.length ? `Reads ${home.formats.join(', ')}` : undefined}
              onClick={() => input.click()}
              onDragEnter={(e) => {
                if (!carriesFiles(e)) return;
                e.preventDefault();
                over.set(true);
              }}
              onDragOver={(e) => {
                if (!carriesFiles(e)) return;
                e.preventDefault();
              }}
              onDragLeave={() => over.set(false)}
              onDrop={(e) => {
                e.preventDefault();
                over.set(false);
                pick([...(e.dataTransfer?.files ?? [])]);
              }}
            >
              <Icon name="folder-open" size={20} />
              <span class={classes.dropText}>
                Drop files here, or <span class={classes.dropLink}>browse</span>
              </span>
            </button>
            {input}
          </div>
        )}
      </Show>
    </div>
  );
}

/** One action: a header that folds, its form, and the apply button (Home's actions, the state tree's). */
export function ActionBlock(props: { action: HomeAction; open: boolean }) {
  const { action } = props;
  const open = observableState(props.open);
  return (
    <div class={classes.action} data-open={open.asReadonly.pipe(map((o) => (o ? '' : undefined)))}>
      <button type="button" class={classes.actionHeader} aria-expanded={open} onClick={() => open.set(!open.value)}>
        <span class={classes.actionTitle}>{action.title}</span>
        <span class={classes.actionDescription}>{action.description}</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      <Show when={open} reveal>
        {() => (
          <div class={classes.actionBody}>
            {action.form ? <ParamsForm form={action.form} /> : null}
            <ActionFooter
              icon={action.icon}
              label={action.applyLabel ?? 'Apply'}
              busy={action.busy}
              onApply={() => action.apply()}
              onReset={action.form ? () => action.form?.reset() : undefined}
            />
          </div>
        )}
      </Show>
    </div>
  );
}

/**
 * What is loaded, while something is: a line saying what, and Clear, which asks once (it removes everything loaded;
 * saved snapshots and files stay).
 */
function Session(props: { session: SessionModel }) {
  const { session } = props;
  const confirming = observableState(false);
  return (
    <Show when={session.isEmpty.pipe(map((e) => !e))}>
      {() => (
        <section class={classes.group} aria-label="Current session">
          <div class={classes.captionRow}>
            <h3 class={classes.caption}>Current session</h3>
          </div>
          <div class={classes.session}>
            <span class={classes.sessionSummary} data-tip-overflow="">
              {session.summary}
            </span>
            <Show
              when={confirming}
              fallback={() => (
                <Button size="sm" variant="ghost" icon={<Icon name="trash-2" />} onClick={() => confirming.set(true)}>
                  Clear
                </Button>
              )}
            >
              {() => (
                <span class={classes.confirm} role="group" aria-label="Clear everything loaded?">
                  <span class={classes.confirmText}>Clear everything?</span>
                  <Button size="sm" variant="ghost" onClick={() => confirming.set(false)}>
                    Cancel
                  </Button>
                  <Button
                    size="sm"
                    variant="danger"
                    onClick={() => {
                      confirming.set(false);
                      session.clear();
                    }}
                  >
                    Clear
                  </Button>
                </span>
              )}
            </Show>
          </div>
        </section>
      )}
    </Show>
  );
}

/**
 * The result line after a load (parity.md §6a): "AF-P04637 has pLDDT" with Show (`annotations.apply`, one step) and a
 * dismiss. Only for what an extension brings; it goes once used or dismissed.
 */
function AnnotationOfferLine(props: { model: AnnotationOfferModel }) {
  const { model } = props;
  return (
    <Show when={model.offer}>
      {(offer) => (
        <section class={classes.offer} aria-label="Loaded">
          <Icon name="sparkles" size={14} class={classes.offerIcon} />
          <span class={classes.offerText} data-tip-overflow="" title={offer.description}>
            {offer.label} has {offer.title}
          </span>
          <Button size="sm" variant="ghost" onClick={() => model.accept()}>
            Show
          </Button>
          <Button size="xs" variant="ghost" icon={<Icon name="x" />} title="Dismiss" onClick={() => model.dismiss()} />
        </section>
      )}
    </Show>
  );
}

/**
 * Sessions on a state server: the list (newest first; a click opens one), refresh, and the server behind a settings
 * button. Home's own group (`caption`), or the body of the `remote-states` feature's section in `home.import`, whose
 * header says the title.
 */
/** A server's host name ("webchem.ncbr.muni.cz"); the address itself when it is not a URL. */
function hostOf(url: string): string {
  try {
    return new URL(url).host;
  } catch {
    return url;
  }
}

export function RemoteStatesList(props: { model: RemoteStatesModel; caption?: boolean }) {
  const { model } = props;
  const settings = observableState(false);
  return (
    <section class={classes.group} aria-label="Remote states">
      <div class={classes.captionRow}>
        {/* Under a header that says the title (the feature's section): the server the list is from, what the buttons
            act on. */}
        {props.caption ? (
          <h3 class={classes.caption}>Remote states</h3>
        ) : (
          <span class={classes.serverHost} data-tip-overflow="">
            {model.serverUrl.pipe(map(hostOf))}
          </span>
        )}
        <span class={classes.captionActions}>
          <Button
            size="xs"
            variant="ghost"
            icon={<Icon name="settings" />}
            title="Server"
            pressed={settings.asReadonly}
            onClick={() => settings.set(!settings.value)}
          />
          <Button
            size="xs"
            variant="ghost"
            icon={
              <Show when={model.loading} fallback={() => <Icon name="refresh-cw" />}>
                {() => <Icon name="loader-circle" class={classes.spin} />}
              </Show>
            }
            title="Refresh"
            onClick={() => model.refresh()}
          />
        </span>
      </div>
      {/* The server is a setting: behind the button, in a well under the caption (always shown in the Remote states card). */}
      <Show when={settings}>
        {() => (
          <label class={classes.server}>
            <span>Server</span>
            <TextInput
              size="sm"
              aria-label="State server URL"
              value={model.serverUrl}
              onSubmit={(v) => {
                model.setServerUrl(v);
                settings.set(false);
              }}
            />
          </label>
        )}
      </Show>
      {/* A long list scrolls on its own. */}
      <ScrollArea class={classes.states} contentClass={classes.statesContent}>
        <For
          each={model.states}
          key={(r) => r.id}
          fallback={() => <p class={classes.none}>No states on this server.</p>}
        >
          {(state) => (
            <button type="button" class={classes.state} onClick={() => model.open(state.id)}>
              <span class={classes.stateLabel} data-tip-overflow="">
                {state.label}
              </span>
              <span class={classes.stateDate}>{state.date}</span>
            </button>
          )}
        </For>
      </ScrollArea>
    </section>
  );
}
