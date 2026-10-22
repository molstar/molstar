/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ContentsSummary } from '../../model/scene.js';
import { classes } from './_generated/contents-summary.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// A structure's contents in its inspector: its counts as figures (entities, chains, residues, atoms), then its entities,
// each a row: a dot for its kind, its name, its chains and atoms; a click opens it in the tree. Many entities (a
// mesoscale structure): the first few, then Show all, as the tree's Contents.

const PageSize = 5;

/** The kinds of entity, as their dots say them (mmCIF's entity types). */
const Kinds: Readonly<Record<string, { label: string; class: string }>> = {
  polymer: { label: 'Polymer', class: classes.polymer },
  'non-polymer': { label: 'Ligand', class: classes.ligand },
  branched: { label: 'Carbohydrate', class: classes.branched },
  water: { label: 'Water', class: classes.water },
};

export function ContentsSummaryView(props: { summary: ContentsSummary }) {
  const { summary } = props;
  const chains = (list: readonly string[]) =>
    list.length > 4 ? `${list.length} chains` : `${list.length === 1 ? 'Chain' : 'Chains'} ${list.join(', ')}`;
  const all = observableState(summary.entities.length <= PageSize + 1);
  const shown = derive(all.asReadonly, (a) => (a ? summary.entities : summary.entities.slice(0, PageSize)));
  return (
    <section class={classes.contents} aria-label="Contents">
      <dl class={classes.counts}>
        {summary.counts.map((c) => (
          <div class={classes.count}>
            <dd class={classes.value}>{c.value}</dd>
            <dt class={classes.label}>{c.label}</dt>
          </div>
        ))}
      </dl>
      <div class={classes.entities}>
        <div class={classes.list} role="list" aria-label="Entities">
          <For each={shown} by={(e) => e.id}>
            {(e) => {
              const kind = Kinds[e.type];
              return (
                <button type="button" role="listitem" class={classes.entity} onClick={() => summary.reveal(e.id)}>
                  <span
                    class={`${classes.dot} ${kind?.class ?? classes.other}`}
                    data-tip={kind?.label ?? e.type}
                    aria-label={kind?.label ?? e.type}
                  />
                  <span class={classes.text}>
                    <span class={classes.name}>{e.name}</span>
                    <span class={classes.meta}>
                      {kind?.label ?? e.type} · {chains(e.chains)} · {e.atoms.toLocaleString('en')} atoms
                    </span>
                  </span>
                  <Icon name="chevron-right" size={14} class={classes.chevron} />
                </button>
              );
            }}
          </For>
        </div>
        <Show when={derive(all.asReadonly, (a) => !a)}>
          {() => (
            <button type="button" class={classes.more} onClick={() => all.set(true)}>
              {`Show all ${summary.entities.length}`}
              <span class={classes.moreCount}>{`${summary.entities.length - PageSize} more`}</span>
            </button>
          )}
        </Show>
      </div>
    </section>
  );
}
