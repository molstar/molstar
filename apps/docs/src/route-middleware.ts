/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineRouteMiddleware } from '@astrojs/starlight/route-data';
import { counterpartUrl } from './pages';
import { DEFAULT_VERSION, relationToDefault, versionOf } from './versions';

type Route = Parameters<Parameters<typeof defineRouteMiddleware>[0]>[0]['locals']['starlightRoute'];
type Entry = Route['sidebar'][number];
type Link = Extract<Entry, { type: 'link' }>;

function links(entries: Entry[]): Link[] {
  return entries.flatMap((e) => (e.type === 'link' ? [e] : links(e.entries)));
}

/**
 * The sidebar config holds one group per version; a page shows only the group of its own version, and the previous and
 * next links stay inside it. Pages outside the versions (the version list) have no sidebar.
 */
export const onRequest = defineRouteMiddleware(async (context) => {
  const route = context.locals.starlightRoute;
  const version = versionOf(route.id);
  if (!version) return;

  const group = route.sidebar.find((e) => e.type === 'group' && e.label === version);
  const entries = group?.type === 'group' ? group.entries : [];
  route.sidebar = entries;
  const all = links(entries);
  const index = all.findIndex((l) => l.isCurrent);
  route.pagination = {
    prev: index > 0 ? all[index - 1] : undefined,
    next: index >= 0 ? all[index + 1] : undefined,
  };

  const relation = relationToDefault(version);
  if (relation === 'default') return;
  const href = await counterpartUrl(route.id, DEFAULT_VERSION);
  route.entry.data.banner = {
    content:
      relation === 'newer'
        ? `This is the <strong>${version}</strong> documentation, still in development. <a href="${href}">Switch to the ${DEFAULT_VERSION} documentation</a>.`
        : `This is the documentation of an older version (<strong>${version}</strong>). <a href="${href}">Switch to the latest, ${DEFAULT_VERSION}</a>.`,
  };
});
