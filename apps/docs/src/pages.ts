/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getCollection } from 'astro:content';
import type { VersionId } from './versions';

/** The URL of a content collection id (a page, or a version's index) under the site's base. */
export function urlOf(id: string): string {
  const base = import.meta.env.BASE_URL.replace(/\/$/, '');
  return `${base}/${id}/`;
}

let ids: Promise<Set<string>> | undefined;
function docIds() {
  return (ids ??= getCollection('docs').then((docs) => new Set(docs.map((d) => d.id))));
}

/** The same page in another version, or that version's start page if it has no such page (yet). */
export async function counterpartUrl(id: string, target: VersionId): Promise<string> {
  const rest = id.split('/').slice(1).join('/');
  const candidate = rest ? `${target}/${rest}` : target;
  return urlOf((await docIds()).has(candidate) ? candidate : target);
}
