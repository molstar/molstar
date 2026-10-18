/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** The documentation versions, newest first. Each one lives in `src/content/docs/<id>/` and is served at `/docs/<id>/`. */
export const VERSIONS = [{ id: 'v6' }, { id: 'v5' }] as const;

export type VersionId = (typeof VERSIONS)[number]['id'];

/**
 * The version `/docs/` redirects to, and the one the other versions point back to. When v6 is released, change this to
 * `'v6'`: the redirect, the version menu and the banners all follow.
 */
export const DEFAULT_VERSION: VersionId = 'v5';

/** The version a content collection id (`v5/plugin/instance`) belongs to; undefined for the pages outside the versions. */
export function versionOf(id: string): VersionId | undefined {
  const head = id.split('/')[0];
  return VERSIONS.find((v) => v.id === head)?.id;
}

/** `newer`: ahead of the default (a preview), `older`: behind it, `default`: the version `/docs/` leads to. */
export function relationToDefault(version: VersionId): 'newer' | 'older' | 'default' {
  if (version === DEFAULT_VERSION) return 'default';
  const index = (id: VersionId) => VERSIONS.findIndex((v) => v.id === id);
  return index(version) < index(DEFAULT_VERSION) ? 'newer' : 'older';
}

export function versionLabel(version: VersionId): string {
  switch (relationToDefault(version)) {
    case 'default':
      return `${version} (latest)`;
    case 'newer':
      return `${version} (preview)`;
    default:
      return version;
  }
}
