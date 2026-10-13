# Mol\* documentation

The documentation site, built with [Astro](https://astro.build) and its [Starlight](https://starlight.astro.build)
theme. It is served at `https://molstar.org/docs/` and holds one folder per documentation version:

| Folder                          | URL               | What                                          |
| ------------------------------- | ----------------- | --------------------------------------------- |
| `src/content/docs/v5/`          | `/docs/v5/`       | the v5 documentation (the former MkDocs site) |
| `src/content/docs/v6/`          | `/docs/v6/`       | the v6 documentation, in progress             |
| `src/content/docs/versions.mdx` | `/docs/versions/` | the list of versions                          |

`/docs/` redirects to the default version, and the old v5 page URLs (`/docs/plugin/instance/`, ...) redirect to their
`/docs/v5/` counterparts.

```sh
pnpm --filter @molstar/docs dev      # http://localhost:4321/docs/
pnpm --filter @molstar/docs build    # static site in dist/
```

Contributions to the documentation are highly welcome! Please make a pull request with your changes.

## Writing pages

- A page is a Markdown file with a `title` in its frontmatter (the page's heading, so the body starts below it). Math is
  written as `$$ ... $$` and rendered with KaTeX.
- Link to other pages with site-absolute URLs that end in a slash (`/docs/v5/plugin/instance/`); relative `.md` links
  are not rewritten.
- A new page also goes into the `sidebar` of its version in `astro.config.ts`.
- Pages that exist in both versions under the same path are linked to each other by the version menu; other pages link
  to the other version's start page.

## Releasing v6

Set `DEFAULT_VERSION` in `src/versions.ts` to `'v6'`. `/docs/` then redirects to v6, the version menu marks v6 as the
latest, and the v5 pages get the "older version" banner. The v5 URLs do not change.
