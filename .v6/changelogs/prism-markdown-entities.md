## Prism: markdown in descriptions

- Character references in snapshot descriptions and hover labels are the characters they name: `&#x2B05;`, `&nbsp;`,
  `&mdash;`, `&Aring;` no longer show as written (code stays literal)
- Bold and italic hold links and code: a story's `**[heme group](!query=…)**` is a command link again, not its source
- Core: `@molstar/core/util/html-entities` (`decodeHtmlEntities`), one DOM-free decoder for numeric and common named
  references, used by core's markdown parser, Prism's markdown and its hover labels
