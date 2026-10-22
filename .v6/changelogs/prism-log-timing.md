## Prism: the Log

- The Log says how long each of a person's runs took (beside the message; the canvas's status line does not), and lists
  the plugin's transforms with their time, as v5's log did, under their own filter (within the Log's one bound since
  perf audit 2: `perf-audit2-prism-ui.md`)
- The Log filters: All, Messages, Transforms, Problems (icons, in one line with the text filter, kept in view while the
  entries scroll), and by text
- `NotificationsManager.post(text, kind, { detail })` and `log(text, kind, source, { detail })`: a detail for the log
  only, and log-only entries by source, written once a task; the Log reads `logged`, and `history` stays the messages
  the status line showed
