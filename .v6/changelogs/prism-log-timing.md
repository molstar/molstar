## Prism: the Log

- The Log says how long each of a person's runs took (beside the message; the canvas's status line does not), and lists
  the plugin's transforms with their time, as v5's log did: the latest 100, under their own filter
- The Log filters: All, Messages, Transforms, Problems, and by text
- `NotificationsManager.post(text, kind, { detail })` and `log(text, kind, source, { detail, keep })`: a detail for the
  log only, and log-only entries by source, each source with its own bound, written once a task
