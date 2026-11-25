## Seq* (next): track headers in the viewer

- Track headers show action slots on hover (hide, color), count chips that show hidden tracks and layers again, and a
  menu on a secondary click, the menu key or a long press; the cursor reaches headers. A host adds its own
  `header-action` entries (namespaced) or places its own DOM parts over the published geometry.
- The headers' late parts (the menu, DOM parts, icons, gestures) load on first use; `viewer.settle()` waits for them.
- The "n filtered" chip clears every filtered layer of its track in one step, with or without the headers behavior.
