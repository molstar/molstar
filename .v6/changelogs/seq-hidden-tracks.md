## Seq* (next): hidden tracks and flyout submenus

- Seq*: hiding a track adds nothing to the view (no section line or chip); a section with every track hidden
  collapses fully; the application lists and restores hidden tracks (`viewer.hiddenTracks`, `track.show` with `all`,
  the menu's "Show hidden", `mountHiddenTracks` with a count badge); a jump to an item only on hidden tracks says so; a
  ramp chosen from a layer's Color menu takes the data's extent as its domain.
- `@molstar/ui`: context menus open submenus as flyouts beside their row (`ContextMenuItem.items`), by click, keys or
  resting.
