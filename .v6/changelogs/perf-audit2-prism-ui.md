## Performance: perf audit 2, Prism's UI core and the assistant UI

- Prism: the Log keeps the last 300 entries of every source together (`NotificationsManager` option `logSize`), warnings
  and errors dropped last, so a burst of transforms or plugin warnings never pushes a problem out while anything else is
  left; it replaces the bound per source (`LogOptions.keep` is gone, and with it the 100 transforms of their own). The
  entries are a virtual list of one-line rows (a long line fades, its tooltip has the rest): a flush of 100 lines with
  the Log shown 4.5 → 0.39 ms, 300 → 54 rows in the DOM; the sources it filters by are `logSources`
- Prism: a running task's pill is made once and follows its task by id (`NotificationsManager.lineTaskIds`,
  `lineTask(id)`): its spinner no longer restarts on every progress event, and Cancel is no longer replaced under a press
- Prism: a pinch (ctrl+wheel) is held over the canvas only; over a panel, a card or the rail it zooms the page. Panels
  scroll on the compositor again while a load or a rebuild holds the main thread
- Prism: `<For key={…}>` is a type error (the JSX runtime never passed `key`, so those lists were keyed by identity); pass
  `by={…}`. The 49 lists in the repository do
- Prism: `Show`, `Switch` and `Dynamic` hold no rxjs chain on a state (4.1 → 2.8 KB a `Show`); a keydown checks only the
  commands on its key (17 → 0.06 µs at 437 commands); `shallowEqual`, the default compare of every record state,
  925 → 120 ns and 656 → 32 B a compare; pointer moves are not delegated (no root walks a move's path); a tree's arrow
  keys find the active row without a scan (9,991 → 5 lookups at the end of 5,000 rows)
- Prism assistant: the Inspector's overview and trace list follow turns and steps by id: a settle while a turn runs makes
  only its new steps (495 → 30 blocks made over a 30-step turn) and reconciles only its rows (1.26 → 0.64 ms at 300
  turns); History says when an item ran as the time of day (a "N min ago" was written once and went stale)
