## Prism

- The empty canvas's splash no longer shows while a slow load is in progress: it waits on the new `state.working`
  (an operation from its start, a lazy extension's load, the plugin updating; not debounced) instead of the plugin's
  `isBusy`, which turns on 750 ms late and misses work that is not the plugin's.
