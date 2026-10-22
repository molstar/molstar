## Fixes: volume streaming between states, the Log from an error

- Volume streaming: a state that streams other maps than the one before (an MVS state with
  `molstar_volume_streaming`, then another that streams an EM map) recreates the streaming behavior instead of failing
  with "Cannot destructure property 'emDefaultContourLevel'"
- Prism: clicking an error in the status line opens the Log at Problems; another message opens it at All. The Log's
  filter is `NotificationsManager.logFilter` (`setLogFilter`), shared by every Log view
