## `@molstar/ui`: a disabled toolbar button reads as disabled

- A disabled toolbar button is dimmed by the toolbar itself (opacity 0.4), so it reads as disabled in both themes; the
  `--msui-text-subtle` color alone was a few percent from an enabled button's. Seq\*'s `mountToolbar` no longer
  carries its own dimming rule.
