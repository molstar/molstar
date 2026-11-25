## Seq* (next): custom navigators

- Seq*: navigator pictures are registry entries: a host or an extension draws its own navigator picture
  (`content: '<id>:<name>'`, its params plain), the built-ins written to the same contract; a reader without it draws
  the axis and says so (`seq.view.navigator.picture`). A host can draw a navigator itself through `NavigatorModel`
  (`xRange`, `positionAt`).
- Headless export loads the navigators' drawing only when a figure includes a navigator.
