## MVS: fixes

- A `clip` plane whose normal points straight down (antiparallel to +y, e.g. `normal: [0, -1, 0]`) clips the lower half
  again, instead of turning into an upward plane and clipping the upper half
