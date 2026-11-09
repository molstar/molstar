## Graphics

- Direct volume `renderMode`: `composite` (as before) blends the samples of a ray front to back by the transfer
  function; `maximum` keeps the brightest sample along the ray, unlit (a maximum intensity projection), as bright as
  the transfer function makes its value. A sample below the value the transfer function first rises above zero at adds
  nothing (its bricks are skipped). On WebGL and WebGPU (#2066)
- Direct volume `renderMode` `attenuated` (the sample brightest by the transfer function once each is dimmed by the density in front of it: depth for the projection) and `additive` (the samples summed), with `attenuation` (the dimming, or the gain); the projection
  presets of the control points (linear, bright, solid) (#2066)
- Fix the direct-volume transfer function going nearly transparent right past a control point of opacity 1 (the
  spline's overshoot wrapped around in its bytes) (#2066)
