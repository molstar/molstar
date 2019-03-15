#if defined(dColorType_uniform)
    vec4 material = vec4(uColor, uAlpha);
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance)
    vec4 material = vec4(vColor.rgb, uAlpha);
#elif defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
    vec4 material = uPickable == 1 ? vColor : vec4(0.0, 0.0, 0.0, 1.0); // set to empty picking id
#endif

// mix material with overpaint
#if defined(dOverpaint) && (defined(dColorType_uniform) || defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance))
    material.rgb = mix(material.rgb, vOverpaint.rgb, vOverpaint.a);
#endif

// apply transparency
#if defined(dTransparency) && (defined(dColorType_uniform) || defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance))
    float ma = material.a * (1.0 - vTransparency);
    ivec2 pixelCoord = ivec2(gl_FragCoord.xy);

    // const mat4 thresholdMatrix = mat4(
    //     1.0 / 17.0,  9.0 / 17.0,  3.0 / 17.0, 11.0 / 17.0,
    //     13.0 / 17.0,  5.0 / 17.0, 15.0 / 17.0,  7.0 / 17.0,
    //     4.0 / 17.0, 12.0 / 17.0,  2.0 / 17.0, 10.0 / 17.0,
    //     16.0 / 17.0,  8.0 / 17.0, 14.0 / 17.0,  6.0 / 17.0
    // );
    // float at = thresholdMatrix[pixelCoord.x % 4][pixelCoord.y % 4];

    // https://research.nvidia.com/publication/hashed-alpha-testing
    // Find the discretized derivatives of our coordinates
    float maxDeriv = max(length(dFdx(vViewPosition)), length(dFdy(vViewPosition)));
    float pixScale = 1.0 / maxDeriv;
    // Find two nearest log-discretized noise scales
    vec2 pixScales = vec2(exp2(floor(log2(pixScale))), exp2(ceil(log2(pixScale))));
    // Compute alpha thresholds at our two noise scales
    vec2 alpha = vec2(hash3d(floor(pixScales.x * vViewPosition)), hash3d(floor(pixScales.y * vViewPosition)));
    // Factor to interpolate lerp with
    float lerpFactor = fract(log2(pixScale));
    // Interpolate alpha threshold from noise at two scales
    float x = (1.0 - lerpFactor) * alpha.x + lerpFactor * alpha.y;
    // Pass into CDF to compute uniformly distrib threshold
    float a = min(lerpFactor, 1.0 - lerpFactor);
    vec3 cases = vec3(
        x * x / (2.0 * a * (1.0 - a)),
        (x - 0.5 * a) / (1.0 - a),
        1.0 - ((1.0 - x) * (1.0 - x) / (2.0 * a * (1.0 - a)))
    );
    // Find our final, uniformly distributed alpha threshold
    float at = (x < (1.0 - a)) ? ((x < a) ? cases.x : cases.y) : cases.z;
    // Avoids ατ == 0. Could also do
    at = clamp(at, 1.0e-6, 1.0);

    if (ma < 0.99 && (ma < 0.01 || ma < at)) discard;
#endif