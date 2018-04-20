#if defined(ATTRIBUTE_COLOR)
    vColor = color;
#elif defined(INSTANCE_COLOR)
    vColor = read_vec3(colorTex, instanceId, colorTexSize);
#elif defined(ELEMENT_COLOR)
    vColor = read_vec3(colorTex, elementId, colorTexSize);
#elif defined(ELEMENT_INSTANCE_COLOR)
    vColor = read_vec3(colorTex, instanceId * float(elementCount) + elementId, colorTexSize);
#endif