#if defined(UNIFORM_COLOR)
    vec3 material = color;
#elif defined(ATTRIBUTE_COLOR) || defined(INSTANCE_COLOR) || defined(ELEMENT_COLOR) || defined(ELEMENT_INSTANCE_COLOR)
    vec3 material = vColor;
#endif