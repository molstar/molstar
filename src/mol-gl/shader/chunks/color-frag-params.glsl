#if defined(UNIFORM_COLOR)
    uniform vec3 color;
#elif defined(ATTRIBUTE_COLOR) || defined(INSTANCE_COLOR) || defined(ELEMENT_COLOR) || defined(ELEMENT_INSTANCE_COLOR)
    varying vec3 vColor;
#endif