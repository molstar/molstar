#if defined(dColorType_attribute)
    vColor = aColor;
#elif defined(dColorType_instance)
    vColor = read_vec3(tColor, aInstanceId, uColorTexSize);
#elif defined(dColorType_element)
    vColor = read_vec3(tColor, aElementId, uColorTexSize);
#elif defined(dColorType_elementInstance)
    vColor = read_vec3(tColor, aInstanceId * float(uElementCount) + aElementId, uColorTexSize);
#endif