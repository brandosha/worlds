struct Sphere {
  vec3 center;
  float radius;
};

const vec3 planetCenter = vec3(0.0, 0.0, 0.0);

const float planetRadius = 100.0;
const float atmosphereRadius = 20.0;

uniform Sphere sun;
const Sphere planet = Sphere(planetCenter, planetRadius);
const Sphere atmosphere = Sphere(planetCenter, planetRadius + atmosphereRadius);

uniform vec3 cameraPos;
uniform mat3 cameraRotation;

// parameters

const float convergedStepSize = 0.001;
const float firstStepSize = 0.5;
const float atmosphereSteps = 15.0;

const float maxDist = 1000.0;

const float scatteringStrength = 0.1;
const float brightness = 1.0;
const float densityFalloff = 1.0;
const float maxDensity = 0.7;

const vec3 wavelengths = vec3(700, 530, 460);
const vec3 scatteringCoefficients = vec3(
  pow(400.0 / wavelengths.r, 4.0),
  pow(400.0 / wavelengths.g, 4.0),
  pow(400.0 / wavelengths.b, 4.0)
) * scatteringStrength;

float sqrSunRadius = sun.radius * sun.radius;
float sqrPlanetRadius = planetRadius * planetRadius;
float sqrAtmosphereRadius = atmosphere.radius * atmosphere.radius;

const int randomFrequencyCount = 5;
uniform vec2 randomFrequencies[randomFrequencyCount];

const float maxFloat = 3.402823466e+38;

struct Ray {
  vec3 origin;
  vec3 direction;
};

// Returns vector (dstToSphere, dstThroughSphere)
// If ray origin is inside sphere, dstToSphere = 0
// If ray misses sphere, dstToSphere = maxValue; dstThroughSphere = 0
vec2 raySphere(Sphere sphere, Ray ray) {
  vec3 offset = ray.origin - sphere.center;
  float a = 1.0; // Set to dot(rayDir, rayDir) if rayDir might not be normalized
  float b = 2.0 * dot(offset, ray.direction);
  float c = dot(offset, offset) - sphere.radius * sphere.radius;
  float d = b * b - 4.0 * a * c; // Discriminant from quadratic formula

  // Number of intersections: 0 when d < 0; 1 when d = 0; 2 when d > 0
  if (d > 0.0) {
    float s = sqrt(d);
    float dstToSphereNear = max(0.0, (-b - s) / (2.0 * a));
    float dstToSphereFar = (-b + s) / (2.0 * a);

    // Ignore intersections that occur behind the ray
    if (dstToSphereFar >= 0.0) {
      return vec2(dstToSphereNear, dstToSphereFar - dstToSphereNear);
    }
  }
  // Ray did not intersect sphere
  return vec2(maxFloat, 0.0);
}

// Simplex 2D noise
//
vec3 permute(vec3 x) { return mod(((x*34.0)+1.0)*x, 289.0); }

float snoise(vec2 v){
  const vec4 C = vec4(0.211324865405187, 0.366025403784439,
           -0.577350269189626, 0.024390243902439);
  vec2 i  = floor(v + dot(v, C.yy) );
  vec2 x0 = v -   i + dot(i, C.xx);
  vec2 i1;
  i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
  vec4 x12 = x0.xyxy + C.xxzz;
  x12.xy -= i1;
  i = mod(i, 289.0);
  vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 ))
  + i.x + vec3(0.0, i1.x, 1.0 ));
  vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy),
    dot(x12.zw,x12.zw)), 0.0);
  m = m*m ;
  m = m*m ;
  vec3 x = 2.0 * fract(p * C.www) - 1.0;
  vec3 h = abs(x) - 0.5;
  vec3 ox = floor(x + 0.5);
  vec3 a0 = x - ox;
  m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );
  vec3 g;
  g.x  = a0.x  * x0.x  + h.x  * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 130.0 * dot(m, g);
}

float sdfSphere(vec3 p, Sphere sphere) {
  return length(p - sphere.center) - sphere.radius;
}
float sdfPlanet(vec3 p, float eyeDist) {
  vec3 surfacePoint = acos(normalize(p));

  float sphereDist = length(p) - planetRadius;
  // float displacement = sin(surfacePoint.x * 20.0) * sin(surfacePoint.y * 20.0) * 3.0;
  // float displacement = snoise(surfacePoint.xy * 8.0) * 0.5;
  float displacement = 0.0;
  if (sphereDist < 3.0) {
    vec2 periodic = sin(surfacePoint.xy * 30.0) + sin(surfacePoint.xy * (20.0 + randomFrequencies[0].x));
    if (periodic.x * periodic.y > 0.7) {
      displacement = -3.0;
    }
  }
  

  return sphereDist + displacement;
}

float sqrDist(vec3 from, vec3 to) {
  vec3 diff = from - to;
  return dot(diff, diff);
}

float atmosphericDensity(float sqrDistance) {
  float height = sqrt(sqrDistance) - planetRadius;
  float height01 = height / atmosphereRadius;
  return exp(-height01 * densityFalloff) * maxDensity;
}
float atmosphericDensityHeight(float height) {
  float height01 = height / atmosphereRadius;
  return exp(-height01 * densityFalloff) * maxDensity;
}

float sunRayOpticalDepthAtPoint0(vec3 origin, float eyeDist) {
  float densitySum = 0.0;
  float res = 0.0;

  vec3 sunDirection = normalize(sun.center - origin);
  vec2 atmosphereIntersection = raySphere(atmosphere, Ray(origin, sunDirection));

  float atmosphereStepSize = atmosphereIntersection.y / atmosphereSteps;
  float rayLength = firstStepSize;
  float nextAtmosphereStep = rayLength;

  for (int i = 0; i < 50; i++) {
    vec3 point = origin + sunDirection * rayLength;

    float sqrAtmosphereDist = sqrDist(point, atmosphere.center);
    if (sqrAtmosphereDist > sqrAtmosphereRadius) {
      break;
    }

    if (rayLength > nextAtmosphereStep - 0.01) {
      densitySum += atmosphericDensity(sqrAtmosphereDist);
      nextAtmosphereStep += atmosphereStepSize;
    }

    float planetDist = sdfPlanet(point, eyeDist + rayLength);
    if (planetDist < convergedStepSize) {
      return maxFloat;
    }
    rayLength += min(planetDist, nextAtmosphereStep - rayLength);
  }

  return densitySum * atmosphereStepSize;
}

float sunRayOpticalDepthAtPoint(vec3 origin, float eyeDist) {
  float opticalDepth = 0.0;

  vec3 sunDirection = normalize(sun.center - origin);
  
  float rayLength = firstStepSize;
  float stepSize = 0.5;

  for (int i = 0; i < 30; i++) {
    vec3 point = origin + sunDirection * rayLength;

    float sqrAtmosphereDist = sqrDist(point, atmosphere.center);
    if (sqrAtmosphereDist > sqrAtmosphereRadius) {
      break;
    }

    if (sdfPlanet(point, eyeDist + rayLength) < 0.0) {
      return maxFloat;
    }

    float density = atmosphericDensity(sqrAtmosphereDist);
    opticalDepth += density * stepSize;

    stepSize *= 1.05;
    rayLength += stepSize;
  }

  return opticalDepth;
}

Ray rayForPixel(vec2 coord) {
  vec2 p = coord * 0.4;
  p.y = -p.y;

  vec3 dir = vec3(p, -1.0);
  dir = cameraRotation * dir;
  
  return Ray(cameraPos, normalize(dir));
}

vec4 pixel1(vec2 coord, vec2 gridCoord) {
  vec3 color = vec3(0.0);
  bool hit = false;

  Ray ray = rayForPixel(gridCoord);

  float rayLength = firstStepSize;
  float stepSize = 0.05;

  float viewRayOpticalDepth = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  for (int i = 0; i < 100; i++) {
    vec3 point = cameraPos + ray.direction * rayLength;

    float planetDist = sdfPlanet(point, rayLength);
    if (planetDist < 0.0) {
      hit = true;

      color = vec3(0, 0.43, 0.04);

      break;
    }

    float sqrAtmosphereDist = sqrDist(point, atmosphere.center);
    if (sqrAtmosphereDist > sqrAtmosphereRadius) {
      break;
    }

    float density = atmosphericDensity(sqrAtmosphereDist);
    viewRayOpticalDepth += density * stepSize;
    float sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

    vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth) * scatteringCoefficients);
    inScatteredLight += transmittance * density;

    rayLength += stepSize;
    stepSize *= 1.06;
  }

  inScatteredLight *= scatteringCoefficients * brightness;

  if (!hit) {
    vec2 sunIntersection = raySphere(sun, ray);
    if (sunIntersection.x < maxFloat) {
      color = vec3(1, 0.85, 0);
    }
  }

  return vec4(color + inScatteredLight, 1.0);
}

vec4 pixel0(vec2 coord, vec2 gridCoord) {
  Ray ray = rayForPixel(gridCoord);

  vec3 color = vec3(0.0);
  bool hit = false;

  vec2 atmoshpereIntersection = raySphere(atmosphere, ray);
  float nextAtmoshpereStep = atmoshpereIntersection.x;
  float atmosphereStepSize = atmoshpereIntersection.y / atmosphereSteps;
  float lastAtmosphereStep = atmoshpereIntersection.x + atmoshpereIntersection.y;

  float viewRayOpticalDepth = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  float rayLength = firstStepSize;

  for (int i = 0; i < 200; i++) {
    vec3 point = cameraPos + ray.direction * rayLength;

    float sqrAtmosphereDist = sqrDist(point, planetCenter);

    float planetDist = sdfPlanet(point, rayLength);
    if (planetDist < convergedStepSize) {
      float density = atmosphericDensity(sqrAtmosphereDist);
      viewRayOpticalDepth += density * (nextAtmoshpereStep - rayLength);
      float sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

      vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth) * scatteringCoefficients);
      inScatteredLight += density * transmittance;

      color = vec3(0, 0.41, 0.06) * transmittance;
      hit = true;
      break;
    }

    float nextStepSize = planetDist;

    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      if (rayLength > nextAtmoshpereStep - 0.1) {
        float density = atmosphericDensity(sqrAtmosphereDist);
        viewRayOpticalDepth += density * atmosphereStepSize;
        float sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

        vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth) * scatteringCoefficients);
        inScatteredLight += density * transmittance;
        
        nextAtmoshpereStep += atmosphereStepSize;
      }

      nextStepSize = min(nextStepSize, nextAtmoshpereStep - rayLength);
    }

    rayLength += nextStepSize;

    if (rayLength > maxDist) {
      break;
    }
  }

  inScatteredLight *= scatteringCoefficients * brightness * atmosphereStepSize;

  if (!hit) {
    vec2 sunIntersection = raySphere(sun, ray);
    if (sunIntersection.x < maxFloat) {
      color = vec3(1, 0.85, 0);
    }
  }

  return vec4(color + inScatteredLight, 1.0);
}

vec4 pixel(vec2 coord, vec2 gridCoord) {
  Ray ray = rayForPixel(gridCoord);

  vec3 color = vec3(0.0);
  bool hit = false;

  vec2 atmoshpereIntersection = raySphere(atmosphere, ray);
  float nextAtmoshpereStep = atmoshpereIntersection.x;
  float atmosphereStepSize = atmoshpereIntersection.y / atmosphereSteps;
  float lastAtmosphereStep = atmoshpereIntersection.x + atmoshpereIntersection.y;

  float viewRayOpticalDepth = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  float rayLength = firstStepSize;
  float stepSize = 0.05;

  for (int i = 0; i < 200; i++) {
    vec3 point = cameraPos + ray.direction * rayLength;

    float sqrAtmosphereDist = sqrDist(point, planetCenter);

    float planetDist = sdfPlanet(point, rayLength);
    if (planetDist < 0.0) {
      float density = atmosphericDensity(sqrAtmosphereDist);
      viewRayOpticalDepth += density * (nextAtmoshpereStep - rayLength);
      float sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

      vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth) * scatteringCoefficients);
      inScatteredLight += density * transmittance;

      color = vec3(0, 0.41, 0.06) * (transmittance + 0.2);
      hit = true;
      break;
    }

    stepSize *= 1.01;
    float nextStepSize = stepSize;

    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      if (rayLength > nextAtmoshpereStep - 0.01) {
        float density = atmosphericDensity(sqrAtmosphereDist);
        viewRayOpticalDepth += density * atmosphereStepSize;
        float sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

        vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth) * scatteringCoefficients);
        inScatteredLight += density * transmittance;
        
        nextAtmoshpereStep += atmosphereStepSize;
      }

      nextStepSize = min(nextStepSize, nextAtmoshpereStep - rayLength);
    }

    rayLength += nextStepSize;

    // if (rayLength > maxDist) {
    //   break;
    // }
  }

  inScatteredLight *= scatteringCoefficients * brightness * atmosphereStepSize;

  if (!hit) {
    vec2 sunIntersection = raySphere(sun, ray);
    if (sunIntersection.x < maxFloat) {
      color = vec3(1, 0.85, 0);
    }
  }

  return vec4(color + inScatteredLight, 1.0);
}