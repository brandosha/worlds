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

uniform sampler2D blueNoise;

// parameters

const float convergedStepSize = 0.001;
const float firstStepSize = 0.05;
const float atmosphereSteps = 10.0;

const float maxDist = 1000.0;

const float scatteringStrength = 0.15;
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

Ray rayForPixel(vec2 coord) {
  vec2 p = coord * 0.4;
  p.y = -p.y;

  vec3 dir = vec3(p, -1.0);
  dir = cameraRotation * dir;
  
  return Ray(cameraPos, normalize(dir));
}

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

float sqrDist(vec3 from, vec3 to) {
  vec3 diff = from - to;
  return dot(diff, diff);
}

float atmosphericDensity(float sqrDistance) {
  float height = sqrt(sqrDistance) - planetRadius;
  float height01 = height / atmosphereRadius;
  return exp(-height01 * densityFalloff) * maxDensity;
}

// https://iquilezles.org/articles/smin
float smin( float a, float b, float k ) {
  float h = max(k-abs(a-b),0.0);
  return min(a, b) - h*h*0.25/k;
}

float sdfSphere(vec3 p, Sphere sphere) {
  return length(p - sphere.center) - sphere.radius;
}

vec2 sphericalCoordinates(vec3 pos) {
  vec2 uv = vec2(acos(pos.y), atan(pos.z, pos.x));
  return uv;
}
vec3 sphericalToCartesian(vec2 uv) {
  float theta = uv.x;
  float phi = uv.y;
  return vec3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
}

float sdfPlanet(vec3 p, float eyeDist) {
  float d = length(p) - planetRadius;
  // vec3 normal = normalize(p);

  const vec3 tilePattern = vec3(1.0 / 8.0);
  vec3 ballPos = (floor(p * tilePattern) + 0.5) / tilePattern;

  float ballSurfaceDist = length(ballPos) - planetRadius;
  if (ballSurfaceDist < 3.0) {
    float radius = (sin(ballPos.x * randomFrequencies[0].x) + sin(ballPos.x * randomFrequencies[1].x) + sin(ballPos.y * randomFrequencies[0].y) + sin(ballPos.y * randomFrequencies[1].y)) / 4.0;
    radius = radius * radius * 4.0;
    // radius = min(max(radius, ballSurfaceDist), 4.0);
    Sphere ball = Sphere(ballPos, radius);

    float ballDist = sdfSphere(p, ball);
    d = smin(d, ballDist, 0.75);
    // if (ballDist < d) {
    //   normal = normalize(p - ballPos);
    // }
  }

  return d;

  /*float sphereDist = length(p) - planetRadius;

  if (sphereDist > 3.1) {
    return sphereDist - 3.0;
  } else {
    // vec3 normal = normalize(p);
    float d = sphereDist;

    const vec3 tilePattern = vec3(1.0 / 8.0);
    vec3 ballPos = (floor(p * tilePattern) + 0.5) / tilePattern;
    float radius = (sin(ballPos.x * randomFrequencies[0].x) + sin(ballPos.x * randomFrequencies[1].x) + sin(ballPos.y * randomFrequencies[0].y) + sin(ballPos.y * randomFrequencies[1].y)) / 4.0;
    Sphere ball = Sphere(ballPos, radius * radius * 4.0);

    float ballDist = sdfSphere(p, ball);
    if (ballDist < sphereDist) {
      d = ballDist;
      // normal = normalize(p - ballPos);
    }

    return d;
  }

  return min(sphereDist, sdfSphere(p, Sphere(vec3(0, 100.5, -4), 0.5)));

  if (sphereDist > 1.0) {
    return sphereDist;
  }

  vec3 surfacePoint = normalize(p);

  vec2 uv = sphericalCoordinates(surfacePoint);
  const float uvGridSize = 3.14159265 / 32.0;
  vec2 uvGrid = floor(uv / uvGridSize) * uvGridSize;

  return min(sphereDist,
              min(
                min(
                  sdfSphere(p, Sphere(planetRadius * sphericalToCartesian(uvGrid), 0.5)),
                  sdfSphere(p, Sphere(planetRadius * sphericalToCartesian(uvGrid + vec2(uvGridSize, 0.0)), 0.5))
                ),
                min(
                  sdfSphere(p, Sphere(planetRadius * sphericalToCartesian(uvGrid + vec2(0.0, uvGridSize)), 0.5)),
                  sdfSphere(p, Sphere(planetRadius * sphericalToCartesian(uvGrid + vec2(uvGridSize)), 0.5))
                )
              )
            );*/
}

vec4 sdfPlanetN(vec3 p, float eyeDist) {
  float d = length(p) - planetRadius;
  vec3 normal = normalize(p);

  const vec3 tilePattern = vec3(1.0 / 8.0);
  vec3 ballPos = (floor(p * tilePattern) + 0.5) / tilePattern;

  float ballSurfaceDist = length(ballPos) - planetRadius;
  if (ballSurfaceDist < 3.0) {
    float radius = (sin(ballPos.x * randomFrequencies[0].x) + sin(ballPos.x * randomFrequencies[1].x) + sin(ballPos.y * randomFrequencies[0].y) + sin(ballPos.y * randomFrequencies[1].y)) / 4.0;
    radius = radius * radius * 4.0;
    // radius = min(max(radius, ballSurfaceDist), 4.0);
    Sphere ball = Sphere(ballPos, radius);

    float ballDist = sdfSphere(p, ball);
    d = smin(d, ballDist, 0.75);
    if (ballDist < d) {
      normal = normalize(p - ballPos);
    }
  }

  return vec4(normal, d);

  /*float sphereDist = length(p) - planetRadius;

  if (sphereDist > 3.1) {
    return vec4(normalize(p), sphereDist - 3.0);
  } else {
    vec3 normal = normalize(p);
    float d = sphereDist;

    const vec3 tilePattern = vec3(1.0 / 8.0);
    vec3 ballPos = (floor(p * tilePattern) + 0.5) / tilePattern;
    float radius = (sin(ballPos.x * randomFrequencies[0].x) + sin(ballPos.x * randomFrequencies[1].x) + sin(ballPos.y * randomFrequencies[0].y) + sin(ballPos.y * randomFrequencies[1].y)) / 4.0;
    Sphere ball = Sphere(ballPos, radius * radius * 4.0);

    float ballDist = sdfSphere(p, ball);
    if (ballDist < sphereDist) {
      d = ballDist;
      normal = normalize(p - ballPos);
    }
    
    return vec4(normal, d);
  }*/
}

// Returns vec2(opticalDepth, shadow)
vec2 sunRayOpticalDepthAtPoint(vec3 origin, float eyeDist) {
  float densitySum = 0.0;
  float res = 1.0;

  vec3 sunDirection = normalize(sun.center - origin);
  Ray sunRay = Ray(origin, sunDirection);
  vec2 atmosphereIntersection = raySphere(atmosphere, sunRay);
  vec2 basePlanetIntersection = raySphere(planet, sunRay);
  if (basePlanetIntersection.x < maxFloat) {
    return vec2(maxFloat, 0.0);
  }

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
      return vec2(maxFloat, 0.0);
    }

    res = min(res, planetDist / rayLength);

    rayLength += min(planetDist, nextAtmosphereStep - rayLength);
  }

  return vec2(densitySum * atmosphereStepSize, clamp(res * 16.0, 0.0, 1.0));
}

vec4 pixel(vec2 coord, vec2 gridCoord) {
  Ray ray = rayForPixel(gridCoord);

  vec3 color = vec3(0.0);
  bool reachedSpace = false;
  bool hit = false;

  vec2 atmoshpereIntersection = raySphere(atmosphere, ray);
  vec2 basePlanetIntersection = raySphere(planet, ray);
  float distanceThroughAtmosphere = min(atmoshpereIntersection.y, basePlanetIntersection.x - atmoshpereIntersection.x);
  float nextAtmoshpereStep = atmoshpereIntersection.x;
  float atmosphereStepSize = distanceThroughAtmosphere / atmosphereSteps;
  float lastAtmosphereStep = atmoshpereIntersection.x + distanceThroughAtmosphere;

  float viewRayOpticalDepth = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  float rayLength = firstStepSize;

  for (int i = 0; i < 300; i++) {
    vec3 point = cameraPos + ray.direction * rayLength;

    float sqrAtmosphereDist = sqrDist(point, planetCenter);

    vec4 planetDist = sdfPlanetN(point, rayLength);
    if (planetDist.w < convergedStepSize) {
      // float density = atmosphericDensity(sqrAtmosphereDist);
      // viewRayOpticalDepth += density * (rayLength - (nextAtmoshpereStep - atmosphereStepSize));
      vec2 sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

      // vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth.x) * scatteringCoefficients);
      // inScatteredLight += density * transmittance;

      color = vec3(0, 0.35, 0.05) /* * max(vec3(0.15), transmittance) */ * sunRayOpticalDepth.y * max(0.15, dot(planetDist.xyz, normalize(sun.center - point)));
      hit = true;
      break;
    }

    float nextStepSize = planetDist.w;

    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      if (rayLength > nextAtmoshpereStep - 0.1) {
        float density = atmosphericDensity(sqrAtmosphereDist);
        viewRayOpticalDepth += density * atmosphereStepSize;
        vec2 sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

        vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth.x) * scatteringCoefficients);
        inScatteredLight += density * transmittance;
        
        nextAtmoshpereStep += atmosphereStepSize;
      }

      nextStepSize = min(nextStepSize, nextAtmoshpereStep - rayLength);
    }

    rayLength += nextStepSize;

    if (rayLength > maxDist) {
      reachedSpace = true;
      break;
    }
  }

  inScatteredLight *= scatteringCoefficients * brightness * atmosphereStepSize;

  if (reachedSpace) {
    vec2 sunIntersection = raySphere(sun, ray);
    if (sunIntersection.x < maxFloat) {
      color = vec3(1, 0.85, 0);
    }
  }

  // if (!hit && !reachedSpace) {
  //   color = vec3(1, 0, 1);
  // }

  // if (!reachedSpace) {
  //   inScatteredLight *= 0.0;
  // }

  inScatteredLight += (texture2D(blueNoise, coord / 64.0).rgb - 0.5) / 16.0;
  return vec4(color + inScatteredLight, 1.0);
}