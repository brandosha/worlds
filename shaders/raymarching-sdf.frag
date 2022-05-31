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

uniform float t;
uniform vec3 jitter;
uniform sampler2D prevFrame;

// parameters

const float convergedStepSize = 0.001;
const float firstStepSize = 0.05;
const float atmosphereSteps = 10.0;
const float opticalDepthSteps = 4.0;

const float maxDist = 1000.0;

const float scatteringStrength = 0.25;
const float brightness = 1.0;
const float densityFalloff = 1.0;
const float maxDensity = 0.25;

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

  vec2 j = jitter.xy / uViewSize;

  vec3 dir = vec3(p + j, -1.0);
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
  return vec2(acos(pos.y), atan(pos.z, pos.x));
}
vec3 sphericalToCartesian(vec2 uv) {
  float theta = uv.x;
  float phi = uv.y;
  return vec3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
}

const float sqareSubdivisions = 8.0;
vec3 discreteNormalize(vec3 p) {
  vec3 cubePoint = vec3(sqareSubdivisions);
  vec3 a = abs(p);
  if (a.x > a.y && a.x > a.z) {
    cubePoint *= p / a.x;
    cubePoint.yz = floor(cubePoint.yz) + 0.5;
  } else if (a.y > a.z) {
    cubePoint *= p / a.y;
    cubePoint.xz = floor(cubePoint.xz) + 0.5;
  } else {
    cubePoint *= p / a.z;
    cubePoint.xy = floor(cubePoint.xy) + 0.5;
  }

  return normalize(cubePoint);
}

mat4 discretNormalize4(vec3 p) {
  const vec2 v10 = vec2(1.0, 0.0);

  vec3 cubePoint = vec3(sqareSubdivisions);
  vec3 a = abs(p);
  if (a.x > a.y && a.x > a.z) {
    cubePoint *= p / a.x;
    vec3 p00 = floor(cubePoint);
    return mat4(
      vec4(p00, 0.0),
      vec4(p00 + v10.yxy, 0.0),
      vec4(p00 + v10.yyx, 0.0),
      vec4(p00 + v10.yxx, 0.0)
    );
  } else if (a.y > a.z) {
    cubePoint *= p / a.y;
    vec3 p00 = floor(cubePoint);
    return mat4(
      vec4(p00, 0.0),
      vec4(p00 + v10.xyy, 0.0),
      vec4(p00 + v10.yyx, 0.0),
      vec4(p00 + v10.xyx, 0.0)
    );
  } else {
    cubePoint *= p / a.z;
    vec3 p00 = floor(cubePoint);
    return mat4(
      vec4(p00, 0.0),
      vec4(p00 + v10.xyy, 0.0),
      vec4(p00 + v10.yxy, 0.0),
      vec4(p00 + v10.xxy, 0.0)
    );
  }
}

Sphere randomSphere(vec3 pos) {
  float radius = 1.0;
  // float height = 0.5;
  float height = 8.0 * abs((sin(pos.x * randomFrequencies[0].x + t) + sin(pos.x * randomFrequencies[1].x + t) + sin(pos.y * randomFrequencies[0].y + t) + sin(pos.y * randomFrequencies[1].y + t)) / 4.0) - 3.0;

  return Sphere(pos * (planetRadius + height), radius);
}

float sdfPlanet(vec3 p, float eyeDist) {
  float distanceFromOrigin = length(p);
  float d = length(distanceFromOrigin) - planetRadius;

  if (d > 10.0) {
    return d - 9.0;
  }

  vec3 ballPos = discreteNormalize(p);
  Sphere ball = randomSphere(ballPos);

  float ballDist = sdfSphere(p, ball);
  d = smin(d, ballDist, 1.75);

  return d;
}

vec4 sdfPlanetN(vec3 p, float eyeDist) {
  float distanceFromOrigin = length(p);
  float d = length(distanceFromOrigin) - planetRadius;
  vec3 normal = normalize(p);

  if (d > 10.0) {
    return vec4(normal, d - 9.0);
  }

  vec3 ballPos = discreteNormalize(p);
  Sphere ball = randomSphere(ballPos);

  float ballDist = sdfSphere(p, ball);
  d = smin(d, ballDist, 1.75);
  if (ballDist < d) {
    normal = normalize(p - ball.center);
  }

  return vec4(normal, d);
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

  float atmosphereStepSize = atmosphereIntersection.y / opticalDepthSteps;
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

vec4 debugSdf(Ray ray) {
  vec2 basePlanetIntersection = raySphere(planet, ray);
  if (basePlanetIntersection.x < maxFloat) {
    vec3 point = ray.origin + ray.direction * basePlanetIntersection.x;
    vec3 g = discreteNormalize(point);
    vec3 col = (g + 1.0) / 2.0;

    return vec4(col, 1.0);

    return vec4(0, 0.35, 0.05, 1);
  }

  vec2 sunIntersection = raySphere(sun, ray);
  if (sunIntersection.x < maxFloat) {
    return vec4(1, 0.85, 0, 1.0);
  }

  const float rayLength = 5.0;
  vec3 p = cameraPos + ray.direction * rayLength;
  float basePlanetDistance = sdfSphere(p, planet);
  vec4 res = sdfPlanetN(p, rayLength);

  if (res.w == basePlanetDistance) {
    return vec4(0.1, 0, 0.1, 1);
  }

  return vec4(vec3(exp(-res.w)), 1.0);
}

/*float opRep( in vec3 p, in vec3 c, in sdf3d primitive )
{
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    return primitive( q );
}
vec3 opRepLim( in vec3 p, in float c, in vec3 l, in sdf3d primitive )
{
    vec3 q = p-c*clamp(round(p/c),-l,l);
    return primitive( q );
}*/

vec3 surfaceShader(float viewRayOpticalDepth, vec2 sunRayOpticalDepth, vec4 planetDist, vec3 point) {
  vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth.x) * scatteringCoefficients);
  vec3 surfaceColor = vec3(0, 0.35, 0.05);
  vec3 color = surfaceColor * transmittance * sunRayOpticalDepth.y * dot(planetDist.xyz, normalize(sun.center - point));
  // ambient light
  color += vec3(0.16, 0.23, 0.54) * surfaceColor;
  return color;
}

vec4 pixel(vec2 coord, vec2 gridCoord, vec2 coord01) {
  Ray ray = rayForPixel(gridCoord);

  // return debugSdf(ray);

  vec3 color = vec3(0.0);
  bool reachedSpace = false;
  bool hit = false;

  vec2 atmoshpereIntersection = raySphere(atmosphere, ray);
  vec2 basePlanetIntersection = raySphere(planet, ray);
  float distanceThroughAtmosphere = min(atmoshpereIntersection.y, basePlanetIntersection.x - atmoshpereIntersection.x);
  float nextAtmosphereStep = atmoshpereIntersection.x;
  // nextAtmoshpereStep = maxFloat; // Uncomment to disable atmosphere
  float atmosphereStepSize = distanceThroughAtmosphere / atmosphereSteps;
  float lastAtmosphereStep = atmoshpereIntersection.x + distanceThroughAtmosphere;

  float viewRayOpticalDepth = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  float rayLength = firstStepSize;

  for (int i = 0; i < 300; i++) {
    vec3 point = ray.origin + ray.direction * rayLength;

    float sqrAtmosphereDist = sqrDist(point, planetCenter);

    vec4 planetDist = sdfPlanetN(point, rayLength);
    if (planetDist.w < convergedStepSize) {
      float density = atmosphericDensity(sqrAtmosphereDist);
      viewRayOpticalDepth += density * (rayLength - (nextAtmosphereStep - atmosphereStepSize));
      vec2 sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

      color = surfaceShader(viewRayOpticalDepth, sunRayOpticalDepth, planetDist, point);
      hit = true;
      break;
    }

    float nextStepSize = planetDist.w;

    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      if (rayLength > nextAtmosphereStep - 0.1) {
        float density = atmosphericDensity(sqrAtmosphereDist);
        viewRayOpticalDepth += density * atmosphereStepSize;
        vec2 sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point, rayLength);

        vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth.x) * scatteringCoefficients);
        inScatteredLight += density * transmittance;
        
        nextAtmosphereStep += atmosphereStepSize;
      }

      nextStepSize = min(nextStepSize, nextAtmosphereStep - rayLength);
    }

    rayLength += nextStepSize;

    if (rayLength > maxDist) {
      reachedSpace = true;
      break;
    }
  }

  // float atmosphereDist = min(rayLength - atmoshpereIntersection.x, distanceThroughAtmosphere);
  // if (atmosphereDist > 0.0) {
  //   inScatteredLight = exp(-atmosphereDist / 64.0 * scatteringCoefficients) * scatteringCoefficients * brightness;
  // }
  
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

  inScatteredLight += (texture2D(blueNoise, coord / 32.0).rgb - 0.5) / 16.0;
  color += inScatteredLight;

  coord01.y = 1.0 - coord01.y;
  vec3 prevPixel = texture2D(prevFrame, coord01).rgb;
  return vec4(mix(color, prevPixel, 0.25), 1.0);

  
  return vec4(color + inScatteredLight, 1.0);
}