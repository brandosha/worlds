vec3 cameraPos = vec3(0.0);
vec3 cameraDir = vec3(0.0, 0.0, -1.0);

const float scatteringStrength = 20.0;
const vec3 wavelengths = vec3(700, 530, 460);
const vec3 scatteringCoefficients = vec3(
  pow(400.0 / wavelengths.r, 4.0),
  pow(400.0 / wavelengths.g, 4.0),
  pow(400.0 / wavelengths.b, 4.0)
) * scatteringStrength;

struct Sphere {
  vec3 center;
  float radius;
};

uniform vec3 planetCenter;

uniform Sphere sun;
Sphere planet = Sphere(planetCenter, 0.8);
Sphere atmosphere = Sphere(planetCenter, 1.0);

float sqrSunRadius = sun.radius * sun.radius;
float sqrPlanetRadius = planet.radius * planet.radius;
float sqrAtmosphereRadius = atmosphere.radius * atmosphere.radius;

struct Ray {
  vec3 origin;
  vec3 direction;
};

Ray rayForPixel(float screenX, float screenY) {
  vec2 p = vec2(screenX, 1.0 - screenY) * 2.0 - 1.0;

  vec3 dir = normalize(cameraDir + vec3(p, 0.0));
  return Ray(cameraPos, dir);
}

float sqrDist(vec3 from, vec3 to) {
  vec3 diff = from - to;
  return diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
}

float densityMultiplier = 1.0;

float atmosphericDensity(float sqrDistance) {
  float d = sqrt(sqrDistance);
  return exp(-d * densityMultiplier) * (1.0 - d);
}

float sunRayOpticalDepth(vec3 origin) {
  float opticalDepthSum = 0.0;

  vec3 rayDirection = normalize(sun.center - origin);

  const float step1 = 0.05;
  for (float m = 0.0; m < 5.0; m += step1) {
    vec3 point = origin + rayDirection * m;

    float sqrPlanetDist = sqrDist(point, planet.center);
    if (sqrPlanetDist < sqrPlanetRadius) {
      return 1000.0;
    }

    float sqrAtmosphereDist = sqrDist(point, atmosphere.center);
    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      opticalDepthSum += atmosphericDensity(sqrAtmosphereDist);
    } else {
      break;
    }
  }

  return opticalDepthSum;
}

vec3 inScatteredLightAtPoint(vec3 point, float localDensity, float viewRayOpticalDepth) {
  vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth(point)) * scatteringCoefficients);
  return localDensity * transmittance;
}

vec4 pixel(float x, float y, float screenX, float screenY) {
  Ray ray = rayForPixel(screenX, screenY);

  vec3 surfaceColor = vec3(0.0, 0.0, 0.0);

  float viewRayOpticalDepthSum = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  const float step = 0.01;
  for (float m = 0.0; m < 5.0; m += step) {
    vec3 point = cameraPos + ray.direction * m;

    float sqrSunDist = sqrDist(point, sun.center);
    if (sqrSunDist < sqrSunRadius) {
      surfaceColor = vec3(1.0, 0.75, 0.0);
      break;
    }

    float sqrPlanetDist = sqrDist(point, planet.center);
    if (sqrPlanetDist < sqrPlanetRadius) {
      surfaceColor = vec3(0, 0.44, 0.1);
      break;
    }

    float sqrAtmosphereDist = sqrDist(point, atmosphere.center);
    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      float localDensity = atmosphericDensity(sqrAtmosphereDist);
      viewRayOpticalDepthSum += localDensity;
      float viewRayOpticalDepth = viewRayOpticalDepthSum * step;

      inScatteredLight += inScatteredLightAtPoint(point, localDensity, viewRayOpticalDepth);
    }
  }

  inScatteredLight *= scatteringCoefficients * step * 15.0;


  return vec4(surfaceColor * (1.0 - inScatteredLight) + inScatteredLight, 1.0);
  // return vec4(surfaceColor + inScatteredLight, 1.0);
  // return vec4(screenX, screenY, 1.0 - screenX, 1.0);
}
