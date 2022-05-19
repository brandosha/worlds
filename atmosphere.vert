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


// const vec3 cameraPos = vec3(0.0);
// const vec3 cameraDir = vec3(0.0, 0.0, -1.0);
uniform vec3 cameraPos;
uniform mat3 cameraRotation;
// uniform mat4 cameraTransform;

// parameters

const float convergedStepSize = 0.001;
const float minStepSize = 0.5;
const float atmosphereStepSize = atmosphereRadius / 5.0;
// const float opticalDepthStepSize = 0.1;

const float scatteringStrength = 1.0;
const float brightness = 1.0;
const float densityFalloff = 3.0;
// const float minDensity = 0.0;
const float maxDensity = 0.2;
// const float opticalDepthMultiplier = 1.0;

const vec3 wavelengths = vec3(700, 530, 460);
const vec3 scatteringCoefficients = vec3(
  pow(400.0 / wavelengths.r, 4.0),
  pow(400.0 / wavelengths.g, 4.0),
  pow(400.0 / wavelengths.b, 4.0)
) * scatteringStrength;

float sqrSunRadius = sun.radius * sun.radius;
float sqrPlanetRadius = planetRadius * planetRadius;
float sqrAtmosphereRadius = atmosphere.radius * atmosphere.radius;

float sdf_sphere(vec3 p, Sphere sphere) {
  return length(p - sphere.center) - sphere.radius;
}

float sdfPlanetSurface(vec3 p) {
  vec3 surfaceVector = normalize(planetCenter - p);
  // float displacement = 3.0 * sin(surfaceVector.z * 100.0) + 2.5 * cos(surfaceVector.x * 197.0) + 0.5 * sin(surfaceVector.y * 221.0) + 0.5 * cos(surfaceVector.z * 231.0);
  // displacement *= 0.4;
  float displacement = sin(surfaceVector.z * 11.0) + cos(surfaceVector.x * 13.0);
  return sdf_sphere(p, planet) + displacement;
}

// vec3 planetSurfaceNormal(vec3 p) {
//   vec3 surfaceVector = normalize(planetCenter - p);
//   float displacement = sin(surfaceVector.z * 80.0) + cos(surfaceVector.x * 81.0) + sin(surfaceVector.y * 82.0);
//   displacement *= 0.1;
//   return normalize(surfaceVector + vec3(displacement, displacement, displacement));
// }

// vec3 planetSurfaceNormal(vec3 p ) 
// {
//   vec2 d = vec2(0.0001,0.0);
//   return normalize( 
//     vec3(
//       sdfPlanetSurface(p+d.xyy)-sdfPlanetSurface(p-d.xyy),
//       sdfPlanetSurface(p+d.yxy)-sdfPlanetSurface(p-d.yxy),
//       sdfPlanetSurface(p+d.yyx)-sdfPlanetSurface(p-d.yyx)) 
//   );
// }

float sdf_scene(vec3 p) {
  float sdf_planet = sdf_sphere(p, planet);
  float sdf_sun = sdf_sphere(p, sun);

  return min(sdf_planet, sdf_sun);
}

struct Ray {
  vec3 origin;
  vec3 direction;
};

uniform sampler2D blueNoise;

Ray rayForPixel(vec2 coord) {
  vec2 p = coord * 0.2;
  p.y = -p.y;

  vec3 dir = vec3(p, -1.0);
  dir = cameraRotation * dir;
  
  return Ray(cameraPos, normalize(dir));
}

float sqrDist(vec3 from, vec3 to) {
  vec3 diff = from - to;
  return diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
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

float sunRayOpticalDepthAtPoint(vec3 origin) {
  float densitySum = 0.0;

  vec3 rayDirection = normalize(sun.center - origin);

  for (float m = 0.0; m < planet.radius; m += atmosphereStepSize) {
    vec3 point = origin + rayDirection * m;

    // float sqrPlanetDist = sqrDist(point, planet.center);
    if (sdfPlanetSurface(point) < 0.0) {
      return 1000.0;
    }

    float sqrAtmosphereDist = sqrDist(point, atmosphere.center);
    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      densitySum += atmosphericDensity(sqrAtmosphereDist);
    } else {
      break;
    }
  }

  return densitySum * atmosphereStepSize;
}

vec4 pixel(vec2 coord, vec2 gridCoord) {
  Ray ray = rayForPixel(gridCoord);

  vec3 surfaceColor = vec3(0.0, 0.0, 0.0);

  float viewRayOpticalDepth = 0.0;
  float sunRayOpticalDepth = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  float stepSize = minStepSize;
  float rayLength = 0.0;
  vec3 point = cameraPos;
  for (int i = 0; i < 300; i++) {
    // vec3 point = cameraPos + ray.direction * m;
    point += ray.direction * stepSize;
    rayLength += stepSize;

    if (rayLength > 1000.0) {
      break;
    }

    float sdfPlanet = sdfPlanetSurface(point); // sdf_sphere(point, planet);
    if (sdfPlanet < convergedStepSize) {
      sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point);
      vec3 sunLight = exp(-sunRayOpticalDepth * scatteringCoefficients);
      surfaceColor = vec3(0, 0.44, 0.1) * (sunLight + 0.2);
      break;
    }

    float sdfSun = sdf_sphere(point, sun);
    if (sdfSun < convergedStepSize) {
      surfaceColor = vec3(1.0, 0.8, 0.0);
      break;
    }

    float nextStepSize = min(sdfSun, sdfPlanet);

    // float sqrSunDist = sqrDist(point, sun.center);
    // if (sqrSunDist < sqrSunRadius) {
    //   surfaceColor = vec3(1.0, 0.75, 0.0);
    //   break;
    // }

    // float sqrPlanetDist = sqrDist(point, planet.center);
    // if (sqrPlanetDist < sqrPlanetRadius) {
    //   vec3 sunLight = exp(-sunRayOpticalDepth * scatteringCoefficients);
    //   surfaceColor = vec3(0, 0.44, 0.1) * (sunLight + 0.2);
    //   break;
    // }

    

    float sqrAtmosphereDist = sqrDist(point, atmosphere.center);
    if (sqrAtmosphereDist < sqrAtmosphereRadius) {
      float localDensity = atmosphericDensityHeight(sdfPlanet);
      viewRayOpticalDepth += localDensity * stepSize;
      sunRayOpticalDepth = sunRayOpticalDepthAtPoint(point);
      // float viewRayOpticalDepth = viewRayOpticalDepthSum * stepSize;
      vec3 transmittance = exp(-(viewRayOpticalDepth + sunRayOpticalDepth) * scatteringCoefficients);
      inScatteredLight += localDensity * transmittance * stepSize;

      nextStepSize = min(nextStepSize, atmosphereStepSize);
    }

    stepSize = nextStepSize;
    // stepSize = min(stepSize, atmosphereStepSize)

    // stepSize *= 1.02;
  }

  inScatteredLight *= scatteringCoefficients * brightness;
  inScatteredLight += (texture2D(blueNoise, coord / 64.0 / 2.0).rgb - 0.5) * 0.05;
  // return vec4(texture2D(blueNoise, coord / 64.0).rgb, 1.0);

  // return vec4(surfaceColor, 1.0);
  // return vec4(vec3(1.0 / rayLength), 1.0);
  return vec4(surfaceColor * (1.0 - inScatteredLight) + inScatteredLight, 1.0);
  // return vec4(surfaceColor + inScatteredLight, 1.0);
  // return vec4(screenX, screenY, 1.0 - screenX, 1.0);
}
