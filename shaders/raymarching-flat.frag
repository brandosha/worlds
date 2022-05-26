uniform vec3 cameraPos;
uniform mat3 cameraRotation;

const int randomFrequencyCount = 3;
uniform vec2 randomFrequencies[randomFrequencyCount];

struct Sphere {
  vec3 center;
  float radius;
};

uniform Sphere sun;

struct Ray {
  vec3 origin;
  vec3 direction;
};

const float maxFloat = 3.402823466e+38;

const float scatteringStrength = 0.15;

const vec3 wavelengths = vec3(700, 530, 460);
const vec3 scatteringCoefficients = vec3(
  pow(400.0 / wavelengths.r, 4.0),
  pow(400.0 / wavelengths.g, 4.0),
  pow(400.0 / wavelengths.b, 4.0)
) * scatteringStrength;

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

Ray rayForPixel(vec2 coord) {
  vec2 p = coord * 0.4;
  p.y = -p.y;

  vec3 dir = vec3(p, -1.0);
  dir = cameraRotation * dir;
  
  return Ray(cameraPos, normalize(dir));
}

float heightMap(vec3 point) {
  float h = 90.0;
  vec2 pos = point.xz / 10.0 + 100.0;

  float multiplier = 1.0;
  float frequency = 1.0;
  
  for (int i = 0; i < 1; i++) {
    for (int i = 0; i < randomFrequencyCount; i++) {
      vec2 s = sin(pos * (randomFrequencies[i] * frequency) + pos.yx) * multiplier;
      h += s.x + s.y;
    }

    multiplier /= 4.0;
    frequency *= 4.0;

    if (abs(point.y - h) > multiplier * 4.0) {
      break;
    }
  }

  return h;
}

vec4 heightMapNormal(vec3 point) {
  float h = 90.0;
  vec2 pos = point.xz / 10.0 + 100.0;

  // vec3 normal = vec3(0.0);
  vec3 dx = vec3(1.0, 0.0, 0.0);
  vec3 dz = vec3(0.0, 0.0, 1.0);

  float multiplier = 1.0;
  float frequency = 1.0;
  
  for (int i = 0; i < 1; i++) {
    for (int i = 0; i < randomFrequencyCount; i++) {
      vec2 f = randomFrequencies[i] * frequency;
      vec2 s = sin(pos * f + pos.yx) * multiplier;
      h += s.x + s.y;

      dx.y += cos(pos.x * f.x + pos.y) * multiplier * f.x;
      dz.y += cos(pos.y * f.y + pos.x) * multiplier * f.y;
    }

    multiplier /= 4.0;
    frequency *= 4.0;

    if (abs(point.y - heightMap(point)) > multiplier * 2.0) {
      break;
    }
  }

  return vec4(normalize(cross(dx, dz)), h);
}

vec3 surfaceMaterial(vec3 normal, vec3 pos) {
  vec3 color = vec3(0.09, 0.43, 0);
  if (dot(normal, vec3(0.0, -1.0, 0.0)) < 0.5) {
    color = vec3(0.45);
  }

  float sunLight = max(-dot(normal, normalize(sun.center - pos)), 0.1);
  if (sun.center.y < 90.0) {
    sunLight = 0.1;
  }

  return color * sunLight;
}

vec4 pixel(vec2 coord, vec2 gridCoord) {
  vec3 color = vec3(0.0);

  Ray ray = rayForPixel(gridCoord);

  float viewRayOpticalDepth = 0.0;
  vec3 inScatteredLight = vec3(0.0);

  float stepSize = 0.01;
  float rayLength = 0.5;

  bool hit = false;

  for (int i = 0; i < 200; i++) {
    vec3 pos = cameraPos + ray.direction * rayLength;

    // float h = heightMap(pos);
    vec4 hn = heightMapNormal(pos);
    if (pos.y < hn.w) {
      color = surfaceMaterial(hn.xyz, pos);
      
      hit = true;
      break;
    }

    viewRayOpticalDepth += stepSize * exp(-pos.y);
    inScatteredLight += stepSize * exp(-viewRayOpticalDepth * scatteringCoefficients);

    rayLength += stepSize;
    stepSize *= 1.03;

    if (rayLength > 1000.0) {
      break;
    }
  }

  if (!hit) {
    vec2 sunIntersection = raySphere(sun, ray);
    if (sunIntersection.x < maxFloat) {
      color = vec3(1, 0.85, 0);
    }
  }

  inScatteredLight *= scatteringCoefficients / 128.0;

  return vec4(color + inScatteredLight, 1.0);
}