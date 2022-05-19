var cos = Math.cos
var sin = Math.sin

function vecAdd(a, b) {
  return [
    a[0] + b[0],
    a[1] + b[1],
    a[2] + b[2]
  ]
}

function vecScale(vec, scale) {
  return [
    vec[0] * scale,
    vec[1] * scale,
    vec[2] * scale
  ]
}

function matTransform3(vec, mat) {
  return [
    vec[0] * mat[0] + vec[1] * mat[3] + vec[2] * mat[6],
    vec[0] * mat[1] + vec[1] * mat[4] + vec[2] * mat[7],
    vec[0] * mat[2] + vec[1] * mat[5] + vec[2] * mat[8]
  ]
}

function rotationMatrix(x, y, z) {
  let cosX = cos(x)
  let sinX = sin(x)

  let cosY = cos(y)
  let sinY = sin(y)

  let cosZ = cos(z)
  let sinZ = sin(z)

  return [
    cosY * cosZ, -cosX * sinZ + sinX * sinY * cosZ, sinX * sinZ + cosX * sinY * cosZ,
    cosY * sinZ, cosX * cosZ + sinX * sinY * sinZ, -sinX * cosZ + cosX * sinY * sinZ,
    -sinY, sinX * cosY, cosX * cosY
  ]
}

function matMul3(a, b) {
  return [
    a[0] * b[0] + a[1] * b[3] + a[2] * b[6],
    a[0] * b[1] + a[1] * b[4] + a[2] * b[7],
    a[0] * b[2] + a[1] * b[5] + a[2] * b[8],
    a[3] * b[0] + a[4] * b[3] + a[5] * b[6],
    a[3] * b[1] + a[4] * b[4] + a[5] * b[7],
    a[3] * b[2] + a[4] * b[5] + a[5] * b[8],
    a[6] * b[0] + a[7] * b[3] + a[8] * b[6],
    a[6] * b[1] + a[7] * b[4] + a[8] * b[7],
    a[6] * b[2] + a[7] * b[5] + a[8] * b[8]
  ]
}
