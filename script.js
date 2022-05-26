const canvas = document.getElementById('canvas')
// canvas.width = window.innerWidth * window.devicePixelRatio
// canvas.height = window.innerHeight * window.devicePixelRatio
canvas.style.width = canvas.width / window.devicePixelRatio + 'px'
canvas.style.height = canvas.height / window.devicePixelRatio + 'px'



const kernel = new GPU2D(canvas)
const { gl } = kernel

for (const input of document.getElementsByTagName('input')) {
  let uniformName = input.id.split('-').join('.')
  uniformName = uniformName.replace('planetcenter', 'planetCenter')
  input.addEventListener('input', e => {
    if (uniformName.endsWith('.x') || uniformName.endsWith('.y') || uniformName.endsWith('.z')) {
      let vecUniformName = uniformName.slice(0, -2)

      let arr = [0, 0, 0]
      let vecId = input.id.slice(0, -2)
      for (let i = 0; i < 3; i++) {
        let id = vecId + '-' + 'xyz'[i]
        let input = document.getElementById(id)
        arr[i] = parseFloat(input.value)
      }

      kernel.setUniformValue(vecUniformName, (gl, l) => gl.uniform3fv(l, arr))
    } else {
      kernel.setUniformValue(uniformName, (gl, l) => gl.uniform1f(l, +e.target.value))
    }
  })
}

const blueNoiseTexture = gl.createTexture()
const blueNoise = new Image()

gl.bindTexture(gl.TEXTURE_2D, blueNoiseTexture)
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT)
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT)
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST)
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST)
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array([255, 255, 255, 255]))

blueNoise.onload = () => {
  gl.bindTexture(gl.TEXTURE_2D, blueNoiseTexture)
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, blueNoise)
}
blueNoise.src = 'textures/blue-noise-rgba.png'

fetch('shaders/raymarching-sdf.frag').then(response => response.text()).then(code => {
  const program = kernel.create2dProgram(code)
  kernel.useProgram(program)

  // kernel.setUniformValue('sun.center', (gl, l) => gl.uniform3f(l, 0, 1.5, -3))
  kernel.setUniformValue('sun.radius', (gl, l) => gl.uniform1f(l, 10))

  const randomFrequencies = Array.from({ length: 2 * 5 }, (_, i) => Math.random() + Math.floor(i / 2))
  console.log(randomFrequencies)
  kernel.setUniformValue('randomFrequencies', (gl, l) => gl.uniform2fv(l, randomFrequencies))
  // kernel.setUniformValue('planetCenter', (gl, l) => gl.uniform3f(l, 0, -0.801, 0))

  // let cameraRotation = rotationMatrix(0, 0.1, 0)
  // const cameraRotation = [
  //   1, 0, 0,
  //   0, 1, 0,
  //   0, 0, 1
  // ]
  // kernel.setUniformValue('cameraRotation', (gl, l) => gl.uniformMatrix3fv(l, false, cameraRotation))

  gl.activeTexture(gl.TEXTURE0)
  gl.bindTexture(gl.TEXTURE_2D, blueNoiseTexture)
  kernel.setUniformValue('blueNoise', (gl, l) => gl.uniform1i(l, 0))

  requestAnimationFrame(render)
})

let keysDown = new Set()
window.addEventListener('keydown', e => {
  keysDown.add(e.key)
})
window.addEventListener('keyup', e => {
  keysDown.delete(e.key)
})

let velocity = [0, 0, 0]
const acceleration = 0.005

let cameraPosition = [0, 101, 0]
// let cameraRotation = [0, 0, 0]
let rotationalVelocities = [0, 0, 0]
let cameraRotation = [
  1, 0, 0,
  0, 1, 0,
  0, 0, 1
]
const rotationalAcceleration = 0.1

let prev = 0
function render(now) {
  const delta = now - prev
  prev = now
  const fps = 1000 / delta
  document.getElementById('fps').innerText = fps.toFixed(2)

  const deltaA = acceleration * delta / 1000
  const deltaR = rotationalAcceleration * delta / 1000

  if (keysDown.has('ArrowLeft')) {
    rotationalVelocities[1] -= deltaR
    // cameraRotation[2] -= deltaA
  } else if (keysDown.has('ArrowRight')) {
    rotationalVelocities[1] += deltaR
    // cameraRotation[2] += deltaA
  }

  if (keysDown.has('ArrowUp')) {
    rotationalVelocities[0] -= deltaR
  } else if (keysDown.has('ArrowDown')) {
    rotationalVelocities[0] += deltaR
  }

  rotationalVelocities = vecScale(rotationalVelocities, 0.9)
  cameraRotation = matMul3(rotationMatrix(rotationalVelocities[0], rotationalVelocities[1], rotationalVelocities[2]), cameraRotation) // vecAdd(cameraRotation, rotationalVelocities)
  // cameraRotation[2] *= 0.9

  if (keysDown.has('w')) {
    velocity[1] += deltaA
  }
  if (keysDown.has('s')) {
    velocity[1] -= deltaA
  }
  if (keysDown.has('a')) {
    velocity[0] -= deltaA
  }
  if (keysDown.has('d')) {
    velocity[0] += deltaA
  }
  if (keysDown.has('x')) {
    velocity[2] -= deltaA
  }
  if (keysDown.has('z')) {
    velocity[2] += deltaA
  }

  velocity = vecScale(velocity, 0.99)

  // let velocityVector = [0, 0, -velocity]
  // let cameraRotationMatrix = rotationMatrix(cameraRotation[0], cameraRotation[1], cameraRotation[2])
  // cameraRotationMatrix = 
  let velocityVector = matTransform3(velocity, cameraRotation)
  cameraPosition = vecAdd(cameraPosition, vecScale(velocityVector, delta))
  
  kernel.setUniformValue('cameraPos', (gl, l) => gl.uniform3fv(l, cameraPosition))
  kernel.setUniformValue('cameraRotation', (gl, l) => gl.uniformMatrix3fv(l, false, cameraRotation))

  const t = now / 1000
  const dayLength = 120
  kernel.setUniformValue('sun.center', (gl, l) => gl.uniform3f(l, 0, 300 * cos(t / dayLength - 1.5), 300 * sin(t / dayLength - 1.5)))
  
  kernel.draw()

  requestAnimationFrame(render)
}